# GTF/FASTA/SQANTI/domain-mapping ingestion logic for Database.
# Extracted out of core/database.py so Database itself stays a thin
# engine/session-lifecycle class; each function here takes the Database
# instance explicitly instead of being a bound method, and is called by a
# same-named thin wrapper method on Database.
import csv
import re
from itertools import chain, groupby
from operator import itemgetter
from typing import Callable, Iterable

from Bio import SeqIO
from biosurfer.core.alignments import CodonAlignment
from biosurfer.core.constants import (APPRIS, SQANTI, STOP_CODONS, AminoAcid,
                                      FeatureType, Strand)
from biosurfer.core.io_helpers import (FastaHeaderFields,
                                       bulk_upsert, count_lines, read_gtf_line)
from biosurfer.core.models.biomolecules import (ORF, Chromosome, Exon,
                                                GencodeTranscript, Gene,
                                                PacBioTranscript, Protein,
                                                Transcript)
from biosurfer.core.models.features import (Domain, Feature, ProjectedFeature,
                                            ProteinFeature)
from more_itertools import chunked
from sqlalchemy import delete, select
from sqlalchemy.dialects.sqlite import insert
from sqlalchemy.orm import contains_eager, joinedload, raiseload
from sqlalchemy.sql.expression import desc
from sqlalchemy.sql.functions import func
from tqdm import tqdm

CHUNK_SIZE = 5000
SQANTI_DICT = {
    'full-splice_match': SQANTI.FSM,
    'incomplete-splice_match': SQANTI.ISM,
    'novel_in_catalog': SQANTI.NIC,
    'novel_not_in_catalog': SQANTI.NNC
}


def _process_exons(transcripts_to_exons, minus_transcripts):
    t = tqdm(
        transcripts_to_exons.items(),
        desc = 'Calculating transcript-relative exon coords',
        total = len(transcripts_to_exons),
        unit = 'transcripts'
    )
    for transcript_id, exon_list in t:
        exon_list.sort(key=itemgetter('start'), reverse=transcript_id in minus_transcripts)
        tx_idx = 1
        for i, exon in enumerate(exon_list, start=1):
            exon['position'] = i
            if 'accession' not in exon:
                exon['accession'] = transcript_id + f':EXON{i}'
            exon_length = exon['stop'] - exon['start'] + 1
            exon['transcript_start'] = tx_idx
            exon['transcript_stop'] = tx_idx + exon_length - 1
            tx_idx += exon_length


def _process_orfs(transcripts_to_cdss, transcripts_to_exons, minus_transcripts):
    orfs_to_upsert = []
    t = tqdm(
        transcripts_to_cdss.items(),
        desc = 'Calculating transcript-relative ORF coords',
        total = len(transcripts_to_cdss),
        unit = 'transcripts'
    )
    for transcript_id, cds_list in t:
        exon_list = transcripts_to_exons[transcript_id]
        cds_list.sort(key=itemgetter(0), reverse=transcript_id in minus_transcripts)
        first_cds, last_cds = cds_list[0], cds_list[-1]
        orf_start = first_cds[0]
        orf_stop = last_cds[1]
        if transcript_id in minus_transcripts:
            orf_start, orf_stop = orf_stop, orf_start
        first_exon = next((exon for exon in exon_list if exon['start'] <= first_cds[0] and first_cds[1] <= exon['stop']), None)
        last_exon = next((exon for exon in reversed(exon_list) if exon['start'] <= last_cds[0] and last_cds[1] <= exon['stop']), None)
        if transcript_id not in minus_transcripts:
            first_offset = first_cds[0] - first_exon['start']
            last_offset = last_exon['stop'] - last_cds[1]
        else:
            first_offset = first_exon['stop'] - first_cds[1]
            last_offset = last_cds[0] - last_exon['start']
        orf_tx_start = first_exon['transcript_start'] + first_offset
        orf_tx_stop = last_exon['transcript_stop'] - last_offset
        orf_length = orf_tx_stop - orf_tx_start + 1
        if orf_length % 3 != 0:
            continue
        orfs_to_upsert.append({
            'transcript_id': transcript_id,
            'position': 1,
            'transcript_start': orf_tx_start,
            'transcript_stop': orf_tx_stop,
            'has_stop_codon': False,
        })
    return orfs_to_upsert


def _upsert_exons_and_orfs(session, transcripts_to_exons, transcripts_to_cdss, minus_transcripts, exons_to_upsert):
    """Shared tail logic for load_gencode_gtf/load_pacbio_gtf: compute transcript-relative
    exon/ORF coordinates and bulk-upsert both tables."""
    _process_exons(transcripts_to_exons, minus_transcripts)
    t = tqdm(
        exons_to_upsert,
        desc = 'Upserting exons',
        total = len(exons_to_upsert),
        unit = 'exons'
    )
    for exons in chunked(t, CHUNK_SIZE):
        bulk_upsert(session, Exon.__table__, exons, primary_keys=('accession', 'transcript_id'))

    orfs_to_upsert = _process_orfs(transcripts_to_cdss, transcripts_to_exons, minus_transcripts)
    t = tqdm(
        orfs_to_upsert,
        desc = 'Upserting ORFs',
        total = len(orfs_to_upsert),
        unit = 'ORFs'
    )
    for orfs in chunked(t, CHUNK_SIZE):
        bulk_upsert(session, ORF.__table__, orfs, primary_keys=('transcript_id', 'position'))


def load_gencode_gtf(db, gtf_file: str, overwrite=False) -> None:
    with db.get_session() as session:
        chromosomes = {}
        genes_to_upsert = []
        transcripts_to_upsert = []
        exons_to_upsert = []
        orfs_to_upsert = []

        minus_transcripts = set()
        transcripts_to_exons = {}
        transcripts_to_cdss = {}

        with open(gtf_file) as gtf:
            lines = count_lines(gtf)
            t = tqdm(gtf, desc=f'Reading GENCODE annotations', total=lines, unit='lines')
            for i, line in enumerate(t, start=1):
                if line.startswith("#"):
                    continue
                chr, _, feature, start, stop, _, strand, _, attributes, tags = read_gtf_line(line)
                if any('_NF' in tag for tag in tags):
                    continue
                gene_id = attributes['gene_id']
                gene_name = attributes['gene_name']
                if attributes['gene_type'] != 'protein_coding':
                    continue

                if feature == 'gene':
                    if chr not in chromosomes:
                        with session.begin():
                            chromosome = session.merge(Chromosome(name=chr))
                        chromosomes[chr] = chromosome
                    else:
                        chromosome = chromosomes[chr]
                    gene = {
                        'accession': gene_id,
                        'name': gene_name,
                        'strand': Strand.from_symbol(strand),
                        'chromosome_id': chr
                    }
                    genes_to_upsert.append(gene)

                elif feature == 'transcript':
                    transcript_id = attributes['transcript_id']
                    transcript_name = attributes['transcript_name']
                    appris = APPRIS.NONE
                    start_nf, end_nf = False, False
                    for tag in tags:
                        if 'appris_principal' in tag:
                            appris = APPRIS.PRINCIPAL
                        if 'appris_alternative' in tag:
                            appris = APPRIS.ALTERNATIVE
                        start_nf, end_nf = False, False
                    transcript = {
                        'accession': transcript_id,
                        'name': transcript_name,
                        'type': 'gencodetranscript',
                        'gene_id': gene_id,
                        'strand': Strand.from_symbol(strand),
                        'appris': appris,
                        'start_nf': start_nf,
                        'end_nf': end_nf
                    }
                    transcripts_to_upsert.append(transcript)
                    if Strand.from_symbol(strand) is Strand.MINUS:
                        minus_transcripts.add(transcript_id)

                elif feature == 'exon':
                    transcript_id = attributes['transcript_id']
                    exon_id = attributes['exon_id']
                    exon = {
                        'accession': exon_id,
                        'start': start,
                        'stop': stop,
                        'transcript_id': transcript_id
                    }
                    if transcript_id in transcripts_to_exons:
                        transcripts_to_exons[transcript_id].append(exon)
                    else:
                        transcripts_to_exons[transcript_id] = [exon,]
                    exons_to_upsert.append(exon)

                elif feature == 'CDS':
                    transcript_id = attributes['transcript_id']
                    protein_id = attributes['protein_id']
                    cds = (start, stop, protein_id)
                    if transcript_id in transcripts_to_cdss:
                        transcripts_to_cdss[transcript_id].append(cds)
                    else:
                        transcripts_to_cdss[transcript_id] = [cds,]

                if i % CHUNK_SIZE == 0:
                    bulk_upsert(session, Gene.__table__, genes_to_upsert)
                    bulk_upsert(session, GencodeTranscript, transcripts_to_upsert)
            bulk_upsert(session, Gene.__table__, genes_to_upsert)
            bulk_upsert(session, GencodeTranscript, transcripts_to_upsert)

        _upsert_exons_and_orfs(session, transcripts_to_exons, transcripts_to_cdss, minus_transcripts, exons_to_upsert)


def load_pacbio_gtf(db, gtf_file: str, overwrite=False) -> None:
    with db.get_session() as session:
        with session.begin():
            existing_genes = {row.name: row.accession for row in session.query(Gene.name, Gene.accession).all()}

        chromosomes = {}
        genes_to_insert = []
        transcripts_to_upsert = []
        exons_to_upsert = []
        orfs_to_upsert = []

        minus_transcripts = set()
        transcripts_to_exons = {}
        transcripts_to_cdss = {}

        with open(gtf_file) as gtf:
            lines = count_lines(gtf)
            t = tqdm(gtf, desc=f'Reading PacBio annotations', total=lines, unit='lines')
            for i, line in enumerate(t, start=1):
                if line.startswith("#"):
                    continue
                chr, _, feature, start, stop, _, strand, _, attributes, _ = read_gtf_line(line)
                if chr not in chromosomes:
                    with session.begin():
                        chromosome = session.merge(Chromosome(name=chr))
                    chromosomes[chr] = chromosome
                gene_name = attributes['gene_id']
                if gene_name not in existing_genes:
                    tqdm.write(f'Could not find Ensembl accession for gene \'{gene_name}\'')
                    existing_genes[gene_name] = gene_name
                    gene = {
                        'accession': gene_name,
                        'name': gene_name,
                        'strand': Strand.from_symbol(strand),
                        'chromosome_id': chr
                    }
                    genes_to_insert.append(gene)
                gene_id = existing_genes[gene_name]

                if feature == 'transcript':
                    transcript_id = attributes['transcript_id'].split('|')[1]
                    transcript_name = gene_name + '|' + transcript_id
                    transcript = {
                        'accession': transcript_id,
                        'name': transcript_name,
                        'type': 'pacbiotranscript',
                        'gene_id': gene_id,
                        'strand': Strand.from_symbol(strand),
                    }
                    transcripts_to_upsert.append(transcript)
                    if Strand.from_symbol(strand) is Strand.MINUS:
                        minus_transcripts.add(transcript_id)

                elif feature == 'exon':
                    transcript_id = attributes['transcript_id'].split('|')[1]
                    exon = {
                        'start': start,
                        'stop': stop,
                        'transcript_id': transcript_id
                    }
                    if transcript_id in transcripts_to_exons:
                        transcripts_to_exons[transcript_id].append(exon)
                    else:
                        transcripts_to_exons[transcript_id] = [exon,]
                    exons_to_upsert.append(exon)

                elif feature == 'CDS':
                    transcript_id = attributes['transcript_id'].split('|')[1]
                    cds = (start, stop)
                    if transcript_id in transcripts_to_cdss:
                        transcripts_to_cdss[transcript_id].append(cds)
                    else:
                        transcripts_to_cdss[transcript_id] = [cds,]

                if i % CHUNK_SIZE == 0:
                    bulk_upsert(session, Gene.__table__, genes_to_insert)
                    bulk_upsert(session, PacBioTranscript, transcripts_to_upsert)
            bulk_upsert(session, Gene.__table__, genes_to_insert)
            bulk_upsert(session, PacBioTranscript, transcripts_to_upsert)

        _upsert_exons_and_orfs(session, transcripts_to_exons, transcripts_to_cdss, minus_transcripts, exons_to_upsert)


def load_transcript_fasta(db, transcript_fasta: str, id_extractor: Callable[[str], 'FastaHeaderFields'], id_filter: Callable[[str], bool] = lambda x: False):
    with db.get_session() as session:
        with session.begin():
            existing_transcripts = {row.accession for row in session.query(Transcript.accession)}
            existing_orfs = {
                row.transcript_id: (row.position, row.transcript_start, row.transcript_stop)
                for row in session.execute(select(ORF.position, ORF.transcript_id, ORF.transcript_start, ORF.transcript_stop))
            }
        transcripts_to_update = []
        orfs_to_update = []
        with open(transcript_fasta) as f:
            records = count_lines(f, only=lambda line: line.startswith('>'))
        t = tqdm(SeqIO.parse(transcript_fasta, 'fasta'), desc='Reading transcripts fasta', total=records, unit='seqs')
        for record in t:
            if id_filter(record.id):
                continue
            ids = id_extractor(record.id)
            transcript_id = ids.transcript_id
            sequence = str(record.seq)
            if transcript_id in existing_transcripts:
                transcript = {'accession': transcript_id, 'sequence': sequence}
                transcripts_to_update.append(transcript)
            if transcript_id in existing_orfs:
                position, tx_start, tx_stop = existing_orfs[transcript_id]
                if sequence[tx_stop:tx_stop+3] in STOP_CODONS:
                    orfs_to_update.append({
                        'transcript_id': transcript_id,
                        'position': position,
                        'transcript_start': tx_start,
                        'transcript_stop': tx_stop + 3,
                        'has_stop_codon': True
                    })

            if len(transcripts_to_update) == CHUNK_SIZE:
                bulk_upsert(session, Transcript.__table__, transcripts_to_update)
                bulk_upsert(session, ORF.__table__, orfs_to_update, primary_keys=('transcript_id', 'position'))
        bulk_upsert(session, Transcript.__table__, transcripts_to_update)
        bulk_upsert(session, ORF.__table__, orfs_to_update, primary_keys=('transcript_id', 'position'))


def load_translation_fasta(db, translation_fasta: str, id_extractor: Callable[[str], 'FastaHeaderFields'], id_filter: Callable[[str], bool] = lambda x: False, overwrite: bool = False):
    with db.get_session() as session:
        with session.begin():
            existing_orfs = {}
            for transcript_id, position, start, stop, has_stop_codon in session.query(ORF.transcript_id, ORF.position, ORF.transcript_start, ORF.transcript_stop, ORF.has_stop_codon):
                existing_orfs.setdefault(transcript_id, []).append((position, start, stop, has_stop_codon))

        orfs_to_update = []
        proteins_to_upsert = []
        with open(translation_fasta) as f:
            records = count_lines(f, only=lambda line: line.startswith('>'))
        t = tqdm(SeqIO.parse(translation_fasta, 'fasta'), desc='Reading translations fasta', total=records, unit='seqs')
        for i, record in enumerate(t, start=1):
            if id_filter(record.id):
                continue
            ids = id_extractor(record.id)
            transcript_id = ids.transcript_id
            protein_id = ids.protein_id
            sequence = str(record.seq)
            seq_length = len(sequence)

            protein = {'accession': protein_id, 'sequence': sequence}
            proteins_to_upsert.append(protein)
            if transcript_id in existing_orfs:
                orfs = existing_orfs[transcript_id]
                for position, start, stop, has_stop_codon in orfs:
                    orf_nt_length = stop - start + 1
                    orf_aa_length = seq_length + int(has_stop_codon)
                    if orf_nt_length == orf_aa_length*3:
                        orfs_to_update.append({
                            'transcript_id': transcript_id,
                            'position': position,
                            'transcript_start': start,
                            'transcript_stop': stop,
                            'protein_id': protein_id
                        })
                        if has_stop_codon:
                            protein['sequence'] = protein['sequence'] + AminoAcid.STOP.value
                        break

            if i % CHUNK_SIZE == 0:
                bulk_upsert(session, Protein.__table__, proteins_to_upsert)
                bulk_upsert(session, ORF.__table__, orfs_to_update, primary_keys=('transcript_id', 'position'))
        bulk_upsert(session, Protein.__table__, proteins_to_upsert)
        bulk_upsert(session, ORF.__table__, orfs_to_update, primary_keys=('transcript_id', 'position'))


def load_sqanti_classifications(db, sqanti_file: str):
    with db.get_session() as session:
        with session.begin():
            existing_gencode_transcripts = {row.accession for row in session.query(GencodeTranscript.accession)}
            existing_pacbio_transcripts = {row.accession for row in session.query(PacBioTranscript.accession)}
        transcripts_to_update = []
        with open(sqanti_file) as f:
            f.readline()  # skip header
            lines = count_lines(f)
            reader = csv.DictReader(f, delimiter='\t')
            for row in tqdm(reader, desc='Reading SQANTI classifications', total=lines, unit='lines'):
                if row['isoform'] in existing_pacbio_transcripts:
                    tx = {
                        'accession': row['isoform'],
                        'sqanti': SQANTI_DICT.get(row['structural_category'], SQANTI.OTHER)
                    }
                    associated_tx = row['associated_transcript']
                    if tx['sqanti'] in {SQANTI.FSM, SQANTI.ISM} and associated_tx in existing_gencode_transcripts:
                        tx['gencode_id'] = associated_tx
                    else:
                        tx['gencode_id'] = None
                    transcripts_to_update.append(tx)
                if len(transcripts_to_update) == CHUNK_SIZE:
                    bulk_upsert(session, PacBioTranscript.__table__, transcripts_to_update)
            bulk_upsert(session, PacBioTranscript.__table__, transcripts_to_update)


def load_domains(db, domain_file: str, overwrite: bool = False):
    with open(domain_file) as f:
        lines = count_lines(f)
        reader = csv.reader(f, delimiter='\t')
        t = tqdm(reader, desc='Reading domain info', total=lines, unit='domains')
        domains_to_upsert = (
            {
                'type': FeatureType.DOMAIN,
                'accession': acc,
                'name': name,
                'description': desc
            } for acc, name, _, desc, *_ in t
        )
        domains_to_upsert = list(domains_to_upsert)
    with db.get_session() as session:
        bulk_upsert(session, Domain.__table__, domains_to_upsert)


def load_patterns(db, pattern_file: str):
    with open(pattern_file) as f:
        name_getter = re.compile(r'(ID   )(\S+)(; PATTERN.)')
        acc_getter = re.compile(r'(AC   )(\S+)(;)')
        desc_getter = re.compile(r'(DE   )(.+)')

        lines = count_lines(f, only=lambda s: name_getter.match(s))
        t = tqdm(None, desc='Reading pattern info', total=lines, unit='patterns')
        patterns_to_upsert = []
        for line in f:
            if (name := name_getter.match(line)):
                pattern = {'name': name.group(2), 'type': FeatureType.NONE}
                t.update()
            elif (acc := acc_getter.match(line)):
                pattern['accession'] = acc.group(2)
            elif (desc := desc_getter.match(line)):
                pattern['description'] = desc.group(2)
                patterns_to_upsert.append(pattern)
    with db.get_session() as session:
        bulk_upsert(session, Feature.__table__, patterns_to_upsert)


def load_feature_mappings(db, domain_mapping_file: str, appris_only: bool = True, overwrite: bool = False):
    with db.get_session() as session:
        feature_types = [
            {
                'type': FeatureType.IDR,
                'accession': 'mobidb-lite',
                'name': 'MobiDB',
                'description': 'intrinsically disordered region (MobiDB)'
            }
        ]
        bulk_upsert(session, Feature.__table__, feature_types)

        with session.begin():
            if overwrite:
                print(f'Clearing table \'{ProteinFeature.__tablename__}\'...')
                session.execute(delete(ProteinFeature.__table__))
            if appris_only:
                protein_length = func.length(Protein.sequence)
                subq = (
                    select(GencodeTranscript.gene_id, Protein.accession, GencodeTranscript.appris, protein_length).
                    select_from(Protein).
                    join(Protein.orf).
                    join(ORF.transcript).
                    order_by(GencodeTranscript.gene_id).
                    order_by(desc(GencodeTranscript.appris)).
                    order_by(desc(protein_length)).
                    order_by(GencodeTranscript.name).
                    subquery()
                )
                proteins_to_map = set(
                    session.execute(
                        select(subq.c.accession).select_from(subq).group_by(subq.c.gene_id)
                    ).scalars()
                )
            else:
                proteins_to_map = set(session.execute(select(Protein.accession)).scalars())
            existing_features = set(session.execute(select(Feature.accession)).scalars())
        with open(domain_mapping_file) as f:
            f.readline()  # skip header
            lines = count_lines(f)
            reader = csv.DictReader(f, delimiter='\t')
            t = tqdm(reader, desc='Reading reference domain mappings', total=lines, unit='mappings')
            domains_to_insert = (
                {
                    'feature_id': row['feature id'],
                    'protein_id': row['Protein stable ID version'],
                    'protein_start': row['start'],
                    'protein_stop': row['end'],
                    'reference': True
                }
                for row in t if row['feature id'] in existing_features and row['Protein stable ID version'] in proteins_to_map
            )
            domains_to_insert = list(domains_to_insert)
        with session.begin():
            session.execute(insert(ProteinFeature.__table__).on_conflict_do_nothing(), domains_to_insert)


def project_feature_mappings(db, gene_ids: Iterable[str] = None, overwrite: bool = False):
    with db.get_session() as session:
        protein_tx = contains_eager(Protein.orf).contains_eager(ORF.transcript)
        q = (
            select(Protein, Transcript.gene_id).
            join(Protein.orf).
            join(ORF.transcript).
            order_by(Transcript.gene_id).
            order_by(desc(Transcript.__table__.c.appris)).
            order_by(desc(ORF.length)).
            order_by(Transcript.name).
            options(
                protein_tx.joinedload(Transcript.orfs),
                protein_tx.joinedload(Transcript.exons).joinedload(Exon.transcript),
                protein_tx.joinedload(Transcript.gene),
                contains_eager(Protein.orf).joinedload(ORF.protein),
                joinedload(Protein.features).joinedload(ProteinFeature.protein),
                joinedload(Protein.features).joinedload(ProteinFeature.feature),
                raiseload('*')
            )
        )
        if not gene_ids:
            gene_ids = list(session.execute(
                select(Transcript.gene_id).distinct().
                select_from(ORF).
                join(ORF.transcript).
                order_by(Transcript.gene_id)
            ).scalars())
        nrows = session.execute(
            select(func.count(Transcript.accession)).
            select_from(ORF).
            join(ORF.transcript).
            where(Transcript.gene_id.in_(gene_ids))
        ).scalars().first()
        rows = chain.from_iterable(
            session.execute(
                q.where(Transcript.gene_id.in_(gene_chunk))
            ).unique()
            for gene_chunk in chunked(gene_ids, 200)
        )
        t = tqdm(None, desc='Projecting domain mappings', total=nrows, unit='proteins', mininterval=0.2)
        domains_to_insert = []
        for gene_chunk in chunked(gene_ids, 200):
            rows = session.execute(
                q.where(Transcript.gene_id.in_(gene_chunk))
            ).unique()
            rows_by_gene = groupby(rows, key=itemgetter(1))
            for gene, group in rows_by_gene:
                proteins = [row[0] for row in group]
                anchor = proteins[0]
                for other in proteins[1:]:
                    if not anchor.features:
                        continue
                    try:
                        aln = CodonAlignment.from_proteins(anchor, other)
                    except Exception as e:
                        tqdm.write(f'{anchor}|{other}\t{e}')
                        continue
                    for feat in anchor.features:
                        proj_feat = aln.project_feature(feat)
                        if proj_feat:
                            record = {
                                k: getattr(proj_feat, k)
                                for k in (
                                    'protein_start',
                                    'protein_stop',
                                    'reference',
                                    '_differences'
                                )
                            }
                            record['feature_id'] = proj_feat.feature.accession
                            record['protein_id'] = other.accession
                            record['anchor_id'] = proj_feat.anchor.id
                            domains_to_insert.append(record)
                t.update(len(proteins))
            if domains_to_insert:
                session.execute(insert(ProteinFeature.__table__).on_conflict_do_nothing(), domains_to_insert)
                session.commit()
                domains_to_insert[:] = []
