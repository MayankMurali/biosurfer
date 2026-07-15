from pathlib import Path
from sqlite3 import Connection as SQLite3Connection
from typing import Callable, Dict, Iterable

from biosurfer.core import data_loaders
from biosurfer.core.io_helpers import FastaHeaderFields
from biosurfer.core.models.base import Base
from biosurfer.core.models.biomolecules import Gene
from sqlalchemy import create_engine, event, select, and_
from sqlalchemy.engine import Engine
from sqlalchemy.orm import scoped_session, sessionmaker


# Adding genetics functionality
import pysam
import pandas as pd
from biosurfer.core.models.genetics import GenomicVariant, GWASStatistic, SampleGenotype


class Database:
    _databases_dir = Path(__file__).parent.parent.parent/'databases'
    registry: Dict[str, 'Database'] = {}

    @staticmethod
    def _get_db_url_from_name(name: str):
        if name:
            db_file = f'{name}.sqlite3'
            return f'sqlite:///{Database._databases_dir/db_file}'
        else:
            return 'sqlite://'

    def __new__(cls, name: str = None, *, url: str = None, **kwargs):
        if url is None:
            url = Database._get_db_url_from_name(name)
        if url in Database.registry:
            return Database.registry[url]
        else:
            obj = super().__new__(cls)
            Database.registry[url] = obj
            return obj

    def __init__(self, name: str = None, *, url: str = None, sessionfactory=None):
        if url is None:
            url = Database._get_db_url_from_name(name)
        self.url = url
        self._engine = create_engine(self.url)
        Base.metadata.create_all(self._engine)
        if sessionfactory is None:
            self._sessionmaker = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=self.engine, future=True))
        else:
            self._sessionmaker = sessionfactory

    def __repr__(self):
        return f'Database(url=\'{self.url}\')'

    @property
    def engine(self):
        return self._engine

    def get_session(self, **kwargs):
        return self._sessionmaker(**kwargs)

    def recreate_tables(self):
        print(f'Recreating tables in {self.url} ...')
        Base.metadata.drop_all(bind=self._engine)
        Base.metadata.create_all(bind=self._engine)

    # --- NEW METHOD: Get variants for plotting ---
    def get_gene_variants(self, session, gene_name: str) -> list['GenomicVariant']:
        """Retrieves genomic variants located within the start/stop boundaries of a gene."""
        gene = Gene.from_name(session, gene_name)
        if not gene:
            return []

        # Query variants on the same chromosome within the gene's range
        query = select(GenomicVariant).where(
            and_(
                GenomicVariant.chromosome == gene.chromosome_id,
                GenomicVariant.position >= gene.start,
                GenomicVariant.position <= gene.stop
            )
        )
        return session.execute(query).scalars().all()

    def load_genetics_data(self, vcf_path: str, gwas_path: str, gene_name: str, trait: str = 'GWAS'):
        """
        Loads VCF and GWAS data specifically for the coordinates of a target gene.
        """
        with self.get_session() as session:
            # 1. Get Gene Coordinates
            gene = Gene.from_name(session, gene_name)
            if not gene:
                print(f"Gene {gene_name} not found in DB. Load GTF first.")
                return

            chrom = gene.chromosome_id
            start = gene.start
            stop = gene.stop

            print(f"Loading genetics for {gene_name} ({chrom}:{start}-{stop})...")

            # 2. Parse VCF (Genotyping) using pysam (requires indexed VCF)
            vcf = pysam.VariantFile(vcf_path)

            # Buffer to bulk insert
            variants_map = {} # Key: (chr, pos, ref, alt) -> Variant Object

            try:
                for record in vcf.fetch(chrom, start, stop):
                    # Basic filtering for SNPs/Indels
                    ref = record.ref
                    for alt in record.alts:
                        key = (chrom, record.pos, ref, alt)

                        if key not in variants_map:
                            var = GenomicVariant(
                                chromosome=chrom,
                                position=record.pos,
                                reference_allele=ref,
                                alternative_allele=alt,
                                rsid=record.id
                            )
                            session.add(var)
                            session.flush() # Get ID
                            variants_map[key] = var

                        # Store Genotypes for your cohort
                        for sample in record.samples:
                            gt = record.samples[sample]['GT'] # Returns (0, 1) etc
                            if gt != (0, 0): # Skip homozygous reference to save space
                                session.add(SampleGenotype(
                                    variant_id=variants_map[key].id,
                                    sample_id=sample,
                                    genotype=f"{gt[0]}/{gt[1]}"
                                ))
            except ValueError as e:
                print(f"Error fetching VCF region: {e}. Ensure VCF is tabix indexed.")

            # 3. Load GWAS Summary Stats
            # Assuming a TAB-separated file with columns: CHR, POS, REF, ALT, P, BETA
            chunk_size = 100000
            for chunk in pd.read_csv(gwas_path, sep='\t', chunksize=chunk_size):
                # Filter for gene region
                region_mask = (chunk['CHR'].astype(str) == chrom) & \
                            (chunk['POS'] >= start) & \
                            (chunk['POS'] <= stop)

                subset = chunk[region_mask]

                for _, row in subset.iterrows():
                    # Check if variant exists from VCF load, or create new
                    key = (str(row['CHR']), row['POS'], row['REF'], row['ALT'])

                    if key not in variants_map:
                        var = GenomicVariant(
                            chromosome=str(row['CHR']),
                            position=row['POS'],
                            reference_allele=row['REF'],
                            alternative_allele=row['ALT']
                        )
                        session.add(var)
                        session.flush()
                        variants_map[key] = var

                    session.add(GWASStatistic(
                        variant_id=variants_map[key].id,
                        p_value=row['P'],
                        beta=row['BETA'],
                        trait=trait
                    ))

            session.commit()

    def load_gencode_gtf(self, gtf_file: str, overwrite=False) -> None:
        data_loaders.load_gencode_gtf(self, gtf_file, overwrite=overwrite)

    def load_pacbio_gtf(self, gtf_file: str, overwrite=False) -> None:
        data_loaders.load_pacbio_gtf(self, gtf_file, overwrite=overwrite)

    def load_transcript_fasta(self, transcript_fasta: str, id_extractor: Callable[[str], 'FastaHeaderFields'], id_filter: Callable[[str], bool] = lambda x: False):
        data_loaders.load_transcript_fasta(self, transcript_fasta, id_extractor, id_filter=id_filter)

    def load_translation_fasta(self, translation_fasta: str, id_extractor: Callable[[str], 'FastaHeaderFields'], id_filter: Callable[[str], bool] = lambda x: False, overwrite: bool = False):
        data_loaders.load_translation_fasta(self, translation_fasta, id_extractor, id_filter=id_filter, overwrite=overwrite)

    def load_sqanti_classifications(self, sqanti_file: str):
        data_loaders.load_sqanti_classifications(self, sqanti_file)

    def load_domains(self, domain_file: str, overwrite: bool = False):
        data_loaders.load_domains(self, domain_file, overwrite=overwrite)

    def load_patterns(self, pattern_file: str):
        data_loaders.load_patterns(self, pattern_file)

    def load_feature_mappings(self, domain_mapping_file: str, appris_only: bool = True, overwrite: bool = False):
        data_loaders.load_feature_mappings(self, domain_mapping_file, appris_only=appris_only, overwrite=overwrite)

    def project_feature_mappings(self, gene_ids: Iterable[str] = None, overwrite: bool = False):
        data_loaders.project_feature_mappings(self, gene_ids=gene_ids, overwrite=overwrite)


if not Database._databases_dir.exists():
    Database._databases_dir.mkdir()


@event.listens_for(Engine, "connect")
def _set_sqlite_pragma(dbapi_connection, connection_record):
    if isinstance(dbapi_connection, SQLite3Connection):
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA foreign_keys=ON;")
        cursor.close()
