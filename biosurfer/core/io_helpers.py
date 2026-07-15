# GTF/FASTA parsing and bulk database-loading helpers
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Optional

from sqlalchemy.dialects.sqlite.dml import insert

if TYPE_CHECKING:
    from io import TextIOBase


@dataclass
class FastaHeaderFields:
    transcript_id: str = None
    protein_id: str = None


def count_lines(file_handle: 'TextIOBase', only: Optional[Callable[..., bool]] = None):
    lines = sum(1 for _ in filter(only, file_handle))
    file_handle.seek(0)
    return lines


def bulk_upsert(session, table, records, primary_keys=('accession',)):
    if records:
        fields = [field for field in records[0] if field not in primary_keys]

        stmt = insert(table)
        upsert_stmt = stmt.on_conflict_do_update(
            index_elements=primary_keys,
            set_={field: stmt.excluded[field] for field in fields}
        )

        # FIX: Check if transaction is active.
        # If active, participate in it. If not, start a new one.
        if session.in_transaction():
            session.execute(upsert_stmt, records)
        else:
            with session.begin():
                session.execute(upsert_stmt, records)

        records[:] = []


def read_gtf_line(line: str) -> list:
    """Read and parse a single gtf line

    Args:
        line (str): unbroken line of a gtf file

    Returns:
        list: gtf attributes
        chromosome : str
        source : str
        feature : str
        start : int
        stop : int
        score : str
        strand : str
        phase : str
        attributes: dict
        tags: list

    """
    chromosome, source, feature, start, stop, score, strand, phase, attributes = line.split('\t')
    start = int(start)
    stop = int(stop)
    attributes = attributes.split(';')[:-1]
    attributes = [att.strip(' ').split(' ') for att in attributes]
    tags = [att[1].strip('"') for att in attributes if att[0] == 'tag']
    attributes = {att[0]: att[1].strip('"') for att in attributes if att[0] != 'tag'}
    return chromosome, source, feature, start, stop, score, strand, phase, attributes, tags


def get_ids_from_gencode_fasta(header: str):
    fields = [field for field in header.split('|') if field and not field.startswith(('UTR', 'CDS'))]
    transcript_id = next((field for field in fields if field.startswith('ENST')))
    protein_id = next((field for field in fields if field.startswith('ENSP')), None)
    return FastaHeaderFields(transcript_id, protein_id)

def get_ids_from_pacbio_fasta(header: str):
    return FastaHeaderFields(header, None)

def get_ids_from_lrp_fasta(header: str):
    fields = header.split('|')
    return FastaHeaderFields(fields[1], fields[1] + ':PROT1')

def skip_par_y(header: str):  # these have duplicate ENSEMBL accessions and that makes SQLAlchemy very sad
    return 'PAR_Y' in header

def skip_gencode(header: str):
    return header.startswith('gc')
