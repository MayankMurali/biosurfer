#%%
from biosurfer.core.database import Database
from biosurfer.core.helpers import get_ids_from_gencode_fasta, skip_par_y
import os

#%%
def create_gencode(gencode_gtf, gencode_tx, gencode_tl, gencode_doms, pfam_dom_info, prosite_patterns, db_name):
    """ Creating new SQLite3 gencode database 
    Args:
        gencode_gtf: Gene annotation file (GTF)
        gencode_tx: Transcript reference sequence file (FASTA)
        gencode_tl: Translation reference sequence file (FASTA)
        gencode_doms: grch38 protein feature file (TSV)
        pfam_dom_info: Protein Family mapping file (TSV)
        prosite_patterns: PROSITE pattern data file 
        db_name: User input database name

    Returns:
      Nothing
    """
    db = Database(db_name)
    db.recreate_tables()
    db.load_gencode_gtf(os.path.abspath(gencode_gtf), overwrite=True)
    db.load_transcript_fasta(os.path.abspath(gencode_tx), get_ids_from_gencode_fasta, skip_par_y)
    db.load_translation_fasta(os.path.abspath(gencode_tl), get_ids_from_gencode_fasta, skip_par_y, overwrite=True)
    db.load_domains(os.path.abspath(pfam_dom_info))
    db.load_patterns(os.path.abspath(prosite_patterns))
    db.load_feature_mappings(os.path.abspath(gencode_doms), overwrite=False)
