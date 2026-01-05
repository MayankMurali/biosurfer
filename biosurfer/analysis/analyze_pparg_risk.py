from biosurfer.core.database import Database
from biosurfer.analysis.genetics_analyzer import analyze_nterm_risk

# Connect to your populated DB
db = Database("biosurfer_db") # or whatever your DB name is

with db.get_session() as session:
    # Run the N-term specific analysis
    analyze_nterm_risk(session, gene_name="PPARG")
