from sqlalchemy import Column, Integer, String, Float, ForeignKey, Enum
from sqlalchemy.orm import relationship
from biosurfer.core.models.base import Base, TablenameMixin
from biosurfer.core.models.nonpersistent import Position
from biosurfer.core.constants import Strand

class VariantType(str):
    SNP = "SNP"
    INDEL = "INDEL"

class GenomicVariant(Base, TablenameMixin):
    """
    Represents a specific genomic change (e.g., chr3:12345 A>G).
    Stored in 1-based genomic coordinates.
    """
    id = Column(Integer, primary_key=True)
    chromosome = Column(String, index=True)
    position = Column(Integer, index=True)
    reference_allele = Column(String)
    alternative_allele = Column(String)
    # Optional: rsID if available
    rsid = Column(String, index=True, nullable=True)
    
    # Relationship to GWAS data
    gwas_stats = relationship("GWASStatistic", back_populates="variant")
    
    def __repr__(self):
        return f"{self.chromosome}:{self.position}_{self.reference_allele}>{self.alternative_allele}"

class GWASStatistic(Base, TablenameMixin):
    """
    Stores summary statistics for a variant.
    """
    id = Column(Integer, primary_key=True)
    variant_id = Column(Integer, ForeignKey("genomicvariant.id"))
    variant = relationship("GenomicVariant", back_populates="gwas_stats")
    
    p_value = Column(Float)
    beta = Column(Float) # Effect size
    se = Column(Float)   # Standard Error
    trait = Column(String) # e.g., "T2D"
    
class SampleGenotype(Base, TablenameMixin):
    """
    Stores individual patient genotypes (from VCF).
    """
    id = Column(Integer, primary_key=True)
    variant_id = Column(Integer, ForeignKey("genomicvariant.id"))
    sample_id = Column(String, index=True)
    genotype = Column(String) # e.g., "0/1", "1/1"
