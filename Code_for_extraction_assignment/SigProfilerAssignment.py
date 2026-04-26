# ===============================================================================
# Mutational Signature Assignment Workflow Using SigProfilerAssignment(172-type)
# ===============================================================================
# We downloaded SigProfilerAssignment v0.2.3 from https://github.com/SigProfilerSuite/SigProfilerAssignment
# Input: Catalog of 172-type-catalogs
# Output: Assignment results and reconstruction plots saved to specified directory
# We ran SigProfilerAssignment for all genomes except: high TMB genomes (TMB>=5000) and ID1/ID2 dominated genomes (Ins:1:T:5+ and Del:1:T:5+ >= 90%)
# Python version 3.10

from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(
  samples="./example.data/assignment.example.data/example.172.type.catalog.txt",
  input_type='matrix',
  output='./output/',
  collapse_to_SBS96=False, 
  cosmic_version=None, 
  exome=False,
  genome_build="GRCh37", 
  signature_database="./example.data/assignment.example.data/example.172.type.signatures.txt",
  exclude_signature_subgroups=None, 
  export_probabilities=True,
  export_probabilities_per_mutation=True, 
  make_plots=True,
  sample_reconstruction_plots=True, 
  verbose=False
)
