# Run SigProfilerExtractor for de novo extraction of mutational signatures from (1) all genomes; (2) genomes of each cancer type; (3) high-TMB genomes.
# We downloaded SigProfilerExtractor v1.2.1 from https://github.com/SigProfilerSuite/SigProfilerExtractor
# The examples shows the 89-type extraction. The same parameters were applied to 83-type and 476-type extraction
# Python version 3.10

import os

# Include 21 cancer types:"Myeloid","Other","Pancreas","Prostate","Skin","Stomach","Thymus","Uterus","Colon",
#  "Esophagus","Head","Kidney","Liver","Lung","Lymphoid","Ovary","Biliary","Bladder",
#  "Bone.SoftTissue","Breast","CNS", "All","TMB" (refers to the genomes with total indel>=5000). 

# Only ran for two example catalogs.
cancer_types = [
  "Bladder", "All"
]


# Please set to your file path
base_dir = "./example.data/sigprofiler.input.89type.catalog/"


for ct in cancer_types:
  from SigProfilerExtractor import sigpro as sig

import multiprocessing
## probably not needed, but to be safe set start method to 'fork'
multiprocessing.set_start_method('fork', force=True)
## Define input and output paths based on the cancer type
output_prefix = os.path.join(base_dir, f"example.89type.{ct}")
input_matrix  = os.path.join(base_dir, f"example.89type.{ct}.txt")

print(f"Running SigProfilerExtractor for {ct}")
print(f"  Input matrix: {input_matrix}")
print(f"  Output base:  {output_prefix}")

## Set different signature number ranges for "All" vs individual cancer types. The rest parameters are the same as default. 
## Please enable GPU support. It is much faster than CPU. 
if ct == "All":
  min_sig, max_sig = 20, 40
else:
  min_sig, max_sig = 2, 20
try:
  sig.sigProfilerExtractor(
    input_type="matrix",
    output=output_prefix,
    input_data=input_matrix,
    reference_genome="GRCh37",
    opportunity_genome="GRCh37",
    context_type="default",
    exome=False,
    minimum_signatures=min_sig,
    maximum_signatures=max_sig,
    nmf_replicates=100,
    resample=True,
    batch_size=1,
    cpu=-1,
    gpu=True,
    nmf_init="random",
    precision="single",
    matrix_normalization="gmm",
    seeds="random",
    min_nmf_iterations=10000,
    max_nmf_iterations=1000000,
    nmf_test_conv=10000,
    nmf_tolerance=1e-15,
    get_all_signature_matrices=False
  )
except Exception as e:
  print(f"[WARN] {ct} failed with error: {e}")
