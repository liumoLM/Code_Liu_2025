
###############################################################################
# Code: mSigHdp Signature Extraction and Analysis
# Description: Extraction of Indel signatures using the mSigHdp algorithm.
###############################################################################

# --- 1. Library Installation & Loading ---
# Check and install hdpx (Dependency for mSigHdp)
if (!requireNamespace("hdpx", quietly = TRUE)) {
  devtools::install_github("steverozen/hdpx")
}

# Check and install mSigHdp
if (!requireNamespace("mSigHdp", quietly = TRUE)) {
  devtools::install_github("steverozen/mSigHdp")
}

library(hdpx)  # version 1.0.6
library(mSigHdp) # version 2.1.2



# --- 2. Data Loading & Preparation ---
# Load mutational catalogs
# Ensure test_data directory exists relative to your working directory
# We used the same parameters for 83-type, 89-type and 476-type. Here only shows the example for running 89-type
Indel89.catalog <- read.table("./mSigHdp.input.89type.catalog/mSigHdp.test.89type.catalog.txt",sep = "\t",stringsAsFactors =T)


# --- 3. Global Parameter Configuration ---
# Setting recommended parameters for high-fidelity extraction
hdp_args <- list(
  seedNumber = 1234,
  K.guess = 10,
  multi.types = FALSE,
  burnin = 1000,
  burnin.multiplier = 20,   
  post.n = 200,             
  post.space = 100,         
  num.child.process = 20,    
  CPU.cores = 20,           
  high.confidence.prop = 0.5,
  gamma.alpha = 1,
  gamma.beta = 50,
  checkpoint = FALSE,
  verbose = TRUE,
  downsample_threshold = 3000 # Minimizes memory footprint for high-count samples, only applied when running for all 6975 samples and high TMB samples
)

# --- 4. Extraction: All Samples ---
message("Starting Indel89 extraction for the example cohort...")
do.call(mSigHdp::RunHdpxParallel, c(list(
  input.catalog = Indel89.catalog, 
  out.dir = "./example.mSigHdp.extraction.output/"
), hdp_args))

