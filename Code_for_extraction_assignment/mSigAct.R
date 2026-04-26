# =================================================================
# Mutational Signature Assignment Workflow Using mSigAct(172-type)
# =================================================================
# Run mSigAct to decompose mutational spectra into known signatures
# We downloaded SigProfilerAssignment v0.2.3 from https://github.com/SigProfilerSuite/SigProfilerAssignment

# We applied mSigAct to genomes with high tumor mutational burden (TMB ≥ 5000) 
# and to those dominated by ID1/ID2 signatures (Ins:1:T:5+ and Del:1:T:5+ ≥ 90%), 
# to mitigate excessively sparse results from SPA in low‑complexity catalogs.


# Load required libraries
library(mSigAct2)
library(data.table)
library(parallel)

# Load catalogs and signatures
catalog <- fread("./example.data/assignment/example.data/example.172.type.catalog.txt")
sigs <- fread("./example.data/example.172.type.signatures.txt")[, -1]
sigs <- apply(as.matrix(sigs), 2, function(x) x / sum(x))

# Quick QP assignment
QP.assign <- sapply(1:ncol(catalog), function(i) {
  test <- mSigTools::find_best_reconstruction_QP(as.matrix(catalog[, i, drop=FALSE]), sig.universe=sigs)
  test$optimized.exposure
})
colnames(QP.assign) <- colnames(catalog)

# Parallel mSigAct assignment
process_sample <- function(i) {
  samplename <- colnames(QP.assign)[i]
  assign.universe <- rownames(QP.assign)[QP.assign[, i] > 0]
  if (length(assign.universe) == 1) assign.universe <- unique(c(assign.universe, rownames(QP.assign)[1:5]))
  
  my_opts <- DefaultManyOpts(); my_opts$nbinom.size <- 20
  retval <- PresenceAttributeSigActivity(
    spectra = as.matrix(catalog[, i, drop=FALSE]),
    sigs = sigs[, assign.universe, drop=FALSE],
    output.dir = "temp", mc.cores.per.sample = 1, num.parallel.samples = 1, seed = 1234, m.opts = my_opts
  )
  
  list(ID172 = retval$proposed.assignment[,1],
       
       sample = samplename)
}
# -------------------------------------------------------------------------------------
# Select samples (high mutation burden + specific cases). here only shows an example.
# -------------------------------------------------------------------------------------
results <- mclapply(1:ncol(catalog), process_sample, mc.cores=40)
saveRDS(results, "new_result_list.rds")
        
     

        