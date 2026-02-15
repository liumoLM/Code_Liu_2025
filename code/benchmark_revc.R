library(microbenchmark)
library(Biostrings)

# 1. Generate 100 sequences of length 30
set.seed(123)
bases <- c("A", "C", "G", "T")
vec_dna <- replicate(
  100,
  paste(sample(bases, 30, replace = TRUE), collapse = "")
)

# 2. Prepare Biostrings input (vectorized)
dna_set <- DNAStringSet(vec_dna)

# 3. Define methods
# Method A: Standard split/paste (wrapped in sapply)
bench_split <- function(v) {
  sapply(v, function(x) {
    chartr("ATGC", "TACG", paste(rev(strsplit(x, NULL)[[1]]), collapse = ""))
  })
}

# Method B: Integer trick (wrapped in sapply)
bench_integer <- function(v) {
  sapply(v, function(x) {
    chartr("ATGC", "TACG", intToUtf8(rev(utf8ToInt(x))))
  })
}

# Method C: Biostrings (Directly vectorized)
bench_bioc <- function(s) {
  reverseComplement(s)
}

# 4. Run Benchmark
results <- microbenchmark(
  split_paste = bench_split(vec_dna),
  integer_trick = bench_integer(vec_dna),
  biostrings = bench_bioc(dna_set),
  times = 100
)

print(results)
