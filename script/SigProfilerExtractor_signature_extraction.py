## Python code to run SigProfilerExtractor
## We used Python3.10 and SigProfilerExtractor v1.2.1
## We used the recommended parameters from https://github.com/AlexandrovLab/SigProfilerExtractor. The parameters here are only for quick testing.

from SigProfilerExtractor import sigpro as sig

sig.sigProfilerExtractor("matrix", "/transfer/LiuMo_data/Working/Code_Liu_2025/test_output/test_ID83_SigProExt_run", "/transfer/LiuMo_data/Working/Code_Liu_2025/test_data/test_ID83_catalog.tsv", reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "default", exome = False, minimum_signatures=2, maximum_signatures=5, nmf_replicates=10, resample = True, batch_size=1, cpu=-1, gpu=False,    nmf_init="random", precision= "single", matrix_normalization= "gmm", seeds= "random", min_nmf_iterations= 100, max_nmf_iterations=1000, nmf_test_conv= 1000, nmf_tolerance= 1e-15, get_all_signature_matrices= False)

sig.sigProfilerExtractor("matrix", "/transfer/LiuMo_data/Working/Code_Liu_2025/test_output/test_ID89_SigProExt_run", "/transfer/LiuMo_data/Working/Code_Liu_2025/test_data/test_ID89_catalog.tsv", reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "default", exome = False, minimum_signatures=2, maximum_signatures=5, nmf_replicates=10, resample = True, batch_size=1, cpu=-1, gpu=False,    nmf_init="random", precision= "single", matrix_normalization= "gmm", seeds= "random", min_nmf_iterations= 100, max_nmf_iterations=1000, nmf_test_conv= 1000, nmf_tolerance= 1e-15, get_all_signature_matrices= False)


sig.sigProfilerExtractor("matrix", "/transfer/LiuMo_data/Working/Code_Liu_2025/test_output/test_ID476_SigProExt_run", "/transfer/LiuMo_data/Working/Code_Liu_2025/test_data/test_ID89_catalog.tsv", reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "default", exome = False, minimum_signatures=2, maximum_signatures=5, nmf_replicates=10, resample = True, batch_size=1, cpu=-1, gpu=False,    nmf_init="random", precision= "single", matrix_normalization= "gmm", seeds= "random", min_nmf_iterations= 100, max_nmf_iterations=1000, nmf_test_conv= 1000, nmf_tolerance= 1e-15, get_all_signature_matrices= False)
