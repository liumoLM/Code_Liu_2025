write a new test script called checkSP.R.  

1. It lists directory ~/MEGA/important_mut_sig_data/pcawg_indel_vcfs and selects the 200 largest files that match "annotated.indel.vcf.gz". 
2. In a loop, for each of the selected file names do:
2.1 Read the files using ICAMS::ReadVCF with filter.status = NULL into variable `vcf`
2.2 Using the logic in @msi_study/msi_study.qmd, code block tract-lengths. In the body of the loop "for (i in seq_along(patterns...". 
keep the logic of lines 169 to 172 and save the counts to 1 of 4 different tibbles that contain the count tables for different values of R / repeat_count.
3. When the loop over ecah selected file name is finished, for each of the 4 tibbes, aggregate the repat count in column n based on the repeat_count grouping column.
4. using the login in msi_study.qmd, same code block, plot the histogram for each of the 4 patterns.
  