For each tumor (column) in Mansusript_data/jinalized_cap9/liu_et_al_89_assignment.tsv do:

calculate the partial credit of each signature to each mutation class in the tumor's 476 spectrum

then sum the partial credits for deletions of T from strectes of 6 T's, 7 T's, 8 T's and 9 T's and provide one vector for the tumor that maps each signature to its contribution to those summed partial credits.

To calculate "partial credit" for one signature, `sig` with contribution (assignment) of `sig_count_mutations` to one channel `channel` multiply the fraction of `sig` that `channel` is responsible for * `sig_count_mutations`

Example channels for deletions of T from stretches of 6 T's, 7 T's etc are 
"A[Del(T):R6]A", the regex is approximately [ACG]\[Del(T):R([6,7,8]|\(9,\))\][ACG].  You can look in ICAMS::catalog.row.order$ID476 for the full list.

So another way of describing this is that we want the partial credit for all signatures of detetions of 1 T from a stretch of 6,7,8, or 9 Ts.

I also want this assoicated with cancer type, so probably a place to start is transpose the assignment file and join with Manuscript_data/sample_info.tsv

