# ~/annovar$ perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp42c humandb/

# Run ithis in ~/annovar/
perl table_annovar.pl ~/github/Code_Liu_2025/ID15_ID16/my_results.avinput humandb/ \
  -buildver hg19 \
  -out ~/github/Code_Liu_2025/ID15_ID16/impact_results \
  -protocol refGene,dbnsfp42c \
  -operation g,f \
  -nastring . \
  -otherinfo