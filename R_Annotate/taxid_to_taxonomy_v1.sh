efetch -db taxonomy -id MY_LIST -format xml | \
xtract -pattern Taxon -first TaxId -element Taxon -block '*/Taxon'  \
-unless Rank -equals 'no rank' -tab ',' -sep '_' -element Rank,ScientificName