Outbreak With VarGenes
Bridenbecker/MalariaSqlReport - commit a8759ad
6/8/2020
DanB

This executable adds the ability to specify the exact set of Var Genes to use instead of having the model
determine them.  The idea is to allow th user (i.e. Jon) to understand how the model works when the var genes
are either very different or very similar.

=== New Config.json Parameter ===

"Malaria_Strain_Model": "FALCIPARUM_FIXED_STRAIN"

=== New Campaign Intervention ===

"class": "OutbreakIndividualMalariaVarGenes",
"MSP_Type": 1,
"Non_Spec_Type": 0,
"IRBC_Type": [ 
    1,74,147,220,293,366,439,512,585,658,731,804,877,950,23,96,169,242,315,388,461,534,
    607,680,753,826,899,972,45,118,191,264,337,410,483,556,629,702,775,848,921,994,67,
    140,213,286,359,432,505,578
],
"Minor_Epitope_Type": [ 
    2,1,1,4,3,0,2,1,0,2,3,3,1,3,4,1,4,
    1,1,4,3,1,4,4,1,4,1,1,1,0,1,2,3,0,
    2,0,0,2,2,2,3,2,3,1,4,4,1,0,1,2
]
