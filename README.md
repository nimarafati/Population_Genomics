# Population_Genomics

In this repostiry you can find a collection of scripts used in studying population genetics differentiation. 

**Muscle_FastTree.sh**. 

By this script you can run multiple sequence alignment by **Muscle** and the generate tree by **FastTree**.  The generated tree can be visualised by tools such as **FigTree**.  
Muscle_FastTree.sh sequence.Fa  
Generated files are: sequence.muscle and sequence.tree. You will use sequence.tree for visualising the tree.  

**extract-GT-V2-1-exclude-miss-called-positions-give-GT-too.pl**. 

By this script you can extract genotypes in a vcf file in three formats:  
"fasta" format 
>ind1-1. 
hap-1. 
>ind1-2. 
hap-2. 
"tab" format makes a matrix for haplotypes:  
ind1-1	h	a	p	-	1.  
ind1-2	h	a	p	-	2. 
"tab_GT" 0:Hom.Ref 1:Het 2:Hom.Alt. 
ind-1 0	1	2  
ind-1 0 2	1. 
Here is how you run it: 
extract-GT-V1.pl -format fasta/tab -vcf SNP.vcf -phased T/F -missingGT T/F. 




