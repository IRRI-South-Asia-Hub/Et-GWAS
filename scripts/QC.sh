#!/bin/bash
#~~~~~~~~~~ Geno Quality Control
b_file="" 	#plink bfiles without extension
new_bfile=""	#new name for the filtered plink bfile
geno_qc() {
	
	plink --bfile $b_file --maf 0.05 --noweb --allow-no-sex --make-bed --out temp1
	plink --bfile temp1 --geno 0.1 --noweb --allow-no-sex --make-bed --out temp2
	plink --bfile temp2 --mind 0.1 --noweb --allow-no-sex --make-bed --out temp3
	plink --bfile temp3 --noweb --allow-no-sex --indep-pairwise 50 10 0.5
	plink --bfile temp3 --noweb --allow-no-sex --extract plink.prune.in --make-bed --out $new_bfile
	
	awk '{print $1,$6}' $new_bfile".fam" > $new_bfile".txt"

}
