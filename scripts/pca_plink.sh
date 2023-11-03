#!/bin/bash

b_file="" 	#plink bfiles without extension
trait="Trait symbol"	#Ex: SPY

#~~~~~~~~~~ PCA through PLINK
pca_plink() {
	
	plink2 --bfile $b_file --pca 5 --nonfounders --out $b_file"-pc"
	
	Rscript pca_plink.R $trait $famfile"-pc.eigenvec" $b_file"-pc.eigenval"

}
