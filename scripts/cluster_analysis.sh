#~~~~~~~~~~ Step 1: Clustering using admixture
cluster_admixture() {

	bedfile="bed file"
	clusterfile=cluster.txt
	
	for ((K=1; K<=10; K++))
	do
		admixture --cv $bedfile $K -j8 | tee log${K}.out >> $clusterfile
	done
	
	echo "cluster_admixture" >> $timefile
}

#Edit the cluster files so that only CV for each K value should be there CV.txt. This information can be taken from the cluster.txt generated from previous step.

#~~~~~~~~~~ Step 2: Clustering plot
cluster_plot() {
	
	trait="Trait symbol"
	pheno_file="Phenotypic file"
	
	clusterfile=CV.txt
	
	Rscript --vanilla cluster_Kplot.R $trait $clusterfile $pheno_file

}
