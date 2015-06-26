CouGaR
==========

This is a tool used to search for complex genomic rearrangements from matched normal/tumour next generation sequencing data.

Requirements
-----------
* picard library
* samtools
* gurobi solver

Running
-----------

cd pipeline

bash 01_prepare_sample.sh TEST_TCGA_5055 /dupa-filer/misko/cougar/test_space/bams/tumor.bam /dupa-filer/misko/cougar/test_space/bams/normal.bam hg19 SKCM

bash 02_cluster_mapability.sh `pwd`/TEST_TCGA_5055

bash 03_graph_ip_and_out.sh TCGA_5055 `pwd`/TEST_TCGA_5055


