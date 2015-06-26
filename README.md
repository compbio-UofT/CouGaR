CouGaR
==========

This is a tool used to search for complex genomic rearrangements from matched normal/tumour next generation sequencing data.

Requirements
-----------
* samtools (tested with 0.1.18 (r982:295))
* gurobi solver
* java
* TCGA credentials (if downloading TCGA data)

Running
-----------
The CouGaR pipeline consists of X main stages
* BAM pre-processing
* HMM prunning
* Flow problem formulation and solving
* IP problem formulation and solving

A quick example is as follows

```cd pipeline```

```bash 01_prepare_sample.sh TEST_TCGA_5055 ~misko/test_space/bams/tumor.bam ~misko/test_space/bams/normal.bam hg19 SKCM```

```bash 02_cluster_mapability.sh `pwd`/TEST_TCGA_5055```

```bash 03_graph_ip_and_out.sh TCGA_5055 `pwd`/TEST_TCGA_5055```

BAM pre-processing
-----------
Because BAM files can be extremely large and are not necessary after discordant clusters and coverage are computed they can be easily preprocessed in this first stage. This is how we were able to run CouGaR on so many TCGA samples with very limited storage space (~6TB). For example you can preprocess the BAM files and then re-run downstream analysis without needing them again (unless you change the way clusters or coverage are computed).

Two scripts have been provided to pre-process BAM files.

```01_grab_tcga_sample_and_prepare.sh``` and ```01_prepare_sample.sh```

The first of these grabs the tumor and normal BAM files for a specified TCGA sample (requires a valid TCGA access key) and pre-processes it. The second of these scripts performs the pre-processing operation on local BAM files. In the local case you will need to specify which reference genome [hg18/hg19] is used and also assign a group label to this sample.




