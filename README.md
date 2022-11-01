# SPOT-RNA-2D
Predicting RNA distance-based contact maps by integrated deep learning on physics-inferred secondary structure and evolutionary-derived mutational coupling

## System Requirments

**Hardware Requirments:**
It is recommended that your system should have 32 GB RAM, 500 GB disk space to support the in-memory operations for RNA sequence length less than 500. Multiple CPU threads are also recommended as the MSA generating process is computationally expensive.

**Software Requirments:**
* [Python3.6](https://docs.python-guide.org/starting/install3/linux/)
* [Perl-5.4 or later](https://www.perl.org/get.html)
* [Anaconda or Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
* [CUDA 10.0](https://developer.nvidia.com/cuda-10.0-download-archive) (Optional if running on GPU)
* [cuDNN (>= 7.4.1)](https://developer.nvidia.com/cudnn) (Optional if running on GPU)

SPOT-RNA-2D has been tested on Ubuntu 14.04, 16.04, and 18.04 operating systems.

## Installation of SPOT-RNA-2D and its dependencies

Clone SPOT-RNA-2D github repo:

1. `git clone https://github.com/jaswindersingh2/SPOT-RNA-2D.git && cd SPOT-RNA-2D`

To create **Conda** virtual environment:

2. `conda create -n venv_spotrna_2d python=3.6 && conda activate venv_spotrna_2d`

To install dependencies:

3. `while read p; do conda install --yes $p; done < requirements.txt`

To install **RNAfold** predictor for base pair probability features:

4. `conda install -c bioconda viennarna`

To install **BLAST-N** and **INFERNAL** tools for mulitple-sequence-alignment search:

5. `conda install -c bioconda blast`
6. `conda install -c bioconda infernal`
7. `conda install -c bioconda easel`

To install **PLMC** for DCA features:

8. `git clone https://github.com/debbiemarkslab/plmc && cd plmc && make all-openmp && cd -`


## Usage


### Run Single-sequence based predictor (SPOT-RNA-2D-Single):

SPOT-RNA-2D-Single only requires RNA sequence in directory mentioned with flag `--input_feats` and list of name of RNA sequences with flag `--list_rna_ids`. Set `--single_seq` flag to 1 for single-sequence-based prediction.

9. `./run.py --list_rna_ids datasets/TS1_ids --input_feats input_features/ --single_seq 1 --outputs outputs/`

To run SPOT-RNA-2D-Single for 1 RNA sequence instead of multiple RNAs:

10. `./run.py --rna_id 6ol3_C --input_feats input_features/ --single_seq 1 --outputs outputs/`


### Run co-evolutionary information based predictor (SPOT-RNA-2D):

SPOT-RNA-2D require RNAcmap based MSA features (PSSM, DCA) as an input. RNAcmap is computationally expensive and require NCBI's database of size around 400 GB. To obtain RNAcmap based features, following command can be used. For first run, it can take several hours to few days as it will download, unzip, and format NCBI's database to use with BLAST-N and INFERNAL.

11. `./run_rnacmap.sh input_features/1eiy_C`

Above command creates a folder `1eiy_C_features` in input file directory (`input_features/` in this case). `1eiy_C_features/` contains all alignment files (MSA-1, MSA-2) and features (PSSM and DCA) generated from RNAcmap pipeline.

### Run SPOT-RNA-2D for single RNA:

12. `./run.py --rna_id 1eiy_C --input_feats input_features/1eiy_C_features/ --outputs outputs/`

### Run SPOT-RNA-2D for multiple RNAs:

If RNAcmap based features (PSSM and DCA) for mutiple RNAs already exists than SPOT-RNA-2D can be run for batch of sequences as follows:

13. `./run.py --list_rna_ids datasets/TS1_ids --input_feats input_features/ --outputs outputs/`


For more options:

14. `./run.py --help`

```
usage: run.py [-h] [--rna_id] [--list_rna_ids] [--input_feats] [--single_seq]
              [--outputs] [--gpu] [--cpu]

optional arguments:
  -h, --help       show this help message and exit
  --rna_id         name of RNA sequence file; default = 6p2h_A
  --list_rna_ids   file consists of list name of RNA sequence files for batch
                   prediction; default = datasets/TS1_ids
  --input_feats    Path to directory consists of input features files with the
                   same name as specified with --rna_id flag; default =
                   inputs_features/
  --single_seq     set equal to 1 for SPOT-RNA-2D-Single; default = 0 for
                   SPOT-RNA-2D
  --outputs        Path to directory to save output files; default = outputs/
  --gpu            To run on GPU, specifiy GPU number. If only one GPU in
                   computer specifiy 0; default = -1 (no GPU)
  --cpu            Specify number of cpu threads that program can use; default
                   = 16
```

## Datasets

Datasets used for training, validation, and testing is available in [datasets](https://github.com/jaswindersingh2/SPOT-RNA-2D/tree/main/datasets) folder of this repo.

## Reproduce results in SPOT-RNA-2D paper

Refer to [benchmarking](https://github.com/jaswindersingh2/SPOT-RNA-2D/tree/main/benchmarking) folder of this repo.

## Reproduce results of tertiary modelling in the paper

[https://github.com/tlitfin/SPOT-RNA-2D-tertiary](https://github.com/tlitfin/SPOT-RNA-2D-tertiary)

## Third party programs

* cmbuild, cmcalibrate, and cmsearch from [INFERNAL tool](http://eddylab.org/infernal) version 1.1.4
* esl-reformat from [easel tool](https://anaconda.org/bioconda/easel) version 0.48
* blastn and makeblastdb from [BLAST tool](https://anaconda.org/bioconda/blast) version 2.11.0
* RNAfold from [ViennaRNA](https://anaconda.org/bioconda/viennarna) version 2.4.18
* utils/reformat.pl from [HHsuite-github-repo](https://github.com/soedinglab/hh-suite/tree/master/scripts)
* utils/getpssm.pl and utils/parse\_blastn\_local.pl from [RNAsol standalone program](https://yanglab.nankai.edu.cn/RNAsol/)
* utils/seqkit from [seqkit toolkit](https://bioinf.shenwei.me/seqkit/)
* PLMC from [plmc-github-repo](https://github.com/debbiemarkslab/plmc)

## Citation guide

**If use SPOT-RNA-2D for your research, please cite the following papers:**

Singh, J., Paliwal, K., Litfin, T., Singh, J. and Zhou, Y., 2022. Predicting RNA distance-based contact maps by integrated deep learning on physics-inferred secondary structure and evolutionary-derived mutational coupling. Bioinformatics, 38(16), pp.3900-3910.

**If use SPOT-RNA-2D input feature pipeline, please consider citing the following papers:**

RNAcmap Pipeline:

[1] Zhang, T., Singh, J., Litfin, T., Zhan, J., Paliwal, K. and Zhou, Y., 2021. RNAcmap: a fully automatic pipeline for predicting contact maps of RNAs by evolutionary coupling analysis. Bioinformatics.

RNAfold features:

[2] Lorenz, R., Bernhart, S.H., Zu Siederdissen, C.H., Tafer, H., Flamm, C., Stadler, P.F. and Hofacker, I.L., 2011. ViennaRNA Package 2.0. Algorithms for molecular biology, 6(1), pp.1-14.

INFERNAL:

[3] Nawrocki, E.P. and Eddy, S.R., 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics, 29(22), pp.2933-2935.

BLAST-N:

[4] Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. and Lipman, D.J., 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic acids research, 25(17), pp.3389-3402.

SeqKit:

[5] Shen, W., Le, S., Li, Y. and Hu, F., 2016. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PloS one, 11(10), p.e0163962.

PLMC:

[6] Hopf, T.A., Ingraham, J.B., Poelwijk, F.J., Schärfe, C.P., Springer, M., Sander, C. and Marks, D.S., 2017. Mutation effects predicted from sequence co-variation. Nature biotechnology, 35(2), pp.128-135.

**If use SPOT-RNA-2D datasets, please consider citing the following papers:**

Protein Data Bank (PDB):

[7] Berman, H.M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T.N., Weissig, H., Shindyalov, I.N. and Bourne, P.E., 2000. The protein data bank. Nucleic acids research, 28(1), pp.235-242.

CD-HIT-EST:

[8] Fu, L., Niu, B., Zhu, Z., Wu, S. and Li, W., 2012. CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics, 28(23), pp.3150-3152.

SPOT-RNA-1D:

[9] Singh, J., Paliwal, K., Singh, J. and Zhou, Y., 2021. RNA backbone torsion and pseudotorsion angle prediction using dilated convolutional neural networks. Journal of Chemical Information and Modeling.


Licence
-----
Mozilla Public License 2.0


Contact
-----
jaswinder.singh3@griffithuni.edu.au, yaoqi.zhou@griffith.edu.au
