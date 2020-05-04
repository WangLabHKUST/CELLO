# CELLO - Cancer EvoLution for LOngitudinal data


<div align=center><img width="300" src="logo.jpg" style="display: block; margin: auto;" ></div>

## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk/)
* **Software Development**: Jiguang Wang, Biaobin Jiang, Dong Song and Quanhua Mu
* **Logo Design**: Shengshuo Huang

## Status
Active Development

## Introduction
Cancer EvoLution for LOngitudinal data (CELLO) is a MATLAB toolbox for comprehensive analysis of longitudinal genomic sequencing data in cancer. It was originally developed by Jiguang Wang [1,2], and the implementation has both MATLAB and R versions:
* [MATLAB](./CELLOM/CELLOM.md)
* [R](./CELLOR/CELLOR.md)

## Docker

To ensure reproducibility and improve usability, we have prepared a [docker](https://www.docker.com/) version of CELLO, based and the R implementation. The docker can be retrieved by
```
docker pull qmu123/cellor
```
Alternatively, one can build the docker image from scratch. The docker file is uploaded to this repository. Download the [Dockerfile](./CELLOR/Dockerfile) and [CELLOR.R](./CELLOR/CELLOR.R) into a directory, then build the image from there by:
```
docker build -t cellor .
```
## Datasets

The input [SAVI report](./input.savi.txt) consists of a list of genetic variants from 90 glioblastoma patients [1].

## Reference

[1] Jiguang Wang, Emanuela Cazzato, Erik Ladewig, Veronique Frattini, Daniel I S Rosenbloom, Sakellarios Zairis, Francesco Abate, Zhaoqi Liu, Oliver Elliott, Yong-Jae Shin, Jin-Ku Lee, In-Hee Lee, Woong-Yang Park, Marica Eoli, Andrew J Blumberg, Anna Lasorella, Do-Hyun Nam, Gaetano Finocchiaro, Antonio Iavarone, Raul Rabadan. (2016). [Clonal evolution of glioblastoma under therapy.](https://www.nature.com/articles/ng.3590) **Nature Genetics**, 48(7), 768-776.

[2] Jiguang Wang, Hossein Khiabanian, Davide Rossi, Giulia Fabbri, Valter Gattei, Francesco Forconi, Luca Laurenti, Roberto Marasca, Giovanni Del Poeta, Robin Foà, Laura Pasqualucci, Gianluca Gaidano, Raul Rabadan. (2014). [Tumor evolutionary directed graphs and the history of chronic lymphocytic leukemia.](https://elifesciences.org/articles/02869) **Elife**, 3, e02869.

## Contact

For any questions, please contact Professor Jiguang Wang via email: jgwang AT ust DOT hk
