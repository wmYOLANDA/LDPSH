# LDPSH
This is the matlab implementation of LDPSH method proposed in our paper "A General Framework for Linear Distance Preserving Hashing". The descriptions of files in this directory are listed below:

- `demo.m`: gives an example to show how to run our LDPSH method.
- `/ITQ`: the referred unsupervised hashing method ITQ, you can replace it with other hashing methods.
- `/TSH` contains the main implementation of the proposed approach `LDPSH`.
- `/NN`: contains the implementation of network architectures.
- `/utils`: contains the implementations for basic manipulations.

Data Preparation
---------------
Please go to http://corpus-texmex.irisa.fr/ to download the experimental datasets ANN_SIFT1M and ANN_GIST1M.

Citation
---------------

@article{wang2018general,
title={A General Framework for Linear Distance Preserving Hashing},
author={Wang, Min and Zhou, Wengang and Tian, Qi and Li, Houqiang},
journal={IEEE Transactions on Image Processing},
volume={27},
number={2},
pages={907--922},
year={2018},
publisher={IEEE}
}

@inproceedings{wang2016linear,
  title={Linear Distance Preserving Pseudo-Supervised and Unsupervised Hashing},
  author={Wang, Min and Zhou, Wengang and Tian, Qi and Zha, Zhengjun and Li, Houqiang},
  booktitle={Proceedings of the 2016 ACM on Multimedia Conference},
  pages={1257--1266},
  year={2016},
  organization={ACM}
}

Contact
---------------
Please feel free to leave suggestions or comments to Wang Min (wm123@mail.ustc.edu.cn).
