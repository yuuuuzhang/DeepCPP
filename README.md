# DeepCPP: A deep neural network based on nucleotide bias information and minimum distribution similarity feature selection for RNA coding potential prediction
The codes and data here are used to predict RNA protein coding potential. It will give researchers usegul guidelines to explore the functions of newly found RNAs.

## PUBLICATION
Please cite this paper if using the codes here: Not available yet

## EXPLANATION
This repository contains four folders, code, rawdata, input_files and outputfiles.

### Code folder:
contains the python codes.  
```
utils_lncRNA.py -- functions will be used.  
DeepCPP.ipynb -- user interface.  
```
### rawdata folder:
This folder contains Two new human test datasets in [1].

### input_files folder:
18 models:\
human normal model: **human1.h5**,**human2.h5**,**human3.h5**\
human sORF model: **hsorf1.h5**,**hsorf2.h5**,**hsorf3.h5**\
vertebrate normal model: **vert1.h5**,**vert2.h5**,**vert3.h5**\
vertebrate sORF model: **vsorf1.h5**,**vsorf2.h5**,**vsorf3.h5**\
insect normal model: **insect1.h5**,**insect2.h5**,**insect3.h5**\
insect sORF model: **isorf1.h5**,**isorf2.h5**,**isorf3.h5**\
6 hexamer score files:\
human: **humanlnc_6mermean.csv**,**humanmrna_6mermean.csv**\
vertebrate: **vertlnc_6mermean.csv**,**vertmrna_6mermean.csv**\
insect: **insectlnc_6mermean.csv**,**insectmrna_6mermean.csv**\
2 feature set files:\
used for normal data: **normal566.csv**\
used for sORF data: **sorf450.csv**\
1 example input file:\
human sORF data: **humansorf.fa**, the first 232 samples are mRNA, the last 232 samples are lncRNA.\



### output_files folder:
This folder contains the predicted results of input file.

## WORKING MECHANISM OF DL-CRISPR
The overview of the mechanism of DeepCPP is illustrated in Figure 1. Firstly, 8 kinds of features are extracted from the raw sequence, then we use the first 566 features for normal data and 450 features for sORF data sorted by our newly proposed feature selection method, MDS, as the input feature subsets to the deep network. Finally, the average scores will be calculated from three prebuilt models and the coding potential prediction label will be given.
![alt text](https://github.com/yuuuuzhang/lncRNA/blob/master/overview.jpg)
figure 1.An overview of the mechanism of DeepCPP.
## USAGE:
Based on python3.  
Python modules:  
```
numpy  
pandas  
csv  
keras
Bio
math
```
will be used. keras is backened on tensorflow.  
Download all the files firstly, open DeepCPP.ipynb, change code:  
```
test_model('.../input_files/','.../output_files/',
           'humansorf.fa',species,type)
```
'humansorf.fa' is the input file name, it must be in .fa or .fasta format.\
species: three options are avalibale: 'human','vert' and 'insect'.\
type: 'normal' or 'sorf'. 'normal' means normaldata, 'sorf'means RNAs with short ORF.\
The predicted results is located in output_files.\

More details can be found from [1]

## REFERANCE
[1] Not available

## CONTACT
If you have any inqueries, please contact YU007@e.ntu.edu.sg
