# PCTpred

PCTpred (**P**TM **C**ross-**T**alk **pred**ictor) is a computational method that can accurately identify PTM cross-talk pairs in a given protein sequence or structure. To this end, we first design a group of novel residue pair- and residue-based features which effectively show the preferences of cross-talk sites from both the sequence and structural perspectives. Then we construct two component predictors using random forest with finely selected sequence- and structure-based features. Further combining these two predictors, PCTpred takes full advantage of the complementarity between residue pair- and residue-based features and that between sequence and structural information. Using both pair- and protein-based evaluations, PCTpred yields consistently better results than the state-of-the-art methods.

## Installation

### Dependencies

* Perl ( >= 5.0 ) and file::copy, statistics::descriptive, pdl, bioperl, moose modules
* Python ( >= 2.7 ) and  scikit-learn, numpy, networkx, biopython modules
* Matlab  ( >= R2014a )
* GCC ( >= 5.3.0 )

### Download

```shell
git clone https://github.com/Liulab-HZAU/PCTpred.git
```

### User installation

To install PCTpred, you need to conduct the installation of some third-party software, including BLAST+, EVfold, DisEMBL, HMMER, SIFT, PolyPhen-2, Fpocket2, DSSP, PSAIA, and HBPLUS.

1. You need to install [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and add it into the environment variables and download [non-redundant database](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for BLAST searches.


2. You can run the following shell script to install EVfold, DisEMBL, HMMER, Fpocket2, DSSP, and PSAIA.

```shell
sh install_software.sh
```

3. SIFT and PolyPhen-2 depend on several large databases, and you need to install them by yourself.

```shell
cd ./software/sift/sift6.2.1/  to install SIFT
```

```shell
cd ./software/polyphen-2/polyphen-2.2.2/  to install PolyPhen-2
```

4. [HBPLUS](https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/) is freely available for academic use. You need to complete the confidentiality agreement and send an email to receive the download instructions.  You can install hbplus in `./software/hbplus/`.  If you want to skip this step, you may need to modify the path in PCTpred mentioned later.

## How to use PCTpred ?

### Check the software and databases path

* You will see a `PCTpred.pl` script in the `./bin/` directory, and you need to check and modify the paths related to specific software and databases involved in this script.

### Input file

* The format of the input file is as follows. You can also check the `prefile.txt` file in the `examples` directory. For proteins without structural information, the PDB_ID and PDB_chain information is denoted by NA.

```tex
P84243	T	12	K	10	2l43	A
P62805	K	21	K	17	5ja4	B
P37840	K	12	S	129	2n0a	B
P46937	S	397	S	400	NA	NA
P04637	S	215	S	392	NA	NA
O14920	S	177	K	163	NA	NA
```
### Run it

* Go to `./bin/` directory, and run the following command to predict the probability score of samples being PTM cross-talk pairs.

  ```shell
  source environment_variable.sh    # add environment variable
  ```

  ```shell
  perl PCTpred.pl -f prediction_file -seq /the/directory/of/sequence/ -str /the/directory/of/PDB/ -out /the/directory/of/output/
  ```

  -f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the path of prediction file, such as ../examples/prefile.txt

  -seq&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the path of sequence file (.fasta)

  -str&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the path of structure file (.pdb)

  -out&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the path of prediction result

* For example

  ```shell
  perl PCTpred.pl -f ../examples/prefile.txt -seq ../sequences/ -str ../PDB/ -out ../examples/
  ```

### Prediction result

* You will obtain three output files in the directory `/the/directory/of/output/`, such as `predicted_results.txt`, `seqfeature.txt`, and `strfeature.txt`. The first file provides the prediction results and the second and third files provide the sequence features and structural features for each PTM pair generated by our method, respectively.
* `predicted_results.txt` contains the expected output as follows

```tex
P84243	T	12	K	10	2l43	A	prediction_score	0.977
P62805	K	21	K	17	5ja4	B	prediction_score	0.923
P37840	K	12	S	129	2n0a	B	prediction_score	0.841
P46937	S	397	S	400	NA	NA	prediction_score	0.992
P04637	S	215	S	392	NA	NA	prediction_score	0.827
O14920	S	177	K	163	NA	NA	prediction_score	0.867
```
## Help and Support

### Contact

* If you have any questions or suggestions about PCTpred, please contact us by email: hfliu@webmail.hzau.edu.cn or liurong116@mail.hzau.edu.cn
* This software is free for academic use. For commercial use, please contact with the authors.

### Citation

* Hui-Fang Liu and Rong Liu. Structure-based prediction of post-translational modification cross-talk within proteins using complementary residue- and residue pair-based features. _Briefings in Bioinformatics_. 2018, in press.

