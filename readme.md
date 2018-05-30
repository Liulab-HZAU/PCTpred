# PCRpred

PCRpred (**P**TM **CR**oss-talk **pred**ictor) is a computational method that can effectively identify PTM cross-talk pairs in a protein based on residue pair- and residue-based features.

## Installation

### Dependencies

* Perl ( >= 5.0 ) and file::copy, statistics::descriptive, pdl, bioperl, moose modules
* Python ( >= 2.7 ) and  scikit-learn, numpy, networkx, biopython modules
* Matlab  ( >= R2014a )
* GCC ( >= 5.3.0 )

### User installation

To run PRCpred software requires the installation of some third-party tools, including BLAST+, EVfold, DisEMBL, HMMER, SIFT, PolyPhen-2, Fpocket2, DSSP, PSAIA, and HBPLUS. 

1. You need to install [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and add it to the environment variables and download [non-redundant database](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) for BLAST searches.


2. You can run the following command line to install EVfold, DisEMBL, HMMER, Fpocket2, DSSP, and PSAIA.

```shell
sh install_software.sh
```

3. SIFT and PolyPhen-2 depend on many databases (too big), so you need to install them yourself.

```shell
cd ./softwares/sift/sift6.2.1/  to install SIFT
```

```shell
cd ./softwares/polyphen-2/polyphen-2.2.2/  to install PolyPhen-2
```

4. [HBPLUS](https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/) is freely available for academic use. You need to complete the confidentiality agreement and send an email to receive the download instructions.  You can install hbplus in ./softwares/hbplus/. Of course, you may not want to do this, but you may need to modify the path in PCRpred, which will be mentioned later.

## How to use PCRpred ?

### Check the software and databases path

* Go to ./bin/ directory, and you will see a PCRpred.pl script, you need to check and modify the paths of the software and databases involved in this script.

### Input file

* The format of the input file is as follows. You can also view the prefile.txt file in the examples directory. For proteins without structural information, PDB_ID and PDB_chain are denoted by NA.

  ```tex
  P84243	T	12	K	10	2l43	A
  P62805	K	21	K	17	5ja4	B
  P37840	K	12	S	129	2n0a	B
  P46937	S	397	S	400	NA	NA
  P04637	S	215	S	392	NA	NA
  O14920	S	177	K	163	NA	NA
  ```


### Run it

* Go to ./bin/ directory, and run the following command to predict the probability score of the sample as PTM cross-talk pair.

  ```shell
  source environment_variable.sh    # add environment variable
  ```

  ```shell
  perl PCRpred.pl -f prediction_file -seq /the/directory/of/sequence/ -str /the/directory/of/PDB/ -out /the/directory/of/output/
  ```

  -f&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the path of prediction file, such as ../examples/prefile.txt

  -seq&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the path of sequence file (.fasta)

  -str&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the path of structure file (.pdb)

  -out&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the path of prediction result

* For example

  ```shell
  perl PCRpred.pl -f ../examples/prefile.txt -seq ../sequence -str ../PDB -out ../examples
  ```

### Prediction result

* You will get three output files under /the/directory/of/output/, namely predicted_results.txt, seqfeature.txt, and strfeature.txt. The former is the prediction result, and the latter two correspond to sequence features and structural features, respectively.

* predicted_results.txt contains the expected output as follows

  ```tex
  P84243	T	12	K	10	2l43	A	prediction_score	0.9778
  P62805	K	21	K	17	5ja4	B	prediction_score	0.9177
  P37840	K	12	S	129	2n0a	B	prediction_score	0.8113
  P46937	S	397	S	400	NA	NA	prediction_score	0.7034
  P04637	S	215	S	392	NA	NA	prediction_score	0.4814
  O14920	S	177	K	163	NA	NA	prediction_score	0.6852
  ```

## Help and Support

### Contact

* If you have any questions or suggestions about PCRpred, please contact us by email: XXXXXXXXXXXXXXX
* This software is free for academic use. For commercial use, please contact with the author first.

### Citation

* XXXXXXXXXXXXXXXXXXXXXXXXXXXX


