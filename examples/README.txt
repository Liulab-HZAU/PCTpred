1. prefile.txt is an example for input to PCTpred.pl

Uniprot_ID	Residue1	Site1	Residue2	Site2	PDB_ID	PDB_Chain
P84243	T	12	K	10	2l43	A
P62805	K	21	K	17	5ja4	B
P37840	K	12	S	129	2n0a	B
P46937	S	397	S	400	NA	NA
P04637	S	215	S	392	NA	NA
O14920	S	177	K	163	NA	NA

For proteins without structural information, the PDB_ID and PDB_chain information is denoted by NA.

2. predicted_results.txt provides the prediction results

Uniprot_ID	Residue1	Site1	Residue2	Site2	PDB_ID	PDB_Chain	Probability
P84243	T	12	K	10	2l43	A	prediction_score	0.977
P62805	K	21	K	17	5ja4	B	prediction_score	0.923
P37840	K	12	S	129	2n0a	B	prediction_score	0.841
P46937	S	397	S	400	NA	NA	prediction_score	0.992
P04637	S	215	S	392	NA	NA	prediction_score	0.827
O14920	S	177	K	163	NA	NA	prediction_score	0.867

3. seqfeature.txt provides the sequence features

Sequence distance (1)
Residue co-evolution (2)
Correlated mutation (3,4,5)
Co-localization within domain (6)
Co-localization within disordered region (7)
SIFT score (8,9,10)
Residue conservation score (11,12,13)
Polyphen-2 score (14,15,16)

4. strfeature.txt provides the structural features

Structural distance (1)
Shortest path distance (2)
Co-localization within pocket (3)
Pairwise secondary structure state (4,5,6,7,8,9)
Depth and protrusion indices (10,11,12,13,14,15)
Topological features (16,17,18,19,20,21,22,23,24,25,26,27)
Laplacian norm (28,29,30,31,32,33,34,35,36,37,38,39,40,41,42)
B-factor (43,44,45)
Number of hydrogen bonds (46,47,48)

