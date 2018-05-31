#!/bin/sh

#This script is used to download and install the software required for PCRpred,
#including calculation of NWalign, co-evolution, disorder, domain, pocket,
#secondary structure, depth and protrusion index.

cd software
#==========================download softwares========================
#NWalign
wget https://zhanglab.ccmb.med.umich.edu/NW-align/NWalign.gz

#co-evolution
wget http://evfold.org/evfold-web/download/calculate_evolutionary_constraints_v2.0.tar.gz

#disorder region
wget http://dis.embl.de/DisEMBL-1.4.tgz
wget https://www.pks.mpg.de/~tisean/TISEAN_2.1-linux.tar.gz

#domain region
#download Pfam-A database/246M
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz
#download PfamScan/28K
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz
#download Hmmer/20M
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz

#pocket region
wget https://sourceforge.net/projects/fpocket/files/latest/download/fpocket2.tar.gz

#secondary structure
wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64

#depth and protrusion index
wget http://complex.zesoi.fer.hr/data/PSAIA-1.0.tar.gz

#sift
wget http://siftdna.org/www/sift/public/sift6.2.1.tar.gz

#polyphen-2
wget http://genetics.bwh.harvard.edu/pph2/dokuwiki/_media/polyphen-2.2.2r405c.tar.gz

#=================================install===================================
#NWalign
cd sequence_alignment
mv ../NWalign.gz ./
gzip -d NWalign.gz
chmod +x NWalign

#co-evolution
cd ..
tar -zxvf calculate_evolutionary_constraints_v2.0.tar.gz
rm calculate_evolutionary_constraints_v2.0.tar.gz

#disorder region
tar -zxvf DisEMBL-1.4.tgz -C disorder
tar -zxvf TISEAN_2.1-linux.tar.gz -C disorder
rm DisEMBL-1.4.tgz
rm TISEAN_2.1-linux.tar.gz
cp disorder/bin-linux/sav_gol disorder/DisEMBL-1.4/
cd disorder/DisEMBL-1.4
gcc -O3 disembl.c -o disembl
disorder_real_path=$(readlink -f ..)
disorder_real_path1=${disorder_real_path//\//\\/}
sed -i "s/\/PATH/$disorder_real_path1/g" DisEMBL.py
sed -i "s/cur_record.seq.data/cur_record.seq/" DisEMBL.py

#domain region
cd ../../
tar -zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz -C domain
rm hmmer-3.1b2-linux-intel-x86_64.tar.gz
cd domain/hmmer-3.1b2-linux-intel-x86_64
./configure
make
make check
domain_real_path=$(readlink -f .)
./configure --prefix $domain_real_path
make install
cd easel
make install
cd ../binaries/
domain_bin_real_path=$(readlink -f .)
export PATH=$PATH:$domain_bin_real_path
cd ../../
mv ../active_site.dat.gz ./
mv ../Pfam-A.hmm.dat.gz ./
mv ../Pfam-A.hmm.gz ./
mv ../PfamScan.tar.gz ./
gzip -d *.gz
hmmpress Pfam-A.hmm
tar -vxf PfamScan.tar
rm PfamScan.tar

#pocket region
cd ../fpocket
mv ../fpocket2.tar.gz ./
tar -zxvf fpocket2.tar.gz
rm fpocket2.tar.gz
cd fpocket2
make
make test

#secondary structure
cd ../../dssp/
mv ../dssp-2.0.4-linux-amd64 ./
chmod +x dssp-2.0.4-linux-amd64

#depth and protrusion index
cd ../PSAIA/
mv ../PSAIA-1.0.tar.gz ./
tar -zxvf PSAIA-1.0.tar.gz
rm PSAIA-1.0.tar.gz
mv psa.cfg ./PSAIA-1.0/
mv libpng12.so.0 ./PSAIA-1.0/
cd PSAIA-1.0/
dp_path=$(readlink -f .)
dp_path1=${dp_path//\//\\/}
sed -i "s/\/PATH/$dp_path1/g" psa.cfg
cd ../../

#==================================end===============================

#You need to install sift and polyphen-2 yourself
tar -zxvf sift6.2.1.tar.gz -C sift
rm sift6.2.1.tar.gz
tar -zxvf polyphen-2.2.2r405c.tar.gz -C polyphen-2
rm polyphen-2.2.2r405c.tar.gz






