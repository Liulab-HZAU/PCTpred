#!/bin/sh

#========================add environment variable====================
PCTpred_path=$(readlink -f ..)

#domain
domain_real_path=${PCTpred_path}/software/domain
export PATH=$PATH:${domain_real_path}/hmmer-3.1b2-linux-intel-x86_64/binaries
export PERL5LIB=${domain_real_path}/PfamScan:$PERL5LIB

#PSAIA

psaia_real_path=${PCTpred_path}/software/PSAIA/PSAIA-1.0
LD_LIBRARY_PATH=${psaia_real_path}
export LD_LIBRARY_PATH

