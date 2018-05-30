#!/bin/sh

#========================add environment variable====================
PCRpred_path=$(readlink -f ..)

#domain
domain_real_path=${PCRpred_path}/software/domain
export PATH=$PATH:${domain_real_path}/hmmer-3.1b2-linux-intel-x86_64/binaries
export PERL5LIB=${domain_real_path}/PfamScan:$PERL5LIB

#dp_indel

