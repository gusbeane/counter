#!/bin/bash

export BASEPATH=/n/holylfs05/LABS/hernquist_lab/Users/abeane/counter/

export AREPO_commit=956184545147294eedcbd80516b7e635faae60cc
export REPLACE=flyV6a
export ICS_DIR=${BASEPATH}/ics/SMUGGLE-flyby-R200-Vphi60-ang16-pro

#source load-modules.sh

# create arepo repo

#git clone git@bitbucket.org:volkerspringel/arepo.git

#cd arepo
#git checkout ${AREPO_commit}
#cp ../../Config-SMUGGLE.sh Config.sh
#cp ../../Makefile.systype .
# make -j # commented out bc will be recompiled at runtime anyways
#cd ../

#for i in 5 4 3 2 1
for i in 4 
do
for ang in 0.0625 0.125 0.1875 0.25 0.3125 0.375 0.4375 0.5 0.5625 0.625 0.6875 0.75 0.8125 0.875 0.9375
#for ang in 0.0
do
    cd lvl${i}-ang${ang}
    sbatch job_lvl${i}.sh
    cd ../
done
done
