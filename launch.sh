
#!/bin/bash
set -euo pipefail
source ../env.source
rm pastis-*.log ./build/sim_mat.mtx || true
rm /home/lukb/git/ELBA/ELBA/bin/codelets/algoipu.gp || true
pushd build_release/

#LOWER_KMER_FREQ=20
#UPPER_KMER_FREQ=30
#FILE=runelba

# LOWER_KMER_FREQ=2
# UPPER_KMER_FREQ=4
# FILE=runelba_100


LOWER_KMER_FREQ=31
UPPER_KMER_FREQ=40
FILE=runelba_celegans

cmake .. -DLOWER_KMER_FREQ=${LOWER_KMER_FREQ} -DUPPER_KMER_FREQ=${UPPER_KMER_FREQ}

make -j$(nproc)
cp /home/lukb/git/ELBA/ELBA/build_release/bin/codelets/algoipu.gp /home/lukb/git/ELBA/ELBA/bin/codelets/algoipu.gp
popd

sbatch $FILE 
