#!/bin/bash

OUTPUT_DIR="/home/pranay/iitd/sem8/btp/minotaur/build" # directory where to store the outputs.
INST_DIR="/home/pranay/iitd/sem8/btp/minotaur/test_instances" # directory which contains the instances.
MAX_PROCS=4 # atleast 4
MPI_HOSTFILE=""
TIME_LIMIT=1800 # in seconds
EXEC="/home/pranay/iitd/sem8/btp/minotaur/build/bin/mqgmpi"

PROCS_LIST=($MAX_PROCS)
i=2
while ((i < MAX_PROCS))
do 
  PROCS_LIST+=($i)
  i=$((i * 2))
done

while read -u 9 instname
do
  for cur_procs in "${PROCS_LIST[@]}"
  do 
    output=${instname}_${cur_procs}.output
    mpirun -np $cur_procs --hostfile ${MPI_HOSTFILE} $EXEC $INST_DIR/${instname}.nl -tree_search=bfs -time_limit=${TIME_LIMIT} > $OUTPUT_DIR/${output}
  done
done 9< $1
