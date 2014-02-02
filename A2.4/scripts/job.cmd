 #!/bin/bash
##
## This script runs all the tests for the Firebenchmark, modifying the
## makefile to change between the optimization flags
##
#@ job_name = Team12_mytest_A2.3
#@ job_type = MPICH
#@ class = test
#@ wall_clock_limit = 01:20:00
#@ island_count = 1
#@ node = 2
#@ total_tasks = 64
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = /home/hpc/h039v/h039vaq/SuperComputing/A2.3/code
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification = always
#@ notify_user = ga46laz@mytum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.intel
module load metis
module load papi
module load scalasca
module load valgrind
export SUBJOB
for i in 1 2 4 8 16 32 64
do
  export SUBJOB=${i}
L1 = $( 1)
L2 = $($i)
#sed -n -e "${L1},${L2}p" $LOADL_HOSTFILE
mpiexec -n ${i} ./gccg ../data/cojack.geo.bin cojack_classic
mpiexec -n ${i} ./gccg ../data/cojack.geo.bin cojack_metis dual
#mpiexec -n ${i} ./gccg ../data/cojack.geo.bin cojack_myclassic_${i} myclassical
mpiexec -n ${i} ./gccg ../data/pent.geo.bin pent_metis dual
mpiexec -n ${i} ./gccg ../data/pent.geo.bin pent_classic
#mpiexec -n ${i} ./gccg ../data/pent.geo.bin pent_myclassic_${i} myclassical
done
#wait

