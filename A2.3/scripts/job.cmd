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
#@ node = 4
#@ total_tasks = 16
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = /home/hpc/h039v/di68goz/A2.3/code
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification = always
#@ notify_user = ga39kax@mytum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.intel
module load metis
module load papi
export SUBJOB
for i in 1 2 4 6 8
do
  export SUBJOB=${i}
L1 = $( 1)
L2 = $(16)
sed -n -e "${L1},${L2}p" $LOADL_HOSTFILE
mpiexec -n ${i} ./gccg ./cojack.geo.bin cojack_classic_${i}
mpiexec -n ${i} ./gccg ./cojack.geo.bin cojack_metis_${i} dual
mpiexec -n ${i} ./gccg ./pent.geo.bin pent_classic_${i}
mpiexec -n ${i} ./gccg ./pent.geo.bin pent_metis_${i} dual
done
wait
