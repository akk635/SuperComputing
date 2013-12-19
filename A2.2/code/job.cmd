 #!/bin/bash
##
## This script runs all the tests for the Firebenchmark, modifying the
## makefile to change between the optimization flags 
##
#@ job_name = Team12_mytest_A2.3
#@ job_type = parallel
#@ class = test
#@ wall_clock_limit = 01:20:00
#@ island_count = 1
#@ node = 4
#@ total_tasks = 8 
#@ node_usage = not_shared
#@ initialdir = /home/hpc/h039v/h039vaq/SuperComputing/A2.2/code 
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification = always
#@ notify_user = ga46laz@mytum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module load metis
module load papi

for i in 2 4 6 8
do
mpiexec -n ${i} ../data/cojack.geo.bin cojack_metis_${i} dual 
done
##mpiexec -n 1 ./gccg ../data/cojack.geo.bin drall_classic_1
##mpiexec -n 2 ./gccg  ../data/cojack.geo.bin drall_metis_2 dual
##mpiexec -n 3 ./gccg  ../data/cojack.geo.bin drall_metis_3 dual
##mpiexec -n 4 ./gccg  ../data/cojack.geo.bin drall_metis_4 dual
##mpiexec -n 5 ./gccg  ../data/cojack.geo.bin drall_metis_5 dual
##mpiexec -n 6 ./gccg  ../data/cojack.geo.bin drall_metis_6 dual
##mpiexec -n 7 ./gccg  ../data/cojack.geo.bin drall_metis_7 dual
##mpiexec -n 8 ./gccg  ../data/cojack.geo.bin drall_metis_8 dual


