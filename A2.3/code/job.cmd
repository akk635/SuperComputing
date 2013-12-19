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
#@ node = 1
#@ total_tasks = 8
#@ node_usage = not_shared
#@ initialdir = /home/hpc/h039v/di68goz/A2.2/code
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification = always
#@ notify_user = ga39kax@mytum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module load metis
module load papi
mpiexec -n 8 ./gccg ../data/cojack.geo.bin drall_classic_8
mpiexec -n 1 ./gccg ../data/cojack.geo.bin drall_classic_1
mpiexec -n 8 ./gccg  ../data/cojack.geo.bin drall_metis_8 dual
mpiexec -n 1 ./gccg  ../data/cojack.geo.bin drall_metis_1 dual
