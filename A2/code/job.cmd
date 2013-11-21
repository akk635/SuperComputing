 #!/bin/bash
##
## This script runs all the tests for the Firebenchmark, modifying the
## makefile to change between the optimization flags 
##
#@ job_name = Team12
#@ job_type = parallel
#@ class = fat
#@ wall_clock_limit = 01:20:00
#@ node = 1
#@ total_tasks = 1
#@ node_usage = not_shared
#@ initialdir = /home/hpc/h039v/di68goz/A2/code/
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification = always
#@ notify_user = ga39kax@mytum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module load metis
mpiexec -n 8 ./gccg  ../dataFiles/cojack.geo.bin cojack_classic_proc3
mpiexec -n 8 ./gccg  ../dataFiles/drall.geo.bin drall_classic_proc3
mpiexec -n 8 ./gccg  ../dataFiles/pent.geo.bin pent_classic_proc3
mpiexec -n 8 ./gccg  ../dataFiles/cojack.geo.bin cojack_dual_proc3 dual
mpiexec -n 8 ./gccg  ../dataFiles/drall.geo.bin drall_dual_proc3 dual
mpiexec -n 8 ./gccg  ../dataFiles/pent.geo.bin pent_dual_proc3 dual
mpiexec -n 8 ./gccg  ../dataFiles/cojack.geo.bin cojack_nodal_proc3 nodal
mpiexec -n 8 ./gccg  ../dataFiles/drall.geo.bin drall_nodal_proc3 nodal
mpiexec -n 8 ./gccg  ../dataFiles/pent.geo.bin pent_nodal_proc3 nodal

