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
#@ initialdir = /home/hpc/h039v/h039vaq/SuperComputing/Assignment1/FireBenchmark/flops/o3-no-prefetch
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification = always
#@ notify_user = ga46laz@mytum.de

#@ queue

module load papi

for i in 1 2 3 
do

./gccg txt cojack.dat run${i}_
./gccg txt drall.dat run${i}_
./gccg txt pent.dat run${i}_
./gccg txt tjunc.dat run${i}_

done
