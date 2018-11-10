#!/bin/bash
#SBATCH --account=m888  	# Your repo goes here
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=stream-ai

# VTune SEP driver is required
#SBATCH --perf=vtune

# for Cori, set appropriately. Not required for Edison
#SBATCH --constraint="knl"
##SBATCH --constraint="haswell"

### start of script configuration parameters

# set to yes or no to select tests to run
run_baseline=yes
run_sde=yes
run_vtune=yes

# use -knl for Cori KNL, -hsw for Cori Haswell or -ivb for Edison
if [ $NERSC_HOST == "cori" ]; then
  if [ $SLURM_JOB_CPUS_PER_NODE == "272" ]; then # KNL node
    myhost=cori_knl
    SDE='sde -knl'

    # determin the number of MPI ranks per node
    n=$SLURM_JOB_NUM_NODES

    # number of threads/rank and virtual cores (needed for srun's -c)
    t=64
    vcores=$(( t * 4 ))
  else 						# Haswell node
    myhost=cori_hsw
    SDE='sde -hsw'
    n=$(( 2 * $SLURM_JOB_NUM_NODES ))
    t=16
    vcores=$(( t * 2 ))
  fi

elif [ $NERSC_HOST == "edison" ]; then
  myhost=edison
  SDE='sde -ivb'
  n=$(( 2 * $SLURM_JOB_NUM_NODES ))
  t=12
  vcores=$(( t * 2 ))
fi

#module load sde	# requires version 8.4.0 or later
#module load vtune	# script setup for Vtune 2017 or later

### End of configuration parameters

echo "Running with $n MPI ranks and $t threads"
export OMP_NUM_THREADS=$t
suffix=${n}p${t}t_${SLURM_JOB_ID}
exe=./stream_mpi.exe

if [ "$run_baseline" == "yes" ]; then
  echo ""
  echo "--------------------------------------------------"
  echo "----->> Running Stream w/o Instrumentation <<-----"
  echo "--------------------------------------------------"
  srun -n $n -c $vcores --cpu_bind=cores $exe
fi

if [ "$run_sde" == "yes" ]; then
  echo ""
  echo "--------------------------------------------------"
  echo "----->> Running w/SDE <<-----"
  echo "--------------------------------------------------"
  srun -n $n -c $vcores --cpu_bind=cores $SDE -d -iform 1 -omix sde_${suffix}.out -i -top_blocks 500 -global_region -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- $exe
  echo "----->> Generating SDE Report <<-----"
  echo "For performance, the SDE report is best done on an external login node"
  echo "Run the following command: "
  echo "\$ ./parse-sde.sh sde_${suffix}.out*"
fi

if [ "$run_vtune" == "yes" ]; then
  echo ""
  echo "--------------------------------------------------"
  echo "----->> Running w/Vtune <<-----"
  echo "--------------------------------------------------"
  srun -n $n -c $vcores --cpu_bind=cores amplxe-cl -start-paused -r vtbw_${suffix} -collect memory-access -finalization-mode=none -trace-mpi -- $exe
  echo "----->> Finalizing VTune and generating report <<-----"
  echo "For performance, the finalize and report are best done on an external login node"
  echo "Run the following commands: "
  echo "Note that if using Vtune version 2017 replace \"-report hw-events -group-by=package\" with \"-report summary\" "
  if [ $myhost == "cori_knl" ]; then
    echo "\$ amplxe-cl -report hw-events -group-by=package -r vtbw_${suffix} -column=UNC_M_CAS_COUNT,UNC_E_RPQ_INSERTS,UNC_E_WPQ_INSERTS -format=csv -csv-delimiter=comma > vtbw_${suffix}.summary"
  else
    echo "\$ amplxe-cl -report hw-events -group-by=package -r vtbw_${suffix} -column=UNC_M_CAS_COUNT -format=csv -csv-delimiter=comma > vtbw_${suffix}.summary"
  fi
  echo "\$ ./parse-vtune2018.sh vtbw_${suffix}.summary"
fi

