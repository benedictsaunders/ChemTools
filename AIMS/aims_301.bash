module purge
module load GCC/11.3.0  OpenMPI/4.1.4 ScaLAPACK/2.2.0-fb CMake/3.24.3

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_OFI_PROVIDER=tcp  # Only if using Intel on dedicated301/302 
ulimit -s unlimited
# Run AIMS
CORES=8
mpirun -np $CORES /path/to/aims/aims.abcdefg.x > aims.out