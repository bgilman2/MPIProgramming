module purge
module load ncarenv/1.3
module load pgi/19.9
module load ncarcompilers/0.5.0
module load netcdf/4.7.3
module load openmpi/3.1.4
module list
pgfortran ComputeInverse.f90 -o ComputeInverse.exe -llapack -lblas
./ComputeInverse.exe
~

