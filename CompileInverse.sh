module purge
module load ncarenv/1.3
module load intel/18.0.5
module load ncarcompilers/0.5.0
module load mpt/2.19
module load netcdf/4.7.3
module list
ifort -o Inverse.exe ComputeInverse.f90
./Inverse.exe
~

