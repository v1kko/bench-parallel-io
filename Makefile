all:
	mpif90 -O3 bench-mpi-io.f90 -o bench-mpi-io
	mpif90 -O3 -I/sw/arch/RHEL9/EB_production/2024/software/HDF5/1.14.5-gompi-2024a/include bench-hdf5-io.f90 -L/opt/hdf5/lib -lhdf5_fortran -o bench-hdf5-io 
	mpif90 -O3 bench-naive-io.f90 -o bench-naive-io
clean:
	rm -rf *.o data/*.* bench-mpi-io bench-hdf5-io bench-naive-io
