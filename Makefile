all:
	mpif90 bench-mpi-io.f90 -o bench-mpi-io
	mpif90 -I/opt/hdf5/include bench-hdf5-io.f90 -L/opt/hdf5/lib -lhdf5_fortran -o bench-hdf5-io
	mpif90 bench-naive-io.f90 -o bench-naive-io
clean:
	rm -rf *.o data/*.* bench-mpi-io bench-hdf5-io bench-naive-io
