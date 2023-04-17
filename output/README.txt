
reference      - reference output that was in original distribution
ifort_original - using ifort with default compiler options and original fortran 77 fixed-form code
gfortran       - latest run with gfortran
ifort          - latest run with ifort
nvfortran      - latest run with nvfortran

Do not have extensive unit tests, but seeing values change with different
compilers and compiler versions that may or may not be significant so
until have better criteria using the ifort_original as a base set of
values for building with ifort with no options.


