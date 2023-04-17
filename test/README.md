
The example programs contained in the comments of the original package
updated to more modern syntax and to use the updated ODEPACK:

 + dlsoda_ex.f90
 + dlsodar_ex.f90
 + dlsode_ex.f90
 + dlsodes_ex.f90
 + dlsodi_ex.f90
 + dlsodis_ex.f90
 + dlsoibt_ex.f90

tested using

 + GNU Fortran (Ubuntu 10.3.0-1ubuntu1~20.04) 10.3.0
 + nvfortran 21.5-0 LLVM 64-bit target on x86-64 Linux -tp nehalem 
 + ifort (IFORT) 2021.3.0 20210609

using the commands

```bash
fpm test --compiler gfortran --profile release -flag -std=legacy|tee gf.
fpm test --compiler nvfortran --profile release |tee nv.
fpm test --compiler ifort --profile release |tee if.
```
