## STATUS

On each push or pull request this repository automatically runs "fpm
test" to run its unit tests on Ubuntu(gfortran,ifort), MacOS(gfortran),
and MSWindows(gfortran) using github actions.

The results are able to be viewed using the **Actions** button at the
top of the github page for the repository, or can be linked to in site
documents:

+ [![Build FORD(1) docs](https://github.com/urbanjost/odepack/actions/workflows/deploy_api_docs.yml/badge.svg)](https://github.com/urbanjost/odepack/actions/workflows/deploy_api_docs.yml)
+ [![run fpm test on ubuntu with intel](https://github.com/urbanjost/odepack/actions/workflows/test_intel_ubuntu.yml/badge.svg)](https://github.com/urbanjost/odepack/actions/workflows/test_intel_ubuntu.yml)
+ [![run fpm test on ubuntu with gfortran](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_ubuntu.yml/badge.svg)](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_ubuntu.yml)
+ [![run fpm test on macos with gfortran](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_macos.yml/badge.svg)](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_macos.yml)
+ [![run fpm test on windows with gfortran](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_windows.yml/badge.svg)](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_windows.yml)
+ [![run fpm test on windows with mingw64 ](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_mingw64_windows.yml/badge.svg)](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_mingw64_windows.yml)
+ [![run fpm test on windows with msys gfortran](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_msys_windows.yml/badge.svg)](https://github.com/urbanjost/odepack/actions/workflows/test_gfortran_msys_windows.yml)
