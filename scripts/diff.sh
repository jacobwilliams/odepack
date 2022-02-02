#!/bin/bash
set -x
exec 2>&1
################################################################################
# fpm bug? segfault
#fpm run --profile release --flag -g0 --compiler ifort --example '*'
#fpm run --profile release --flag -g0 --compiler nvfortran --example '*'
################################################################################
SKIP(){
fpm run --profile release --flag -g0 --compiler ifort 
fpm run --profile release --flag -g  --compiler nvfortran 
}
################################################################################
IFORT(){
   [ -f ccout ] && rm -f ccout 
   #fpm run --profile release --flag -g0 --compiler ifort     --verbose --example $NAME >output/ifort/$NAME.txt
   fpm run                  --flag -cpp  --compiler ifort               --example $NAME >output/ifort/$NAME.txt
   [ -f ccount ] && mv ccout output/ifort/ccout.txt
}
################################################################################
NVFORTRAN(){
   [ -f ccout ] && rm -f ccout 
   fpm run --profile release --flag -g  --compiler nvfortran --verbose --example $NAME >output/nvfortran/$NAME.txt
   [ -f ccount ] && mv ccout output/nvfortran/ccout.txt
}
################################################################################
GFORTRAN(){
#      --flag -Wno-argument-mismatch \
#      --flag -Wno-implicit-interface \
#      --flag -Wno-aliasing \
#      --flag -Wno-all \
#      --flag -Wno-array-bounds \
#      --flag -Wno-c-binding-type \
#      --flag -Wno-character-truncation \
#      --flag -Wno-conversion \
#      --flag -Wno-do-subscript \
#      --flag -Wno-function-elimination \
#      --flag -Wno-implicit-interface \
#      --flag -Wno-implicit-procedure \
#      --flag -Wno-intrinsic-shadow \
#      --flag -Wno-intrinsics-std \
#      --flag -Wno-line-truncation \
#      --flag -Wno-real-q-constant \
#      --flag -Wno-surprising \
#      --flag -Wno-underflow \
#      --flag -Wno-unused-parameter \
#      --flag -Wno-use-without-only \
   [ -f ccout ] && rm -f ccout 
   fpm run \
      --flag -g  \
      --flag -std=legacy \
      --compiler gfortran  --verbose --example $NAME >output/gfortran/$NAME.txt
   [ -f ccount ] && mv ccout output/gfortran/ccout.txt
}
################################################################################

for NAME in lsoda lsodar lsode lsodes lsodi lsodis lsodkr lsodpk lsoibt
do
   time IFORT
   time NVFORTRAN
   time GFORTRAN
done
diff -r output/ifort output/ifort_original
exit
