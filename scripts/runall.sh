#!/bin/bash
export GFORTRAN_FLAGS='-fallow-argument-mismatch -std=f2018 -finit-real=snan -Wall -fbacktrace -fcheck=bounds'
export IFORT_FLAGS=' -check bounds  -check uninit'
(
exec 2>&1
 rm ccout
 fpm run                     --compiler nvfortran --example '*' --verbose

 rm ccout
 fpm run --profile release  --flag "$IFORT_FLAGS"  --compiler ifort     --example '*' --verbose

 rm ccout
 fpm run --profile debug --flag "$GFORTRAN_FLAGS"  --compiler gfortran  --example '*' --verbose

 rm ccout
) |tee example.out
(
exec 2>&1
 fpm test                     --compiler nvfortran  --verbose
 fpm test --profile release  --flag "$IFORT_FLAGS"  --compiler ifort      --verbose
 fpm test --flag "$GFORTRAN_FLAGS"  --compiler gfortran   --verbose
)|tee test.out
exit

mkdir -p /tmp/NEW
cd /tmp/NEW
git clone  https://github.com/jacobwilliams/odepack
cd odepack
quiet bash scripts/diff.sh >/tmp/x.out
exit
