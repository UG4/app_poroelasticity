#!/bin/sh

NUMREFS=5
NUMPROCS=1
LIMEX_TOL=0.001

function biot_exec {

mkdir $1
cd $1
ugsubmit $NUMPROCS --- ugshell -ex poroelasticity/scripts/biot_driver.lua --limex-num-stages $2 --limex-num-stages $2 --limex-tol $LIMEX_TOL --num-refs $NUMREFS
cd ..

}


biot_exec "nstages2" 2
biot_exec "nstages3" 3
biot_exec "nstages4" 4