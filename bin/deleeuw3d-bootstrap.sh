#!/bin/sh

NUMREFS=1
NUMPROCS=1
LIMEX_TOL=0.001
PROBLEMID=deleeuw3dTet


function biot_exec {

CWD=$(pwd)

mkdir $1
cd $1
ugsubmit $NUMPROCS --- ugshell -ex poroelasticity/scripts/biot_driver.lua --limex-num-stages $2 --limex-tol $LIMEX_TOL --num-refs $NUMREFS --mg-smoother-type $3 --problem-id $PROBLEMID
cd $CWD

}

SMOOTHER=vanka-ssc
SMOOTHER=uzawa3

mkdir "deLeeuw3D"
NUMREFS=1
biot_exec "deLeeuw3D/ref1" 2 $SMOOTHER
NUMREFS=2
biot_exec "deLeeuw3D/ref2" 2 $SMOOTHER

