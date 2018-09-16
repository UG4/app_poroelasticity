#!/bin/sh

NUMPROCS=1
LIMEX_TOL=0.001
LIMEX_NSTAGES=4
PROBLEMID=bm2D_tri
SMOOTHER=vanka-ssc

function biot_exec {

CWD=$(pwd)

mkdir $1
cd $1
ugsubmit $NUMPROCS --- ugshell -ex poroelasticity/scripts/biot_driver.lua --problem-id bm2D_tri $PROBLEMID --limex-tol $LIMEX_TOL --limex-num-stages $LIMEX_NSTAGES --num-refs $3 --mg-smoother-type $SMOOTHER --bm-napprox $BM_TAYLOR 
# --with-debug-iter 
cd $CWD

}

NUMREFS=2
BM_TAYLOR=256
biot_exec "vanka-nstages2-refs$NUMREFS" 2 $NUMREFS

NUMREFS=3
#BM_TAYLOR=512
biot_exec "vanka-nstages2-refs$NUMREFS" 2 $NUMREFS
# biot_exec "vanka-nstages3-refs$NUMREFS" 3 $NUMREFS
# biot_exec "vanka-nstages4-refs$NUMREFS" 4 $NUMREFS

NUMREFS=4
#BM_TAYLOR=512
biot_exec "vanka-nstages2-refs$NUMREFS" 2 $NUMREFS
# biot_exec "vanka-nstages3-refs$NUMREFS" 3 $NUMREFS
# biot_exec "vanka-nstages4-refs$NUMREFS" 4 $NUMREFS

#NUMREFS=5
#biot_exec "vanka-nstages2-refs$NUMREFS" 2 $NUMREFS