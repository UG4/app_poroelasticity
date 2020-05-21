#!/bin/sh

NUMPROCS=1
LIMEX_TOL=0.001
LIMEX_NSTAGES=4
PROBLEMID=bm2D_tri
SMOOTHER=vanka-ssc
SMOOTHER=uzawa3
SOLVERID=GMGKrylov

function biot_exec {

CWD=$(pwd)

mkdir $1
cd $1
ugsubmit $NUMPROCS --- ugshell -ex poroelasticity/scripts/biot_driver.lua \
--problem-id $PROBLEMID --num-refs $3 \
--solver-id $SOLVERID --mg-smoother-type $SMOOTHER \
--limex-tol $LIMEX_TOL --limex-num-stages $LIMEX_NSTAGES \
--bm-napprox $BM_TAYLOR \
--orderU 1 --stab 0.5
# --with-debug-iter 
cd $CWD

}

DIRNAME="P1P1stab-uzawa3-nu04"
#mkdir $DIRNAME
NUMREFS=2
BM_TAYLOR=256
# biot_exec "$DIRNAME/refs$NUMREFS" 2 $NUMREFS

NUMREFS=3 #BM_TAYLOR=512
# biot_exec "$DIRNAME/refs$NUMREFS" 2 $NUMREFS
# biot_exec "vanka-nstages3-refs$NUMREFS" 3 $NUMREFS
# biot_exec "vanka-nstages4-refs$NUMREFS" 4 $NUMREFS

NUMREFS=4 #BM_TAYLOR=512
# biot_exec "$DIRNAME/refs$NUMREFS" 2 $NUMREFS
# biot_exec "vanka-nstages3-refs$NUMREFS" 3 $NUMREFS
# biot_exec "vanka-nstages4-refs$NUMREFS" 4 $NUMREFS

#NUMREFS=5
#biot_exec "vanka-nstages2-refs$NUMREFS" 2 $NUMREFS
NUMREFS=5 #BM_TAYLOR=512
biot_exec "$DIRNAME/refs$NUMREFS" 2 $NUMREFS
