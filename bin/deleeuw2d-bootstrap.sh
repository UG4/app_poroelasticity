#!/bin/sh

NUMREFS=1
NUMPROCS=1
LIMEX_TOL=0.01
PROBLEMID=deleeuw2d
SOLVERID=GMGKrylov # SuperLU

function biot_exec {

CWD=$(pwd)

mkdir $1
cd $1
ugsubmit $NUMPROCS --- ugshell -ex poroelasticity/scripts/biot_driver.lua\
 --num-refs $NUMREFS --solver-id $SOLVERID\
 --limex-num-stages $2 --limex-tol $LIMEX_TOL\
 --mg-smoother-type $3 --problem-id $PROBLEMID\
  --orderU 1 --stab 0.5 # --use-rap
 # --orderU 1 --stab 0.0833333  # --use-rap


 
cd $CWD

}

SMOOTHER="vanka-ssc"
SMOOTHER="uzawa3"
# SMOOTHER="sgs"


DIRNAME="P1P1stab-tol1e-2"
mkdir $DIRNAME
NUMREFS=2
#biot_exec "$DIRNAME/ref2" 2 $SMOOTHER
NUMREFS=3
biot_exec "$DIRNAME/ref3" 2 $SMOOTHER
NUMREFS=4
biot_exec "$DIRNAME/ref4" 2 $SMOOTHER
NUMREFS=5
biot_exec "$DIRNAME/ref5" 2 $SMOOTHER
NUMREFS=6
biot_exec "$DIRNAME/ref6" 2 $SMOOTHER
#NUMREFS=7
#biot_exec "$DIRNAME/ref7" 1 $SMOOTHER