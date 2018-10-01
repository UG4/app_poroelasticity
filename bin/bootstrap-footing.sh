#!/bin/sh

NUMREFS=1
NUMPROCS=1
MG_CYCLE=W
MG_NSMOOTH=3
LIMEX_TOL=0.001
LIMEX_NSTAGES=1  				# 1: solver evaluation
PROBLEMID=footing3D   # 2D: 3-6 ref
SOLVERID=GMG #GMG # SuperLU

function biot_exec {

CWD=$(pwd)

mkdir $1
cd $1
ugsubmit $NUMPROCS --- ugshell -ex poroelasticity/scripts/biot_driver.lua\
 --num-refs $NUMREFS --problem-id $PROBLEMID \
 --limex-num-stages $LIMEX_NSTAGES --limex-tol $LIMEX_TOL \
 --solver-id $SOLVERID --mg-cycle-type $MG_CYCLE --mg-smoother-type $3 --mg-num-smooth $MG_NSMOOTH\
 $UG4_EXT_PARAMS
cd $CWD

}

SMOOTHER="vanka-ssc"
#SMOOTHER="uzawa3"
# SMOOTHER="sgs"




function HDependenceTest {
DIRNAME=$1
echo "Creating $DIRNAME"
mkdir $DIRNAME
NUMREFS=1
while [ $NUMREFS -le 3 ]
do
	biot_exec "$DIRNAME/ref$NUMREFS" 2 $SMOOTHER
	((NUMREFS++))
done
}




#OPTIONAL:
# --with-transient-mechanics; 	
# --use-rap
# --with-vtk

MG_CYCLE=W; 
UG4_EXT_PARAMS="--orderU 1 --stab 0.0833333"; HDependenceTest "P1P1stab_ass_gmg_cycleW"
UG4_EXT_PARAMS="--orderU 2"; HDependenceTest "P1P2_ass_gmg_cycleW" 

#MG_CYCLE=F; 
#UG4_EXT_PARAMS="--orderU 1 --stab 0.0833333"; HDependenceTest "P1P1stab_ass_gmg_cycleF"
#UG4_EXT_PARAMS="--orderU 2 "; HDependenceTest "P1P2_ass_gmg_cycleF" 

