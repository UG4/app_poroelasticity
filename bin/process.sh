
#!/bin/sh

NUMREFS=5
NUMPROCS=1
LIMEX_TOL=0.001

function biot_eval {

cd $1


cd lastRun

grep "p0:" job.output > p0.txt
grep "LIMEX-ACCEPTING" job.output > accepting.txt
grep "LIMEX-REJECTING" job.output > rejecting.txt

# grep q=3

cd ../..

}


biot_eval "nstages2" 2
biot_eval "nstages3" 3
biot_eval "nstages4" 4