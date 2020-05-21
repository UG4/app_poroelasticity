
#!/bin/sh

NUMREFS=5
NUMPROCS=1
LIMEX_TOL=0.001

function biot_eval {

cd $1
cd lastRun

# grep "p0:" job.output > p0.txt
grep "LIMEX-ACCEPTING" job.output > ../accepting.txt
grep "LIMEX-REJECTING" job.output > ../rejecting.txt

grep "q=\t2" ../accepting.txt > ../accepting2.txt
grep "q=\t3" ../accepting.txt > ../accepting3.txt
grep "q=\t4" ../accepting.txt > ../accepting4.txt

# h1-semi
grep "deltaU1A" job.output > ../deltaU1A.txt # h1-semi
grep "deltaU2A" job.output > ../deltaU2A.txt

# l2-norm
grep "deltaP" job.output > ../deltaP.txt
grep "deltaU1B" job.output > ../deltaU1B.txt
grep "deltaU2B" job.output > ../deltaU2B.txt

cd ../..

}


biot_eval "refs2" 2
biot_eval "refs3" 3
biot_eval "refs4" 4 
biot_eval "refs5" 5 

gnuplot /Users/anaegel/Software/ug4-git/apps/poroelasticity/bin/bm2D.gnuplot 