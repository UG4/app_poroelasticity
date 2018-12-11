
#!/bin/sh

function postprocess {

cd $1/lastRun

# grep "p0:" job.output > p0.txt
grep "LIMEX-ACCEPTING" job.output > ../accepting.txt
grep "LIMEX-REJECTING" job.output > ../rejecting.txt

grep "q=\t2" ../accepting.txt > ../accepting2.txt
grep "q=\t3" ../accepting.txt > ../accepting3.txt
grep "q=\t4" ../accepting.txt > ../accepting4.txt

grep "deltaU1" job.output > ../deltaU1.txt

grep "p0:" job.output > ../p0.txt
grep "errP:" job.output > ../errP.txt

cd ../..

}


postprocess "ref3" 
postprocess "ref4"
postprocess "ref5"
postprocess "ref6"

gnuplot "/Users/anaegel/Software/ug4-git/apps/poroelasticity/bin/deleeuw2d.gnuplot"