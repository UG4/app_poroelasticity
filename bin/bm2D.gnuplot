set terminal pdfcairo 

#TCHAR=1.661e-03
#UCHAR=0.073016760655994
PCHAR=1.0
UCHAR=1.0

#set output 'bm2D-accepted-timesteps.pdf'
#plot "accepting4.txt" using ($2/TCHAR):($3/TCHAR) t "q=4" with points pointtype 4, "accepting3.txt" using ($2/TCHAR):($3/TCHAR) t "q=3" with points pointtype 3, "accepting2.txt" using ($2/TCHAR):($3/TCHAR) t "q=2" with points pointtype 2


#set output 'bm2D-accuracyU.pdf'
# plot "deltaU1A.txt" using ($2/TCHAR):($6/UCHAR) with lines t "Reference", "deltaU1A.txt" using ($2/TCHAR):($5/UCHAR) with points t "Numerical solution"


set grid 

### Displacement
set output 'bm2D-L2errorP.pdf'
set xlabel 'Time'
set ylabel 'L2 error: ||p-p_h|| 
plot "refs2/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=", "refs3/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=", "refs4/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h="

### Pressure
set output 'bm2D-L2errorU.pdf'
set xlabel 'Time'
set ylabel 'L2 error: ||u-u_h||'
plot "refs2/deltaU1B.txt" using ($3):($4/UCHAR) with linespoints t "h=", "refs3/deltaU1B.txt" using ($3):($4/UCHAR) with linespoints t "h=", "refs4/deltaU1B.txt" using ($3):($4/UCHAR) with linespoints t "h="