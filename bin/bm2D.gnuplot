set terminal pdfcairo 

TCHAR=0.0014660765716752
UCHAR=0.018746620403829
PCHAR=142857+2.0*35714.3
UCHAR=1.0

set grid 

### Pressure
set output 'bm2D-L2errorP.pdf'
set xlabel 'Normalized time' # t/({/Symbol p} t^*)'
set ylabel 'L2 error: ||p-p_h|| 
# plot "refs2/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=1/16", "refs3/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=1/32", "refs4/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=1/64"
plot "refs2/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=1/16", "refs3/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=1/32", "refs4/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=1/64", "refs5/deltaP.txt" using ($3):($4/PCHAR) with linespoints t "h=1/128"

### Displacement
set output 'bm2D-L2errorU.pdf'
set xlabel 'Normalized time'
set ylabel 'L2 error: ||u-u_h||'
plot "refs2/deltaU1B.txt" using ($3):($4/UCHAR) with linespoints t "h=1/16", "refs3/deltaU1B.txt" using ($3):($4/UCHAR) with linespoints t "h=1/32", "refs4/deltaU1B.txt" using ($3):($4/UCHAR) with linespoints t "h=1/64"

set output 'bm2D-H1errorU.pdf'
set xlabel 'Normalized time'
set ylabel 'H1 error: ||u-u_h||'
plot "refs2/deltaU1B.txt" using ($3):($4/UCHAR) with linespoints t "h=1/16", "refs3/deltaU1A.txt" using ($3):($4/UCHAR) with linespoints t "h=1/32", "refs4/deltaU1A.txt" using ($3):($4/UCHAR) with linespoints t "h=1/64"

set output 'bm2D-steps.pdf'
set xlabel 'Normalized time'
set ylabel 'Step size'
plot "refs2/accepting.txt" using ($2/TCHAR):($4/TCHAR) t "h=1/16", "refs3/accepting.txt" using ($2/TCHAR):($4/TCHAR) t "h=1/32", "refs4/accepting.txt" using ($2/TCHAR):($4/TCHAR)  t "h=1/64"