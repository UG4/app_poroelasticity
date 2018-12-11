set terminal pdfcairo color

T3="Refinement level 3"
T4="4"
T5="5"
T6="6"

set key outside

set ylabel 'Normalized error norm'

# Control: p0 vs. time
set output 'deleeuw2d-p0control.pdf'
set xlabel 'Time'
set ylabel 'Normalized pressure p(x_0)'

set logscale x
set xrange [1e-5:2]
plot "ref3/p0.txt" using 3:4 with lines t "Reference", "ref3/p0.txt" using 3:4 with points t "Simulation"

################
# Punktweiser Fehler im Zentrum.
##############
set output 'deleeuw2d-p0error-pointwise.pdf'
set yrange [-0.02:0.02]
set grid
set ylabel 'Pointwise error p_h(x_0)-p(x_0)'
plot "ref3/p0.txt" using 3:6 t T3, "ref4/p0.txt" using 3:6 t T4, "ref5/p0.txt" using 3:6 t T5, "ref6/p0.txt" using 3:6 t '6',0 notitle

##################
# L2-Fehler
##################
set output 'deleeuw2d-l2error-log.pdf'
set grid

set xlabel 'Time'
set logscale x
set xrange [1e-5:2]

set yrange [0:0.04]
set ylabel 'L2 error ||p-p_h||'

plot "ref3/errP.txt" using 3:4, "ref4/errP.txt" using 3:4, "ref5/errP.txt" using 3:4, "ref6/errP.txt" using 3:4

# L2 Fehler, logx liny
set output 'deleeuw2d-l2error-linear.pdf'
set xlabel 'Time'
unset logscale x
set xrange [0.001:0.5]

set ylabel 'L2 error ||p-p_h||'
plot "ref3/errP.txt" using 3:4, "ref4/errP.txt" using 3:4, "ref5/errP.txt" using 3:4, "ref6/errP.txt" using 3:4


# L2 Fehler, loglog
set output 'deleeuw2d-l2error-loglog.pdf'
set xrange [1e-5:2]
unset yrange
set xlabel 'Time'
set logscale x
set logscale y


set ylabel 'L2 error ||p-p_h||'
plot "ref3/errP.txt" using 3:4, "ref4/errP.txt" using 3:4, "ref5/errP.txt" using 3:4, "ref6/errP.txt" using 3:4

##################
# Accepted steps
##################
CHAR_TIME=0.0010416666666667
set output 'deleeuw2d-accepting2-time.pdf'
set xlabel 'Time'
unset logscale x
set ylabel 'Time step size'
set logscale y
unset xrange
plot "ref3/accepting2.txt" using ($2/CHAR_TIME):($3/CHAR_TIME) t T3, "ref4/accepting2.txt" using ($2/CHAR_TIME):($3/CHAR_TIME) t T4, "ref5/accepting2.txt" using ($2/CHAR_TIME):($3/CHAR_TIME) t T5, "ref6/accepting2.txt" using ($2/CHAR_TIME):($3/CHAR_TIME) t T6

set output 'deleeuw2d-accepting2-steps.pdf'
set xlabel 'Step'
unset logscale x
set ylabel 'Time step size'
set logscale y
unset xrange
plot "ref3/accepting2.txt" using ($3/CHAR_TIME) t T3, "ref4/accepting2.txt" using ($3/CHAR_TIME) t T4, "ref5/accepting2.txt" using ($3/CHAR_TIME) t T5, "ref6/accepting2.txt" using ($3/CHAR_TIME) t T6



