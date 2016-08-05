#/bin/sh
#
# This script helps utilizing gnuplot to plot cross-section
# through the pipe. Give it 4 or 5 arguments
# $1 is the axis that one looks along
# $2 ViewNum is the plane you want to look at
# $3 the number of the species to look at
# $4 time frame at which to look at. -> t_${3}.dat
# $5 plot type, can be pm3d or points
# $6 extra gnuplot commands to be added

# assing variables from shell to more helpful names
ViewDir=$1
ViewNum=$2
Species=$3
ViewTime=$4
PlotType=$5
GnuPlotAdd=$6

case "$ViewDir" in
  "x") ViewDirCol=1
       OtherCol1=3
       OtherCol2=2
    ;;
  "y") ViewDirCol=2
       OtherCol1=3
       OtherCol2=1
    ;;
  "z") ViewDirCol=3
       OtherCol1=1
       OtherCol2=2
    ;;
  *) echo "ERROR: This axis does not exit, we're in 3 dimensions"
    ;;
esac

SpeciesCol=$(echo "$Species + 3" | bc)

# setting the Gnuplot view
echo "# GnuPlot script for plotting concentration values" > pp.gnuplt
if [ "$ViewDir" == "z" ]
  then 
    echo "set view equal" >> pp.gnuplt
fi
case "$PlotType" in
  "pm3d") 	echo "set dgrid3d 80,80 qnorm 100" >> pp.gnuplt
		PlotLook=" "
    ;;
  "points")	PlotLook="pointtype 7 pointsize 2.5"
    ;;
  "lines")	echo "set dgrid3d 80,80 qnorm 100" >> pp.gnuplt
    ;;
  *)		PlotLook="pointtype 7 pointsize 2.5"
esac
echo "set view 0,0" >> pp.gnuplt
echo "unset ztics" >> pp.gnuplt
echo "unset key" >> pp.gnuplt
echo "set palette rgbformulae 33,13,10" >> pp.gnuplt
echo "set terminal gif size 1000,1000" >> pp.gnuplt
echo "set output \"${ViewDir}${ViewNum}_Subs${Species}_t${ViewTime}.gif\" " >> pp.gnuplt
printf "$GnuPlotAdd" >> pp.gnuplt
echo " " >> pp.gnuplt

#echo "splot \"< awk '\$${ViewDirCol}==${ViewNum} { print \$${OtherCol1}, \$${OtherCol2}, \$${SpeciesCol} }' t_${ViewTime}.dat\" using 1:2:3 with points pointtype 18 pointsize 2.5 palette" >> pp.gnuplt
echo "splot \"< awk '\$${ViewDirCol}==${ViewNum} { print \$${OtherCol1}, \$${OtherCol2}, \$${SpeciesCol} }' t_${ViewTime}.dat | sed 's/D/E/g'\" using 1:2:3 with $PlotType $PlotLook palette" >> pp.gnuplt

gnuplot pp.gnuplt
