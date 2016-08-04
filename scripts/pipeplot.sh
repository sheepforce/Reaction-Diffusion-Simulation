#/bin/sh
#
# This script helps utilizing gnuplot to plot cross-section
# through the pipe. Give it three arguments
# $1 is the axis that one looks along
# $2 ViewNum is the plane you want to look at
# $3 the number of the species to look at
# $4 time frame at which to look at. -> t_${3}.dat

# assing variables from shell to more helpful names
ViewDir=$1
ViewNum=$2
Species=$3
ViewTime=$4

case "$ViewDir" in
  "x") ViewDirCol=1
       OtherCol1=2
       OtherCol2=3
    ;;
  "y") ViewDirCol=2
       OtherCol1=1
       OtherCol2=3
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
echo "set view equal" > pp.gnuplt
echo "set view 0,0" >> pp.gnuplt
echo "unset ztics" >> pp.gnuplt
echo "unset key" >> pp.gnuplt
echo "set palette rgbformulae 33,13,10" >> pp.gnuplt
echo "set dgrid3d 80,80 qnorm 100" >> pp.gnuplt
#echo "set terminal png transparent size 1000,1000" >> pp.gnuplt
echo "set terminal gif size 1000,1000" >> pp.gnuplt
echo "set output \"${ViewDir}${ViewNum}_Subs${Species}_t${ViewTime}.gif\" " >> pp.gnuplt
echo " " >> pp.gnuplt

#echo "splot \"< awk '\$${ViewDirCol}==${ViewNum} { print \$${OtherCol1}, \$${OtherCol2}, \$${SpeciesCol} }' t_${ViewTime}.dat\" using 1:2:3 with points pointtype 18 pointsize 2.5 palette" >> pp.gnuplt
echo "splot \"< awk '\$${ViewDirCol}==${ViewNum} { print \$${OtherCol1}, \$${OtherCol2}, \$${SpeciesCol} }' t_${ViewTime}.dat | sed 's/D/E/g'\" using 1:2:3 with pm3d palette" >> pp.gnuplt

gnuplot pp.gnuplt
