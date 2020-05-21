#!/bin/sh


function CreateIterXX {
DIRNAME=$1
echo "Changing to: $DIRNAME"
cd $DIRNAME
NUMREFS=$2
while [ $NUMREFS -le $3 ]
do
	if [ -d "ref$NUMREFS" ]; then
		cd "ref$NUMREFS"
		
		# This line finds all lines with %
		eval $(grep "%" lastRun/job.output > resultsTmp.txt)
		
		# Remove '%'
		eval $(tr -d '%' < resultsTmp.txt | awk '{$1=$1};1' > resultsTmp2.txt)
		
		# We can now split this
		eval $(csplit -k -f iter resultsTmp2.txt '/Iteration/' '{99}')
		for f in iter??; do mv $f "$f.txt"; 
		done;
		
		cd ..
	else 
		echo "ERROR: directory ref$NUMREFS not found!"
	fi
	((NUMREFS++))
done

cd ..
}


function CreateTable {
DIRNAME=$1
echo "Changing to: $DIRNAME"
cd $DIRNAME


# ls;
NUMREFS=$2
while [ $NUMREFS -le $3 ]
do
	cd "ref$NUMREFS"
	echo "ref$NUMREFS"
	for f in iter??.txt ; do 
	
		#echo $f
		MYVAR1="$(grep DTFACTOR $f)"
		MYVAR1=${MYVAR1//DTFACTOR= /}
		
		MYVAR1=${MYVAR1// index= /:}
		#echo $MYVAR1
		
		# Parse convergence history
		MYVAR2="$(grep -E '^[0-9]+:\W' $f)"
		#echo $MYVAR
		MYVAR2=${MYVAR2//[ ]/:}
		echo "$MYVAR2" > "$f.tmp"
		MYVAR3="$(tail -n 1 "$f.tmp")"
		
	
		# Output variable (for German Excel)
		MYOUT="$MYVAR1:$MYVAR3"
		MYOUT=${MYOUT//./,}
		echo $MYOUT
 	done;
	cd ..;
	((NUMREFS++))
done

cd ..
}




#OPTIONAL:
# --with-transient-mechanics; 	
# --use-rap
# --with-vtk

# creates iterXX files
FROM=4
TO=4

#CreateIterXX "P1P2_ass_gmg_cycleW" $FROM $TO
#CreateIterXX "P1P2_ass_gmg_cycleF" $FROM $TO
CreateIterXX "P1P1stab_ass_gmg_cycleW" $FROM $TO
CreateIterXX "P1P1stab_ass_gmg_cycleF" $FROM $TO
#CreateIterXX "P1P1stab_ass_gmg_cycleV" $FROM $TO


#CreateTable "P1P2_ass_gmg_cycleW" $FROM $TO
#CreateTable "P1P2_ass_gmg_cycleF" $FROM $TO
CreateTable "P1P1stab_ass_gmg_cycleW" $FROM $TO
CreateTable "P1P1stab_ass_gmg_cycleF" $FROM $TO
#CreateTable "P1P1stab_ass_gmg_cycleV" $FROM $TO


