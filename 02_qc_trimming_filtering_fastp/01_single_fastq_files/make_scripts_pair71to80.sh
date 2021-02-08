UNIT=PAIR
SCRIPT=fastp #exclude ".sh"
UNITS_PER_SCRIPT=2
#N_UNITS=80
FROM=71
TO=80

#for i in `seq 1 $UNITS_PER_SCRIPT $N_UNITS` 
for i in `seq $FROM $UNITS_PER_SCRIPT $TO`
do
	#echo $i
	#echo "from $i"
	#echo "to $(( i + $(( $UNITS_PER_SCRIPT -1 )) ))"
	#printf "\n"
	from=$i
	to=$(( i + $(( $UNITS_PER_SCRIPT -1 )) ))
	cp ${SCRIPT}.sh $SCRIPT_$UNIT${from}_to_$UNIT${to}_$SCRIPT.sh
	nohup $SCRIPT_$UNIT${from}_to_$UNIT${to}_$SCRIPT.sh $from $to &> $SCRIPT_$UNIT${from}_to_$UNIT${to}_$SCRIPT.out &
done
