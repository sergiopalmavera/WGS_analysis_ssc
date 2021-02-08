UNIT=pair
SCRIPT=BWAmem #exclude ".sh"
UNITS_PER_SCRIPT=5
N_UNITS=80
#FROM=1
#TO=10

for i in `seq 1 $UNITS_PER_SCRIPT $N_UNITS` 
#for i in `seq $FROM $UNITS_PER_SCRIPT $TO`
do
	#echo $i
	#echo "from $i"
	#echo "to $(( i + $(( $UNITS_PER_SCRIPT -1 )) ))"
	#printf "\n"
	from=$i
	to=$(( i + $(( $UNITS_PER_SCRIPT -1 )) ))
	NEW_NAME=$UNIT${from}_to_$UNIT${to}_$SCRIPT
	#echo $NEW_NAME
	#printf "\n"
	#cp $SCRIPT.sh $NEW_NAME.sh 
	sed 's/#FROM/FROM='$from'/;s/#TO/TO='$to'/' $SCRIPT.sh > $NEW_NAME.sh
	chmod +x $NEW_NAME.sh
	#nohup $NEW_NAME.sh &> $NEW_NAME.out &
done
