#!/bin/sh

###########################
# Variables used for training are put in the file VAR_FILE. Each time one variable will be hidden.
VAR_FILE=TMVA_variables.txt
CODE=trainingBDT_FCNC_tZ.C
OUTPUTFILE=trainingBDT_FCNC_zut.root
###########################

`eval scramv1 ru -sh`
NVAR=`wc -l $VAR_FILE | awk '{print $1}'`
NVAR=15
OUTPUT=`echo $OUTPUTFILE | awk -F. '{print $1}'`
mv $CODE $CODE.save

# LOOP OVER VARIABLES
for ((I=1; I<= $NVAR; I++))
do
 # REMOVE THE VARIABLE
 VARTOREMOVE=`cat $VAR_FILE | sed -n ${I}p`
 echo $VARTOREMOVE
 cat ${CODE}.save | sed  "/AddVariable(\"$VARTOREMOVE\"/d" >> $CODE
 # RUN ROOT
 root -l -b -q $CODE
 mv weights weights_$VARTOREMOVE 
 echo "mv $OUTPUTFILE ${OUTPUT}_${VARTOREMOVE}.root"
 mv $OUTPUTFILE ${OUTPUT}_${VARTOREMOVE}.root
 rm -rf $CODE
done

mv -f $CODE.save $CODE

