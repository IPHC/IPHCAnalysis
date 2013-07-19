#!/bin/sh

###########################
# Variables used for training are put in this file :
VAR_FILE=TMVA_variables.txt
CODE=ReaderBDT.C
OUTPUTROOT=TMVApp.root
THEVERTEX="zut"
###########################

mkdir BDTControlPlots

`eval scramv1 ru -sh`
NVAR=`wc -l $VAR_FILE | awk '{print $1}'`
mv $CODE $CODE.save
OUTPUT=`echo $OUTPUTROOT | awk -F. '{print $1}'`

# LOOP OVER VARIABLES
for ((I=1; I<= $NVAR; I++))
do
 # REMOVE THE VARIABLE
 VARTOREMOVE=`cat $VAR_FILE | sed -n ${I}p`
 VARTOREMOVEROOT=`echo $VARTOREMOVE | awk -F_ '{print $2}'`
 echo $VARTOREMOVE
 rm -rf $CODE
# cat ${CODE}.save | sed  -e "/AddVariable(\"$VARTOREMOVE\"/d" -e "s/weights/weights_$VARTOREMOVE/" -e "s/$OUTPUT/${OUTPUT}_${VARTOREMOVE}/" >> $CODE
 cat ${CODE}.save | sed  -e "/AddVariable(\"$VARTOREMOVE\"/d" -e "s/weights/weights_$VARTOREMOVE/" -e "s/$OUTPUT/${OUTPUT}_${VARTOREMOVE}/" >> $CODE
 #cp $CODE ${CODE}_${VARTOREMOVE} # for debugging
 # RUN ROOT
 root -l -b -q $CODE+
 echo "root -l -b -q PlotBDToutput.C\(\"$THEVERTEX\",\"${VARTOREMOVEROOT}_BDTcut\",\"HistoBDToutput/TMVApp_${VARTOREMOVE}_${THEVERTEX}_bdtcutnom.root\"\)"
 root -l -b -q PlotBDToutput.C\(\"$THEVERTEX\",\"${VARTOREMOVEROOT}_BDTcut\",\"HistoBDToutput/TMVApp_${VARTOREMOVE}_${THEVERTEX}_bdtcutnom.root\"\)
done

mv -f $CODE.save $CODE

