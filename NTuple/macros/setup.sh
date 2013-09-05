#!/bin/bash

echo "#################################" 

dir=`pwd`
export NTUPLEANA_PATH=$dir
export NTUPLEANA=$dir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib:$dir/../../../IPHCDataFormat/NTFormat/src/.:$dir/../.lib/
export NTUPLEDATAFORMAT_PATH=$dir/../../../IPHCDataFormat/NTFormat/


echo NTUPLEANA_PATH=$NTUPLEANA_PATH
echo NTUPLEANA=$NTUPLEANA_PATH

export LIBRARY_PATH=$NTUPLEANA_PATH/../LHAPDF/lib
export LD_LIBRARY_PATH=$NTUPLEANA_PATH/../LHAPDF/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$NTUPLEANA_PATH/LHAPDF/../include:$CPLUS_INCLUDE_PATH
export LHAPATH=$NTUPLEANA_PATH/../LHAPDF/share/lhapdf/

if [ -d ../.lib ] ; then 
	echo "List of librairies"
	ls ../.lib
else 
	mkdir ../.lib
fi

echo " Setup is DONE"
echo "#################################" 
