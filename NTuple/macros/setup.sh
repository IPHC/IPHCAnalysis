#!/bin/bash

echo "#################################" 
cmsenv
dir=`pwd`
export NTUPLEANA_PATH=$dir
export NTUPLEANA=$dir
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib:$dir/../../../IPHCDataFormat/NTFormat/src/.:$dir/../.lib/
export NTUPLEDATAFORMAT_PATH=$dir/../../../IPHCDataFormat/NTFormat/

# PoD
export SBGPODPATH=$dir/../../../PoD/
export PATH=$PATH:$SBGPODPATH

echo NTUPLEANA_PATH=$NTUPLEANA_PATH
echo NTUPLEANA=$NTUPLEANA_PATH

LHAPDF_PATH=$NTUPLEANA_PATH/../LHAPDF/
JETCORRECTIONS_PATH=$NTUPLEANA_PATH/../JetCorrections/

export LIBRARY_PATH=$LHAPDF_PATH
export LD_LIBRARY_PATH=$LHAPDF_PATH/lib/:$JETCORRECTIONS_PATH:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$LHAPDF_PATH/include:/opt/exp_soft/cms/slc5_amd64_gcc462/external/boost/1.50.0-cms/include/:$CPLUS_INCLUDE_PATH
export LHAPATH=$LHAPDF_PATH/share/lhapdf/

if [ -d ../.lib ] ; then 
	echo "List of librairies"
	ls ../.lib
else 
	mkdir ../.lib
fi

echo " Setup is DONE"
echo "#################################" 


export CMSSW_BASE=/opt/sbg/data/safe1/cms/jandrea/UpdateFramework_2013_08_22/CMSSW_5_3_11/
