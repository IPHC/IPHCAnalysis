#!/bin/bash

SIGNAL=zut
WZ_PRIVATE=false
WZ_SYS=true
TZq_SYS

mkdir Output

mkdir Output/withoutJES
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutJES/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doJES=false'

mkdir Output/withoutJER
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutJER/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doJER=false'

mkdir Output/withoutBTag
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutBTag/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doBTag=false'

mkdir Output/withoutPU
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutPU/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doPU=false'

mkdir Output/withoutLept
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutLept/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doLept=false'

mkdir Output/withoutTopMass
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutTopMass/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doTopMass=false'

mkdir Output/withoutPDF
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutPDF/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doPDF=false'

mkdir Output/withoutScale
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutScale/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doScale=false'

mkdir Output/withoutMatch
./prodtemplate 'signal='$SIGNAL 'WZ_PRIVATE='$WZ_PRIVATE 'WZ_SYS='$WZ_SYS 'TZq_SYS='$TZq_SYS  'output=Output/withoutMatch/NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root' 'doMatch=false'
