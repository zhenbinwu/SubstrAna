#!/bin/csh -v

set TAG  = $1
set GEN  = $2
set SelZ = $3
set MET  = $4
set List = $5

set DEL       = DELDIR
set EXE       = ${DEL}/bin/slc6_amd64_gcc472/DELEXE

#============================================================================#
#-----------------------------   Setup the env   ----------------------------#
#============================================================================#
echo "============  Running on" $HOST " ============"
set LOCAL=`dirname $0`
cd $LOCAL
set LOCAL=`pwd`  
source /uscmst1/prod/sw/cms/cshrc prod
cd ${DEL}/src
cmsenv
cd $LOCAL
pwd
ls

#============================================================================#
#--------------------------   To Run the Process   --------------------------#
#============================================================================#
echo $EXE -1 $TAG $GEN $SelZ $MET ${List}
$EXE -1 $TAG $GEN $SelZ $MET ${List}
