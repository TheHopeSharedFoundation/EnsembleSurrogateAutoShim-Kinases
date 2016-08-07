#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -N ICM_AutoShimModelBuilding
#$ -V





#DEFINITIONS
HME=`pwd`
cd $HME
OUT=$HME/TaskArrayOutput.txt
MDL=$SGE_TASK_ID
DEC=$(($MDL-1)) 
MOD=$HME/14_Models_SET02
ACT=$HME/13_Assays
FTR=$HME/12_Features
CTF=$FTR/CatMinOut.txt

#CREATE AN ARRAY OF EXISTING ASSAY TEXT FILES
cd $ACT
array=($(ls *.txt))
CNT=${array[$DEC]}

 
#CREATE MODEL DIRECTORIES
cd $MOD
BSE=`echo $CNT |  awk -F'/' '{print $NF}' | awk -F'.' '{print $1}'`
mkdir $MOD/$BSE
cd $MOD/$BSE

#SPLIT THE ASSAY ACTIVITES 25:75 INTO TEST AND TRAINING SETS, RESPECTIVELY
less $ACT/$CNT | grep CHIRONID > AssaysColumnHeader.txt
less $ACT/$CNT | grep -v CHIRONID | awk '{if(rand()<0.25)print $1,$2 > "TestTail.txt"     ; else print $1,$2 > "TrainTail.txt"}'
/bin/sh /home/emv-compchem/protos/archive/sullida/Programs/WinToUnix.scr TestTail.txt
/bin/sh /home/emv-compchem/protos/archive/sullida/Programs/WinToUnix.scr TrainTail.txt

#CREATE THE TARGET/MODEL SPECIFIC ACTIVITY FILE
cat TrainTail.txt TestTail.txt | sort -k1,1 | uniq > t
#cat TrainSet.txt TestSet.txt | sort -k1,1 | uniq | grep CHIRONID > h
mv t activity.txt
txtfile=activity.txt
rm t 

#SORT THE TRAIN AND TEST SETS
cat TrainTail.txt | awk '{print $1}' | sort -k1,1 > temp
mv temp TrainingCHIRONIDonlyTail.txt
cat TestTail.txt | awk '{print $1}' | sort -k1,1 > temp
mv temp TestCHIRONIDonlyTail.txt
rm temp

#JOIN TARGET/MODEL SPECIFIC ACTIVITY FILE WITH SORTED TEST AND TRAIN FILES
cat $txtfile | join - TrainingCHIRONIDonlyTail.txt > TrainingJoinedActivityNoTop.txt
cat $txtfile | join -j 1 - TestCHIRONIDonlyTail.txt > TestJoinedActivityNoTop.txt

#ADD TOP TO JOINED TEST AND TRAIN FILES
cat AssaysColumnHeader.txt TrainingJoinedActivityNoTop.txt > TrainingJoinedActivity.txt
cat AssaysColumnHeader.txt TestJoinedActivityNoTop.txt > TestJoinedActivity.txt

#JOIN TO CatMinOut.txt FILE
cat TrainingJoinedActivity.txt | join - $CTF | awk '{printf"%sx%s %s %s %s", $3,$5,$1,$2,$4;for(i=6;i<=NF;i++)printf" %s",$i;printf"\n"}' \
| sed 's/ /\t/g' | sed 's/_//g' >> ./temptrain
cat TestJoinedActivity.txt | join -j 1 - $CTF | awk '{printf"%sx%s %s %s %s", $3,$5,$1,$2,$4;for(i=6;i<=NF;i++)printf" %s",$i;printf"\n"}' \
| sed 's/ /\t/g' | sed 's/_//g' >> ./temptest

#CREATE R INPUT FILES
grep CHIRONID ./temptrain | head -1 > ./temphead
# In the two following lines, -k4,4nr was replaced by -k4,4n because ICM minimization scores are typically NEGATIVE.  More NEGATIVE scores are more favorable.
grep -v CHIRONID ./temptrain | sort  -k2,2 -k4,4n | cat ./temphead - > ./Rinput.txt
grep -v CHIRONID ./temptest  | sort  -k2,2 -k4,4n | cat ./temphead - >  ./Rinput_test.txt
gzip ./Rinput.txt ./Rinput_test.txt

#RUN R SCRIPTS
export R_LIBS=/usr/prog/R/2.14.2_intel/lib64/R/library
/usr/prog/R/2.14.2_intel/bin/R --no-save --args < /home/emv-compchem/protos/archive/mukhepr3/autoshim_scripts/skywalker/PlsOnRp_newR.r $txtfile


cd $HME

