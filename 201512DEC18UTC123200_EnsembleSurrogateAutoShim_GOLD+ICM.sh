#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -V 
#$ -N AutoShim_UKER
#$ -e /dev/null
#$ -o /dev/null

#************************************************DESCRIPTION*****************************************************************
# 2014JUL08 16:43:00 PDT
# ICM AutoShim Job Array Script 

# To run this script, do the following;

#	1) Create a home directory and cd into this directory.
#	2) Create a directory named Features.
#	3) Copy all necessary ICM AutoShim directories (listed below) into the home directory.
#		a) EnsembleDefinitions
#			i) This directory contains text files which contain the names of the protein members of each ensemble.
#		b) FeatureGeneration
#			i) This directory contains a Python script, and supporting files, which calculates AutoShim features.
#		c) KinseHingeFilter
#			i) This directory contains a Python script to calculate hydrogen bonding to the kinase hinge region.
#		d) ScoreRanking
#			i) This directory contains Perl scripts which order compound poses by ICM docking score.
#		e) VanDerWaalsFilter
#			i) This directory contains a Perl script which eliminates all poses with a postiive van der Waals energy.
#		f) EnsembleStructures
#			i) This directory contains sub-directories for each member of the ensemble.  These sub-directories contain
#			   files which are necessary for ICM docking.
#		g) RedundantPoseFilter
#			i) This directory contains a C program which eliminates redundant poses.
#	4) Create a directory which contains all of the ocmpounds which are to be docked by ICM (called "input sdf file directory" here).
#	   Every compound must have its own sdf file.  Each sdf file must have the following naming convention: InputSdfCompound_*.sdf 
#	   where * is a placeholder for the compound number/index.  Compounds must be numbered sequentially.  The input sdf file has a 
#          special format. An example of how to make such a file can be found in the Pipeline Pilot 8.0 protocol library of samudbe1 with 
#          the protocol name: 201301JAN23PST0845_Make74AssayFullDataSetSdfFile
#	5) Execute this job array 

# Useage: qsub -t 1-maxNumberOfCompounds:1  fullPathToThisJobArrayScript  fullPathToInputSdfFileDirectory  fullPathToOutputFile 

# The input sdf directory specified by the argument "fullPathToInputSdfFileDirectory" contains an sdf file for each compound that is to be docked.  OMIT THE FORWARD SLASH IN THE PATH 
# DEFINITION.  You can name this directory whatever you would like. 

# The ouput file specified by the argument "fullPathToOutputFile" is a space separated text file which contains information such as, for each compound,  how long it took each segment 
# of this job array to execute and how many compounds and poses passed each segment, etc.  You can name this file whatever you would like. 

# Process: 

# A) For each compound contained in the SdfFilePieces directory and for each ensemble and for each ensemble member specified in the directory
#    EnsembleDefinitions, the following is done; 
#	1) ICM docking (Utilizing a call to ICM, referring to an input sdf file directory, and the directories: EnsembleDefinitions and EnsembleStructures)
#		a) generate poses
#		b) dock
#		c) minimize
#		d) score
#	2) Kinase hinge hydrogen bonding filter (Utilizing the directory: KinaseHingeFilter)
#		a) Keep each pose based on whether it forms at least one hydrogen bond to the kinase hinge region
# B) For each ocmpound that successfully docks to at least one ensemble member and has at least one pose which passes the kinase hinge hydrogen
#    bonding filter, the following is done;
#	3) The poses are collected into a list (Utilizing the directory: ScoreRanking).
#	4) The poses in the pose list are ranked.  The poses having the most negative ICM docking scores (favorable) are placed at the top of the list 
#          (Utilizing the directory: ScoreRanking).
#	5) The poses in the pose list are filtered for steric clashes (Utilizing the directory: VanDerWaalsFilter)
#	6) The pose list is filtered to eliminate redundant pose structures (Utilizing the directory: RedundantPoseFilter).
#	7) AutoShim features are generated for each pose (Utilizing the directories: FeatureGeneration and outputting the generated features to the directory: Features)

#	After the generated features files have been deposited into the Features directory, they will have the following naming convention:
#	GeneratedFeaturesCompound_*.txt where * is a placeholder for the compound number/index.  You must manually append all of these files into one file and change
#	the name of this file to: CatMinOut.txt.  You must also ensure that there is only one features header line at the top of the file.   The CatMinOut.txt file must be sorted such that
#       the smallest CHIRONID is at the top of the file and that, for each CHIRONID, the energies are sorted with the most favorable energy towards the top of the file.  The CatMinOut.txt 
#	file can now be used for model building or prediciton.

 
# Benjamin M. Samudio
# Novartis Institutes for Biomedical Research
# Computational Chemistry Group
# Emeryville, California, U.S.A.
# benjamin.samudio@novartis.com
# bmsamudio@ucdavis.edu
# Directed by Eric Martin

# "Toward Alleviating Suffering"
#******************************START OF JOB ARRAY SCRIPT********************************************************************** 
START_Script=$(date +%s)
#****************************************************PATH DEFINITIONS***************************************************
# FULL PATH TO HOME DIRECTORY 
HME=$PWD
OUT=$HME/00_Benchmarking
# PATH TO LIGAND DIRECTORY
LIG=$HME/01_StartLigands
# PATH TO DOCK POSE GENERATION DIRECTORY
DCK=$HME/02_GoldPoseGeneration
# PATH TO KINASE HINGE FILTER STRUCTURE
HNG=$HME/03_KinaseHingeFilter
#PATH TO MINIMIZATION DIRECTORY
MNZ=$HME/04_IcmMinimization
# PATH TO PROTEIN ENSEMBLE DEFINITIONS
DEF=$HME/05_EnsembleDefinitions
# PATH TO PROTEIN ENSEMBLE STRUCTURES
ENS=$HME/06_EnsembleStructures
# PATH TO SDF POST-FORMATTING SCRIPTS
FMT=$HME/07_ReFormatting
# PATH TO SCORE RANKING
SCR=$HME/08_ScoreRanking
# PATH TO Van der Waals (VDW) FILTER
VDW=$HME/09_VanDerWaalsFilter
# PATH TO REDUNDANT POSE FILTER
PSE=$HME/10_RedundantPoseFilter
# PATH TO FEATURE GENERATION
FTR=$HME/11_FeatureGeneration
# PATH TO GENERATED FEATURES DIRECTORY
GEN=$HME/12_Features


CMP=$SGE_TASK_ID


#*************************************************************************************************************************
INI="BlockIndex: $CMP"
THR=1
CFM=150
MAX=5
BNC=21 #Batch size (The maximum batch size is 21)
STR=$((1 + ($BNC * ($CMP - 1))))
FSH=$(($BNC + ($BNC * ($CMP - 1))))
#*****************************************MKER ENSEMBLE MEMBERS********************************************************
for NME in `less $DEF/MKER_8F18.txt | awk '{print $(NF-1)}'`
do
# DOCKING
#  (cat $AUTOSCRIPTS/cmd.cex $ENZDIR/StructurePrep/protein+plptypes.cex \
#  $ENZDIR/StructurePrep/spheres.cex \
#  $LIGANDDIR/ligands+types.cex | dock | gzip > dockit.cex.gz)>& errordockit.txt
#  cat $DEF/DockPoseGenerationCommands.cex $ENS/$NME/protein+plptypes.cex $ENS/$NME/spheres.cex $LIG/ligands+types.cex | dock

#  /home/samudbe1/201303MAR07PST1612_UKERensembleIndividualeTdtToSdfFiles/P38chan/P38chan_10145.sdf

# CONVERT FROM STRING ENCODING 
START_ConvertString=$(date +%s)
less $LIG/*.sdf | obabel -isdf -osdf -f $STR -l $FSH --errorlevel 1 --filter "MW <= 700" -h -p 7 --title --append "InChIKey" | tee $HME/"$NME"_PostString_PreDock_$CMP.sdf  | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 1_CnvStgStr: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 1_CnvStgStr: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_ConvertString=$(date +%s)
DELTA_ConvertString=$(($FINISH_ConvertString - $START_ConvertString))

REPLIG=$HME/"$NME"_PostString_PreDock_$CMP.sdf
REPDIR=$DCK/"$NME"_$CMP
rm -r $REPDIR
mkdir $REPDIR
REPOUT=$HME/"$NME"_PostDock_PreFormat_$CMP.sdf
REPPRO=$DCK/"$NME".mol2


# less /home/samudbe1/201411NOV24UTC08NEG050700_TestGOLDconformationGenerationInICMautoShim/02_GoldPoseGeneration/GoldDockingTemplate.conf

# less /home/samudbe1/201411NOV24UTC08NEG050700_TestGOLDconformationGenerationInICMautoShim/02_GoldPoseGeneration/GoldDockingTemplate.conf | perl -ne 'if($_ =~ m/REPLACE_LIGAND/){$NewStr = $_; $NewStr =~ s/REPLACE_LIGAND/TEST/;print "$NewStr";} elsif($_ =~ m/REPLACE_DIRECTORY/){$NewStr = $_; $NewStr =~ s/REPLACE_DIRECTORY/TEST/;print "$NewStr";}elsif($_ =~ m/REPLACE_OUTPUT/){$NewStr = $_; $NewStr =~ s/REPLACE_OUTPUT/TEST/;print "$NewStr";} elsif($_ =~ m/REPLACE_PROTEIN/){$NewStr = $_; $NewStr =~ s/REPLACE_PROTEIN/TEST/; print "$NewStr";}else{print "$_";}' > /home/samudbe1/201411NOV24UTC08NEG050700_TestGOLDconformationGenerationInICMautoShim/02_GoldPoseGeneration/aurora_0_Gold.conf



#less $DCK/GoldDockingTemplate.conf | perl -ne 'if($_ =~ m/REPLACE_LIGAND/){$NewStr = $_; $NewStr =~ s/REPLACE_LIGAND/'$REPLIG'/;print "$NewStr";}if($_ =~ m/REPLACE_DIRECTORY/){$NewStr = $_; $NewStr =~ s/REPLACE_DIRECTORY/'$REPDIR'/;print "$NewStr";}if($_ =~ m/REPLACE_OUTPUT/){$NewStr = $_; $NewStr =~ s/REPLACE_OUTPUT/'$REPOUT'/;print "$NewStr";}if($_ =~ m/REPLACE_PROTEIN/){$NewStr = $_; $NewStr =~ s/REPLACE_PROTEIN/'$REPPRO'/; print "$NewStr";}else{print "$_";}' >> $REPDIR/"$NME"_"$CMP"_Gold.conf

cat $DCK/"$NME"_head.conf $DCK/All_Body_DiverseSampling.template > $REPDIR/"$NME"_"$CMP"_combo.template 

less $REPDIR/"$NME"_"$CMP"_combo.template | perl -ne 'while(<>){$NewStr = $_; $Var1 = "'$REPLIG'"; $Var2 = "'$REPDIR'"; $Var3 = "'$REPOUT'"; $Var4 = "'$REPPRO'"; $NewStr =~ s/REPLACE_LIGAND/$Var1/;$NewStr =~ s/REPLACE_DIRECTORY/$Var2/;$NewStr =~ s/REPLACE_OUTPUT/$Var3/;$NewStr =~ s/REPLACE_PROTEIN/$Var4/;print "$NewStr";}' > $REPDIR/"$NME"_"$CMP"_Gold.conf


 /usr/prog/ccdc/gold/522/bin/gold_auto $REPDIR/"$NME"_"$CMP"_Gold.conf
#/usr/prog/ccdc/gold/530/bin/gold_auto $REPDIR/"$NME"_"$CMP"_Gold.conf

less $HME/"$NME"_PostDock_PreFormat_$CMP.sdf | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$Inchi = $1; print "$Inchi\n";}else{print "$_";}' > $HME/"$NME"_PostDock_PreHinge_$CMP.sdf



# FIRST INSTANCE OF HINGE BINDING FILTER
START_HingeFirst=$(date +%s)
python $HNG/filterPose.py --receptorFile="$HNG/$NME/$NME.mol2" --ph4FiltFile="$HNG/$NME/hingeph4Const.dat" $HME/"$NME"_PostDock_PreHinge_$CMP.sdf | tee $HME/"$NME"_PostHingeFirst_PreMinimization_$CMP.sdf | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 2_FstHngFlt: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 2_FstHngFlt: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_HingeFirst=$(date +%s)
DELTA_HingeFirst=$(($FINISH_HingeFirst - $START_HingeFirst))
# MINIMIZATION
START_Min=$(date +%s)
less  $HME/"$NME"_PostHingeFirst_PreMinimization_$CMP.sdf | $ICMHOME/icm64 -T $MNZ/201408AUG19PDT053500_novartis.icm $ENS/$NME/MKER_8_F_1_8_$NME.pdb - | tee $HME/"$NME"_PostMinimization_PreHingeSecond_$CMP.sdf | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 3_IcmMinPse: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 3_IcmMinPse: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_Min=$(date +%s)
DELTA_Min=$(($FINISH_Min - $START_Min))
# SECOND INSTANCE OF HINGE BINDING FILTER
START_HingeSecond=$(date +%s)
python $HNG/filterPose.py --receptorFile="$HNG/$NME/$NME.mol2" --ph4FiltFile="$HNG/$NME/hingeph4Const.dat" $HME/"$NME"_PostMinimization_PreHingeSecond_$CMP.sdf | tee $HME/"$NME"_PostHingeSecond_PreReFormat_$CMP.sdf | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 4_SndHngFlt: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 4_SndHngFlt: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_HingeSecond=$(date +%s)
DELTA_HingeSecond=$(($FINISH_HingeSecond - $START_HingeSecond))

# VDW FILTER AND REFORMAT SDF FILE
START_ReFormat=$(date +%s)
less $HME/"$NME"_PostHingeSecond_PreReFormat_$CMP.sdf | babel -isdf -osdf --filter "Steric < 0" --addinindex --append ScoreAfterMin | obabel -isdf -osdf --addtotitle " $NME" --add title --delete projPath | tee  $HME/"$NME"_PostReFormat_PreConsolidate_$CMP.sdf | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; $DupPse = $PseCnt / 2; print "'$CMP' $CmpCnt $OldIdn '$NME' 5_VdwFltRft: $DupPse\n"; $DupPse = 0; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; $DupPse = $PseCnt / 2; print "'$CMP' $CmpCnt $OldIdn '$NME' 5_VdwFltRft: $DupPse\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_ReFormat=$(date +%s)
DELTA_ReFormat=$(($FINISH_ReFormat - $START_ReFormat))


rm $HME/"$NME"_PostHingeFirst_PreMinimization_$CMP.sdf $HME/"$NME"_PostMinimization_PreHingeSecond_$CMP.sdf $HME/"$NME"_PostHingeSecond_PreReFormat_$CMP.sdf $HME/"$NME"_PostDock_PreHinge_$CMP.sdf $HME/"$NME"_PostString_PreDock_$CMP.sdf $HME/"$NME"_PostDock_PreFormat_$CMP.sdf 
rm -r $REPDIR  
done
#************************************************************************************************************************
NME=""
#******************************************UKER-MKER ENSEMBLE MEMBERS*****************************************************
for NME in `less $DEF/UKER.txt | awk '{print $(NF-1)}'`
do

# CONVERT FROM STRING ENCODING
START_ConvertString=$(date +%s)
less $LIG/*.sdf | obabel -isdf -osdf -f $STR -l $FSH --errorlevel 1 --filter "MW <= 700" -h -p 7 --title --append "InChIKey" | tee $HME/"$NME"_PostString_PreDock_$CMP.sdf  | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 1_CnvStgStr: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 1_CnvStgStr: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_ConvertString=$(date +%s)
DELTA_ConvertString=$(($FINISH_ConvertString - $START_ConvertString))

REPLIG=$HME/"$NME"_PostString_PreDock_$CMP.sdf
REPDIR=$DCK/"$NME"_$CMP
rm -r $REPDIR
mkdir $REPDIR
REPOUT=$HME/"$NME"_PostDock_PreFormat_$CMP.sdf
REPPRO=$DCK/"$NME".mol2

cat $DCK/"$NME"_head.conf $DCK/All_Body_DiverseSampling.template > $REPDIR/"$NME"_"$CMP"_combo.template

less $REPDIR/"$NME"_"$CMP"_combo.template | perl -ne 'while(<>){$NewStr = $_; $Var1 = "'$REPLIG'"; $Var2 = "'$REPDIR'"; $Var3 = "'$REPOUT'"; $Var4 = "'$REPPRO'"; $NewStr =~ s/REPLACE_LIGAND/$Var1/;$NewStr =~ s/REPLACE_DIRECTORY/$Var2/;$NewStr =~ s/REPLACE_OUTPUT/$Var3/;$NewStr =~ s/REPLACE_PROTEIN/$Var4/;print "$NewStr";}' > $REPDIR/"$NME"_"$CMP"_Gold.conf

 /usr/prog/ccdc/gold/522/bin/gold_auto $REPDIR/"$NME"_"$CMP"_Gold.conf
#/usr/prog/ccdc/gold/530/bin/gold_auto $REPDIR/"$NME"_"$CMP"_Gold.conf

less $HME/"$NME"_PostDock_PreFormat_$CMP.sdf | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$Inchi = $1; print "$Inchi\n";}else{print "$_";}' > $HME/"$NME"_PostDock_PreHinge_$CMP.sdf

# FIRST INSTANCE OF HINGE BINDING FILTER
START_HingeFirst=$(date +%s)
python $HNG/filterPose.py --receptorFile="$HNG/$NME/$NME.mol2" --ph4FiltFile="$HNG/$NME/hingeph4Const.dat" $HME/"$NME"_PostDock_PreHinge_$CMP.sdf | tee $HME/"$NME"_PostHingeFirst_PreMinimization_$CMP.sdf | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 2_FstHngFlt: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 2_FstHngFlt: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_HingeFirst=$(date +%s)
DELTA_HingeFirst=$(($FINISH_HingeFirst - $START_HingeFirst))
# MINIMIZATION
START_Min=$(date +%s)
less  $HME/"$NME"_PostHingeFirst_PreMinimization_$CMP.sdf | $ICMHOME/icm64 -T $MNZ/201408AUG19PDT053500_novartis.icm $ENS/$NME/UKER.pdb - | tee $HME/"$NME"_PostMinimization_PreHingeSecond_$CMP.sdf | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 3_IcmMinPse: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 3_IcmMinPse: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_Min=$(date +%s)
DELTA_Min=$(($FINISH_Min - $START_Min))
# SECOND INSTANCE OF HINGE BINDING FILTER
START_HingeSecond=$(date +%s)
python $HNG/filterPose.py --receptorFile="$HNG/$NME/$NME.mol2" --ph4FiltFile="$HNG/$NME/hingeph4Const.dat" $HME/"$NME"_PostMinimization_PreHingeSecond_$CMP.sdf  | tee $HME/"$NME"_PostHingeSecond_PreReFormat_$CMP.sdf |  perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 4_SndHngFlt: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; print "'$CMP' $CmpCnt $OldIdn '$NME' 4_SndHngFlt: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_HingeSecond=$(date +%s)
DELTA_HingeSecond=$(($FINISH_HingeSecond - $START_HingeSecond))
# VDW FILTER AND REFORMAT SDF FILE
START_ReFormat=$(date +%s)
less $HME/"$NME"_PostHingeSecond_PreReFormat_$CMP.sdf | obabel -isdf -osdf --filter "Steric < 0" --addinindex --append ScoreAfterMin | obabel -isdf -osdf --addtotitle " $NME" --add title --delete projPath | tee $HME/"$NME"_PostReFormat_PreConsolidate_$CMP.sdf | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $CmpCnt++; $DupPse = $PseCnt / 2; print "'$CMP' $CmpCnt $OldIdn '$NME' 5_VdwFltRft: $DupPse\n"; $DupPse = 0; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $CmpCnt++; $DupPse = $PseCnt / 2; print "'$CMP' $CmpCnt $OldIdn '$NME' 5_VdwFltRft: $DupPse\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_ReFormat=$(date +%s)
DELTA_ReFormat=$(($FINISH_ReFormat - $START_ReFormat))
rm $HME/"$NME"_PostHingeFirst_PreMinimization_$CMP.sdf $HME/"$NME"_PostMinimization_PreHingeSecond_$CMP.sdf $HME/"$NME"_PostHingeSecond_PreReFormat_$CMP.sdf $HME/"$NME"_PostDock_PreHinge_$CMP.sdf $HME/"$NME"_PostString_PreDock_$CMP.sdf $HME/"$NME"_PostDock_PreFormat_$CMP.sdf
rm -r $REPDIR
done
#**************************************************CONSOLIDATE FILES********************************************************
START_Consolidate=$(date +%s)
cat $HME/*_PostReFormat_PreConsolidate_$CMP.sdf > $HME/All_PostConsolidate_PreScoreRanking_$CMP.sdf
FINISH_Consolidate=$(date +%s)
DELTA_Consolidate=$(($FINISH_Consolidate - $START_Consolidate))
DOCPSE="$DOCPSE ConsolidateTime: $DELTA_Consolidate"
rm $HME/*_PostReFormat_PreConsolidate_$CMP.sdf
#******************************************FLATTEN AND SORT BY SCORE********************************************************
START_Sort=$(date +%s)
less $HME/All_PostConsolidate_PreScoreRanking_$CMP.sdf | perl -ne 'push(@Ary,$_);if($_ =~ m/\$\$\$\$/){chomp(@Ary);$String=join("%",@Ary);print "$String\n";@Ary=();}' > $HME/PreSortFile_$CMP.txt
less $HME/PreSortFile_$CMP.txt | perl -ne 'if($_ !~ m/^%?\$\$\$\$/){print "$_"}' > $HME/PrePim_$CMP.txt 
cat $HME/PrePim_$CMP.txt $HME/PostSortFile_"$CMP"_*.txt | sort -k1,1 -k3n,3 > $HME/AllSort_$CMP.txt
perl $SCR/201410OCT16UTC07NEG113800_SeparateIndividualCompounds.pl $CMP $HME $HME/AllSort_$CMP.txt
FINISH_Sort=$(date +%s)
DELTA_Sort=$(($FINISH_Sort - $START_Sort))
DOCPSE="$DOCPSE SortTime: $DELTA_Sort"
rm $HME/All_PostConsolidate_PreScoreRanking_$CMP.sdf $HME/PreSortFile_$CMP.txt $HME/AllSort_$CMP.txt $HME/PrePim_$CMP.txt
#**************************************************************************************************************************

IND=0

for TMP in `ls $HME/PostSortFile_"$CMP"_*.txt`
   do
IND=$((IND+1))

#*****************************************************EXPAND***************************************************************
START_Expand=$(date +%s)
less $HME/PostSortFile_"$CMP"_$IND.txt | perl -ne '@Ary=split("%",$_);foreach $val (@Ary){print "$val\n";}' >  $HME/ScoreRankedCompound_"$CMP"_$IND.sdf
FINISH_Expand=$(date +%s)
DELTA_Expand=$(($FINISH_Expand - $START_Expand))
#rm $HME/PostSortFile_"$CMP"_$IND.txt  This line was commented out on 201502FEB15UTC08NEG132300 by SamudBe1 for GOLD_ICM Surrogate AutoShim improvement testing
#****************************************************FILTER REDUNDANT POSES*************************************************
START_Redundant=$(date +%s)
$PSE/201503MAR23UTC07NEG112700_CorrectionToTitlePropertyPoseNumberBug $HME/ScoreRankedCompound_"$CMP"_$IND.sdf | tee $HME/ConformationFilteredCompound_"$CMP"_$IND.sdf | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; $DupPse = $PseCnt / 2; print "'$CMP' '$IND' $OldIdn 6_RedPseFlt: $DupPse\n"; $DupPse = 0; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; $DupPse = $PseCnt / 2; print "'$CMP' '$IND' $OldIdn 6_RedPseFlt: $DupPse\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_Redundant=$(date +%s)
DELTA_Redundant=$(($FINISH_Redundant - $START_Redundant))
rm $HME/ScoreRankedCompound_"$CMP"_$IND.sdf
#****************************************************GENERATE AUTOSHIM FEATURES*********************************************
START_Features=$(date +%s)
python $FTR/genDescriptor.py --ph4desc="$FTR/autoshim_ph4_pts.dat" --ph4DescOutFile="$GEN/GeneratedFeaturesCompound_"$CMP"_$IND.txt" --ph4FiltFile="$FTR/AutoShimBaseFeatures.fdef" $HME/ConformationFilteredCompound_"$CMP"_$IND.sdf
less $GEN/GeneratedFeaturesCompound_"$CMP"_$IND.txt | perl -ne  'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){print "$1\n";}' | sort | perl -ne 'if($_ =~ m/(\w{14}-\w{10}-\w{1}).*/){$NewIdn = $1; $LneCnt++;  if($OldIdn eq $NewIdn){$PseCnt++;} if($LneCnt == 1){}elsif($OldIdn ne $NewIdn){$PseCnt++; print "'$CMP' '$IND' $OldIdn 7_PseFtrGen: $PseCnt\n"; $PseCnt = 0;} $OldIdn = $NewIdn; if(eof){$PseCnt++; print "'$CMP' '$IND' $OldIdn 7_PseFtrGen: $PseCnt\n";}}' >> $OUT/Benchmarking_$CMP.txt
FINISH_Features=$(date +%s)
DELTA_Features=$(($FINISH_Features - $START_Features))

cat $GEN/GeneratedFeaturesCompound_"$CMP"_$IND.txt >> $GEN/GeneratedFeaturesCompound_$CMP.txt

rm $GEN/GeneratedFeaturesCompound_"$CMP"_$IND.txt
rm $HME/ConformationFilteredCompound_"$CMP"_$IND.sdf

done

#****************************************************************************************************************************


FINISH_Script=$(date +%s)
DIFF_Script=$(($FINISH_Script - $START_Script))
DIFF_Script=$(($FINISH_Script - $START_Script))


PSTPRC="StartScript: $START_Script FinishScript: $FINISH_Script TotalScriptTime: $DIFF_Script"

echo "$PSTPRC" >> $OUT/Benchmarking_$CMP.txt

# Key to terms found in the output benchmarking file;

# Num = Number
# Pse = Pose
# Rdt = Redundant
# Ftr = Feature
# Gen = Generation
# Pas = Pass
# Cpd = Compound
# VDW = van der Waals
# Doc = Dock
# Hng = Hinge
# Tme = Time
# All times are in units of seconds



