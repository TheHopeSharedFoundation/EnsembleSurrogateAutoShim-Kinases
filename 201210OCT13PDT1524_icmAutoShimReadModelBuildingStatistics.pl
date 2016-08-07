chomp($HmeDir = `pwd`);
$RotDir = $HmeDir . "/14_Models";
#print "RootDir = $RotDir\n";

$MdlAry[0] = "01";
$MdlAry[1] = "02";
$MdlAry[2] = "03";
$MdlAry[3] = "04";
$MdlAry[4] = "05";
$MdlAry[5] = "06";
$MdlAry[6] = "07";
$MdlAry[7] = "08";
$MdlAry[8] = "09";
$MdlAry[9] = "10";

open(Output,">>AllR2_extStatistics.txt");

foreach $Val (@MdlAry)
       {
$MdlDir = $RotDir . "_SET$Val";

#print "ModelAry = $MdlDir\n";

chdir($MdlDir);
chomp(@DirLst = `ls -d */`);
@StaAry = ();

 foreach $DirStr  (@DirLst)
       {
#	   print "DirStr = $DirStr \n";	   
       if($DirStr =~ m/^(.*)\//)
         {
         $DirNme = $1;          
           $DirCnt++;
#           print "$DirCnt  $DirNme\n";
           $AsyDir = $MdlDir . "/" . $DirNme;
           chdir($AsyDir);
           open(Stats,"R2_test.txt");
           (@StaAry = <Stats>);
           $CrtSta = pop(@StaAry);
          close(Stats);
		   if($CrtSta =~ m/^\d+\.\d+\s+.*/)
			 {  
           print "SET$Val $DirCnt  $DirNme $CrtSta";
		   print Output "SET$Val $DirCnt  $DirNme $CrtSta";
             }
		   @StaAry = ();
         }
       }
$MdlDir = "";
$Val = "";
$DirNme = "";
$DirCnt = 0;
	 }

close(Output);	 
