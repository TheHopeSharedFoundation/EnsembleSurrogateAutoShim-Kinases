$CmpNum = $ARGV[0];
$HmeDir = $ARGV[1];
$InpFle = $ARGV[2];
$NewChi = "";
$OldChi = "";
@PseAry = ();
$PseStr = "";
$LneCnt = 0;


open(InpFle,"$InpFle");
while(<InpFle>)
     {
	 $PseStr = $_;
	 if($PseStr =~ m/^(\w{14}-\w{10}-\w{1})\s+\d+\s+\D?\d+\.\d+\s+.*/)
		{
		$NewChi = $1;
		$LneCnt++;
		 if($LneCnt == 1)
		   {
 #          push(@PseAry,$PseStr);
		   $OldChi = $NewChi;
		   }
         if($NewChi eq $OldChi)
           {
		   push(@PseAry,$PseStr);	   
		   }
         if($NewChi ne $OldChi) 
           {	   
		   $PseNum++;	   
		   #$HME/PostSortFile_"$CMP"_$IND.txt	   
		   $PntOut = "PostSortFile_" . "$CmpNum" . "_$PseNum.txt"; 
		   open(OutFle,">>$HmeDir/$PntOut");
		   foreach $PseVal (@PseAry)	
		         {
				 print OutFle "$PseVal";	 
				 }
		   close(OutFle);
		   @PseAry = ();
		    push(@PseAry,$PseStr);
		    $OldChi = $NewChi;	
		   }
		$PntOut = "";
		$PseVal = "";
	    } 
		$PseStr = "";  
	 }

$PseNum++;
$PntOut = "PostSortFile_" . "$CmpNum" . "_$PseNum.txt";
open(OutFle,">>$HmeDir/$PntOut");
foreach $PseVal (@PseAry)
	   {
	   print OutFle "$PseVal";
	   }
close(OutFle);

	  



close(OutFle);
close(InpFle);
