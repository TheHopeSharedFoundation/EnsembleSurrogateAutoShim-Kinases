


$pthInpSdfFle = $ARGV[0];
$pthOutSdfFle = $ARGV[1];

chdir($pthInpSdfFle);


if($pthInpSdfFle =~ m/.*PostSortFile_\d+_\d+.txt/)
  {
  $OutSdfFle = $pthOutSdfFle;
  open(InpFle,"$pthInpSdfFle");
  open(OutOne,">>$OutSdfFle");
      while(<InpFle>)
           {
           $InpSdfStr = $_;
           @PntAry = split("%",$InpSdfStr);   
           foreach $PntAryStr (@PntAry)
                 {
                 print OutOne "$PntAryStr\n";
                 }
           }
  close(OutOne);
  close(InpFle); 
  }



