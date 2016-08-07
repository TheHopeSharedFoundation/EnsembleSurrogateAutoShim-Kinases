chomp($InpPth = $ARGV[0]);
chomp($InpFle = $ARGV[1]);


print "Input path: $InpPth  Input file: $InpFle\n";

chdir($InpPth);

if($InpFle =~ m/^(.*)_(\d+)_out\.sdf$/)
  {
  $EnsMmb = $1;
  $CmpNum = $2;
  

print "Working on file: $InpFle   Ensemble member: $EnsMmb  Compound number: $CmpNum\n";

$InpNme = $InpPth . "/" . $InpFle;

open(Input,"$InpNme");
   @InpAry = <Input>;
close(Input);

$ArySze = @InpAry;

$OutNme = $InpPth . "/" . $EnsMmb . "_PreHinge_" . $CmpNum . ".sdf";

open(Output,">$OutNme");

for($AryCnt = 0; $AryCnt < $ArySze; $AryCnt++)
    {
    $InpAry[$AryCnt] =~ s/> <NAM>/> <NAME>/;
    $InpAry[$AryCnt] =~ s/> <ScoreAfterMin>/> <Score>/;
    $InpAry[$AryCnt] =~ s/> <Steric>/> <VwInt>/;
    if($InpAry[$AryCnt] =~ /^>\s+<PCN>$/)
      {
      $PseCnt++;
      $InpAry[$AryCnt] =~ s/> <PCN>/> <IX>/;
      $InpAry[$AryCnt + 1] = "$PseCnt\n"    
      }    
    if($InpAry[$AryCnt] =~ /^>\s+<PC>$/)
      {
      $InpAry[$AryCnt] =~ s/> <PC>/> <projPath>/;
      $InpAry[$AryCnt + 1] = "/NewAutoShim/$EnsMmb/\n"
      }
    if($InpAry[$AryCnt] =~ /^(\d{6}_\d+)$/)
      {
      $OldChiIdn = $1;
      if($InpAry[$AryCnt] =~ /^(\d{6})_\d+$/)
        {
        $NewChiIdn = $1;
        $InpAry[$AryCnt] =~ s/$OldChiIdn/$NewChiIdn/;
        if($InpAry[$AryCnt + 1] =~ m/^\s+MOLSOFT.*/)
          {
          $InpAry[$AryCnt] =~ s/$NewChiIdn/$NewChiIdn PSE:null SCR:null PRO:null/;
          }
        }  
      }
   print Output "$InpAry[$AryCnt]";
    }
close(Output);
  }

