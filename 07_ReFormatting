$EnsMmb = $ARGV[0];
@InpAry = <STDIN>;

$ArySze = @InpAry;



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
        if($InpAry[$AryCnt + 1] =~ m/^\s+MOLSOFT.*/ || $InpAry[$AryCnt + 1] =~ m/^\s+SMI2MOL*/)
          {
          $InpAry[$AryCnt] =~ s/$NewChiIdn/$NewChiIdn PSE:null SCR:null PRO:null/;
          }
        }  
      }
     print STDOUT "$InpAry[$AryCnt]";
    }

