


$pthInpSdfFle = $ARGV[0];
$pthOutSdfFle = $ARGV[1];

chdir($pthInpSdfFle);

#@PstHngAry = <AllPostHinge_*.sdf>;

#foreach  $PstHngStr  (@PstHngAry)
#{


if($pthInpSdfFle =~ m/.*All_PostConsolidate_PreScoreRanking_(\d+).sdf/)
  { #open 1
  chomp($CmpNum = $1);
  $OutSdfFle = $pthOutSdfFle;
  
  open(InpFle,"$pthInpSdfFle");
  open(OutOne,">$OutSdfFle");
      while(<InpFle>)
           { #open 2
           $InpSdfStr = $_;
           $ELS_PseLneCnt++;
           push(@PseAry,$InpSdfStr);
             #Extract and reformat score
             if($InpSdfStr =~ m/>\s+<Score>/)
               { #open 3
               $ELS_flgScrFnd = 1;
               } #close 3
             if($ELS_flgScrFnd == 1 && $InpSdfStr =~ m/(\D?\d+\.\d+)/)
               { #open 4
               $ELS_NonTrcScr = $1;
               $ELS_NewTotScr = sprintf("%.1f",$ELS_NonTrcScr);
               $ELS_flgScrFnd = 0;
               } #close 4
             if($InpSdfStr =~ m/>\s+<IX>/)
               { #open 5
               $ELS_PseCnt++;
               $ELS_IdxFndLne = $ELS_PseLneCnt;
               $ELS_flgIdxFnd = 1;
               $ELS_IdxFndVal = 0;  
               print "INDEX FOUND\n"; 
               } #close 5
             if($ELS_flgIdxFnd == 1 && $InpSdfStr =~ m/^(\d+)\s?$/)
               { #open 6
               chomp($ELS_IdxFndVal = $1);
               $ELS_flgIdxFnd = 0;
               print "THE INDEX VALUE IS: $ELS_IdxFndVal\n";
               } #close 6
               #Extract protein ensemble name
             if($InpSdfStr =~ m/>\s+<projPath>/)
               { #open 7
               $ELS_flgProFnd = 1;
               $ELS_PrjPrpLne = $ELS_PseLneCnt - 1;
               } #close 7   
             if($ELS_flgProFnd == 1 && $InpSdfStr =~ m/^.*\/(.*)\/\s?$/)
               { #open 8
               chomp($ELS_ProEnsNme = $1);
               print "PROTEIN: $ELS_ProEnsNme\n";
               $ELS_PrjPthLne = $ELS_PseLneCnt - 1;
               $ELS_flgProFnd = 0;
               } #close 8

 
           if($InpSdfStr =~ m/^\$\$\$\$/)
             { #open 9
             $PseArySze = @PseAry;
             print "IDXFNDLNE: $ELS_IdxFndLne\n";
             $Before = $PseAry[$ELS_IdxFndLne];
             $PrjPrpBefore = $PseAry[$ELS_PrjPrpLne];
             $PrjPthBefore = $PseAry[$ELS_PrjPthLne];
             print "BEFORE PRP: $PrjPrpBefore\n";
             print "BEFORE PTH: $PrjPthBefore\n";
             print "BEFORE SUB: $Before\n";
             $PseAry[$ELS_IdxFndLne] =~ s/$ELS_IdxFndVal/$ELS_PseCnt/;
             $PseAry[$ELS_PrjPrpLne] =~ s/projPath/EnsembleMember/;
             $PseAry[$ELS_PrjPthLne] =~ s/.*$ELS_ProEnsNme.*/$ELS_ProEnsNme/;
             $After = $PseAry[$ELS_IdxFndLne];
             $PrjPrpAfter = $PseAry[$ELS_PrjPrpLne];
             $PrjPthAfter = $PseAry[$ELS_PrjPthLne];
             print "AFTER SUB: $After\n";
             print "AFTER PRP: $PrjPrpAfter";
             print "AFTER PTH: $PrjPthAfter";
             $ELS_PseLneCnt = 0;
             $ELS_IdxFndLne = 0;  
             for($PseCnt = 0; $PseCnt < $PseArySze; $PseCnt++)
                { #open 10
                if($PseAry[$PseCnt] =~ m/.*\s+PSE:null\s+SCR:null\s+PRO:null/)
                  { #open 11
                  #chomp($NewChiIdn = $1);
				  #if($NewChiIdn != $OldChiIdn)
					#{ #open 12
					#$ELS_PseCnt = 1;	
				    #} #close 12 
                  $PseAry[$PseCnt] =~ s/PRO:null/$ELS_ProEnsNme/;
                  $PseAry[$PseCnt] =~ s/PSE:null/$ELS_PseCnt/;
                  $PseAry[$PseCnt] =~ s/SCR:null/$ELS_NewTotScr/;
                  $ELS_NewPseHdr = $PseAry[$PseCnt];
                  #$OldChiIdn = $NewChiIdn;
                  } #close 11
                if($PseAry[$PseCnt] =~ s/\s+SMI2MOL\s+.*/$ELS_NewPseHdr/)
                  { #open 13
                  $PseAry[$PseCnt] = $ELS_NewPseHdr;
                  } #close 13
                } #close 10
             chomp(@PseAry);
             $OneLneStr = join("%",@PseAry);
             print OutOne "$OneLneStr\n";
             @PseAry = ();
             } #close 9
           } #close 2
  close(OutOne);
  close(InpFle); 
  } #close 1
#}


