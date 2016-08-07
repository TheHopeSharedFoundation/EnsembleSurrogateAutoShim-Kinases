


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
  $EndFnd = 0;
      while(<InpFle>)
           { #open 2
           $InpSdfStr = $_;
           
           if($InpSdfStr =~ m/^M\s+.*/)
			  { 
            #  $EndFnd = 1;
			 if($InpSdfStr =~ m/^M\s+END/) 
			    {
				push(@PseAry,$InpSdfStr);	
				}  	 
			  }  	  
		   else	  
           # if($EndFnd == 0)
			  {  
              push(@PseAry,$InpSdfStr);
		      }
             
			 #EXTRACT AMD REFORMAT SCORE
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
            
			  #EXTRACT VDW SCORE
              if($InpSdfStr =~ m/>\s+<VwInt>/)
                {
				$ELS_flgVdwFnd = 1;	
			    }
               if($ELS_flgVdwFnd == 1 && $InpSdfStr =~ m/(\D?\d+\.\d+)/)
                {
				$ELS_VdwScore = $1;
				$ELS_flgVdwFnd = 0;
				}



             #EXTRACT POSE INDEX
             if($InpSdfStr =~ m/>\s+<IX>/)
               { #open 5
               $ELS_PseCnt++;
               $ELS_flgIdxFnd = 1;
               } #close 5
             if($ELS_flgIdxFnd == 1 && $InpSdfStr =~ m/^(\d+)\s?$/)
               { #open 6
               chomp($ELS_IdxFndVal = $1);
               $ELS_flgIdxFnd = 0;
               } #close 6

               #EXTRACT PROTEIN ENSEMBLE NAME
             if($InpSdfStr =~ m/>\s+<projPath>/)
               { #open 7
               $ELS_flgProFnd = 1;
               } #close 7   
             if($ELS_flgProFnd == 1 && $InpSdfStr =~ m/^.*\/(.*)\/\s?$/)
               { #open 8
               chomp($ELS_ProEnsNme = $1);
               $ELS_flgProFnd = 0;
               } #close 8

               #EXTRACT PROTEIN ENSEMBLE NAME
		      if($InpSdfStr =~ m/>\s+<NAME>/)
			    { #open 7
			     $ELS_flgNmeFnd = 1;
			     } #close 7
			  if($ELS_flgNmeFnd == 1 && $InpSdfStr =~ m/^(\d+)$/)
				 { #open 8
				 chomp($ELS_Chironid = $1);
				 $ELS_flgNmeFnd = 0;
		         } #close 8



 
           if($InpSdfStr =~ m/^\$\$\$\$/)
             { #open 9
			 $EndFnd = 0;
#			 $Space = "";
#			 $Ending = "\$\$\$\$";
#             push(@PseAry,$Space);
 #            push(@PseAry,$Ending);
			 $PseArySze = @PseAry;
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
                  $PseAry[$PseCnt] =~ s/PSE:null/$ELS_IdxFndVal/;
                  $PseAry[$PseCnt] =~ s/SCR:null/$ELS_NewTotScr/;
                  $ELS_NewPseHdr = $PseAry[$PseCnt];
                  #$OldChiIdn = $NewChiIdn;
                  } #close 11
				$BtmHdr = $ELS_Chironid . "_" . $ELS_IdxFndVal;   
                $PseAry[$PseCnt] =~ s/\s+MOLSOFT\s+.*/$BtmHdr/;
                $PseAry[$PseCnt] =~ s/<projPath>/<EnsembleMember>/;
				$PseAry[$PseCnt] =~ s/\/NewAutoShim\/.*\//$ELS_ProEnsNme/;
				$BtmHdr = "";
                } #close 10
			if($ELS_VdwScore < 0)
			  {	
			  chomp(@PseAry);	
              $OneLneStr = join("%",@PseAry);
              print OutOne "$OneLneStr\n";
			  }
              @PseAry = ();
             } #close 9
           } #close 2
  close(OutOne);
  close(InpFle); 
  } #close 1
#}


