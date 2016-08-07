$InpSdfFle = $ARGV[0];
$OutSdfFle = $ARGV[1];

open(OutFle,">$OutSdfFle");
open(InpFle,"$InpSdfFle");
    while(<InpFle>)
         {
         chomp($InpStr = $_);
         if($InpStr =~ m/.*<VwInt>%-\d+\.\d+%%.*/)
           {
           print OutFle "$InpStr\n";
           }
         }  
close(InpFle);
close(OutFle);

