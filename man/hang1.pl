#!/usr/bin/perl

$src = "fres29";
$targ = "fres29-h";

print "Copy $src to $targ with hanging indents\n";

mkdir($targ,0755);

opendir(SRC,$src) || die "no src $src\n";
opendir(TARG,$targ) || die "no target $targ\n";

while ($file = readdir(SRC)) {
 if($file ne "." && $file ne "..") {

#  print " ... read file $file\n";
  
if($file=~/htm/) {
  print " modify html file $file\n";
  
  open(IN,"$src/$file");
  open(OUT,">$targ/$file");

  $pf = 0 ; 
  while(<IN>) {
   $in = $_;
   chop($in);
     if($in =~ /^<P>$/) { 
       $pf=1;
       }
     else {
      if($pf==1) {
      if($in =~ /(\d+ex)/)  {
       $s = $1;
       $in =~ s/$s//;
       print OUT "<p  style='margin:0mm;margin-bottom:0mm;margin-left:$s;text-indent:-$s;'>\n";
#       print OUT "<p  style>\n";
        } 
       else {
         print OUT "<P>\n";
       }
       }
      print OUT "$in\n";
      $pf=0;
      }           
  }
  close(IN); close(OUT);
}
else
{
#  print "cp to $targ/$file\n";
# system("cp -pf $src/$file $targ/$file");
}



}
}
