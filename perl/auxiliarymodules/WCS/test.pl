use blib;
use PDL;
use Test;
        
BEGIN{
  eval " use WCS; ";
  unless ($@){
    plan tests => 2;
  }
  else {
    plan tests => 1;
    print "ok 1 # Skipped: WCS not installed\n";
    exit;
  }
}


$ra =  167.481458333333;
$dec = 88.94575;

$filename = 'nvss_polaris.fits';

#w$filename = 'iris12_polaris_J2000.fits';

#print "FILE PASS \n";
#($x,$y) = &WCS::wcs2pix($filename,$ra,$dec);
#print "BACK TO PERL $ra $dec ----> $x $y \n";
#
#ok($x > 0);
#

my @h = rfits($filename,{Data => 0});

#$h = rfits($filename,{Data => 0});

print "HEADER PASS \n";
($x,$y) = &WCS::wcs2pix(\@h,$ra,$dec);
#($x,$y) = &WCS::wcs2pix($h,$ra,$dec);
print "BACK TO PERL $ra $dec ----> $x $y \n";


($ra2,$dec2) = &WCS::pix2wcs(\@h,$x,$y);
print "cross check $ra2 $dec2 \n";


ok($x > 0);





