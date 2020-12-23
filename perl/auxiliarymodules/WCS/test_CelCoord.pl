use blib;
use PDL;
use Test;
        
BEGIN{
  eval " use lib './'; use CelCoord; ";
  unless ($@){
    plan tests => 1;
  }
  else {
    plan tests => 1;
    print "ok 1 # Skipped: CelCoord not installed\n";
    exit;
  }
}



#$file_ref = 'iris12_polaris_J2000.fits';
#$file_ref = 'nvss_polaris.fits';
$file_ref = 'iris12_polaris_J2000.fits';
#$file_raw = 'iris12_polaris_B1950.fits';
$file_raw = 'iris12_polaris_GAL.fits';

print "matching file_raw to file_ref, involving a change of coordinate
system, FK5 to GAL \n";


#$out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs', CoordSys => ['GAL','FK5']});

$h_ref = rfits($file_ref,{Data=>0});

#$out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs', CoordSys => ['GAL','FK5']});

$out = &CelCoord::match($h_ref,$file_raw,{Algo => 'wcs', CoordSys => ['GAL','FK5']});



ok(ref($out) =~ /PDL/);





