=head1 

Example script to resample a FITS image to match another, using CelCoord.pm

=cut
use PDL;
use PDL::NiceSlice;
use CelCoord;
use strict;
use PDL::IO::Dumper;
use Vtools;

my $debug = 1;

#END GLOBALS

{

#REFERENCE IMAGE
    my $file_ref = 'iris12_polaris_J2000.fits';


    my $file_raw = 'iris12_polaris_GAL.fits';  # cannot be gzipped if using libwcs

    my $out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs', CoordSys => ['GAL','FK5']});
#    my $out = &CelCoord::match($file_ref,$file_raw,{Algo => '', CoordSys => ['GAL','FK5']});
    




    
#    my $file_raw = './OphA_Extn2MASS_F_Eq.fits';  # cannot be gzipped if using libwcs
#    my $out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs'});
    
    &Vtools::view($out);

    wfits $out,'test_match.fits';

}
