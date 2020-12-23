
=head1 NAME

    CelCoord - Celestial coordinates from FITS standard astrometry 

=head1 SYNOPSIS

    use CelCoord;

    $h = rfits($file,{Data => 0}); 

    $t = &CelCoord::wcs($file);  # transforms with Mark Calabretta's WCS routines, using the interfaces written by Doug Mink (libwcs)

    $out = $t -> apply($in); # input: pixel coordinates, $in is a 2-elt pdl or is the ndcoords  ouput for a 2D image

    $out = $t -> invert($in); # input: WCS coords in deg

    $out = &CelCoord::match($file_ref,$file_raw);



=head1 DESCRIPTION


CelCoord returns an object transform to convert pixel coordinates to
celestial coordinates. See PDL::Transform for the usage of $t. 

Pixel coordinates are 1-offset, as in the FITS standard (not the PDL standard). 


=head2 FUNCTIONS


=head3 pos 

input ra and dec in string format ('16:20:06.1')

convenient to convert string RA DEC into pixel coord \n"



=head3 match

resamples a FITS image to the WCS of another reference image. Inputs
are the filename of the reference FITS image, and that of the image to
be resampled. A hash reference to a FITS header can be passed instead
of $file_ref, but this alternative involves writing a dummy fits file
for use with libwcs (removed at the end of match).

     $out = &CelCoord::match($file_ref,$file_raw);

$out will be resampled on the WCS info of $file_ref.  $file_raw cannot be gzipped (libwcs does not support gzip)


=head4 options: 

    Nx_out => 100, Ny_out => 100, size of output image

    CoordSys => ['Gal','FK4'] , acronyms for input/output coordinate systems 
                              , case insensitive

    Algo => 'wcs', if set to 'wcs', use libwcs - in this case
    $file_ref cannot be gzipped. If Algo is not set to 'wcs', then the
    only supported projections are ARC, TAN, or SIN (in the CD
    formalism - no CROTA).



=head3 pdl_cd

Implementation of the latest FITS standard described by Calabretta &
Griesen (astro-ph/0207413), without the use of libwcs.  Only ARC, SIN,
or TAN. Not thoroughly tested, but useful when only a header is
available.

$t = &CelCoord::pdl_cd($h);     


=head3 wcs 

Requires libwcs, which can be found as part of the wcstools package at http://tdc-www.harvard.edu/software/wcstools/

libwcs is the interface written by Doug Mink to access the projections
routines by Mark Calabretta, based on the latest FITS standard.


=head3 EXAMPLES

[lio:~/ROPH/Spitzer_ROPH/IRS/]% grep -r 'CelCoord::match' ~/ROPH/*pl  
/Users/simon/ROPH/image_match_Ext2MASS.pl:    my $out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs', CoordSys => ['GAL','FK5']});
/Users/simon/ROPH/image_match_call.pl:    my $out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs', CoordSys => ['GAL','FK5']});
/Users/simon/ROPH/image_match_call.pl:#    my $out = &CelCoord::match($file_ref,$file_raw,{Algo => '', CoordSys => ['GAL','FK5']});
/Users/simon/ROPH/image_match_call.pl:#    my $out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs'});
/Users/simon/ROPH/proc_rgb_dev.pl:	my $out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs'});
/Users/simon/ROPH/proc_rgb_dev.pl:	my $out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs'});
/Users/simon/ROPH/proc_rgb_dev.pl:	my $out = &CelCoord::match($file_ref,$file_raw,{Algo => 'wcs'});

=head1 OUTPUT

STDOUT

=head1 AUTHOR

simon

=cut
package CelCoord;
use PDL;
use PDL::Transform;
use PDL::Complex;
use PDL::NiceSlice;
use PDL::Slatec;
use Astro::Time;
use Astro::Coord;
use Data::Dumper;
use WCS;
use strict;
use warnings;


my $Pi= 2*acos(0) ;
my $pi = $Pi->sclr;

my $dodumfits = 1;

sub match {



#REFERENCE IMAGE
    my ($file_ref, $file_raw, $opts) = @_;


    my $h;
    if (ref($file_ref) =~ /HASH/ ) {
	$h = $file_ref;
    } else {
	$h = rfits($file_ref,{Data => 0});
    }

    my %args = (Algo => '', Nx_out => $$h{NAXIS1}, Ny_out => $$h{NAXIS1}, CoordSys => '');

    if (defined($opts)) {
	%args = (%args,  %$opts);
    }

    
    my $out = zeroes($args{Nx_out},$args{Ny_out});

    my $tf_out; 

    my $h_out =  PDL::_hdr_copy($h) ; # in case we need to edit the reference header

    my $doclean = 0;
    if ( $args{Algo} =~ /wcs/i ) {
	if (ref($file_ref) =~ /HASH/ ) {
	    my $im = zeroes($$h{NAXIS1},$$h{NAXIS2});
	    $im->sethdr($h);
	    print "writing a reference fits file for use with libwcs \n";
	    unlink 'dum.fits' if -e 'dum.fits';

	    wfits $im,'dum.fits';
	    $doclean = 1;
	    $file_ref = 'dum.fits';
	} 

	print "using file_ref $file_ref \n";
	$tf_out = &CelCoord::wcs($file_ref); # libwcs 
    } else {
	$tf_out = &CelCoord::pdl_cd($h_out); # official FITS standard  
    
    }


    
    my $domap = 1;
    my $dogzip = 0;
    if ($domap) {
	if ($file_raw =~ /\.gz/ ) {
	    $dogzip = 1;
	    print "GUNZIP enforced \n";
	    system("gunzip $file_raw \n");
	    $file_raw  =~ s/\.gz//; 
	}
	my $im1 = rfits($file_raw);
	my $h = $im1->gethdr; 
	my $tf = &CelCoord::wcs($file_raw); 

	my $in = pdl($$h{CRPIX1}, $$h{CRPIX2}) ;

	if ($args{CoordSys} ne '' ) {

	    my $tsys = &sys($args{CoordSys});

	    $tf = $tsys -> compose($tf); 
	}
	
	my $t = $tf -> inverse   -> compose($tf_out);

	print "MAPPING \n";
	$out = &mapit($im1, $out, $t->inverse, {Boundary => 'extend'});
	$out = $out -> setbadtoval(0);
	$out = float($out);
	$out -> sethdr($h_out);
	if ($dogzip) {
	    print "GZIP enforced \n";
	    system("gzip $file_raw\n");
	}
    } 

    if ($doclean) {
	unlink 'dum.fits';
    }
    return $out;
}

sub sys {
# coordinate system transforms - based on SLA, as implemented in Astro::Coord

    use Astro::Coord;
    use Astro::Time;
    
    my $acronym = shift;
    
    	
    my $t = t_code(\&sys_forward,\&sys_reverse,{p=> $acronym, name=>'Proj', idim=>2, odim=>2});


    return $t;

}

sub sys_forward {
    return &proj_sys(1,@_);
} 

sub sys_reverse {
    return &proj_sys(-1,@_);
} 


sub pos {
    # convenient to convert string RA DEC into pixel coord \n"
    my ($t, $ra_str,$dec_str) = @_;
    my $ra = str2deg($ra_str,'H');
    my $dec = str2deg($dec_str,'D');
    my $in = pdl[$ra,$dec];
    my $pix = $t->invert($in);
    return $pix;
}

sub proj_sys {

    my $dir = shift;
    my $in = shift;

    my $acronym = shift;

    my $out = $in -> copy;

    my $sysname;

    if ($dir == 1) {
	
	$sysname = "$acronym->[0]$acronym->[1]";
	
    } elsif ($dir == -1) {
	
	$sysname = "$acronym->[1]$acronym->[0]";

    }



    if ($sysname =~ /GALFK5/i ) {
	my $x  = $out((0),:,:)/360. ;  my $y  = $out((1),:,:)/360.;
	
	my ($raB, $decB) = galfk4($x,$y);
	
	my ($ra, $dec) = fk4fk5($raB,$decB);
	
	$out((0),:,:) .= $ra * 360.;
	$out((1),:,:) .= $dec * 360.;
 	
    } elsif ($sysname =~ /FK5GAL/i ) {
	my $x  = $out((0),:,:)/360. ;  my $y  = $out((1),:,:)/360.;
	
	my ($raB, $decB) = fk5fk4($x,$y);
	
	my ($l, $b) = fk4gal($raB,$decB);
	
	$out((0),:,:) .= $l * 360.;
	$out((1),:,:) .= $b * 360.;
 	
    } elsif ($sysname =~ /FK5FK4/i ) {
	my $x  = $out((0),:,:)/360. ;  my $y  = $out((1),:,:)/360.;
	
	my ($raB, $decB) = fk5fk4($x,$y);
	
	$out((0),:,:) .= $raB * 360.;
	$out((1),:,:) .= $decB * 360.;
 	
    } elsif ($sysname =~ /FK4FK5/i ) {
	my $x  = $out((0),:,:)/360. ;  my $y  = $out((1),:,:)/360.;
	
	my ($raJ, $decJ) = fk4fk5($x,$y);
	
	$out((0),:,:) .= $raJ * 360.;
	$out((1),:,:) .= $decJ * 360.;
 	
    }




    return $out;



}




sub wcs  {
# projection transform
# input: filename, cannot be gziped.
# output: transform object
# convention: apply, (i,j) -> (a, d).


    my $file_ref = shift;

    my $doclean = 0;
    if ( (ref($file_ref) =~ /HASH/ )  ) {
	if ($dodumfits) {
	    my $h = $file_ref;
	    my $im = zeroes($$h{NAXIS1},$$h{NAXIS2});
	    $im -> sethdr($file_ref);
	    unlink 'dum.fits' if -e 'dum.fits';
	    print "writing temp fits file\n";
	    wfits $im,'dum.fits';
	    $doclean = 1;
	} 
	$file_ref = 'dum.fits';
	$dodumfits = 0;
    } 



    my $t = t_code(\&wcs_forward,\&wcs_reverse,{p=> $file_ref, name=>'Proj', idim=>2, odim=>2});

#    if ($doclean) {
#	unlink 'dum.fits';
#    }

    return $t;

}

sub wcs_forward {
    return &proj_wcs(1,@_);
} 

sub wcs_reverse {
    return &proj_wcs(-1,@_);
} 


sub proj_wcs {

    my $dir = shift;
    
    # ij -> ad for dir = 1;

    my $in = shift;
    my $file = shift;
    
    my $out = $in -> copy;



    if ($dir == 1) {
	
	my $x  = $out((0),:,:) ;  my $y  = $out((1),:,:);

	my ($ra, $dec) = &WCS::pix2wcs($file,$x,$y);


	$out((0),:,:) .= $ra;
	$out((1),:,:) .= $dec;
	
    } elsif ($dir == -1) {
	
	my $ra  = $out((0),:,:) ;  my $dec  = $out((1),:,:);

	my ($x, $y) = &WCS::wcs2pix($file,$ra,$dec);


	$out((0),:,:) .= $x;
	$out((1),:,:) .= $y;
 	
    }

    return $out;


}

sub pdl_cd  {
# projection transform
# input: FITS header.
# output: transform object
# convention: apply, (i,j) -> (a, d).
# only zenithal projections, Calabretta & Griesen astro-ph/0207413
# no external lib

    my $h0 = shift;

    my $h = PDL::_hdr_copy($h0);


#   OJO! NOT ALL PROJ MATCH NEW CD MATRIX 

    my $docrot = 1 if (defined($$h{CROTA2}) && ($$h{CROTA2} != 0));

    die "cannot do crota without libwcs\n" if ($docrot);

    my $t = t_code(\&proj_forward,\&proj_reverse,{p=> $h, name=>'Proj', idim=>2, odim=>2});
    
    return $t;

}


sub proj_forward {
    # ij -> ad

    my $in = shift;
    my $h = shift;
    
    my $out = $in -> copy;

#    print "out 0 $out \n";
    my $crval1=$$h{CRVAL1};
    my $crval2=$$h{CRVAL2};
    my $crpix1=$$h{CRPIX1};
    my $crpix2=$$h{CRPIX2};


    my $cd;

    if (defined($$h{'CD1_1'})) {
	my $cd1_1 = $$h{cd1_1}; ###CDELT1
	my $cd1_2 = $$h{cd1_2};
	my $cd2_1 = $$h{cd2_1};
	my $cd2_2 = $$h{cd2_2};
	$cd = pdl [ [$cd1_1, $cd1_2],[$cd2_1, $cd2_2]] ;
    } else {
	my $cdelt1=$$h{CDELT1};
	my $cdelt2=$$h{CDELT2};
	my $crota=$$h{CROTA2};
	$crota=$$h{CROTA1} unless ($crota);  
	$crota=0 unless ($crota); 
	my $rho=deg2rad($crota);
	$cd = pdl [  [$cdelt1* cos($rho) , -1*$cdelt2 * sin($rho)] ,
		     [$cdelt1* sin($rho) ,  $cdelt2 * cos ($rho) ] ] ;
    }


    my $ctype1=$$h{CTYPE1};
    my $ctype2=$$h{CTYPE2};




    my $dp0 = $out((0),:,:) - $crpix1; 
    my $dp1 = $out((1),:,:) - $crpix1; 


    $out((0),:,:) .= $cd((0),(0)) * $dp0 +   $cd((1),(0)) * $dp1;

    $out((1),:,:) .= $cd((1),(0)) * $dp0 +   $cd((1),(1)) * $dp1;

    #native longitude and latitude:
    $out = deg2rad($out);

    my $x  = $out((0),:,:) ;  my $y  = $out((1),:,:);
    my $r = 1 * sqrt($x**2+$y**2);

    my $c = $out -> copy;

    $c((0),:,:) .= -1 * $out((1),:,:);
    $c((1),:,:) .= $out((0),:,:);

    my $phi =  Carg($c);  

    my $theta;
    $ctype1 = 'UNKNOWN' unless ($ctype1);
    if ($ctype1 =~ /TAN/ ) { 
	#GNOMONIC PROJECTION 
	$theta = atan ( 1/$r ) ;
    } elsif  ($ctype1 =~ /ARC/ ) { 
        #ZENITHAL EQUIDISTANT
	$theta = $pi/2 - $r ; 
    } elsif ($ctype1 =~ /SIN/ ) { 
	#SLANT ORTHOGRAPHIC 
	$theta = acos($r) ;
    } else {
	print "NOT A RECOGNISED PROJECTION, assuming TAN\n";
	$theta = atan ( 1/$r ) ;
    }

    #alpha delta :
    $crval2= deg2rad($crval2) ; 


    my $alphap = deg2rad($crval1);
    my $deltap = $crval2;
    
    #define projection center: Zenithal projections:
    my $phip=-$pi;
    $c((0),:,:) .= sin($theta) * cos($deltap) - cos($theta) * sin($deltap)  * cos($phi - $phip) ;
    $c((1),:,:) .= -1 * cos($theta) * sin($phi - $phip);
    my $alpha = $alphap + Carg($c);



    my $delta = asin(sin($theta) * sin($crval2)-cos($theta) * cos($phi) * cos($crval2)) ;
    
    #    my $delta = asin(sin($theta) * sin($crval2)-cos($theta) * cos($phi) * cos($crval2)) ;
    #    my $alpha = deg2rad($crval1) + asin(cos($theta)*sin($phi)/cos($delta)) ;

    

    my $a=rad2deg($alpha);
    my $d=rad2deg($delta);


    $a -> where($a < 0) .= $a -> where($a <0) + 360.;


    $out((0),:,:) .= $a; 
    $out((1),:,:) .= $d; 

    return $out;


}


sub proj_reverse {

    # ad -> ij

#    print " ad -> ij \n";

    my $in = shift;
    my $h = shift;
    
    my $out = $in -> copy;



    my $crval1=$$h{CRVAL1};
    my $crval2=$$h{CRVAL2};
    my $crpix1=$$h{CRPIX1};
    my $crpix2=$$h{CRPIX2};


    my $cd;

    if (defined($$h{'CD1_1'})) {
	my $cd1_1 = $$h{cd1_1}; ###CDELT1
	my $cd1_2 = $$h{cd1_2};
	my $cd2_1 = $$h{cd2_1};
	my $cd2_2 = $$h{cd2_2};
	$cd = pdl [ [$cd1_1, $cd1_2],[$cd2_1, $cd2_2]] ;
    } else {
	my $cdelt1=$$h{CDELT1};
	my $cdelt2=$$h{CDELT2};
	my $crota=$$h{CROTA2};
	$crota=$$h{CROTA1} unless ($crota);  
	$crota=0 unless ($crota);
	my $rho=deg2rad($crota);
	$cd = pdl [  [$cdelt1* cos($rho) , -1*$cdelt2 * sin($rho)] ,
		     [$cdelt1* sin($rho) ,  $cdelt2 * cos ($rho) ] ] ;


    }


    my $ctype1=$$h{CTYPE1};



    my $a = $out((0),:,:);
    my $d = $out((1),:,:);




    my $alpha=deg2rad($a);
    my $delta=deg2rad($d);


    #define projection center: Zenithal projections:
    my $phip=-$pi;

    #native longitude and latitude:
    my $deltap=deg2rad($crval2);
    my $alphap=deg2rad($crval1);

    my $c = $out -> copy;

    $c((0),:,:) .= sin($delta) * cos($deltap) -cos($delta) * sin($deltap) * cos($alpha-$alphap); 
    $c((1),:,:) .= -cos($delta) * sin($alpha-$alphap);

#    my $phi = $phip + Carg(pdl [ sin($delta) * cos($deltap) -cos($delta) * sin($deltap) * cos($alpha-$alphap), -cos($delta) * sin($alpha-$alphap)] ) ;

    
    my $phi = $phip + Carg( $c);


    my $theta =  asin( sin($delta)*sin($deltap) + cos($delta)*cos($deltap) * 
		cos($alpha-$alphap)) ;

    #intermediate world coordinates
    my $r;
    $ctype1 = 'UNKNOWN' unless ($ctype1);
    if ($ctype1 =~ /TAN/ ) { 
	# GNOMONIC PROJECTION 
	$r = 1/ tan($theta) ;
    } elsif  ($ctype1 =~ /ARC/ ) { 
#	print "Zenithal Equidistant : ARC \n"; 
	# Zenithal Equidistant : ARC 
	$r = $pi/2 - $theta ; 
    } elsif ($ctype1 =~ /SIN/ ) { 
	#GNOMONIC PROJECTION 
	$r = cos ($theta) ;
    } else {
#	print "NOT A RECOGNISED PROJECTION, assuming TAN";
	$r = 1/ tan($theta) ;
    }

    my $x = rad2deg($r * sin($phi));
    my $y = rad2deg(-$r * cos($phi));


#    my $rho=deg2rad($crota);

#    $cd = pdl [  [$cdelt1* cos($rho) , -1*$cdelt2 * sin($rho)] ,
#	     [$cdelt1* sin($rho) ,  $cdelt2 * cos ($rho) ] ] ;

    my $det = det $cd; 

#    my $cdinv = (1/$det) * pdl [  [ $cdelt2 * cos ($rho) , $cdelt2 * sin($rho)] ,
#	     [-1 * $cdelt1* sin($rho) , $cdelt1* cos($rho)  ] ] ;


##also OK, same as matinv    my $cdinv = (1/$det) * pdl [  [ $cd2_2 , -$cd1_2 ] ,
#				  [-$cd2_1  , $cd1_1  ] ] ;


    my $cdinv = matinv($cd);


    $out((0),:,:) .= $cdinv((0),(0)) * $x + 
	$cdinv((0),(1)) * $y + $crpix1;

    $out((1),:,:) .= $cdinv((1),(0)) * $x + 
	$cdinv((1),(1)) * $y + $crpix2;


    return $out;


}



sub mapit {
    my $in = shift;
    my $out = shift;
    my $t = shift;
    my $opts = shift;
    my $l = ndcoords($out);
    my $inds = $t->invert(double($l));
    my $a = $in -> interpND($inds,$opts);
    return $a;
}





1;
