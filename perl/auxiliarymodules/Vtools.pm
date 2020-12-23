=head1 NAME

    Vtools - library for easy FITS display

=head1 SYNOPSIS


    &Vtools::view($im);


    $im2 = &Vtools::view($im1);


    &Vtools::view($file);


    &Vtools::spec($file);



=head1 DESCRIPTION


Vtools.pm is an application of PDL::Graphics::PGPLOT::Window for
interactive display of fits images or spectra. 



=head1 FUNCTIONS



=head2 view

=head3 options: 

=head4 code 

         eval block, as in      

         &Vtools::view($im,{code=>"\$win->hold; \$win->ellipse($mux,$muy,$lmax,$lmin,$alpha,{color => 'black',filltype => 'outline'})"});

=head4 File 

         add filename if $im is associated to a file - important for a succesful application of CelCoord (i.e. WCS info)

         &Vtools::view($im,{File=>$filename});



=head3 keystrokes

=head4 g  :           lut                     change colour look-up-table 

=head4 a  :   aspectratio                       change image aspect ratio 

=head4 u  :  boxgreyscale change grey scale extrema to max/min within box 

=head4 q  :      killquit                                             die 

=head4 d  :      distance                    distances between two points 

=head4 s  :    statistics                                     basic stats 

=head4 e  :   extractcuts                                        x,y cuts 

=head4 f  :      gaussfit                                             fun 

=head4 2  :             2         quit Vtools, returns grey scale extrema 

=head4 t  :     transpose                                    change x / y 

=head4 o  :      printgif                        print to gif, pgplot.gif 

=head4 p  :       printps                   print to postscript pgplot.ps 

=head4 1  :             1                          quit Vtools, returns 1 

=head4 C-H:         reset                              return to defaults 

=head4 h  :          help       list keystrokes - see also perldoc Vtools 

=head4 O  : varcentrezoom  select a square subimage about cursor location 

=head4 >  :        unzoom                      undo zoom - better use 'r' 

=head4 P  :    astrometry      change to world coords - could be spectral 

=head4 -  :        reduce                           shrink the image size 

=head4 H  :         hours                sexagesimal axis - does it work? 

=head4 +  :       magnify                         increase the image size 

=head4 X  :          quit                       quit Vtools, close window 

=head4 z  :          zoom                                 select subimage 

=head4 Z  :    centrezoom       select square subimage about image CRVALs 

=head4 w  :     writefits           write current (sub)image to view.fits 

=head4 r  :       reverse                                            undo 

=head4 0  :    removeaxis     get rid of those annoying ticks and numbers 

=head4 J  :       justify       justify pixels - useful for spectral data 

=head4 K  :         retim                   quit Vtools, return image pdl 

=head4 k  :        setbad                      set to bad selected region





=head2 spec


=head3 options: 

=head4 ERRORB

       plot error bars

       &Vtools::spec($x,$y, {ERRORB => $sigma}); 


=head4 OVERPLOT

       &Vtools::spec($callbda,$subdata->sumover,{Overplot => [ [$callbda, $refspec] ] });


=head4 Colors

       color for overplots:

       &Vtools::spec($callbda,$subdata->sumover,{Overplot => [ [$callbda, $refspec], [$callbda, $refspec2] ], Color => [ 2 ,  4 ] });

	&Vtools::spec($lbda_reference,$calspec/max($calspec),{Overplot => [ [$lbda_reference, $refspec/max($refspec)] ]});


=head4 Points

       use points rather than a solid line

       &Vtools::spec($xxs,$spectrum,{Points => 1});

=head4 REVXAXIS 

       reverse x-axis (as when plotting RA)

       &Vtools::spec($xxs,$spectrum,{REVXAXIS => 1});


=head3 keystrokes


=head4 q  :   killquit                       die, quit 

=head4 g  :  gaussfit0        no-baseline gaussian fit 

=head4 D  :     points              change line/points 

=head4 1  :          1           quit, return 'select' 

=head4 p  :    printps   print to postscript pgplot.ps 

=head4 o  :   printgif        print to gif, pgplot.gif 

=head4 ^H :      reset               starting defaults 

=head4 f  :     varfit                    various fits 

=head4 h  :       help                 list keystrokes 

=head4 -  :     invert                     take 1/spec 

=head4 x  :     xrange          cursor adjust  x-scale 

=head4 y  :     yrange          cursor adjust  y-scale 

=head4 Y  :     Yrange           STDIN adjust  y-scale 

=head4 X  :       quit        right-click, quit Vtools 

=head4 s  :      stats                     rough stats 

=head4 r  :    reverse                            undo 

=head4 L  :    autolog                 toggle log axis 

=head4 W  :   toascii2              like w, with noise 

=head4 w  :    toascii                   dump spec.dat


=head2 cube

Collapse datacube into an image, interacting with the cube through
'view' and 'spec'.

=head3 keystrokes from 'view'


=head4 v    :         cubespec      extract spec, change depth range retoptions=>0 

#=head4 v    :         scan          collapse contiguous spectral ranges 


=head3 keystrokes from 'spec'

=head4 0    :       retoptions                       quit, return display settings 

=head4 m    :       medover                         like retoptions with medover 



=head1 PENDING 

cube: switch to two windows for spec and view, so as to pass settings from spec to view with keystrokes without having to return from spec and kill $win.

spec: var fits 

wfits - view.fits with P 



=head1 INSTALLATION REQUIREMENTS


apart from Gauss, CelCoord (and libwcs), available from Simon, also
requires Astro::Time and Astro::Coord from CPAN.

=head1 SEE ALSO

PDL::Graphics::PGPLOT::Window

=head1 AUTHOR

simon casassus, simon@das.uchile.cl

=cut
package Vtools;
use PGPLOT;
use PDL;
use PDL::Graphics::PGPLOT::Window;
#use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');;
use PDL::Fit::Gaussian;
use strict;
use Data::Dumper;
use PDL::Fit::LM;
use PDL::Math;
use PDL::NiceSlice;
use PDL::Fit::Polynomial;
use PDL::Graphics::LUT;
use PDL::IO::Dumper;
use Stat;
use Astro::Time;
use Gauss;
#use Gaussf77;
#use Ftools;
use Gauss1D;
use CelCoord;
use Storable qw(dclone);

$PDL::BIGPDL = 1;

my $cube;
my $stack_spec;
my @view_back_spec;

#my $retsigs;

my $pi = 2*acos(0);


my $win_im;
my $win_spec;
my $h_ref;

my $wcs = 0;

my $nz;
my $cbz;
my $kstep = 0;

my $log = 0;

sub cube {
    my @ori=@_;
    $cube = $ori[0];
    $h_ref = $cube->gethdr;

    $cbz = $cube-> copy;
    $nz = $cbz->getdim(2);

    my $ndims = $cube -> getndims;
    if ($ndims == 4) {
#	my $nx = $cube ->getdim(0);
#	my $ny = $cube ->getdim(1);
#	my $nz = $cube ->getdim(2);
	$cube = $cube(:,:,:,(0));
	$$h_ref{NAXIS} = 3;
	$cube-> sethdr($h_ref);
    }
    

    my $opts = deep_copy($ori[1]);
    my @view = ($cube, $h_ref, $opts);

#    my $im_c = $cube->xchg(0,2)->sumover->xchg(0,1);

    my $im_c = $cube->xchg(0,2)->medover->xchg(0,1);

    my $h_im = deep_copy($h_ref);
    $im_c -> sethdr($h_im);

    
    my $dumopts = $opts;
#    %$dumopts =  %$opts;

    my $kstep = 0;
    while(1) {

	$dumopts = &Vtools::view($im_c, $dumopts);

	my %args;
	

	if (ref($dumopts) =~ /HASH/ ) {
	    
	    $args{Medover} = 0;

	    if (defined($dumopts)) {
		%args = (%args, %$dumopts) ;
	    }

	    my $retsigs =\%args;

	    my $k_1 = $retsigs->{k_1};
	    my $k_2 = $retsigs->{k_2};
#	    print Dumper($retsigs);
#	    print "k_1 $k_1 k_2 $k_2 \n";
	    $cbz = $cube(:,:,$k_1:$k_2);
	    $nz = $cbz->getdim(2);

#	    $im_c = $cbz->xchg(0,2)->sumover->xchg(0,1);

	    if ($retsigs->{Medover}) {
		$im_c = $cbz->xchg(0,2)->medover->xchg(0,1);
	    } elsif ($retsigs->{Sumover}) {
		$im_c = $cbz->xchg(0,2)->sumover->xchg(0,1);
	    } 

#            elsif ($retsigs->{Plane} == +1) {
#		$kstep++;
#		$kstep = $nz-1 if ($kstep >= $nz);
#		print "plane $kstep\n";
#		$im_c = $cbz(:,:,($kstep));
#	    } elsif ($retsigs->{Plane} == -1) {
#		$kstep--;
#		$kstep = 0 if ($kstep <= 0);
#		print "plane $kstep\n";
#		$im_c = $cbz(:,:,($kstep));
#	    }
	    

	    $im_c -> sethdr($h_im);
	    

	} elsif ($dumopts == -1) {
	    last;

	}
	

    }

    # modify im_c with retsigs

    
    
}



sub view {
    my @ori=@_;
    my $im_ref=$ori[0]->copy;
    my $h2=$ori[0]->gethdr;
    my $h3;
    if (defined($h2)) {
	$h3 = deep_copy($h2);
	$im_ref -> sethdr($h3);
    }
    my $opts = deep_copy($ori[1]);
    my @view = ($im_ref, $opts, $h3);


#    my @view=@ori;
    my $redraw=1;
    my $dev='/xwin';
    my @view_back;
    my $stack=1; #stack switch  
    my $options_0=$view[1];
    my $code=$options_0->{'code'};
    if (defined($code)) {
#	print "EVAL CODE $code \n";
    }
    my $w = $options_0 -> {'add'};

    my $genwin;

    if (defined($opts->{ViewWin})) {
	$genwin = 0;
	$win_im = $opts->{ViewWin};
    } else {
	$genwin = 1;
    }

    my $hdr_ref=deep_copy($h2); 

    
    my $im_ref= $view[0] -> copy;
    my $im_ref2= $view[0] -> copy;  # only used in `unzoom', no good.


#    my $hdr = dclone($hdr_ref) if (ref($hdr_ref) =~ /HASH/);
    my $hdr = deep_copy($hdr_ref) if (ref($hdr_ref) =~ /HASH/);
    
    my $file=$options_0->{'File'};
    if (!defined($file)) {
	$file = $hdr;
    } 

    my $hours = 0; 
    my $arcsec = 0; 
    my $dowedge = 0;
    my $ctbl = 'ramp';
    my $inv = 0;
    my $ramp = 'ramp';
    my $pi = 2*acos(0);
    my $nojustify=1;
    my $winscale = 1;
    while ($redraw) {
	if ($stack) {
	    my $im2=$view[0]->copy;
	    my $h2=$view[0]->gethdr;
	    my $h3;
	    if (defined($h2)) {
		$h3 = deep_copy($h2);
	    }
	    my $opts = deep_copy($view[1]);
	    push @view_back, [$im2,$opts,$h3]; 
	    if ((scalar @view_back) > 50) {shift @view_back}
	}
	$stack=1;
	my $im= $view[0] ;
	my $icube=$im -> getndims;
	if ($icube > 2) {
	    my $nx = $im ->getdim(0);
	    my $ny = $im ->getdim(1);
	    $im-> reshape($nx,$ny);
	    print "Image is larger than 2D \n";
	}
	$hdr = $im -> gethdr;
	
	my $options=$view[1];
	$options->{Device}=$dev;
	$options->{WindowWidth}=5;
	$options->{Justify}= $nojustify;
	$options->{PIX}=1 if ($hours);
	my $naxis1 = $im -> getdim(0);
	my $naxis2 = $im -> getdim(1);
	my $aspect = 1;
	$options->{ASPECT} = $aspect;
       
	$options->{Size} = [$winscale * 100, $winscale * 100];
        $options->{Unit} = 'mm';

	if ($genwin) {
	    $win_im = PDL::Graphics::PGPLOT::Window->new($options);

       	}


	$genwin=0;
	my $gsmin=$options->{'MIN'};
	my $gsmax=$options->{'MAX'};
	if (defined($options->{Transform})) { 

	    my $tr = &gettransform($hdr, $hours, $win_im, {Arcsec => $arcsec});
	    $options->{Transform} = $tr;
	    my $scale = 1;
	    $scale = 3600 if ($arcsec);
	
	    my $i = 0; my $j =0;
	    my ($a0, $d0) = &ijad_linear($tr,$i,$j);
#	    $a0 *= $scale;
#	    $d0 *= $scale;

	    $i = ($im->getdim(0))-1; $j = ($im->getdim(1))-1;
	    my ($a1, $d1) = &ijad_linear($tr,$i,$j);
#	    $a1 *= $scale;
#	    $d1 *= $scale;
	    
#	    die "a1 $a1 d1 $d1 a0 $a0 \n";
	    print $options->{Xtitle},"<<<< \n";
	    $win_im -> env($a0,$a1,$d0,$d1,$options) ;

	} else {

#	    $win_im -> env(map {$_--} ($im->dims),$options) ;

	}
	if (($dev =~ /ps/) or ($dev eq '/gif')) {
#	    $win-> ctab('Gray'); #bug fix - Jarle Brinchman
#	     if (!(defined $gsmin) & !(defined $gsmax)) { 
#		 $win-> cont($im);
#	     } else {
#		 my $contours = sequence(10)*($gsmax-$gsmin)/10.;
#		 $win-> cont($im, $contours);
#	     }
#	     $win -> close; $dev = '/xwin'; next;
#	 }
#	    print Dumper($options);
	    unlink 'pgplot.ps' if -e 'pgplot.ps';
	    my $win2 = PDL::Graphics::PGPLOT::Window->new($options);
	    $win2 -> ctab( lut_data($ctbl,$inv,$ramp) );
	    if (defined($options->{Transform})) { 
#		my $tr = &gettransform($hdr, $hours, $win);
#		$options->{Transform} = $tr;
		my $tr=$options->{Transform};

		my $i = 0; my $j =0;
		my ($a0, $d0) = &ijad_linear($tr,$i,$j);
		$i = ($im->getdim(0))-1; $j = ($im->getdim(1))-1;
		my ($a1, $d1) = &ijad_linear($tr,$i,$j);
		$win2 -> env($a0,$a1,$d0,$d1,$options) ;
	    }
	    $win2 -> imag($im,$options);
	    $win2 -> close; $dev = '/xserve'; next;
	} else {


	    $win_im -> ctab( lut_data($ctbl,$inv,$ramp) );
	    $win_im -> imag($im,$options);

	}
	if (defined($code)) {eval $code; print "error in code: $@ \n"  if ($@) ; }
	my $xref; my $yref;
	while(1) {

	    (my $x, my $y, my $ch,  $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref,'Type' => 'CrossHair'});
	    $xref=$x;
	    $yref=$y;

#	    print "x $x y $y \n";

#	    my $key = {killquit => "q", magnify => "+", reduce => "-",
#	    quit => 'X', astrometry => "P", removeaxis => "0", reverse
#	    => "r", statistics => "s", boxgreyscale => "u", setbad =>
#	    "k", retim => "K", distance => "d", writefits => "w",
#	    transpose => "t", gaussfit => "f", justify => "J", unzoom => ">" , zoom => "z", centrezoom => "Z"};
##
#

	    my $key = {
		'killquit' => ['q', 'die'],
		'quit' => ['X', 'right-click, quit Vtools, close window'] ,
		'ret' => ['V', 'cube: return'] ,
		'writefits' => ['w', 'write current (sub)image to view.fits'],
		'retim' => ['K', 'quit Vtools, return image pdl'],
		'reset' => ['^H', 'return to defaults'], 
		'help' => ['h', 'list keystrokes - see also perldoc Vtools'],
		'1' => ['1', 'quit Vtools, returns 1'], 
		'2' => ['2', 'quit Vtools, returns grey scale extrema'],

		'reverse' => ['r', 'undo'],

		'setbad' => ['k', 'set to bad selected region'],

		'boxgreyscale' => ['u', 'change grey scale extrema to max/min within box'],

		'distance' => ['d', 'distances between two points'],
		'statistics' => ['s', 'basic stats'],
		'gaussfit' => ['f', 'fun - try it out'],

		'astrometry' => ['P', 'change to world coords - could be spectral'],

		'transpose' => ['t', 'change x / y'],
		'reduce' => ['-', 'shrink the image size'],
		'magnify' => ['+', 'increase the image size'],
		'justify' => ['J', 'justify pixels - useful for spectral data'],
 		'aspectratio' => ['a', 'change image aspect ratio'],
		'extractcuts' => ['e', 'x,y cuts'], 
		'unzoom' => ['>', 'undo zoom - better use \'r\''],
		'zoom' => ['z', 'select subimage'],
		'cubespec' => ['v', "extract spec, change depth range retoptions=>0"],
		'cubepointspec' => ['b', "extract spec, change depth range retoptions=>0"],
		'cubespecmed' => ['m', "extract spec, medover, change depth range retoptions=>0"],
#		'scan' => ['b', 'collapse contiguous spectral ranges'],
		'centrezoom' => ['Z', 'select square subimage about image CRVALs'],
		'varcentrezoom' => ['O', 'select a square subimage about cursor location'],
		'removeaxis' => ['0', 'get rid of those annoying ticks and numbers'],
		'printps' => ['p', 'print to postscript pgplot.ps'],
		'printgif' => ['o', 'print to gif, pgplot.gif'],
 
		'lut' => ['g', 'change colour look-up-table'],
		'hours' => ['H', 'sexagesimal axis - does it work?'],
		'arcsec' => ['"', 'arsec labels'],
		'dowedge' => ['T', 'draw/undraw wedge'],
		'addcont' => ['c', 'add contours'],
		'zstepneg'   =>  ['[','cube: decrement plane z-coord'],
		'zsteppos'   =>  [']','cube: increment plane z-coord']
	    };

#	    print Dumper($key);
#	    die;


	    if ($ch eq $key->{killquit}->[0]) {
		$win_im -> close;
		die;
	    }
	    if ($ch eq $key->{magnify}->[0]) {
		$winscale *= 1.2; 
		$win_im -> close;
		$genwin=1;
		@view=($im,$options);
		last;
	    }
	    if ($ch eq $key->{reduce}->[0]) {
		$winscale *= 0.8; 
		$win_im -> close;
		$genwin=1;
		@view=($im,$options);
		last;
	    }

	    if (defined($options->{Transform})) {
		my $tr=$options->{Transform};


		my ($i, $j) = &adij_lin($tr,$x,$y);

		(((0<$i) && ($i<$naxis1-1)) && ((0<$j) && ($j<$naxis2-1))) ? printf "$ch \( %.1f, %.1f\) = %.2e \n",  $i, $j, $im->at($i,$j) : 1;

		my $offsets = $options->{offsets};

		my ($a,$d);

		if (defined($cube)) {
		    
		    ($a,$d) = &towcs($cube->gethdr, $offsets, $i, $j);
		    
		} else {
		    

		    ($a,$d) = &towcs($file, $offsets, $i, $j);
		}


		$i = sprintf("%d",$i);
		$j = sprintf("%d",$j);
		
		print "$ch: ",deg2str($a,'H',4),"  ",deg2str($d,'D',3),"  DEG: ",$a,"  ",$d," \n";  # EXACT WCS RA DEC

		

	    } else {
		(((0<$x) && ($x<$naxis1-1)) && ((0<$y) && ($y<$naxis2-1))) ?  printf "$ch \( %.1f, %.1f\) = %.2e \n",  $x, $y, $im->at($x,$y) : 1;


		$x=0 if ($x<0);
		$x=$naxis1-1 if ($x>$naxis1-1);
		$y=0 if ($y<0);
		$y=$naxis2-1 if ($y>$naxis2-1);
	    }
	    if ($ch eq $key->{quit}->[0]) {		
		$win_im -> close;
		$redraw=0; 
		return -1;
	    }
	    if ($ch eq $key->{ret}->[0]) {		
#		$win_im -> close;
#		undef($win_im);
		return 1;
	    }


	    if ($ch eq $key->{astrometry}->[0]) {
		my $dowedge = 1;
		my $hdr=$im->gethdr;
		if (!defined($hdr->{CRVAL1})) {
		    print "NO ASTROMETRY\n";
		    last;
		}

		my $transform = &gettransform($hdr, $hours, $win_im);

		$wcs = 1;

		$options->{Transform} = $transform;
		$options->{Axis} = ["BCINST","BCINST"];
		$options->{Axis} = ["ZYHBCINST","ZYDBCINST"] if ($hours);
		$options->{DrawWedge}=$dowedge;
		$win_im -> close;
		$genwin = 1;
		@view=($im,$options);
		last;
	    }


	    if ($ch eq $key->{dowedge}->[0]) {
		$dowedge = $dowedge?0:1;
		$options->{DrawWedge}=$dowedge;
		$win_im -> close;
		$genwin = 1;
		@view=($im,$options);
		last;
	    }


	    if ($ch eq $key->{arcsec}->[0]) {
		$dowedge = 1;
		my $hdr=$im->gethdr;
		if (!defined($hdr->{CRVAL1})) {
		    print "NO ASTROMETRY\n";
		    last;
		}
		$arcsec= 1;
		my $transform = &gettransform($hdr, $hours, $win_im, {Arcsec=> 1});

		$wcs = 1;

		$options->{Transform} = $transform;
		$options->{Axis} = ["BCINST","BCINST"];
		$options->{Axis} = ["ZYHBCINST","ZYDBCINST"] if ($hours);
		$options->{DrawWedge}=$dowedge;
		$win_im -> close;
		$genwin = 1;
		@view=($im,$options);
		last;
	    }
	    if ($ch eq $key->{removeaxis}->[0]) {
		$options->{Axis} = [" "," "];
		$win_im -> close;
		$genwin = 0;
		@view=($im,$options);
		last;
	    }


	    if ($ch eq $key->{reverse}->[0]) { #reverse
		$win_im -> close;
		if ((scalar @view_back) > 0) {
		    my $temp = pop @view_back;
		    $view[0]=$temp->[0];
		    $view[1]=$temp->[1];
		    my $hdr2=$temp->[2];
		    $view[0] -> sethdr($hdr2);
		}
		$stack=0;
		$genwin = 1;



		last;
	    }
	    if ($ch eq $key->{statistics}->[0]) {
		#statistics
		my $iref; my $jref;
		my $dim1=$im->getdim(0);
		my $dim2=$im->getdim(1);
		my $beam; my $pix;
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};

#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $iref=($xref-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $jref=($yref-$t4-$iref*$t5)/$t6;

		    ($iref, $jref) = &adij_lin($tr,$xref,$yref);
		    
		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		    #Jy/beam fluxes
		    my $h = $im->gethdr();
		    if ($$h{BUNIT}=~/Jy\/beam/i) {
			$beam = $pi*$$h{BMAJ}*$$h{BMIN}/(4*log(2));
			$pix = abs($$h{CDELT1}*$$h{CDELT2});
		    }
		} else {
		    $iref=$xref;
		    $jref=$yref;
		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		}
		($x, $y, my $ch, $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref, 'Type' => 'Rectangle'});
		my $i; my $j;
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};

#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $i=($x-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $j=($y-$t4-$i*$t5)/$t6;

		    ($i, $j) = &adij_lin($tr,$x,$y);

		    $i=0 if ($i<0);
		    $i=$naxis1-1 if ($i>$naxis1-1);
		    $j=0 if ($j<0);
		    $j=$naxis2-1 if ($j>$naxis2-1);
		} else {
		    $x=0 if ($x<0);
		    $x=$naxis1-1 if ($x>$naxis1-1);
		    $y=0 if ($y<0);
		    $y=$naxis2-1 if ($y>$naxis2-1);
		    $i=$x;
		    $j=$y;
		}
		my $x1=int(min(pdl [$i,$iref]));
		my $x2=int(max(pdl [$i,$iref]));
		my $y1=int(min(pdl [$j,$jref]));
		my $y2=int(max(pdl [$j,$jref]));
		my $subim  = $im($x1:$x2,$y1:$y2);
		my $min_region = min ($subim);
		print "Min:  $min_region\n";
		my $max_region = max ($subim);
		print "Max:  $max_region\n";
		my $total_flux=sum($subim);
		my $npix=$subim -> nelem;
		print "Total flux: $total_flux   in  $npix pixels\n"; 
		if (defined($beam)) {
		    my $h = $im->gethdr();
		    print "Which is ",$total_flux*$pix/$beam," Jy\n";
		    printf "for a beam of %.4f x %.4f or  %e  deg2 \n ",$$h{BMAJ}*3600,$$h{BMIN}*3600., $beam;
		    print "and pixels with $pix deg2 solid angle \n";
		}
		$hdr=$im->gethdr;
		my $bunits=$$hdr{BUNIT};
		if (defined($bunits)) {
		    if ($bunits =~ /MJY\/sr/i) {
			print "Units are MJy/sr, so the total flux",sum($subim)," is \n";
			my $dx=$$hdr{CDELT1};
			$dx = defined($dx)?$dx:$$hdr{CD2_2};
			print "dx = $dx \n";
			my $domega=($dx*$pi/180.)**2;
			print (sum($subim)*$domega*1E6);
			print "   Jy  domega $domega \n";
		    }
		}
		my $mean=$total_flux/$npix;
		print "Mean: $mean  ";
		my $vari=sum(($subim -$mean )**2)/float(abs($npix));
		my $rms=sqrt($vari);
		print "rms: $rms ";
		my $median=median($subim);
		print "median: $median \n ";
		print "strike middle mouse button to attempt a noise estimate, or u to specify units for a flux measurement (MJy/sr), any other key  to proceed \n";
		($x, $y, $ch, $xref, $yref) = $win_im -> cursor({'XRef' => $x,'YRef' => $y, 'Type' => 'Crosshair'});

		if ($ch eq "D") { 
		    my @noisearr=(&Stat::noise($subim)); 
		    print "noise: ", shift @noisearr,"\n";
		}
		if ($ch eq "u") {
		    print "If units are MJy/sr,  the total flux is \n";
		    my $dx=$$hdr{cdelt1};
		    my $do=($dx*$pi/180.)**2;
		    print (sum($subim)*$do*1E6);
		    print "with $do sterad per pixel \n";
		};
	    }
	    if ($ch eq $key->{boxgreyscale}->[0]) {
		#new color scale, from the min and max within a box

		my $iref; my $jref;
		my $dim1=$im->getdim(0);
		my $dim2=$im->getdim(1);
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};

		    ($iref, $jref) = &adij_lin($tr,$xref,$yref);

#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $iref=($xref-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $jref=($yref-$t4-$iref*$t5)/$t6;

		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		} else {
		    $iref=$xref;
		    $jref=$yref;
		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		}
		($x, $y, my $ch, $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref, 'Type' => 'Rectangle'});
		my $i; my $j;
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};

		    ($i, $j) = &adij_lin($tr,$x,$y);

#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $i=($x-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $j=($y-$t4-$i*$t5)/$t6;
		    
		    
		    
		    $i=0 if ($i<0);
		    $i=$naxis1-1 if ($i>$naxis1-1);
		    $j=0 if ($j<0);
		    $j=$naxis2-1 if ($j>$naxis2-1);
		} else {
		    $x=0 if ($x<0);
		    $x=$naxis1-1 if ($x>$naxis1-1);
		    $y=0 if ($y<0);
		    $y=$naxis2-1 if ($y>$naxis2-1);
		    $i=$x;
		    $j=$y;
		}
		my $x1=int(min(pdl [$i,$iref]));
		my $x2=int(max(pdl [$i,$iref]));
		my $y1=int(min(pdl [$j,$jref]));
		my $y2=int(max(pdl [$j,$jref]));
		my $min_region = min ($im($x1:$x2,$y1:$y2));
		printf "Min in box is: %.2e \n",  $min_region;
		my $max_region = max ($im($x1:$x2,$y1:$y2)); 
		printf "Max in box is: %.2e \n", $max_region;
		$win_im -> erase;
		$genwin=0;
 		$options->{"MIN"} = $min_region;
 		$options->{"MAX"} = $max_region;
		@view=($im,$options);
		last;
	    }

######################################################################
# IT'S A CUBE!
######################################################################

#
#	    if ($ch eq $key->{scan}->[0]) {  
#		print "enter number of contiguous channels\n";
#		my $retsigs;
#		$retsigs->{scan} = <>;
#		return $retsigs;
#	    }
#

	    
	    if ($ch eq $key->{zsteppos}->[0]) {
		$kstep++;
		$kstep = $nz-1 if ($kstep >= $nz);
		print "plane ",$kstep+1," nz $nz \n";
		$im = $cbz(:,:,($kstep));
		$im -> sethdr($h_ref);
		$genwin = 0;
		$options->{Xtitle} = "plane $kstep";
		@view=($im,$options);
		last;
	    }


	    if ($ch eq $key->{zstepneg}->[0]) {
		$kstep--;
		$kstep = 0 if ($kstep <= 0);
		print "plane ",$kstep+1," nz $nz \n";
		$im = $cbz(:,:,($kstep));
		$im -> sethdr($h_ref);
		$genwin = 0;
		$options->{Xtitle} = "plane $kstep";
		@view=($im,$options);
		last;
	    }


	    if ($ch eq $key->{cubepointspec}->[0]) {  # $key->{boxgreyscale}->[0]) {

		#new color scale, from the min and max within a box
		my $iref; my $jref;
		my $dim1=$im->getdim(0);
		my $dim2=$im->getdim(1);
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};

		    ($iref, $jref) = &adij_lin($tr,$xref,$yref);

#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $iref=($xref-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $jref=($yref-$t4-$iref*$t5)/$t6;

		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		} else {
		    $iref=$xref;
		    $jref=$yref;
		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		}
		my $x1=int($iref);
		my $y1=int($jref);

#		my $min_region = min ($im($x1:$x2,$y1:$y2));
#		print "Min in box is:  $min_region";
#		my $max_region = max ($im($x1:$x2,$y1:$y2)); 
#		print "Max in box is:  $max_region";
#		
		
		if (defined($cube)) {
		    my $hcube = $cube->gethdr;

		    
#		    print "CUBE info ",$cube->info,"\n";
#		    print  "x1 $x1: x2 $x2,, y1 $y1:y2 $y2 \n";
		    my $spec  = $cube(($x1),($y1),:);

		    print "spec info ",$spec->info,"\n";
		    
		    my $cdelt3 = $$hcube{CDELT3};
		    my $crval3 = $$hcube{CRVAL3};

		    $cdelt3 /= 1E9 if ($$hcube{CTYPE3} =~ /FREQ/ );
		    $crval3 /= 1E9 if ($$hcube{CTYPE3} =~ /FREQ/ );

		    my $lbda = $cdelt3 *  (sequence($$hcube{NAXIS3}) - $$hcube{CRPIX3}) + $crval3;


		    my $specopts  = &Vtools::spec($lbda,$spec,{ViewWin => $win_spec} );
		    
		    
		    if (ref($specopts) =~ /HASH/ ) {
			print "z axis depth>>\n";
			printf "%.2e %.2e  \n ",$specopts->{XMIN},  $specopts->{XMAX} ;

			my $k_1 = minimum_ind(abs($lbda-$specopts->{XMIN}));
			my $k_2 = minimum_ind(abs($lbda-$specopts->{XMAX}));
			my $retsigs;
			$retsigs->{k_1} = $k_1;
			$retsigs->{k_2} = $k_2;
			return $retsigs;
		    }
		    
		    
		}


#		$win_im -> close;
		$genwin=0;

# 		$options->{"MIN"} = $min_region;
# 		$options->{"MAX"} = $max_region;
		@view=($im,$options);
		last;
	    }




	    if (($ch eq $key->{cubespec}->[0]) || ($ch eq $key->{cubespecmed}->[0]) ) {  # $key->{boxgreyscale}->[0]) {

		#new color scale, from the min and max within a box
		my $iref; my $jref;
		my $dim1=$im->getdim(0);
		my $dim2=$im->getdim(1);
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};

		    ($iref, $jref) = &adij_lin($tr,$xref,$yref);

#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $iref=($xref-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $jref=($yref-$t4-$iref*$t5)/$t6;

		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		} else {
		    $iref=$xref;
		    $jref=$yref;
		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		}
		($x, $y, my $ch, $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref, 'Type' => 'Rectangle'});
		my $i; my $j;
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};

		    ($i, $j) = &adij_lin($tr,$x,$y);

#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $i=($x-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $j=($y-$t4-$i*$t5)/$t6;
		    
		    
		    
		    $i=0 if ($i<0);
		    $i=$naxis1-1 if ($i>$naxis1-1);
		    $j=0 if ($j<0);
		    $j=$naxis2-1 if ($j>$naxis2-1);
		} else {
		    $x=0 if ($x<0);
		    $x=$naxis1-1 if ($x>$naxis1-1);
		    $y=0 if ($y<0);
		    $y=$naxis2-1 if ($y>$naxis2-1);
		    $i=$x;
		    $j=$y;
		}
		my $x1=int(min(pdl [$i,$iref]));
		my $x2=int(max(pdl [$i,$iref]));
		my $y1=int(min(pdl [$j,$jref]));
		my $y2=int(max(pdl [$j,$jref]));

#		my $min_region = min ($im($x1:$x2,$y1:$y2));
#		print "Min in box is:  $min_region";
#		my $max_region = max ($im($x1:$x2,$y1:$y2)); 
#		print "Max in box is:  $max_region";
#		
		
		if (defined($cube)) {
		    my $hcube = $cube->gethdr;

		    
#		    print "CUBE info ",$cube->info,"\n";
#		    print  "x1 $x1: x2 $x2,, y1 $y1:y2 $y2 \n";
		    my $minicube  = $cube($x1:$x2,$y1:$y2,:);

		    my $spec;
		    if (($ch eq $key->{cubespecmed}->[0])) {

			
			$spec =  $minicube -> medover -> medover;
		     
		    } else {
			
			$spec =  $minicube -> sumover -> sumover;

		    }
		    print "spec info ",$spec->info,"\n";
		    
		    my $cdelt3 = $$hcube{CDELT3};
		    my $crval3 = $$hcube{CRVAL3};

		    $cdelt3 /= 1E9 if ($$hcube{CTYPE3} =~ /FREQ/ );
		    $crval3 /= 1E9 if ($$hcube{CTYPE3} =~ /FREQ/ );

		    my $lbda = $cdelt3 *  (sequence($$hcube{NAXIS3}) - $$hcube{CRPIX3}) + $crval3;


		    my $specopts  = &Vtools::spec($lbda,$spec,{ViewWin => $win_spec} );
		    
		    
		    if (ref($specopts) =~ /HASH/ ) {
			print "z axis depth>>\n";
			printf "%.2e %.2e  \n ",$specopts->{XMIN},  $specopts->{XMAX} ;

			my $k_1 = minimum_ind(abs($lbda-$specopts->{XMIN}));
			my $k_2 = minimum_ind(abs($lbda-$specopts->{XMAX}));
			my $retsigs;
			$retsigs->{k_1} = $k_1;
			$retsigs->{k_2} = $k_2;
			return $retsigs;
		    }
		    
		    
		}


#		$win_im -> close;
		$genwin=0;

# 		$options->{"MIN"} = $min_region;
# 		$options->{"MAX"} = $max_region;
		@view=($im,$options);
		last;
	    }



	    if ($ch eq $key->{setbad}->[0]) {
		#set to bad a region of the image 
		my $iref; my $jref;
		my $dim1=$im->getdim(0);
		my $dim2=$im->getdim(1);
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};
#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $iref=($xref-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $jref=($yref-$t4-$iref*$t5)/$t6;
		    
		    ($iref, $jref) = &adij_lin($tr,$xref,$yref);


		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		} else {
		    $iref=$xref;
		    $jref=$yref;
		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		}
		($x, $y, my $ch, $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref, 'Type' => 'Rectangle'});
		my $i; my $j;
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};

#		    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		    $i=($x-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		    $j=($y-$t4-$i*$t5)/$t6;

		    ($i, $j) = &adij_lin($tr,$x,$y);

		    $i=0 if ($i<0);
		    $i=$naxis1-1 if ($i>($naxis1-1));
		    $j=0 if ($j<0);
		    $j=$naxis2-1 if ($j>($naxis2-1));
		} else {
		    $x=0 if ($x<0);
		    $x=$naxis1-1 if ($x>($naxis1-1));
		    $y=0 if ($y<0);
		    $y=$naxis2-1 if ($y>($naxis2-1));
		    $i=$x;
		    $j=$y;
		}
		my $x1=int(min(pdl [$i,$iref]));
		my $x2=int(max(pdl [$i,$iref]));
		my $y1=int(min(pdl [$j,$jref]));
		my $y2=int(max(pdl [$j,$jref]));
		my $min_region = min ($im($x1:$x2,$y1:$y2));
		print "Min in box is: %.2e \n", $min_region;
		my $max_region = max ($im($x1:$x2,$y1:$y2)); 
		print "Max in box is: %.2e \n", $max_region;

		$im($x1:$x2,$y1:$y2) .= -1E30;
		$im = $im ->setbadif($im<-1E20);

		$im ->sethdr($hdr);

		$win_im -> erase;
		$genwin=0;
 		$options->{"MIN"} = $min_region;
 		$options->{"MAX"} = $max_region;
		@view=($im,$options);
		last;
	    }

	    if ($ch eq $key->{retim}->[0]) {
		return $im;
	    }

	    if ($ch eq $key->{distance}->[0]) {
		($x, $y, $ch,  $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref, 'Type' => 'RadialLine'});
		if (defined($options->{Transform})) {
		    my $tr=$options->{Transform};

		    my ($i, $j) = &adij_lin($tr,$x,$y);
		    my ($iref, $jref) = &adij_lin($tr,$xref,$yref);
		    my $h = $im->gethdr;

#		    print Dumper($h);
##  		    my $t = &CelCoord::proj($h);
##  		    my $in = pdl($i,$j);
##  		    my $out = $t->apply($in);
##  		    my ($a,$d) = list($out);
##  		    
##  #		    my ($a,$d)=&Ftools::ijad_linear($im,$i,$j);
##  		    my $inref = pdl($iref,$jref);
##  		    my $outref = $t->apply($inref);
##  		    my ($aref,$dref) = list($outref);
##  
##  #		    my ($aref,$dref)=&Ftools::ijad_linear($im,$iref,$jref);
##  
##  		    my $dist=sqrt((($a-$aref)*cos(deg2rad($d)))**2+($d-$dref)**2);

		    my $dist_pix = sqrt( ($i-$iref)**2 + ($j-$jref)**2);

#		    print "scale : $$h{CRVAL2} \n";

		    my $pixscale;
		    if (defined($$h{CDELT1})) {
			$pixscale = $h->{CDELT2} ;
#			$cdelt2 = $h->{CDELT2}  * 3600. if ($hours) ;
#			$cdelt1 = $h->{CDELT1} ;
#			$cdelt1 = $h->{CDELT1} * 3600.0 /(cos(deg2rad($crval2/3600)) * 15) if ($hours) ;
			
		    } elsif (defined($$h{CD1_1})) {

#			$cdelt1 = -1 * sqrt( ($h ->{CD1_1})**2 +  ($h ->{CD2_1})**2 );
			$pixscale = sqrt( ($h ->{CD2_2})**2 +  ($h ->{CD1_2})**2 );
		    }
		    
		    my $dist = $dist_pix * $pixscale;


		    printf "distance: %d pix  %e  deg  %e arcmin %e arcsec\n",$dist_pix, $dist, $dist*60, $dist*3600;
		} else {
#		    (((0<$x) && ($x<$naxis1-1)) && ((0<$y) && ($y<$naxis2-1))) ? print "$ch image( $x, $y) = ", $im->at($x,$y),"\n" : 1;
		    $x=0 if ($x<0);
		    $x=$naxis1-1 if ($x>$naxis1-1);
		    $y=0 if ($y<0);
		    $y=$naxis2-1 if ($y>$naxis2-1);
		    my $dist=sqrt(($x-$xref)**2+($y-$yref)**2);
		    print "angular distance $dist pixels\n";


		}
	    }
	    if ($ch eq $key->{writefits}->[0]) {
		#		my $dum = shift @view;
		#		my $h = $dum2->gethdr;
		#		print Dumper($h);
		#		my $dum = $dum2->copy;
		#		my $dum = zeroes($dum2->dims);
		#		$dum(:,:) .= float($dum2(:,:));
		#		$dum->sethdr($h);
		#		print $dum(100,100);
		my $saveim=$view[0];
		my $h = $saveim->gethdr;
		$saveim -> sethdr($h);
		wfits $saveim,'view.fits';
		print "wrote output file view.fits \n";
	    }
	    if ($ch eq $key->{transpose}->[0]) {
		@view=($im->transpose,$options);
		$genwin = 0;
		last;
	    }


	    if ($ch eq $key->{gaussfit}->[0]) {

		print "Fitting elliptical gaussian \n"; 

#M: Minuit with err bars, no baseline; 
#e: pdl LM; 

		print "d: Minuit, with Leg baseline;  
p: pdl LM with parabolic baseline; 
o: pdl LM with Leg baseline; 
0: pdl LM no baseline \n";

		($x, $y, $ch, $xref, $yref) = $win_im -> cursor({'XRef' =>  $xref,'YRef' => $yref,Type=>'Default'});

		$win_im -> hold; 
		my $imfit=$im->copy; 
		my $nx = $im->getdim(0); 
		my $ny = $im->getdim(1); 

		my $params = &Gauss::adjustellipse($win_im,$im);

#		print ">>>>>>>>>>>>>>>>>>>>>>>>  ",$win->{PlotOptions}->{CURRENT}->{Transform},"\n";
#		print ">>>>>>>>>>>>>>>>>>>>>>>>  ",$win->{PlotOptions}->{DEFAULTS}->{Transform},"\n";



#		my $bestfit = &Gauss::cholesky($imfit,{Pars=>$params,Mode=>'view'});
		my $bestfit;

		if ($ch =~ /d/i ) {
		    print "enter order for Legendre 2D baseline (2 for quadratic, 0 for constant)  \n";
		    my $ans =  <>;
		    $bestfit = &Gauss::cholesky($imfit,{Pars=>$params, Minimizer => 'Minuit', DoErrBars => 0, Baseline => $ans, UseGrad => 0, Mode => ''}); # viewmod
		} elsif ($ch =~ /o/i) {
		    print "strike a digit for the order of a  2D  Legendre baseline (2 for quadratic, 0 for constant)  \n";

		    my ($x, $y, $ans, $xref, $yref) = $win_im -> cursor({'XRef' =>  $xref,'YRef' => $yref,Type=>'Default'});
		    
		    ($bestfit, my $bkg) = &Gauss::cholesky($imfit,{Pars=>$params,LegPoly=>$ans,Mode=> '', ReturnBkg => 1}) if ($ch =~ /o/i); # LegPoly carries order for baseline

		    wfits $bkg, 'bkg.fits';


		} elsif ($ch =~ /0/i) {
		    $bestfit = &Gauss::cholesky($imfit,{Pars=>$params,NoBaseline=>1,Mode=> ''}) if ($ch =~ /0/i); # LegPoly carries order for baseline

		} else {
		    
		    $bestfit = &Gauss::cholesky($imfit,{Pars=>$params,Poly=>1, Mode=> ''}) if ($ch =~ /p/i); # viewmod

		    $bestfit = &Gauss::cholesky($imfit,{Pars=>$params, Minimizer => 'Minuit', DoErrBars => 1}) if ($ch eq 'M');
		    
		    $bestfit = &Gauss::cholesky($imfit,{Pars=>$params, Minimizer => 'PDL_LM'}) if ($ch eq 'e');
		    

		}


		(my $mux, my $muy, my $lmax, my $lmin, my $alpha, my $flux) = &Gauss::translatebestfit($bestfit);

		print "Centroid $mux $muy; Bmaj/2 $lmax; Bmin/2 $lmin;\nPA $alpha rad North of East; peak intensity $flux \n"; 
		my $influx = ($pi/(4*log(2))) * $lmax*2 * $lmin*2 *$flux;
		print "INTEGRATED FLUX $influx \n";

		my $wintr = $win_im->{PlotOptions}->{CURRENT}->{Transform};
		if (defined($wintr)) {
		    print "scaling to wcs coords \n";
		    my $pixscale = abs($wintr(1)->sclr);
		    $lmax *= $pixscale;
		    $lmin *= $pixscale;


		    my $offsets = $options->{offsets};


		    my ($a,$d) = &towcs($file, $offsets, $mux, $muy);

		    print "RA $a DEC $d \n";
		    print "Exact WCS Centroid: ",deg2str($a,'H',4),"  ",deg2str($d,'D',4)," \n";

		    ($mux, $muy) = &Vtools::ijad_linear($wintr, $mux, $muy); 

		    print "Linear WCS Centroid $mux $muy\n";

		    print "Bmaj/2 $lmax; Bmin/2 $lmin;\nPA", $alpha * 180 / $pi, "deg \(South of East\);\n";
		    print "BMAJ=",$lmax*2,"\n";
		    print "BMIN=",$lmin*2,"\n";
		    print "PA=",90.+$alpha * 180 / $pi,"deg (East of Nort)\n";
		    print "inc=", (180.0/ $pi)*acos($lmin/$lmax)," deg\n";

		    #Jy/beam fluxes
		    my $h = $im->gethdr();
		    if ($$h{BUNIT}=~/Jy\/beam/i) {
			my $beam = $pi*$$h{BMAJ}*$$h{BMIN}/(4*log(2));
			print "integrated flux density in Jy (should be equal to intensity for correct `beam' value:\n";
			print ($influx *  ($$h{CDELT2}**2 / $beam ));
			print "\n";

		    }



		}


		my $revert_orient = $wcs?-1:1;

		$win_im->hold; 
		$win_im->ellipse($mux,$muy,$lmax,$lmin, $revert_orient * $alpha,{color => 'red',filltype => 'outline'});

		print "strike r to subtract gaussian, any other key to proceed  \n";

		($x, $y, $ch, $xref, $yref) = $win_im -> cursor({'XRef' =>  $xref,'YRef' => $yref,Type=>'Default'});


		if ( $ch eq 'r')  {
		    print "subtracting model \n";
		    my $im_mod = &Gauss::mgauss(zeroes($im->dims),$bestfit);
		    $im -= $im_mod;

#		    $genwin = 0;

		    my $offsets = $options->{offsets};
		    my $x1 = 0; my $y1 = 0;
		    if (defined($offsets)) {
			my @offs = @$offsets;
			foreach (@offs) {
			    $x1 += $_->[0];
			    $y1 += $_->[2];
			}
		    }

		    my $im_ref = $ori[0] ;
		    $im_ref($x1 : $x1 + $im -> getdim(0) -1 ,  $y1 : $y1 + $im -> getdim(1) -1 ) .= $im;

#			    if ((scalar @view_back) > 0) {
#				my $temp = pop @view_back;
#				$view[0]=$temp->[0];
#				$view[1]=$temp->[1];
#				my $hdr2=$temp->[2];
#				$view[0] -> sethdr($hdr2);
#
#				push @view_back, $temp;
#			    }
#			    	    push @view_back, [$im2,$opts,$h3]; 
#
		    

		    last;
		}
		



		
#		my $a; 
#		if ($ch eq 's') { 
#		    $a = &Gauss::paramstoa($params); 
#		    $a(0:1) .= $a(0:1) +1;
##		    &Gaussf77::choleskylm($imfit,$a);
#		    $a(0:1) .= $a(0:1) -1;
#		    $bestfit = &Gauss::atoparams($a);
#		}


		last;


	    }

	    if ($ch eq $key->{justify}->[0]) {
		if ($nojustify) {
		    print "justfiying pixels ";
		    $nojustify=0;
		} else {
		    print "not justfiying pixels ";
		    $nojustify=1;
		}

		last;
	    }
	    if ($ch eq $key->{unzoom}->[0]) {
		print "unzoom \n";
		my $dumim = $im_ref2->copy;
		my $offsets =  $options->{"offsets"};


		my @offs = @$offsets;

		my $off = shift @offs;

		$options->{"offsets"} = \@offs;
		
#		my @shifts = list $offsets;
		
		my $x1 = $off->[0];
		my $x2 = $off->[1];
		my $y1 = $off->[2];
		my $y2 = $off->[3];

		$dumim($x1:$x2,$y1:$y2) .= $im;
		$im = $dumim;
		$redraw = 1;
		@view=($im,$options);
		$im_ref2 = $im->copy;
#		$options->{"offsets"}  = [$x1,$x2,$y1,$y2];

		last;
	    }

	    if ($ch eq $key->{zoom}->[0]) {
		#zooming
		my $hdr=$im->gethdr;
		my $iref; my $jref;
		my $dim1=$im->getdim(0);
		my $dim2=$im->getdim(1);
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};


		    ($iref, $jref) = &adij_lin($tr,$xref,$yref);


		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		} else {
		    $iref=$xref;
		    $jref=$yref;
		    $iref=0 if ($iref<0);
		    $iref=$naxis1-1 if ($iref>$naxis1-1);
		    $jref=0 if ($jref<0);
		    $jref=$naxis2-1 if ($jref>$naxis2-1);
		}
		($x, $y, my $ch, $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref, 'Type' => 'Rectangle'});
		my $i; my $j;
		if (defined $options->{Transform}) {
		    my $tr=$options->{Transform};


		    ($i, $j) = &adij_lin($tr,$x,$y);


		    $i=0 if ($i<0);
		    $i=$naxis1-1 if ($i>$naxis1-1);
		    $j=0 if ($j<0);
		    $j=$naxis2-1 if ($j>$naxis2-1);
		} else {
		    $x=0 if ($x<0);
		    $x=$naxis1-1 if ($x>$naxis1-1);
		    $y=0 if ($y<0);
		    $y=$naxis2-1 if ($y>$naxis2-1);
		    $i=$x;
		    $j=$y;
		}
		my $x1=int(min(pdl [$i,$iref]));
		my $x2=int(max(pdl [$i,$iref]));
		my $y1=int(min(pdl [$j,$jref]));
		my $y2=int(max(pdl [$j,$jref]));
		my $min_region = min ($im($x1:$x2,$y1:$y2));
		my $max_region = max ($im($x1:$x2,$y1:$y2)); 
		$win_im -> erase;

		$genwin=0;

 		$options->{"MIN"} = $min_region;
 		$options->{"MAX"} = $max_region;

		if (defined($options->{"offsets"})) {
		    my $off = $options->{"offsets"};
		    my @offs = @$off;
		    push @offs, [$x1,$x2,$y1,$y2];
		    $options->{"offsets"} = \@offs;
		} else {
		    my @offs = ( [$x1,$x2,$y1,$y2] );
		    $options->{"offsets"} = \@offs;
		}

		print "Offsets are $x1 $x2 $y1 $y2 \n";
		my $sub_im=$im($x1:$x2,$y1:$y2);
		$$hdr{NAXIS1}=$sub_im->getdim(0);
		$$hdr{NAXIS2}=$sub_im->getdim(1);
#		if (defined $options->{Transform}) {
		if (defined $$hdr{CRPIX1}) {
		    $$hdr{CRPIX1}=$$hdr{CRPIX1}-$x1;
		    $$hdr{CRPIX2}=$$hdr{CRPIX2}-$y1;
		}
		$sub_im->sethdr($hdr);
		@view=($sub_im,$options);
		last;
	    }
	    if ($ch eq $key->{centrezoom}->[0]) {
		my $hdr=$im->gethdr;
		print "zooming about CRVAL .... \n";
		if (!defined $options->{Transform}) {
		    print "Strike P first \n";
		    last;
		}
		print  "Select lower-left corner \n";
		my $h=$im->gethdr;
		my $crval1=$$h{'CRVAL1'};
		my $crval2=$$h{'CRVAL2'};
		my $crpix1=$$h{'CRPIX1'};
		my $crpix2=$$h{'CRPIX2'};
		$win_im -> hold;
		print "crval1 $crval1  crval2 $crval2 \n"; 
		$win_im -> points($crval1,$crval2,{Color => 'Blue'});
		(my $x, my $y, my $ch, $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref, 'Type' => 'CrossHair'});
		my $i; my $j;
		my $tr=$options->{Transform};

#		my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;
#		$i=($x-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
#		$j=($y-$t4-$i*$t5)/$t6;


		my ($i, $j) = &adij_lin($tr,$x,$y);

		$i=0 if ($i<0);
		$i=$naxis1-1 if ($i>$naxis1-1);
		$j=0 if ($j<0);
		$j=$naxis2-1 if ($j>$naxis2-1);
		my $x1=int($i);
		my $x2=int($i+2*($crpix1-$i));
		#for arbitrary rectangular zoom
		my $y1=int($j);
		my $y2=int($j+2*($crpix2-$j));
		#for a square zoom
		$y1=int($crpix2-($crpix1-$i));
		$y2=int($crpix2-($crpix1-$i)+2*($crpix2-$crpix2+($crpix1-$i)));
		my $min_region = min ($im($x1:$x2,$y1:$y2));
		my $max_region = max ($im($x1:$x2,$y1:$y2)); 
		$win_im -> erase;
		$genwin=0;
 		$options->{"MIN"} = $min_region;
 		$options->{"MAX"} = $max_region;


		if (defined($options->{"offsets"})) {
		    my $off = $options->{"offsets"};
		    my @offs = @$off;
		    push @offs, [$x1,$x2,$y1,$y2];
		    $options->{"offsets"} = \@offs;
		} else {
		    my @offs = ( [$x1,$x2,$y1,$y2] );
		    $options->{"offsets"} = \@offs;
		}


#		$options->{"offsets"}  = [$x1,$x2,$y1,$y2];



		my $sub_im=$im($x1:$x2,$y1:$y2);
		$$hdr{NAXIS1}=$x2-$x1+1;
		$$hdr{NAXIS2}=$y2-$y1+1;
#		if (defined $options->{Transform}) {
		if (defined $$hdr{CRPIX1}) {
		    $$hdr{CRPIX1}=$$hdr{CRPIX1}-$x1;
		    $$hdr{CRPIX2}=$$hdr{CRPIX2}-$y1;
		}
		$sub_im->sethdr($hdr);
		@view=($sub_im,$options);
		last;
	    }
	    if ($ch eq $key->{"varcentrezoom"}->[0]) {
		my $hdr=$im->gethdr;
		print "zooming about user supplied center .... switch to pixel coords\n";
		my ($ix,  $iy,  $ch, $xref, $yref) = $win_im -> cursor({'Type' => 'CrossHair'});
		print "using center coords $ix $iy \n";
		print  "Select lower-left corner \n";
		$win_im -> hold;
		$win_im -> points($ix,$iy,{Color => 'Blue'});
		(my $i, my $j,  $ch, $xref, $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref, 'Type' => 'CrossHair'});
		
		$i=0 if ($i<0);
		$i=$naxis1-1 if ($i>$naxis1-1);
		$j=0 if ($j<0);
		$j=$naxis2-1 if ($j>$naxis2-1);
		#for a square zoom
		my $side = abs($ix-$i);
		my $x1=int($i);
		my $x2=int($i+2*$side);
		my $y1=int($iy- $side);
		my $y2=int($iy+ $side);
		print "i $i j $j ix $ix iy $iy im($x1:$x2,$y1:$y2)) \n";
		my $min_region = min ($im($x1:$x2,$y1:$y2));
		my $max_region = max ($im($x1:$x2,$y1:$y2)); 
		$win_im -> erase;
		$genwin=1;
 		$options->{"MIN"} = $min_region;
 		$options->{"MAX"} = $max_region;
		$options->{"offsets"}  = [$x1,$x2,$y1,$y2];
		my $sub_im=$im($x1:$x2,$y1:$y2);
		$$hdr{NAXIS1}=$sub_im->getdim(0);
		$$hdr{NAXIS2}=$sub_im->getdim(1);
		if (defined $options->{Transform}) {
		    $$hdr{CRPIX1}=$$hdr{CRPIX1}-$ix;
		    $$hdr{CRPIX2}=$$hdr{CRPIX2}-$iy;
		}
		$sub_im->sethdr($hdr);
		@view=($sub_im,$options);
		last;
	    }
#	    if ($ch eq "y") {
#		print "Left button to confirm selection, Right button refreshes \n";
#		my $xmin= min $im -> xvals;
#		my $xmax= max $im -> xvals;
#		my $xs= pdl [$xmin,$xmax];
#		my $ys= pdl [$y, $y];
#		$win -> hold;
#		$win -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
#		(my $x2,my  $y2, $ch) = $win -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'HorizontalLine'});
#		if ($ch eq 'A') {$win -> close; return $y}
#		if ($ch eq 'X') {last;}
#	    }

	    if ($ch eq $key->{extractcuts}->[0]) {
		print "Extract a row/column  and plot it \n";
		print "Strike x for row, y for column, s for x-spectrum, u for y-spectrum \n";

		(my $x, my $y, my $ch, my $xref,my  $yref) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref,'Type' => 'CrossHair'});


		my $tr; my $dotransf = 0;
		if (defined $options->{Transform}) {
		    $tr=$options->{Transform};
		    $dotransf = 1;
		}

		if ($ch eq "s") {
		    my $xmin= min $im -> xvals;
		    my $xmax= max $im -> xvals;

		    ($xmin, my $dumy) = &ijad_linear($tr, $xmin, 0.) if ($dotransf);
		    ($xmax, $dumy) = &ijad_linear($tr, $xmax, 0.) if ($dotransf);


		    my $xs= pdl [$xmin,$xmax];
		    my $ys= pdl [$y, $y];
		    $win_im -> hold;
		    $win_im -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    print "Strike s again to select other extraction row, Right button refreshes, \n";
		    (my $x2,my  $y2, $ch) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'HorizontalLine'});
		    $ys= pdl [$y2, $y2];
		    if ($ch eq 's') {


			$win_im -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    }
		    if ($ch eq 'X') {last;}
		    print "Left button to confirm selection, Right button refreshes, \n";
		    ($x2, my  $y3, $ch) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'HorizontalLine'});


		    if ($ch eq 'A') {
			print "EXTRACTING A SPECTRUM FROM $y to $y2 \n";

			(my $dumx, $y) = &adij_lin($tr, 0., $y) if ($dotransf);
			($dumx, $y2) = &adij_lin($tr, 0., $y2) if ($dotransf);

			my $spectrum=sumover($im(:,$y:$y2)->xchg(0,1));

			my $indxs = sequence(nelem($im(:,$y)));

			if ($dotransf) {
			    
			    my ($andxs, $dumx) = &ijad_linear($tr, $indxs,0.);
			    &Vtools::spec($andxs,$spectrum->flat,{REVXAXIS => ($xmin > $xmax)});
			    
			} else {		
			    
			    &Vtools::spec($indxs,$spectrum->flat);
			}
			last;
		    }

		    if ($ch eq 'X') {last;}
		}
		if ($ch eq "u") {
		    my $ymin= min $im -> yvals;
		    my $ymax= max $im -> yvals;

		    (my $dumx, $ymin) = &ijad_linear($tr, 0., $ymin) if ($dotransf);
		    ($dumx, $ymax) = &ijad_linear($tr, 0., $ymax) if ($dotransf);


		    my $xs= pdl [$x, $x];
		    my $ys= pdl [$ymin,$ymax];
		    $win_im -> hold;
		    $win_im -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    print "Strike u again to select other extraction row, Right button refreshes, \n";
		    (my $x2,my  $y2, $ch) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});

		    $xs= pdl [$x2, $x2];

		    if ($ch eq 'u') {
			$win_im -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    }
		    if ($ch eq 'X') {last;}
		    print "Left button to confirm selection, Right button refreshes, \n";
		    (my $x3,my  $y3, $ch) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});

#		    ($x2, $y3) = &adij_lin($tr, $x2, $y3) if ($dotransf);

		    if ($ch eq 'A') {
			print "EXTRACTING A SPECTRUM FROM $x to $x2 \n";

			($x,  $y) = &adij_lin($tr,$x,$y) if ($dotransf);
			($x2, $y2) = &adij_lin($tr,$x2,$y2) if ($dotransf);


			my $indys = sequence(nelem($im($x,:)));
			my $spectrum=sumover($im($x:$x2,:));


			if ($dotransf) {
			    my ($dumxs, $andys) = &ijad_linear($tr, 0.,$indys);
			    &Vtools::spec($andys,$spectrum,{REVXAXIS => ($ymin > $ymax)});

			} else {
			    &Vtools::spec($indys,$spectrum);
			}
			last;


		    }
		    if ($ch eq 'X') {last;}
		}

		if ($ch eq "x") {
		    print "Left button to confirm selection, Right button refreshes, x again to select a range of lines \n";

		    my $xmin= min $im -> xvals;
		    my $xmax= max $im -> xvals;
		    ($xmin, my $dumy) = &ijad_linear($tr, $xmin, 0.) if ($dotransf);
		    ($xmax, $dumy) = &ijad_linear($tr, $xmax, 0.) if ($dotransf);

		    my $xs= pdl [$xmin,$xmax];
		    my $ys= pdl [$y, $y];


		    $win_im -> hold;
		    $win_im -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});



		    (my $x2, my  $y2, $ch) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'HorizontalLine'});

		    if ($ch eq 'x') {
			$win_im -> hold;

			my $ys= pdl [$y2, $y2];

			$win_im -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
			
			(my $dumx, $y) = &adij_lin($tr, 0., $y) if ($dotransf);
			($dumx, $y2) = &adij_lin($tr, 0., $y2) if ($dotransf);

			my $spectrum=sumover($im(:,$y:$y2)->xchg(0,1));
			
			my $indxs = sequence(nelem($im(:,$y)));
			
			if ($dotransf) {
			    
			    my ($andxs, $dumx) = &ijad_linear($tr, $indxs,0.);
			    &Vtools::spec($andxs,$spectrum->flat,{REVXAXIS => ($xmin > $xmax)});
			    
			} else {		
			    &Vtools::spec($indxs,$spectrum->flat);
			}

		    }
		    if ($ch eq 'A') {


			($x, $y) = &adij_lin($tr,$x,$y) if ($dotransf);

			print "extracting cut at column j=$y \n";

			my $indxs = sequence(nelem($im(:,(int($y)))));
			my $spectrum = $im(:,($y))->flat; # 			    &Vtools::spec($indxs,$spectrum->flat);

			if ($dotransf) {
			    
			    my ($andxs, $dumx) = &ijad_linear($tr, $indxs,0.);
			    &Vtools::spec($andxs,$spectrum,{REVXAXIS => ($xmin > $xmax)});
			    
			} else {		
			    
			    &Vtools::spec($indxs,$spectrum);
			}


			$genwin = 0; 

			last;
			
		    }
		    
		    if ($ch eq 'X') {
			
			last;
			
			
		    }
		} elsif ($ch eq "y") {
		    print "Left button to confirm selection, Right button refreshes \n";
		    my $ymin= min $im -> yvals;
		    my $ymax= max $im -> yvals;

		    (my $dumx, $ymin) = &ijad_linear($tr, 0., $ymin) if ($dotransf);
		    ($dumx, $ymax) = &ijad_linear($tr, 0., $ymax) if ($dotransf);

		    my $ys= pdl [$ymin,$ymax];
		    my $xs= pdl [$x, $x];
		    $win_im -> hold;
		    $win_im -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    (my $x2,my  $y2, $ch) = $win_im -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});

		    ($x2, $y2) = &adij_lin($tr, $x2, $y2) if ($dotransf);

		    if ($ch eq 'A') {


			($x, $y) = &adij_lin($tr,$x,$y) if ($dotransf);

			print "extracting cut at line i=$x \n";

			my $indys = sequence(nelem($im((int($x)),:)));
			my $spectrum = $im((int($x)),:)->flat; # 			    &Vtools::spec($indxs,$spectrum->flat);

			if ($dotransf) {
			    my ($dumxs, $andys) = &ijad_linear($tr, 0.,$indys);
			    &Vtools::spec($andys,$spectrum,{REVXAXIS => ($ymin > $ymax)});

			} else {
			    &Vtools::spec($indys,$spectrum);
			}
			last;


		    }
		    if ($ch eq 'X') {last; }
		}
	    }

	    if ($ch eq $key->{'addcont'}->[0]) {
		print "Enter contour levels, separate by spaces \n";
		chomp (my $levelstring = <>);
		my @levels = split(/\s+/,$levelstring);
		my $levs = pdl @levels;
		print "levels: $levs \n";

		$win_im->hold;
#		$win_im -> cont($im,{Levels => $levs, LineWidth => 3, Colour => 'black'});
		$win_im -> cont($view[0],{contours => $levs, LineWidth => 3, Colour => 'black'});



		$win_im->hold;

				
	    }
	    if ($ch eq $key->{lut}->[0]) {
		my $test = 1;
		my @tables = lut_names();
		print "AVAILABLE COLOUR TABLES :\n"; 
		map {print "$_ "} @tables;
		print "\n";
		print "Enter colour table : \n";
		print "(C-D for defaults) \n";
		chomp ($ctbl = <>);
		map {$test = 0 if ($ctbl eq $_)} @tables;
		if (!defined($ctbl) or ($test == 1)) {
		    print "unknown colour table, using default grey\n";
		    $ctbl='ramp';
		    $test = 0;
		}
		$ctbl='ramp' if (!defined($ctbl));
		print "AVAILABLE RAMPS  :\n";
		my @ramps = lut_ramps();
		map {print "$_ "} @ramps;
		print "\n";
		print "Enter ramp : \n";
		print "(C-D for defaults) \n";
		chomp ($ramp = <>);
		$test=1;
		map {$test = 0 if ($ramp eq $_)} @ramps;
		if (!defined($ramp) or ($test == 1)) {
		    print "unknown ramp, using default linear\n";
		    $ramp='ramp';
		}

		print "Enter 0 or 1 (invert): \n";
		print "(C-D for defaults) \n";
		chomp ($inv = <>);
		$inv=0 if (!defined($inv));

		print "Enter the range for the linear color scale : \n";
		print "(C-D for defaults) \n";
		print "min: "; 
		chomp (my $min = <>);
		print "max: ";
		chomp (my $max = <>);
		$options->{MIN}=$min;
		$options->{MAX}=$max;
		$win_im->erase;
		$genwin=0;
		last;
	    }



	    

	    if ($ch eq $key->{hours}->[0]) {
		if ($hours == 0) {
		    $hours = 1 ;
		} else {
		    $hours = 0;
		}
		$win_im->erase;
		$genwin=0;
		last;
	    }
	    if ($ch eq $key->{aspectratio}->[0]) {
		print "Enter new aspect ratio:  \n";
		chomp ($aspect = <>);
		$win_im->close;
		$options->{ASPECT}=$aspect;
		last;
	    }
	    if ($ch eq $key->{printps}->[0]) {
		$win_im -> close;
		$dev='/vcps';
		$genwin = 0; 
		last;
	    }
	    if ($ch eq $key->{printgif}->[0]) {
		$win_im -> close;
		$genwin = 0; 
		$dev='/gif';
		last;
	    }
	    if ($ch eq $key ->{'1'}->[0]) {$win_im->close;$redraw=0;return "1"}
	    if ($ch eq $key ->{'2'}->[0]) {
		$win_im->close;$redraw=0;
		return (1,$options->{"MIN"}, $options->{"MAX"},$ramp) 
		}
	    if ($ch eq $key->{help}->[0]) {
		print "########################KEYSTROKES####################################\n";
		foreach my $fn (keys %$key) {
		    printf "%-5s:  %15s  %50s \n", $key->{$fn}->[0], $fn,    $key->{$fn}->[1] ;
		}

		print "######################################################################\n";
#		print "View keystrokes:\n a: change aspect ratio z: zoom g: set colorscale u: linear greyscale from min to max within a box\n s: stats e:various spectrum extractions P: alpha delta plot  ^H: refresh\n   R: reverse last step X: quit\n k: set to 0 region w: write output view.fits d: distance f: gauss fit \n  ";
	    }
	    if ($ch eq "") {

#		undef $view[1] -> {'offsets'};
#		undef $view[1] -> {'Transform'};
		undef $options;
#
#		$hdr = dclone($hdr_ref);
#		my $im = $im_ref -> copy;
#		$im -> sethdr($hdr);
#		$view[0] = $im;
#		$view[1] .= $ori[1];
#
		my $im_ref=$ori[0]->copy;
		my $h2=$ori[0]->gethdr;
		my $h3;
		if (defined($h2)) {
		    $h3 = deep_copy($h2);
		    $im_ref -> sethdr($h3);
		}
		my $opts = deep_copy($ori[1]);
		@view = ($im_ref, $opts, $h3);
		undef $aspect; $ctbl='ramp';  $ramp='ramp';  $inv =0;last;}



	}
	undef $options;
#	$win->close;
    }
}



sub spec {

    my @ori=@_; 
    if ($#ori < 1) {
	unshift @ori, sequence($ori[0]->nelem);
    }

    my @view=@ori;
    my $redraw=1;

    my $dev='/xserve';
    $stack_spec=1; #stack switch  
    my $genwin_spec=1;
    my $x; my $y;
    my %argsref = (REVXAXIS => 0, Dev=> '/xserve', XTitle => '', YTitle => '', autolog => 0, Pointsoverplot => 0, Points => 0); 

    
    my $specref= $ori[1]->copy;
    my $lbdaref= $ori[0]->copy;
    $specref = $specref -> setbadif(!isfinite($specref));

    my $options0=$ori[2];
    if (defined($options0)) {
	%argsref = (%argsref, %$options0) ;
    }
    $dev = $argsref{Dev};
    
    $options0 = \%argsref;  ### DEV SEP 2011


    my @view=($lbdaref,$specref,deep_copy(\%argsref));

    if (defined($options0->{ViewWin})) {
	$win_spec = $options0->{ViewWin};
	$genwin_spec = 0;
    }
    while ($redraw) {
#	print "BEGIN LOOP STACK $stack_spec \n ";
#	print Dumper(@view);
	if ($stack_spec) {
	    if (defined($view[2])) {
		my $opts_bis = deep_copy($view[2]);

#	    my $dum=$view[2];#
#	    my %options_bis=%$dum if defined($dum);
#	    push @view_back_spec,\%options_bis;

		
		push @view_back_spec,$opts_bis;

	    }
	}
	$stack_spec=1;
	my $spec= $view[1];
	my $lbda= $view[0];
	$spec = $spec -> setbadif(!isfinite($spec));
	my $options0=$view[2];
#	if (defined($options0)) {#
#	    %argsref = (%args, %$options0) ;
#	}

#	my $options = \%args;

	my $options = $options0;

	my $code=$options->{'code'};
	my $sig = $options->{'ERRORB'};
	my $overplot = $options->{'Overplot'};
	my $colors = $options->{'Colors'};

	my $xtitle = $options ->{'XTitle'};
	my $ytitle = $options ->{'YTitle'};

#	print Dumper($options0);
#	die "xtitle $xtitle \n";

	my $aspect=$options->{'ASPECT'};
	my $charsize=$options->{'CHARSIZE'};
	my $points=$options->{'Points'};

	my $pointsoverplot = $options->{'Pointsoverplot'};

	my $log=$options->{'autolog'};

	my $naxis1 = $spec -> getdim(0);
	if (!(defined $aspect)) {$aspect=0.5;};

	if (!(defined $charsize)) {$charsize = 1.5;};

	

	if ($genwin_spec) {
	    if ($log) {
		$win_spec = PDL::Graphics::PGPLOT::Window->new(Device => $dev, WindowWidth => 8, AspectRatio => $aspect, AxisColour => 'BLACK', Colour => 'BLACK', LINEWIDTH => 1, Charsize => $charsize, XTitle => $xtitle, YTitle => $ytitle); #, axis => 'LOGXY'
		$win_spec->autolog(1);
	    } else {
		$win_spec = PDL::Graphics::PGPLOT::Window->new(Device => $dev, WindowWidth => 8, AspectRatio => $aspect, AxisColour => 'BLACK', Colour => 'BLACK', LINEWIDTH => 1, Charsize => $charsize, XTitle => $xtitle, YTitle => $ytitle);

	    }
	}
	$genwin_spec=0;
	my $ymin=$options->{'YMIN'};
	my $ymax=$options->{'YMAX'};
	my $xmin=$options->{'XMIN'};
	my $xmax=$options->{'XMAX'};
	if (!defined($ymin)){$ymin=min($spec)}
	if (!defined($ymax)){$ymax=max($spec)}
	if (!defined($xmin)){$xmin=min($lbda)}
	if (!defined($xmax)){$xmax=max($lbda)}


	my $axis = ['BCNTS','BCNTS'] ;
	if ($log) {
	    $axis = ['LBCNTS','LBCNTS'] ;
	}
	if ($options -> {'REVXAXIS'}) {
	    $win_spec -> env($xmax,$xmin,$ymin,$ymax,{Axis => $axis});
	} else {
	    $win_spec -> env($xmin,$xmax,$ymin,$ymax,{Axis => $axis});
	}


	if (defined($sig)) {
	    $win_spec -> errb($lbda,$spec,$sig);
	} elsif ($points) {
	    $win_spec -> points($lbda,$spec,{Symbol => 1}); #, {COLOUR => 'BLACK'});
	} else {
	    $win_spec -> line($lbda,$spec); #, {COLOUR => 'BLACK'});
	} 
	my $style=1;
	if (!defined($colors)) {
	    $colors = pdl [2,3,4,5,6,7,8,14,16,1] ;
	}
	if (defined($overplot)) {
	    my $icolor=0;
	    foreach  (@$overplot) {

		$win_spec -> hold;
		my $color = $style+1;
		if (defined($colors)) {
		    my $idum  = $icolor % nelem($colors) ;
		    $color = $colors($idum)->sclr;
		} 

		if ($pointsoverplot) {
		    $win_spec -> points($_->[0],$_->[1],{Color=>$color, Symbol => 1});
		} else {
#		    print $_->[0]->info,"\n";
#		    print $_->[1]->info,"\n";
		    $win_spec -> line($_->[0],$_->[1],{LineStyle=>$style,Color=>$color});
		}
		$style = 1+$style % 4;
		$icolor++;
	    }
	}
	if (defined($code)) {eval $code; print "error in code: $@ \n"  if ($@) ; }
	if ( ($dev =~ 'ps') && ($argsref{Dev}=~/ps/)) {$win_spec -> close; return}
	if ($dev =~ 'ps') {$win_spec -> close; $dev = '/xserve'; $genwin_spec = 1; next}
	if ($dev =~ 'gif') {$win_spec -> close; $dev = '/xserve'; next}
	my $xref; my $yref;
	while(1){
	    ($x, $y, my $ch,  $xref,  $yref) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,'Type' => 'CrossHair'});
	    $x=min($lbda) if ($x<min($lbda));
	    $x=max($lbda)  if ($x>max($lbda));
	    $xref=$x;
	    $yref=$y;

	    my $key = {
		'quit'       =>  ['X', 'right-click, quit Vtools'], 
		'killquit'   =>  ['q', 'die, quit'], 
		'retoptions' =>  ['0', 'quit, return display settings'],
		'medover'    =>  ['m', 'like retoptions with medover'],
		'1'          =>  ['1', 'quit, return \'select\''],
		'cubesum'    =>  ['d', 'sumover cube projection depth'],
		'cubemed'    =>  ['e', 'medover cube projection depth'],
		'reverse'    =>  ['r', 'undo'],
		'points'     =>  ['D', 'change line/points'],
		'printps'    =>  ['p', 'print to postscript pgplot.ps'],
		'printgif'   =>  ['o', 'print to gif, pgplot.gif'],
		'autolog'    =>  ['L', 'toggle log axis'],
		'Yrange'     =>  ['Y', 'STDIN adjust  y-scale'],
		'yrange'     =>  ['y', 'cursor adjust  y-scale'],
		'arange'     =>  ['R', 'auto  adjust  y-scale'],
		'xrange'     =>  ['x', 'cursor adjust  x-scale'],
		'gaussfit0'  =>  [';', 'builtin no-baseline gaussian fit'],
		'gaussfit'   =>  ['g', 'baseline gaussian fit'],
		'mgaussfit'  =>  ['z', 'multi  gaussian fit'],
		'gaussnobas' =>  ['G', 'gaussian fit'],
		'varfit'     =>  ['f', 'various fits'],
		'stats'      =>  ['s', 'rough stats'],
		'toascii'    =>  ['w', 'dump spec.dat'],
		'toascii2'   =>  ['W', 'like w, with noise'],
		'invert'     =>  ['-', 'take 1/spec'],
		'help'       =>  ['h', 'list keystrokes'],
		'reset'      =>  ['^H', 'starting defaults'],
 		'aspectratio' => ['a', 'change image aspect ratio'],
 		'charsize'   => ['c', 'change charsize']
	    };


	    if ($ch eq 'A') {
		print "x = $x, y = $y \n" ;
		
	    }
	    if ($ch eq $key->{aspectratio}->[0]) {
		print "Enter new aspect ratio:  \n";
		chomp ($aspect = <>);
		$win_spec->close;
		$genwin_spec = 1;
		$options->{ASPECT}=$aspect;
		last;
	    }
	    if ($ch eq $key->{charsize}->[0]) {
		print "Enter new char size:  \n";
		chomp ($charsize = <>);
		$win_spec->close;
		$genwin_spec = 1;
		$options->{CHARSIZE}=$charsize;
		last;
	    }

	    if ($ch eq $key->{quit}->[0]) {
		

		if (defined($cube)) {
		    $win_spec->close;
		    $genwin_spec = 1;
		    return;
		} else {
		    $redraw=0; $win_spec->close; last;
		}
	    }


	    if ($ch eq $key->{1}->[0]) {return 'select';}
	    if ($ch eq $key->{reverse}->[0]) { #reverse
		$win_spec -> erase; $genwin_spec = 0;
		$view[2] = pop @view_back_spec;
		print "number of elts in view_back_spec is ", scalar @view_back_spec,"\n";
		if (!defined($view[2])) {
		    @view=($lbdaref,$specref,\%argsref);
		}

		$stack_spec=0;
		last;
	    }

	    if ($ch eq $key->{retoptions}->[0]) {
		$win_spec -> close;
		return $options;
	    }

 	    if ($ch eq $key->{cubesum}->[0]) {
		
		print "\'d\' again to select other projection boundary, Right button refreshes \n";
		my $ymin= min $spec;
		my $ymax= max $spec;
		my $xs= pdl [$x,$x];
		my $ys= pdl [$ymin,$ymax];
		$win_spec -> hold;
		$win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		(my $x2, my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		my $allxs;

		if ($ch eq 'd') {
		    my $xs= pdl [$x2, $x2];
		    $allxs= pdl [$x,$x2];
		    $xs= pdl [$x2, $x2];
		    $win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    print "Left button to confirm selection,  Right button refreshes \n";
		    (my $x2,my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		    if ($ch eq 'A') {

			
			my $x_1 = min($allxs);
			my $x_2 = max($allxs);

			my $k_1 = minimum_ind(abs($lbdaref - $x_1)) ;
			my $k_2 = minimum_ind(abs($lbdaref - $x_2)) ;

			my $cbz = $cube(:,:,$k_1:$k_2);
			my $im_c = $cbz->xchg(0,2)->sumover->xchg(0,1);
			$im_c -> sethdr($h_ref);
			
			print "k_1 $k_1  k_2 $k_2 \n";
			
#			$win_im -> close;
			

			&Vtools::view($im_c,{ViewWin => $win_im});

			$genwin_spec = 0;
			last;
		    } 		    
		    if ($ch eq 'X') {last;}
		}
		if ($ch eq 'X') {last;}
		
#		$win_spec -> close;#
#		return $options;
	    }

 	    if ($ch eq $key->{cubemed}->[0]) {
		
		print "\'e\' again to select other projection boundary, Right button refreshes \n";
		my $ymin= min $spec;
		my $ymax= max $spec;
		my $xs= pdl [$x,$x];
		my $ys= pdl [$ymin,$ymax];
		$win_spec -> hold;
		$win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		(my $x2, my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		my $allxs;

		if ($ch eq 'e') {
		    my $xs= pdl [$x2, $x2];
		    $allxs= pdl [$x,$x2];
		    $xs= pdl [$x2, $x2];
		    $win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    print "Left button to confirm selection,  Right button refreshes \n";
		    (my $x2,my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		    if ($ch eq 'A') {

			
			my $x_1 = min($allxs);
			my $x_2 = max($allxs);

			my $k_1 = minimum_ind(abs($lbdaref - $x_1)) ;
			my $k_2 = minimum_ind(abs($lbdaref - $x_2)) ;

			my $cbz = $cube(:,:,$k_1:$k_2);
			my $im_c = $cbz->xchg(0,2)->medover->xchg(0,1);
			$im_c -> sethdr($h_ref);
			
			print "x_1 $x_1 k_1 $k_1 x_2 $x_2 k_2 $k_2 \n";
			
#			$win_im -> close;
			

			&Vtools::view($im_c,{ViewWin => $win_im});

			$genwin_spec = 0;

			last;
		    } 		    
		    if ($ch eq 'X') {last;}
		}
		if ($ch eq 'X') {last;}
		
#		$win_spec -> close;#
#		return $options;
	    }


	    if ($ch eq $key->{medover}->[0]) {
		$genwin_spec = 0;
		$win_spec -> close;
		$options->{Medover} = 1;
		return $options;
	    }

	    
	    if ($ch eq $key->{points}->[0]) {
		$win_spec -> points($lbda,$spec);
	    }
	    if ($ch eq $key->{killquit}->[0]) {
		$win_spec -> close;
		die;
	    }
	    if ($ch eq $key->{printps}->[0]) {
		$win_spec -> close;
		$genwin_spec = 1;
		$dev='/vcps';
		last;
	    }
	    if ($ch eq $key->{printgif}->[0]) {
		$win_spec -> close;
#		print "Enter device, choose /gif \n";
		$dev='/gif';
#		print "using $dev \n";
		last;
	    }
	    
	    if ($ch eq $key->{autolog}->[0]) {
		$log = $win_spec -> autolog;
		if ($log) {
		    $win_spec->autolog(0);
		}  else {
		    $win_spec->autolog(1);
		}
#		$win_spec -> close;
		last;
		
	    }

#	    if ($ch eq "L") {
#		$win_spec -> autolog(1);
#	    }


	    if ($ch eq  $key->{arange}->[0]) {
		$win_spec -> erase; $genwin_spec = 0 ; 

		my $i1 = minimum_ind(abs($lbda - $options->{XMIN})); 
		my $i2 = minimum_ind(abs($lbda - $options->{XMAX})); 
		$options->{YMIN}= min($spec($i1:$i2));
		$options->{YMAX}= max($spec($i1:$i2));
		@view=($lbda,$spec,$options);
		last;
	    }

	    if ($ch eq  $key->{Yrange}->[0]) {
		print "Enter Ymin: ";
		$ymin=<>;
		print "Enter Ymax: ";
		$ymax=<>;
		$win_spec -> erase; $genwin_spec = 0 ; 
		$options->{YMIN}= $ymin;
		$options->{YMAX}= $ymax;		
		@view=($lbda,$spec,$options);
		last;
	    }
	    if ($ch eq $key->{yrange}->[0]) {
		print "\'y\' again to select other zooming boundary, Middle button to select bottom, Right button refreshes \n";
		my $xmin= min $lbda;
		my $xmax= max $lbda;
		my $xs= pdl [$xmin,$xmax];
		my $ys= pdl [$y, $y];
		$win_spec -> hold;
		$win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		(my $x2,my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'HorizontalLine'});
		my $allys= pdl [$y,$y2];
		$ys= pdl [$y2, $y2];
		if ($ch eq 'y') {
		    $win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    print "Left button to confirm selection,  Right button refreshes \n";
		    (my $x2,my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'HorizontalLine'});
		    if ($ch eq 'A') {
			$win_spec -> erase; 
			$genwin_spec=0; 
			$options->{YMIN}= min($allys);
			$options->{YMAX}= max($allys);
			@view=($lbda,$spec,$options); last;} 
		    if ($ch eq 'X') {last;}
		}
		if ($ch eq 'D') {
		    $win_spec -> erase; $genwin_spec =0 ; 
		    $options->{YMIN}= min($spec);
		    $options->{YMAX}= max($allys);
		    @view=($lbda,$spec,$options);
		    last;
		}
		if ($ch eq 'X') {last;}
	    }
	    if ($ch eq $key->{xrange}->[0]) {
		print "\'x\' again to select other zooming boundary, Right button refreshes \n";
		my $ymin= min $spec;
		my $ymax= max $spec;
		my $xs= pdl [$x,$x];
		my $ys= pdl [$ymin,$ymax];
		$win_spec -> hold;
		$win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		(my $x2, my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		if ($ch eq 'x') {
		    my $xs= pdl [$x2, $x2];
		    my $allxs= pdl [$x,$x2];
		    $xs= pdl [$x2, $x2];
		    $win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    print "Left button to confirm selection,  Right button refreshes \n";
		    (my $x2,my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		    if ($ch eq 'A') {$win_spec -> erase; $genwin_spec=0; 
				     $options->{XMIN} =min($allxs);
				     $options->{XMAX} =max($allxs);
				     @view=($lbda,$spec,$options);
#				     print "STACK $stack_spec\n";
				     last;} 		    

		    if ($ch eq 'X') {last;}
		}
		if ($ch eq 'X') {last;}
	    }

	    if ($ch eq $key->{gaussfit}->[0]) {
		print "gaussian fit: select wavelength range left button \n";
		print "gaussian fit: init sigma is 1/10 wavelength range \n";
		my $ymin= min $spec;
		my $ymax= max $spec;
		my $ys=pdl [$ymin,$ymax];
		my $x1; my $x2;

		$ch = 0;
		until  ($ch eq 'A') {
		    ($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		}
		my $xs= pdl [$x, $x];
		$win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		$x1=$x;

		$ch = 0;
		until  ($ch eq 'A') {
		    ($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		}
		my $xs= pdl [$x, $x];
		$win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		$x2=$x;

		my $ns = $lbda -> xvals;
		my $i1 = minimum_ind(abs($lbda-$x1));
		my $i2 = minimum_ind(abs($lbda-$x2));

		my $allxs= pdl [$x1,$x2];
		my $is= pdl [$i1,$i2];
		my $n1=int(min($is));
		my $n2=int(max($is));
		
		my $subspec = $spec($n1:$n2);
		my $sublbda = $lbda($n1:$n2);
		my $dlbda = (max($sublbda)-min($sublbda)) / nelem($sublbda);

		print "gaussian fit: place centroid \n";
		$ch = 0;
		until  ($ch eq 'A') {
		    ($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		} 
		my $xs= pdl [$x, $x];
		$win_spec -> line ($xs, $ys, {COLOR => 'blue', LINESTYLE => 3});

		my $icentroid = minimum_ind(abs($sublbda-$x));
		my $amp = $subspec($icentroid)-median($subspec);
		print "gaussian fit: strike a digit for the order of 1D Legendre baseline (2 for quadratic, 0 for constant)  \n";
		($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'Default'});
		my $nord_base = $ch;


		my $opts;
		$opts -> {Init} = {Centroid => $icentroid , Sigma => (abs($n2-$n1)/10.), Amplitude => $amp} ; 
		$opts -> {Baseline} = $nord_base;
		
		my @retvals = &Gauss1D::fit($subspec,$opts);
		
		my $a = $retvals[0];
		print "bestfit  $a \n";
		my $specfit = $retvals[1];
		my $continuum = $retvals[2];

		
		my ($l0, $err) = interpolate($a(0)->sclr,sequence(nelem($sublbda)),$sublbda);
		my ($Fc0, $err) = interpolate($a(0)->sclr,sequence(nelem($sublbda)),$continuum); # continuum at line centre

		my $sigma = $a(1)->sclr * $dlbda;
		my $FWHM = $sigma*2*sqrt(2*log(2));
		my $amplitude = $a(2)->sclr;
		
		my $linearea = sum($continuum-$specfit)*$dlbda;
		my $EQW = $linearea/$Fc0;
		

		printf "centroid: %.6e  FWHM: %.2e amplitude %.2e EQW %.2e  \n", $l0, $FWHM, $amplitude, $EQW;
		$win_spec -> line ($sublbda, $specfit, {COLOR => 'green', LINESTYLE => 3});
		$win_spec -> line ($sublbda, $continuum, {COLOR => 'blue', LINESTYLE => 3});

	    }		

	    if ($ch eq $key->{mgaussfit}->[0]) {
		print "multi-gaussian fit: enter number of Gaussians  \n";
		my $ngauss = <>;

		print "multi-gaussian fit: enter the order of 1D Legendre baseline (2 for quadratic, 0 for constant)  \n";
		my $nord_base = <>;
		
		print "multi-gaussian fit: select wavelength range left button \n";
		my $ymin= min $spec;
		my $ymax= max $spec;
		my $ys=pdl [$ymin,$ymax];
		my $x1; my $x2;
		
		$ch = 0;
		until  ($ch eq 'A') {
		    ($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		}
		my $xs= pdl [$x, $x];
		$win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		$x1=$x;
		
		$ch = 0;
		until  ($ch eq 'A') {
		    ($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		}
		my $xs= pdl [$x, $x];
		$win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		$x2=$x;
		
		my $ns = $lbda -> xvals;
		my $i1 = minimum_ind(abs($lbda-$x1));
		my $i2 = minimum_ind(abs($lbda-$x2));
		
		my $allxs= pdl [$x1,$x2];
		my $is= pdl [$i1,$i2];
		my $n1=int(min($is));
		my $n2=int(max($is));
		
		my $subspec = $spec($n1:$n2);
		my $sublbda = $lbda($n1:$n2);
		my $dlbda = (max($sublbda)-min($sublbda)) / nelem($sublbda);
		
		my @gaussinits;
		foreach my $ng (0..$ngauss-1) {
		    
		    print "gaussian fit: place centroid \n";
		    $ch = 0;
		    until  ($ch eq 'A') {
			($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		    } 
		    my $xs= pdl [$x, $x];
		    $win_spec -> line ($xs, $ys, {COLOR => 'blue', LINESTYLE => 3});

		    my $icentroid = minimum_ind(abs($sublbda-$x));
		    my $amp = $subspec($icentroid)-median($subspec);


		    print "gaussian fit: select HWHM  from centroid \n";

		    $ch = 0;
		    until  ($ch eq 'A') {
			($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		    }
		    my $xs= pdl [$x, $x];
		    $win_spec -> line ($xs, $ys, {COLOR => 'blue', LINESTYLE => 3});
		    $x1=$x;
		    my $iwidth = minimum_ind(abs($lbda-$x1));

		    my $sigma = abs($iwidth-$icentroid) / (sqrt(2*log(2))); 

		    
		    print Dumper({Centroid => $icentroid->sclr , Sigma => $sigma->sclr, Amplitude => $amp->sclr});

		    push @gaussinits, {Centroid => $icentroid->sclr , Sigma => $sigma->sclr, Amplitude => $amp->sclr}; 

		}

		my $opts;
		$opts -> {Init} = \@gaussinits; 
		$opts -> {Baseline} = $nord_base;
		$opts -> {NGauss} = $ngauss;

		my @retvals = &Gauss1D::mfit($subspec,$opts);
		
		my $a = $retvals[0];
		print "bestfit  $a \n";
		my $specfit = $retvals[1];
		my $continuum = $retvals[2];

		
		
		foreach my $ng (0..$ngauss-1) {
		    print "ng = $ng \n";

		    my ($l0, $err) = interpolate($a(0+3*$ng)->sclr,sequence(nelem($sublbda)),$sublbda);
		    my $sigma = $a(1+3*$ng)->sclr * $dlbda;
		    my $FWHM = $sigma*2*sqrt(2*log(2));
		    my $amplitude = $a(2+3*$ng)->sclr;
		    
		    printf "{ centroid=> %.6e,  FWHM => %.2e, amplitude => %.2e } \n", $l0, $FWHM, $amplitude;
		    
		}

		#my ($Fc0, $err) = interpolate($a(0)->sclr,sequence(nelem($sublbda)),$continuum); # continuum at line centre
		#my $linearea = sum($continuum-$specfit)*$dlbda;
		#my $EQW = $linearea/$Fc0;
		#print "EQW = $EQW \n";

		$win_spec -> line ($sublbda, $specfit, {COLOR => 'green', LINESTYLE => 3});
		$win_spec -> line ($sublbda, $continuum, {COLOR => 'blue', LINESTYLE => 3});

		wcols $sublbda,$specfit,$continuum,'mgaussfit.dat'

	    }		

	    if ($ch eq $key->{gaussnobas}->[0]) {
		print "gaussian fit: select wavelength range left button \n";
		print "gaussian fit: init sigma is 1/10 wavelength range \n";
		print "gaussian fit: strike s to skip and use displayed range \n";
		my $ymin= min $spec;
		my $ymax= max $spec;
		my $ys=pdl [$ymin,$ymax];
		my $x1; my $x2;

		$ch = 0;
		until  (($ch eq 'A') || ($ch eq 's')) {
		    ($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		}
		if ($ch eq 's') {
		    $x1 = min($lbda);
		    $x2 = max($lbda);
		} else {
		    my $xs= pdl [$x, $x];
		    $win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    $x1=$x;
		    
		    $ch = 0;
		    until  ($ch eq 'A') {
			($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		    }
		    my $xs= pdl [$x, $x];
		    $win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    $x2=$x;
		}

		my $ns = $lbda -> xvals;
		my $i1 = minimum_ind(abs($lbda-$x1));
		my $i2 = minimum_ind(abs($lbda-$x2));

		my $allxs= pdl [$x1,$x2];
		my $is= pdl [$i1,$i2];
		my $n1=int(min($is));
		my $n2=int(max($is));
		
		my $subspec = $spec($n1:$n2);
		my $sublbda = $lbda($n1:$n2);
		my $dlbda = (max($sublbda)-min($sublbda)) / nelem($sublbda);

		print "gaussian fit: place centroid with left click, or \n";
		print "gaussian fit: strike  \'a\' to place centroid and store the resulting fit into file gaussfit.dat \n";
		$ch = 0;
		my $storefit = 0;
		until  ( ($ch eq 'A') || ($ch eq 'a') ) {
		    ($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		} 
		$storefit = 1 if ($ch eq 'a');
		my $xs= pdl [$x, $x];
		$win_spec -> line ($xs, $ys, {COLOR => 'blue', LINESTYLE => 3});

		

		my $icentroid = minimum_ind(abs($sublbda-$x));
		my $amp = $subspec($icentroid)-median($subspec);


		my $init_sigma = (abs($n2-$n1)/10.);

		my $isubs = sequence(nelem($subspec));
		my $mask =  ( ($isubs < ($icentroid - $init_sigma)) | ($isubs > ($icentroid + $init_sigma)) );
		my ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($subspec->where($mask));

		print "noise outside +-sigma from centroid: ",$rms,"\n";
		my $noise = ones(nelem($subspec)) * $rms; #(1/$rms**2); 

		my $opts;
		$opts -> {Init} = {Centroid => $icentroid , Sigma => $init_sigma , Amplitude => $amp} ; 
#		$opts -> {Baseline} = $nord_base;
		
		$opts ->{Noise} =  $noise;
		$opts ->{Minimizer} = 'Minuit';
		$opts ->{ConservativeErrors} = 1; # set to 0 for errors on amplitude only, keeping the rest at best fit
		my @retvals = &Gauss1D::fit_nobase($subspec,$opts);
		
		my $a = $retvals[0];
		print "bestfit  $a \n";
		my $specfit = $retvals[1];
#		my $continuum = $retvals[2];
		
		
		my $residual = ($subspec - $specfit);
		my ($mean2,$prms2,$median2,$min2,$max2,$adev2,$rms2) = stats($residual);

		print "COMPARE INIT NOISE VALUE: $rms WITH RESIDUAL NOISE: $rms2 \n";

		my $gsig;

		my $downsig;
		my $upsig;

		if ($opts ->{Minimizer} eq 'LM') {

		    my $covar = $retvals[2];
		    $gsig = pdl( sqrt($covar(0,0))->sclr,  sqrt($covar(1,1))->sclr , sqrt($covar(2,2))->sclr);
		    print " covar = $covar \n";

		} else {
		    $downsig = $retvals[2];
		    $upsig  = $retvals[3];
		    $gsig = $retvals[4];

#		    wcols $downsig, $upsig, $gsig;


		}

#		my $dlbda = (max($sublbda)- min($sublbda))/nelem($sublbda);
		
		my ($l0, $err) = interpolate($a(0)->sclr,sequence(nelem($sublbda)),$sublbda);
#		my ($Fc0, $err) = interpolate($a(0)->sclr,sequence(nelem($sublbda)),$continuum); # continuum at line centre
		
		my $dl0 = $gsig(0)->sclr * $dlbda;

		my $sigma = $a(1)->sclr * $dlbda;
		my $dsigma = $gsig(1)->sclr;

		my $FWHM = $sigma*2*sqrt(2*log(2));
		my $dFWHM = $dsigma*2*sqrt(2*log(2));


		my $amplitude = $a(2)->sclr;
		my $damplitude = $gsig(2)->sclr;

#		my $linearea = sum($continuum-$specfit)*$dlbda;
#		my $EQW = $linearea/$Fc0;
		

#		printf "centroid: %.6e  FWHM: %.2e amplitude %.2e EQW %.2e  \n", $l0, $FWHM, $amplitude, $EQW;


		print "input noise: $noise \n";
		printf "centroid: %.6e  FWHM: %.2e amplitude %.2e  \n", $l0, $FWHM, $amplitude;
		printf "errors  : %.6e  FWHM: %.2e amplitude %.2e  \n", $dl0, $dFWHM, $damplitude;
		
		my $header = sprintf "\#input noise  $noise  \n";
	        $header .= sprintf "\#centroid: %.6e  FWHM: %.2e amplitude %.2e \n", $l0, $FWHM, $amplitude;
		$header .= sprintf "\#std errors  : %.6e  FWHM: %.2e amplitude %.2e  \n", $dl0, $dFWHM, $damplitude;
		if ($opts->{Minimizer} eq 'Minuit') {
		    
		    my ($upl0, $err) = interpolate($upsig(0)->sclr,sequence(nelem($sublbda)),$sublbda);
		    my ($dol0, $err) = interpolate($downsig(0)->sclr,sequence(nelem($sublbda)),$sublbda);

		    my $upsigma =  $upsig(1)->sclr * $dlbda;
		    my $dosigma =  $downsig(1)->sclr * $dlbda;
		    my $upFWHM = $upsigma * 2*sqrt(2*log(2));
		    my $doFWHM = $dosigma * 2*sqrt(2*log(2));

		    my $upamplitude = $upsig(2)->sclr;
		    my $doamplitude = $downsig(2)->sclr;

		    printf "formal errors:\n $dol0<centroid<$upl0\n $doFWHM<FWHM<$upFWHM\n $doamplitude<amp<$upamplitude\n";
		    $header .= sprintf "\#formal errors:\n\# $dol0<centroid<$upl0\n\# $doFWHM<FWHM<$upFWHM\n\# $doamplitude<amp<$upamplitude";

		    $header .= sprintf "\#input rms \n";

		    print "compare std and formal errors:";
		    wcols abs($upsig - $downsig)/2, $gsig,{HEADER => "formal std"};
			


		}


		
		$win_spec -> line ($sublbda, $specfit, {COLOR => 'green', LINESTYLE => 3});

		wcols $sublbda, $specfit, 'gaussfit.dat', {HEADER => $header} if $storefit;

#		$win_spec -> line ($sublbda, $continuum, {COLOR => 'blue', LINESTYLE => 3});

		

	    }		


	    if ($ch eq $key->{gaussfit0}->[0]) {

		print "gaussian fit: select wavelength range left button \n";
		($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		my $ymin= min $spec;
		my $ymax= max $spec;
		my $ys=pdl [$ymin,$ymax];
		my $x1; my $x2;
		if ($ch eq 'A') {
		    my $xs= pdl [$x, $x];
		    $win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    $x1=$x;
		}
		($x, $y, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'VerticalLine'});
		if ($ch eq 'A') {
		    my $xs= pdl [$x, $x];
		    $win_spec -> line ($xs, $ys, {COLOR => 'RED', LINESTYLE => 3});
		    $x2=$x;
		}
		my $ns = $lbda -> xvals;
		my $nn1 = max($ns -> flat -> index(which( ($lbda-$x1) <  0)));
		my $nn2 = max($ns -> flat -> index(which( ($lbda-$x2) <  0)));

		$nn1 = $ns -> where(abs($lbda-$x1) == min(abs($lbda-$x1))) 
		    ->sclr;
		$nn2 = $ns -> where(abs($lbda-$x2) == min(abs($lbda-$x2))) 
		    ->sclr;

		my $allxs= pdl [$x1,$x2];
		my $nns= pdl [$nn1,$nn2];
		my $n1=int(min($nns));
		my $n2=int(max($nns));
		print "x1 $x1, x2 $x2 \n";
		print "n1 $n1, n2 $n2 \n";
		(my $cen,my  $pk,my $fwhm,my $back,my $err,my $fit) = fitgauss1d($lbda->slice("$n1:$n2"), $spec->slice("$n1:$n2"));
		print "cen $cen, pk $pk,fwhm $fwhm,back $back,err $err \n";
		$win_spec->points($lbda->slice("$n1:$n2"),$fit); 
	    }		

	    if ($ch eq $key->{varfit}->[0]) {
		my $i1=max(which($lbda < $xmin));
		my $i2=min(which($lbda > $xmax));
		my $x=$lbda($i1:$i2);
		my $y=$spec($i1:$i2);
		my $sig=$x/$x;
		print "Fitting options: \n 
                       g for built-in PDL gaussian fitter \n
                       l for levenberg-marquardt poli+gauss fit \n
                       n for poli fit, and noise estimate \n";
		(my $x2,my  $y2, $ch) = $win_spec -> cursor({'XRef' => $xref,'YRef' => $yref,Type => 'Default'});
		
		my $ym;
		if ($ch eq "l") {
		    my $a=zeroes(5);
		    $a(0) .= max($y)/100;
		    $a(1) .= 0.;
		    $a(2) .= max($y)/2.;
		    $a(3) .= (max($x)-min($x))/2.+min($x);
		    $a(4) .= (max($x)-min($x))/10;
		    
		    print "a is ",$a,"\n";
		    ($ym, $a,my $covar,my $iters) =
			lmfit $x, $y, $sig, \&lmfunc, $a, {Maxiter => 2000, Eps => 1e-4};
		    
		}
		if ($ch eq "n") {
		    print "Type in the order for the polinomial (2 for straight line) \n";
		    my $nord=<>;
		    $ym=fitpoly1d $y,$nord;
		    print "1 sigma: ",&Stat::rms($ym-$y),"\n";
		}
		
		$win_spec->line ($x,$ym, {COLOR => 'GREEN'});
	    }




	    if ($ch eq $key->{stats}->[0]) {
		#statistics
#		my $i1=max(which($lbda <= $xmin));
#		my $i2=min(which($lbda >= $xmax));

		my $i1=minimum_ind(abs($lbda  - $xmin));
		my $i2=minimum_ind(abs($lbda  - $xmax));

		my $x=$lbda($i1:$i2)->sever;
		my $y=$spec($i1:$i2)->sever;
		my $min_region = min ($y);
		print "Min:  $min_region\n";
		my $max_region = max ($y);
		print "Max:  $max_region\n";
		my $npix=nelem($x);
		my $dlba = (max($x)-min($x))/($npix-1);
		my $total_flux = sum($y)*$dlba;
		print "Total flux: $total_flux in $npix pixels from $xmin to $xmax\n"; 
		my $mean=average($y);
		print "Mean: ", average($y),"\n";
		my $vari=sum(($y -$mean )**2)/float(abs($npix));
		my $rms=sqrt($vari);
		print "rms: $rms ";
		my $median=median($y);
		print "median: $median \n ";

#		print "noise: ",(&Stat::noise1D($y))[0],"\n";
#		print "done \n";

#		my $nord = 3;
#		my $w=ones(nelem($y));
#		$w->flat->index(which( ($y > $median+$rms) | ($y < $median-$rms))) .= 0;
#		my $mod = fitpoly1d($y,$nord,{Weights => $w});
#		$rms = &Stat::rms(($mod-$y)->setbadif($w == 0));
#		$nord = 4;
#		$w=ones(nelem($y));
#		$w->flat->index(which(abs($y-$mod) > 5*$rms)) .= 0;
#		$mod = fitpoly1d($y,$nord,{Weights => $w});
#		$rms = &Stat::rms(($mod-$y)->setbadif($w == 0));
#		$nord = 5;
#		$w=ones(nelem($y));
#		$w->flat->index(which(abs($y-$mod) > 4*$rms)) .= 0;
#		$mod = fitpoly1d($y,$nord,{Weights => $w});
#		$rms = &Stat::rms(($mod-$y)->setbadif($w == 0));
#		$nord = 5;
#		$w=ones(nelem($y));
#		$w->flat->index(which(abs($y-$mod) > 2*$rms)) .= 0;
#		$mod = fitpoly1d($y,$nord,{Weights => $w});
#		$rms = &Stat::rms(($mod-$y)->setbadif($w == 0));
#		$nord = 6;
#		$w=ones(nelem($y));
#		$w->flat->index(which(abs($y-$mod) > 3*$rms)) .= 0;
#		$mod = fitpoly1d($y,$nord,{Weights => $w});
#		$rms = &Stat::rms(($mod-$y)->setbadif($w == 0));
#		$nord = 6;
#		$w=ones(nelem($y));
#		$w->flat->index(which(abs($y-$mod) > 2*$rms)) .= 0;
#		$mod = fitpoly1d($y,$nord,{Weights => $w});
#		$rms = &Stat::rms(($mod-$y)->setbadif($w == 0));
#		$nord = 7;
#		$w=ones(nelem($y));
#		$w->flat->index(which(abs($y-$mod) > 2*$rms)) .= 0;
#		$mod = fitpoly1d($y,$nord,{Weights => $w});
#		$rms = &Stat::rms(($mod-$y)->setbadif($w == 0));
#		$nord = 7;
#		$w=ones(nelem($y));
#		$w->flat->index(which(abs($y-$mod) > 3*$rms)) .= 0;
#		$mod = fitpoly1d($y,$nord,{Weights => $w});
#		$rms = &Stat::rms(($mod-$y)->setbadif($w == 0));
#
#		my $badys=$y->flat->index(which(abs($mod-$y) > (5*$rms)));
#		my $badxs=$x->flat->index(which(abs($mod-$y) > (5*$rms)));
#		my $goodx = $x -> setbadif(abs($mod-$y) > (5*$rms));
#		my $goody = $y -> setbadif(abs($mod-$y) > (5*$rms));
#		my $goodmod = $mod -> setbadif(abs($mod-$y) > (5*$rms));
#		$rms = &Stat::rms($goodmod-$goody);
#
#		if ($rms > (&Stat::rms($y))) {print "SOMETHING WRONG WITH &NOISE"};
#		print "noise: ",$rms,"\n";
#		$win->line ($x,$mod, {COLOR => 'GREEN'});
#		$win->points ($badxs,$badys, {COLOR => 'RED'});
	    }
	    if ($ch eq "W") {
		my $i1=which( abs($lbda - $xmin) == min(abs($lbda-$xmin))     );
		my $i2=which( abs($lbda - $xmax) == min(abs($lbda-$xmax))     );
		my $iis = pdl [ $i1->sclr,$i2->sclr];
		$i1 = min($iis);
		$i2 = max($iis);

		my $x=$lbda($i1:$i2)->sever;
		my $y=$spec($i1:$i2)->sever;
		my $sig=$y->copy;
		(my $noise, my $mod) = &Stat::noise($y);
		$sig(:) .= $noise;
		print "noise: ",$noise,"\n";
		$win_spec->line ($x,$mod, {COLOR => 'GREEN'});
		print "writing ASCII file spec.dat \n";
		wcols $x,$y,$sig,'spec.dat';
		print "done.\n";
	    }
	    if ($ch eq $key->{toascii}->[0]) {
		my $i1=which( abs($lbda - $xmin) == min(abs($lbda-$xmin))     );
		my $i2=which( abs($lbda - $xmax) == min(abs($lbda-$xmax))     );
		my $iis = pdl [ $i1->sclr,$i2->sclr];
		$i1 = min($iis);
		$i2 = max($iis);
#		my $i1=max(which($lbda <= $xmin));
#		my $i2=min(which($lbda >= $xmax));
		my $x=$lbda($i1:$i2)->sever;
		my $y=$spec($i1:$i2)->sever;
		print "writing ASCII file spec.dat \n";
		wcols $x,$y,'spec.dat';
	    }
#	    if ($ch eq $key->{toascii}->[0]) {
#		print "inverting spec \n";
#		$spec = 1/$spec;
#		@view=($lbda,$spec,$options);
#		last;
#	    }
	    if ($ch eq $key->{help}->[0]) {
		print "########################KEYSTROKES####################################\n";

		foreach my $fn (keys %$key) {
		    printf "%-5s:  %15s  %50s \n", $key->{$fn}->[0], $fn,    $key->{$fn}->[1] ;
		}

		print "######################################################################\n";
#		print "View Keystrokes:\n a: change aspect ratio z: zoom g: set linear greyscale s: stats X: quit ^H: refresh\n  y: selects a horizontal line";
	    }
	    if ($ch eq $key->{reset}->[0]) {@view=@ori; undef $aspect; last;}
	}
#	$win->close;
    }
    return ($x,$y);
}

sub lmfunc {
   my ($x,$par,$ym,$dyda) = @_;
   my $npol=2; #degree of the fitting polinomial
   $npol=$npol-1; # to use it as index
   my $ngauss=$npol+1;
   my $ngauss2=$npol+3;
   my $indim=nelem($par);
   die if ($indim != ($ngauss2+1));
#The user supplied sub routine reference should accept 4 arguments
#    a vector of independent values $x 
#    a vector of fitting parameters 
#    a vector of dependent variables that will be assigned upon return 
#    a matrix of partial derivatives with respect to the fitting parameters that will be assigned upon return
       
   my @pol = map {$par->slice("($_)")} (0..$npol);
   my $polinome=0;
   foreach my $n (0..$npol) {
       $polinome=$polinome+$pol[$n]*($x**$n);
   }
   my @gauss = map {$par->slice("($_)")} ($ngauss..$ngauss2);
   my $gaussian=$gauss[0]*exp(-0.5*(($x-$gauss[1])**2)/($gauss[2]**2));
   $ym .= $polinome+$gaussian;
   my (@dy) = map {$dyda->slice(",($_)")} (0..$ngauss2);
   foreach my $n (0..$npol) {
       $dy[$n] .= $x**$n;
   }
   my $ngauss1=$ngauss+1;
   $dy[$ngauss] .= exp(-0.5*(($x-$gauss[1])**2)/($gauss[2]**2));
   $dy[$ngauss1] .= (($x-$gauss[1])*$gauss[1]/($gauss[2]**2))*$gaussian; # WRONG!!!
   $dy[$ngauss2] .= (($x-$gauss[1])**2/($gauss[2]**3))*$gaussian;
   
}

#sub lmfunc3 {
#   my ($x,$par,$ym,$dyda) = @_;
#   my $npol=3; #degree of the fitting polinomial
#   $npol=$npol-1; # to use it as index
#   my @pol = map {$par->slice("($_)")} (0..$npol);
#   my $polinome=0;
#   my $indim=nelem($par);
#   die if ($indim != ($npol+1));
#   foreach my $n (0..$npol) {
#	$polinome=$polinome+$pol[$n]*($x**$n);
#   }
#   $ym .= $polinome;
#   my (@dy) = map {$dyda->slice(",($_)")} (0..$npol);
#   foreach my $n (0..$npol) {
#	$dy[$n] .= $x**$n;
#   }
#}
#

sub gettransform {
    my ($hdr,$hours,$win,$args) = @_;

    my $opts={Arcsec => 0}; 
    $opts = {%$opts, %$args} if (ref($args) =~ /HASH/);


    my $scale = 1;
    $scale = 3600 if ($opts->{Arcsec});
    my $cdelt1;
    my $cdelt2;

    my $naxis2 = $hdr->{NAXIS2} ;		
    my $crpix2 = $hdr->{CRPIX2} ;
    my $crval2 = $hdr->{CRVAL2} ;
    $crval2 = $hdr->{CRVAL2} * 3600. if ($hours) ;


    my $naxis1 = $hdr->{NAXIS1} ;		
    my $crpix1 = $hdr->{CRPIX1} ;
    my $crval1 = $hdr->{CRVAL1} ;
    $crval1 = $hdr->{CRVAL1} * 3600.0 / 15  if ($hours) ;


    if (defined($$hdr{CDELT1})) {
	$cdelt2 = $scale  *  $hdr->{CDELT2} ;
	$cdelt2 = $hdr->{CDELT2}  * 3600. if ($hours) ;
	$cdelt1 = $scale  *  $hdr->{CDELT1} ;
	$cdelt1 = $hdr->{CDELT1} * 3600.0 /(cos(deg2rad($crval2/3600)) * 15) if ($hours) ;

    } elsif (defined($$hdr{CD1_1})) {
	$cdelt1 = $scale  * -1 * sqrt( ($hdr ->{CD1_1})**2 +  ($hdr ->{CD2_1})**2 );
	$cdelt2 = $scale  * sqrt( ($hdr ->{CD2_2})**2 +  ($hdr ->{CD1_2})**2 );
    }

#    my($transform) = $win->transform({ImageDimensions=>[$naxis1,$naxis2], Angle=>0 * 3.14159265358979323846264338/180, Pixinc=> [ $cdelt1, $cdelt2],RefPos =>[ [ $crpix1,$crpix2 ], [$crval1,$crval2]]});



    my($transform) = $win->transform({ImageDimensions=>[$naxis1,$naxis2], Angle=>0 * 3.14159265358979323846264338/180, Pixinc=> [ $cdelt1, $cdelt2],RefPos =>[ [ $crpix1-1,$crpix2-1 ], [0,0]]}); # FITS IS ONE-OFFSET!

    return $transform;
}



	    
sub ijad_linear {
    my ( $tr, $i, $j ) = @_;
    my ( $t1, $t2, $t3, $t4, $t5, $t6 ) = $tr -> list ;
    my $a0 = $t1 + $t2 * $i + $t3 * $j ;
    my $d0 = $t4 + $t5 * $i + $t6 * $j ;
    return ( $a0, $d0 );
}

sub adij_lin {
    my ($tr, $x, $y) = @_;
    my ($t1,$t2,$t3,$t4,$t5,$t6) = $tr -> list;

    my $i=($x-$t1-($y*$t3/$t6)+($t3*$t4/$t6))/($t2-$t5*$t3/$t6);
    my $j=($y-$t4-$i*$t5)/$t6;
    return ($i, $j)
}

sub towcs {
    my ($file, $offsets, $i, $j) = @_;
    
#		my $h = $im->gethdr;
    
    my $t;

    if ($file =~ /HASH/) {
	my $hdr_ref = $file;
#	print "CAUTION: using pdl_cd because passing hdr ref to CelCoord\n";
#	$t = &CelCoord::pdl_cd($hdr_ref);
	$t = &CelCoord::wcs($hdr_ref);
    } elsif ($file ne '') {
	$t = &CelCoord::wcs($file);
    } else {
	die "cannot run CelCoord with file $file \n";

    }
    
    my $x1 = 0; my $y1 = 0;
    
    if (defined($offsets)) {
	my @offs = @$offsets;
	foreach (@offs) {
	    $x1 += $_->[0];
	    $y1 += $_->[2];
	}
    }
    my $in = pdl($i+$x1,$j+$y1);
    
    my $out = $t->apply($in);
#    my ($a,$d) = list($out);

    return list($out);

}
		    


1;











