use PDL;
use PDL::NiceSlice;
use GSLMIN;
#use NRMIN;
use Vtools;
use PDL::ImageND;
use PDL::Image2D;
use PDL::FFT;
use Data::Dumper;
use Astro::Time;
use Astro::Coord;
use CelCoord;
#use PDL::FFT;
use strict;
no strict 'refs';

#perl multimem.pl out_clip_39157_fix3/  39157 dochi2 
#perl multimem.pl out_clip_33157_IRAC_fix3/  33157 dochi2  _IRAC 
#perl multimem.pl out_mem7_IRAC_H169_20160_l0_fix2/  20160 _IRAC _H169
#perl multimem.pl out_clip_8800_fix3/ 8800 dochi2 

my $prior; my $data; my $kernel;  # SIGMA IS RMS NOISE 
my $omega; my $fprior; my @residuals_global;
my $globaldebug = 0;
my $viewinit = 0;
my $viewout = 0;

my $noise_mosaic;
my $PBmosaic;
my $datapass;
my @bestpass;

my $dochi2 = 0;
my $doclip = 1;
my $entropy_switch = 0;  # START ENTROPY 
my $entropy_switch0 = 0;  # TARGET ENTROPY - REVERT TO THIS WHEN ITER > ITER_THRESHOLD
my $target_lambda = 0;
my $lambda = 0;


my $freq = '17481';

#my $outdir = './out_mem7_ATCA_H75_39157_l0_noPS/';
my $outdir = './out_mem7_H75_17481_l0_fix2/';
my $simultag = '';
#$simultag = '_IRAC';
#$simultag = '_SPIKES';

my $array = '';
if (@ARGV) {
    $outdir = shift @ARGV;
    $freq = shift @ARGV;
# perl multimem.pl out_mem7_H75_17481_l0_fix2  17481 
    if (@ARGV) {
	my $dum =  shift @ARGV;
	if ($dum =~ /dochi2/) {
	    $dochi2 = 1;
	    $dum = shift @ARGV;
	} 
	if ($dum =~ /H/) {
	    $array = $dum;
	} else {
	    $simultag = $dum;
	}
	if (@ARGV) {
	    my $array =  shift @ARGV;
	}
# perl multimem.pl out_mem7_H75_17481_l0_fix2  17481 IRAC H169
    }
}

my $skipoptim = 0;  # WATCH OUT FOR DELETION OF $outdir 

# Exact WCS Centroid: 16:25:56.26  -24:20:50.4  from gaussit no DC offset
my @pointsources; # = ({ Pos => [ '16:25:56.26', '-24:20:50.4']    , Flux => 0.000486344157484239} );
#if ($freq == '39157') {
#    print "ADDING POINT SOURCES \n";
#    @pointsources =  ({ Pos => [ '16:25:56.26', '-24:20:50.4']    , Flux => 0.000486344157484239} ) ;
#}
#

#if ( ($freq eq '5500')  && ($simultag eq '') && (!$dochi2) )  {
#    @pointsources = (
#	{ Pos => [ '16:26:22.71',  '-24:22:55.4']    , Flux => 6.56e-03} ,
#	{ Pos => [ '16:26:02.95',  '-24:23:36.8']    , Flux => 3.37e-03} 
#	) ;
#}
##Max in box is: 6.56e-03 
##Max in box is: 3.37e-03 
#
#if ( ($freq eq '8800')  &&  ($simultag eq '')  &&  (!$dochi2) )  {   # 
#    @pointsources = ({ Pos => [ '16:26:19.79', '-24:26:56.6']    , Flux => 0.00120} ) ;
#}
#

#Exact WCS Centroid: 16:26:19.79  -24:26:56.6 
#Linear WCS Centroid 0.0864453644461169 -0.101859805933138
#Bmax/2 0.00612013611717531; Bmin/2 0.00289595329081443;
#PA-4.62673851583121deg (South of East);
#peak intensity 0.00120172277617356 
#integrated flux density in Jy (should be equal to intensity for correct `beam' value:
#0.00188312682831355



unless ($dochi2) {
    $entropy_switch = 7;  # START ENTROPY 
    $entropy_switch0 = 7;  # TARGET ENTROPY - REVERT TO THIS WHEN ITER > ITER_THRESHOLD
    $target_lambda = 0; ###OK for s=6
    $lambda = 0; ###OK for s=5
#    $target_lambda = 1; ###OK for s=6
#    $lambda = 1; ###OK for s=5
}


# no further user input below this line

my $pi = 2*acos(0);

my $iterdf =0; my $norm;
my $N; my $Nk; my $clips = 0; 
my $true; my $h;
#my $boundary = 'extend'; # boundary condition for convolutions
#my $boundary = 'extend'; # boundary condition for convolutions
my $boundary = 'periodic'; # boundary condition for convolutions
$boundary = 'truncate';
my $progressive_clip = 0; 
my $minpix; my $target_minpix;
my $aux = 'aux_exp' ;
my $daux = 'aux_exp' ;
my $aux_inv = 'aux_inv_exp';

$aux = 'identity' ;
$daux = 'didentity' ;
$aux_inv = 'identity';


my $iter_threshold = -1;


system("\\rm -r $outdir; mkdir $outdir\n") unless $skipoptim;

system("rsync -va multimem.pl  $outdir/multimem_record.pl\n") unless $skipoptim;

my $reflabel = '1_2';

my $fieldref =     {image => "data_128$array/oph_$reflabel\.$freq$simultag\.fits", psf => "data_128$array/oph_$reflabel\.$freq$simultag\.beam.fits", 
		    sens =>  "data_128$array/oph_$reflabel\.$freq\.sensitivity.fits", rstor => "data_128$array/oph_$reflabel\.$freq$simultag\.rstor.fits"};

if ($freq < 1E4) {
    $reflabel = '0_0';
    $fieldref =     {image => "data_256$array/oph_$reflabel\.$freq$simultag\.fits", psf => "data_256$array/oph_$reflabel\.$freq$simultag\.beam.fits", 
		     sens =>  "data_256$array/oph_$reflabel\.$freq\.sensitivity.fits", rstor => "data_256$array/oph_$reflabel\.$freq$simultag\.rstor.fits"};
} 

if ($array ne '') {
   $fieldref =     {image => "data_256$array/oph_$reflabel\.$freq$simultag\.fits", psf => "data_256$array/oph_$reflabel\.$freq$simultag\.beam.fits", 
		     sens =>  "data_256$array/oph_$reflabel\.$freq\.sensitivity.fits", rstor => "data_256$array/oph_$reflabel\.$freq$simultag\.rstor.fits"};
}

#my $fieldref =     {image => "data_128$array/oph_1_2.$freq$simultag\.fits", psf => "data_128$array/oph_1_2.$freq$simultag\.beam.fits", 
#		    sens =>  "data_128$array/oph_1_2.$freq\.sensitivity.fits", rstor => "data_128$array/oph_1_2.$freq$simultag\.rstor.fits"};

# 128 is enough for model image, but need to resample 256 images so as not to loose the edges

#my $fieldref =     {image => 'data/oph_1_2.17481.fits', psf => 'data/oph_1_2.17481.beam.fits', 
#		    sens =>  'data/oph_1_2.17481.sensitivity.fits', rstor => 'data/oph_1_2.17481.rstor.fits'};


#linmos in=oph_1_1.17481.map,oph_1_2.17481.map,oph_1_3.17481.map,oph_2_1.17481.map,oph_2_2.17481.map,oph_2_3.17481.map out=oph.17481.sensitivity.mos options=sensitivity
#my @datapass = (
#    {image => 'data/oph_1_1.17481.fits', psf => 'data/oph_1_1.17481.beam.fits', 
#     sens =>  'data/oph_1_1.17481.sensitivity.fits', rstor => 'data/oph_1_1.17481.rstor.fits'},
#    );
#

my @labels = ('1_1',  '1_2', '1_3',  '2_1', '2_2', '2_3');


my @datapass = map {    
    {image => "data_256$array/oph_$_\.$freq$simultag\.fits", psf => "data_128$array/oph_$_\.$freq$simultag\.beam.fits", 
     sens =>  "data_256$array/oph_$_\.$freq\.sensitivity.fits", rstor => "data_256$array/oph_$_\.$freq$simultag\.rstor.fits"}
} @labels; 

if ($freq < 1E4) {
    @labels = ('0_0');
	@datapass = map {    
	    {image => "data_256$array/oph_$_\.$freq$simultag\.fits", psf => "data_256$array/oph_$_\.$freq$simultag\.beam.fits", 
	     sens =>  "data_256$array/oph_$_\.$freq\.sensitivity.fits", rstor => "data_256$array/oph_$_\.$freq$simultag\.rstor.fits"}
    } @labels; 
    
};


my $tref;

#die Dumper(@datapass);
my $pixsolidangle_inbeams;
{
    my @PBs;
    my @sigmas;
    my $referencemosaic = rfits($fieldref->{image});
    my $href = $referencemosaic->gethdr;


    my $hrstor = rfits( $fieldref->{rstor},{Data => 0});
    $tref = &CelCoord::wcs($fieldref->{rstor});

    my $dum = 0;
#    my $dumjunk = 0;
    my $norm = 0;

    foreach my $field (@datapass) {
	my $datafile = $field -> {image};
	$data = rfits($datafile); #
	$h = $data->gethdr;
	$data = $data(:,:,(0),(0));
	my $kernel0 = rfits($field->{psf}); #produced by resamp.pl
	my @kdims = $kernel0 -> dims;
#    $norm = max($data);
	$norm = 1;
	$data /= $norm;


	$data->sethdr($h);
#	print "data before matching: \n";
#	&Vtools::view($data); 
#	$PB = &CelCoord::match($fieldref->{image}, $kernel);
	$data = &CelCoord::match($href, $data);
#	print "data after matching: \n";
#	&Vtools::view($data); 






	my $cdelt = $$h{CDELT2};
	

#	&Vtools::view($kernel);
    
	my $sens0 = rfits($field->{sens});
	my $hsens = $sens0->gethdr;
	print "building PB \n";
	print "sens0 info ",$sens0 -> info,"\n";
	my $sens = $sens0(:,:,(0),(0));
	print "sens info ",$sens -> info,"\n";
	$sens =  $sens->setnantobad;
	

	my $sigma = min($sens -> where($sens > 0) ); #1E-4; # noise in diry map 
	
	push @sigmas, $sigma;

	$field->{sigma} = $sigma;

	$sens -> where( isbad($sens) ) .= -1;
	$sens -> where($sens <= 0) .= -1.;
	my $PB = double(1/$sens);
	$PB -> where($PB <= 0) .= 0;

	$PB->sethdr($h);
#	print "PB before matching: \n";
#	&Vtools::view($PB); 
#	$PB = &CelCoord::match($fieldref->{image}, $kernel);
	my $PB_ori = $PB->sever;
	$PB = &CelCoord::match($href, $PB);
#	print "PB after matching: \n";
#	&Vtools::view($PB); 
	$PB /= max($PB);


#    $kernel  = $kernel0(:,:,(0),(0));
#    my $boxk =  254;
	my $boxk =  128; #  $data -> getdim(0); # assumes $data is square, and its size is half that of the beam image
#	$boxk = 128; #trying small kernel

	my $x1 = int($kdims[0]/2-$boxk/2+1); # not sure about the +1
	my $x2 = int($kdims[0]/2+$boxk/2);
	my $y1 = int($kdims[1]/2-$boxk/2+1);
	my $y2 = int($kdims[1]/2+$boxk/2);
	
#    my $x1 = int(255-$boxk/2);
#    my $x2 = int(255+$boxk/2);
#    my $y1 = int(255-$boxk/2);
#    my $y2 = int(255+$boxk/2);
#    
#
#    my $x1 = int(1023-$boxk/2);
#    my $x2 = int(1023+$boxk/2);
#    my $y1 = int(1023-$boxk/2);
#    my $y2 = int(1023+$boxk/2);
#    
# beam centred on 254.999987688296 255.000015071733;    
	
#	my $kernel1  = $kernel0($x1:$x2,
#				$y1:$y2,(0));    

	my $kernel1  = $kernel0(:,:,(0));    


#	$kernel1 = $kernel1 * $PB_ori;



#	my $kernel1  = $kernel0(:,:,(0));    

#    my $kernel1  = $kernel0(:,:,(0));
	
 	
#	&Vtools::view($kernel1);
#    $kernel = kernctr($data,$kernel1);
	$kernel = $kernel1 -> copy;
#	print "kernel norm ", sum($kernel),"\n";

	my $beam_deg2 = ($pi/(4*log(2)))*($$hrstor{BMAJ})*($$hrstor{BMIN});
	
	$pixsolidangle_inbeams = ($$hrstor{CDELT2}**2/$beam_deg2);

#	my $trial_intensity = zeroes($data->dims);
#	$trial_intensity = exp(-0.5*($trial_intensity->rvals)**2/5.**2);
#	$trial_intensity($data->getdim(0)/2,$data->getdim(1)/2)  .= 1.0 / ($$href{CDELT2}**2 / $beam_deg2);  
#	print "trial input \n";
#	&Vtools::view($trial_intensity) if $viewinit;
#	print "kernel \n";
#	&Vtools::view($kernel) if $viewinit;
#	my $test = convolveND( $trial_intensity, $kernel, {Boundary => 'truncate'});
#	print "trial output \n";
#	&Vtools::view($test) if $viewinit;
#	print "normfactor : ",sum($trial_intensity),"\n";
##	$kernel /= sum($test)/sum($trial_intensity);
#	$kernel /= sum($trial_intensity);

	my $knorm = ($beam_deg2/$$href{CDELT2}**2);  
	$kernel /= $knorm;

#	my $test = convolveND( $trial_intensity, $kernel, {Boundary => 'truncate'});
#	print "trial output \n";
#	&Vtools::view($test) if $viewinit;
#	my $testrstored = &rstor($trial_intensity,zeroes($test->dims), $hrstor); # need to pass header for Miriad restored file, which contains BMAJ and BMIN.
#	print "trial output \n";
#	&Vtools::view($testrstored) if $viewinit;
#	die;
#Try with theoretical assess of 1 spike at phase centre. -> rstor should give peak Jy/beam of 1.





	$sens->sethdr($h);
	my $sensmos = &CelCoord::match($href, $sens);
	
#    $PB = exp(-0.5 * (ones(128,128)->rvals**2)/10**2);
	
	
	print "data \n";
	&Vtools::view($data) if $viewinit;
	print "primary beam\n";
	&Vtools::view($PB) if $viewinit;
	print "kernel \n";
	print $kernel -> info,"\n";
	&Vtools::view($kernel) if $viewinit;

	$field->{PB} = $PB;
	$field->{kernel} = $kernel;
	$field->{data} = $data;
	
	push @PBs, $PB;
	

	my $dum1 = 1/$sensmos**2;
	$dum1 -> where($dum1 <= 0) .= 0.;

	$dum += $dum1;
#	$dumjunk += sqrt($dum1);

#	$norm += $PBs[$ifield]**2;
    } 	

    $noise_mosaic = sqrt(1/$dum);
#    my $noise_mosaicjunk = 1/$dumjunk;

    $noise_mosaic -> sethdr($href);
    $noise_mosaic -> where($noise_mosaic  > 1E-1) .= -1;
    wfits $noise_mosaic,$outdir.'noise_mosaic.fits';
    $PBmosaic = 1/$noise_mosaic;
    $PBmosaic->where($noise_mosaic < 0) .= 0;
    $PBmosaic /= max($PBmosaic);

#    wfits $noise_mosaicjunk,'noise_mosaicjunk.fits';
#    print "mosaic noise \n";
#    &Vtools::view($noise_mosaic);
#    die;

#    $PB = 1;
	
    
    $omega = nelem($referencemosaic);
    
    $iter_threshold = -1;
    
    $target_minpix = 0;

    if ($entropy_switch0 > 0) {
	my $sens = rfits($fieldref->{sens});
	my $sigma = min($sens -> where($sens > 0) ); #1E-4; # noise in diry map 
	$target_minpix =  $sigma/1E2;
	$minpix =  $target_minpix;
	if ($progressive_clip) {
	    $minpix =  median($referencemosaic); #$sigma/1E4;
	}    
    }
    
    my $xfree = zeroes($data->dims) -> flat  + $minpix;

    if ($entropy_switch0 == 7) {
	
#	this is IRAC simu - tied to a fixed input prior. Should not require  cross-correlations.
#>> DELETE	  fIxed: must use prior mosaiced from 
#>> DELETEROPH_IRAC8um_ascale_Jy_pixel_atten_oph_2_3.33157
#>> DELETEfor MEM maybe switch to a regularizing term per field - so as to keep track of PB attenuation. 
#>> DELETEfor chi2: model is compared to each field after PB attenuation. 
#
#prior is pixelized according to $filesens, not fileinp!!! so resamp.pl OK in the end. 
#rawi_resamp = &CelCoord::match($filesens, $fileinp); <<<< from simul.pl
#

	my $file_raw = 'ROPH_IRAC8um_ascale_Jy_pixel.fits';
	my $hraw = rfits($file_raw,{Data=>0});
#	my $pix2 = ($$hrstor{CDELT2}**2); # WATCH OUT - hrstor for mistaken simul.pl
	my $pix2 = ($$hraw{CDELT2}**2); # WATCH OUT - hrstor for mistaken simul.pl

	$prior = &CelCoord::match($fieldref->{rstor}, $file_raw);
	print "prior \n";
#	&Vtools::view($prior);


	my $hrstor  = $prior->gethdr;
	my $bmaj = $$hrstor{BMAJ};
	my $bmin = $$hrstor{BMIN};
	my $beam = ($pi / (4 * log(2)) ) * $bmaj  *  $bmin; #beam in deg
#	$prior *= (6.46e-03/ 1.18e-03) * $beam  / $pix2;  fix for 33157 mistake
	$prior *=  $beam  / $pix2;

######################################################################
#     xcorr with IRAC

	my @alphas;
	open ALPHA, ">$outdir"."alpha.txt";
	foreach my $ifield (0..$#labels) {
	    my $fieldref = $datapass[$ifield];
	    my $ATCA = rfits($fieldref->{image});
	    my $tagIRAC = '_IRAC';
	    my $label = $labels[$ifield];
	    
	    my $fieldrefIRAC =     {image => "data_256/oph_$label.$freq$tagIRAC\.fits", psf => "data_256/oph_$label.$freq$tagIRAC\.beam.fits", 
				sens =>  "data_256/oph_$label.$freq\.sensitivity.fits", rstor => "data_256/oph_$label.$freq$tagIRAC\.rstor.fits"};
	    my $IRAC = rfits($fieldrefIRAC ->{image});
	    my $weights0 = 1/rfits($fieldref ->{sens});
	    $weights0 /= max($weights0);
	    my $weights = $weights0 -> flat;
	    print  "ANCHOR \n";
	    print $weights -> info,"\n";
	    my $vIRAC = $IRAC->flat;
	    my $vATCA = $ATCA->flat;
	    print $vATCA->info,"\n";
	    print $vIRAC->info,"\n";

	    my $alpha = sum($vIRAC * $vATCA * $weights) / sum( ($vIRAC)**2 * $weights);     # radio = alpha IR

	    push @alphas, $alpha;
	    print ALPHA "xcorr $label  $alpha \n";
	}
	my $alphas =  pdl @alphas;
#	my $alpha = average($alphas);
	my ($mean,$prms,$median,$min,$max,$adev,$rms) = statsover($alphas);
	my $alpha = $mean;
	my $alpharms = $rms; 
	print ALPHA "alpha = $alpha +- $alpharms \n";
	close ALPHA;
	$prior *= $alpha unless ($simultag ne '');

#	print "prior -> info,",$prior->info,"\n";
	$prior =  $prior ->setnantobad;
	$prior -> where( isbad($prior) ) .= -1;
	$prior -> where($prior <= 0) .= -1.;
#	$prior -> where($prior <= $minpix) .= $minpix;   DEV !!! 

	
	my $fileout_fr = $outdir.'mod_in_prior.fits';
	$prior -> sethdr($hrstor);
	wfits $prior,$fileout_fr;

	print "initial conditions \n";
#	&Vtools::view($prior);
	$fprior = $prior -> flat;
	$xfree(:) .= $fprior(:); # ->copy;


    }

    $xfree = &$aux_inv($xfree);


    my $out;
    if ($skipoptim) {
	$out = rfits($outdir.'mod_out.fits');

	foreach my $ires (0..$#labels) {
	    my $label = $labels[$ires];
	    my $res = rfits($outdir."res_$ires\.fits");
	    push @residuals_global, $res;
	}

#	$residualglobal = rfits('res_out.fits');

    } else {
	print "START MINIMISATION \n";
	
	
#	&NRMIN::fr($xfree,\&f,\&df,{ftol => 1E-6}),"\n";
#	print "----------------------------------------------------------------------\n";
#	&f($xfree);
	
	my $stepsize = 1.; #1E-6 ; # 1.
	my $linmintol = 1E-2;#  1E-2
	my $gradtol = 1E-2;  
	print "init max: ",max($xfree)," \n";
	my $gslmin_retval = &GSLMIN::fr($xfree,\&f,\&df,{InitialStep=> $stepsize, LinMinTol => $linmintol, GradTol => $gradtol}),"\n";
	print " gslmin_retval $gslmin_retval \n";


	print "output max: ",max($xfree)," \n";
	$globaldebug = 0;
	print "----------------------------------------------------------------------\n";
	print " entropy_switch  ",$entropy_switch," \n";
	my $bestL = &f($xfree); # update $global_residual
	open BST, ">$outdir"."best.txt";
	print BST "MINUMUM  L:  $bestL  Chi2: $bestpass[1]   redChi2: ",$bestpass[1]/$bestpass[2]," dofs:  $bestpass[2]    \n";
	close BST;

	
	$out = restorematrix($xfree);
	$out = &$aux($out);
	$out *= $norm;
	$out->sethdr($h);

	print "DONE, clips $clips \n";
	
	wfits $out,$outdir.'mod_out.fits';

	foreach my $ires (0..$#residuals_global) {
	    my $res = $residuals_global[$ires];
	    my $label = $labels[$ires];
	    $res->sethdr($h);
	    wfits $res,$outdir."res_$ires\.fits";
	}

	
    }
    

# MOSAIC residuals;

    my $residual_mosaic = &mos(\@residuals_global, \@PBs, \@sigmas, $href);

    wfits $residual_mosaic,$outdir.'residual_mosaic.fits';

    my $rstored = &rstor($out,$residual_mosaic, $hrstor); # need to pass header for Miriad restored file, which contains BMAJ and BMIN.
    wfits $rstored,$outdir.'rstor_out.fits';
    
    my $rstored_PB = $rstored * $PBmosaic;
    $rstored_PB->sethdr($href);
    wfits $rstored_PB, $outdir.'rstor_out_PBmult.fits';
    
	
}

sub mos {
    my ($reffields, $refPBs, my $refsigmas, my $href ) = @_;

    my @fields = @$reffields;
    my @PBs = @$refPBs;
    my @sigmas = @$refsigmas;
    my $dum = 0;
    my $norm = 0;
    foreach my $ifield (0..$#fields) {
	$dum += $fields[$ifield] * $PBs[$ifield] / ($sigmas[$ifield]**2);
	$norm += $PBs[$ifield]**2 / ($sigmas[$ifield]**2);
    }
    my $mos = $dum / $norm;
    $mos -> sethdr($href);
    return  $mos;
}


sub rstor {
    my ($out, $residuals, $hrstor) = @_;
    
    use GSL_Cholesky;
    use Gauss;

    my $bmaj; my $bmin; my $alpha; my $bpa;
    my $dx;

    my $nkernel=56;
    my $kernel = zeroes($nkernel,$nkernel);
    my $mu  = pdl [$nkernel/2,$nkernel/2];

# mgauss best fit
#Centroid 14.996921325627 12.9844326916882; Bmax/2 4.25505441809659; Bmin/2 3.44810439708965;
#PA -0.181242802032051; peak intensity 1.15788695013428 
    $dx = $$hrstor{CDELT2};
    $bmaj = $$hrstor{BMAJ}/(2*$dx);
    $bmin = $$hrstor{BMIN}/(2*$dx);

    $alpha = $$hrstor{BPA};

    print "BMAJ/2 $bmaj pixels BMIN/2 $bmin pixels PA $alpha  \n";

    $alpha=deg2rad($alpha-90);
#$alpha=deg2rad($alpha);
    print "alpha $alpha or ",rad2deg($alpha),"\n";
    print "pi + alpha ",$pi+$alpha," \n";
    
    my $rot = pdl [[cos($alpha),sin($alpha)],[-sin($alpha),cos($alpha)]];
    my $dm= stretcher(pdl [($bmaj/(sqrt(2*log(2))))**2,($bmin/(sqrt(2*log(2))))**2]);
    my $C =  transpose($rot)  x $dm x ($rot) ;
    my $Cinv = inv($C);
    my $L = &GSL_Cholesky::Choldc($Cinv);
    my $cis =  pdl [$L(0,0)->sclr, $L(1,1)->sclr, $L(0,1)->sclr];
    my $flux=1;
    my $params = {Centroid => $mu, Cinv => $cis, Flux => $flux};
    my $kernel = &Gauss::mgauss($kernel,$params);


    $kernel = $kernel/sum($kernel);  


    print " kernel \n";
#    &Vtools::view($kernel);

    print "model \n";
#    &Vtools::view($out);

    
    my $beam = ($pi/(4*log(2)))*($$hrstor{BMAJ}*$pi/180.)*($$hrstor{BMIN}*$pi/180);
    print "beam solid angle is $beam sr \n";
#    $im1 = conv2d $im1, $kernel, {Boundary => 'Reflect'}; 
    print $out->info,"\n";

    print "model \n";
    &Vtools::view($out) if $viewout;

    print "gaussian beam \n";
    &Vtools::view($kernel) if $viewout;

    my $rstor = convolveND( $out, $kernel, {Boundary => 'truncate'});

    print "convolved model \n";
    &Vtools::view($rstor)  if $viewout;

    print "residuals \n";
    &Vtools::view($residuals)  if $viewout;

    print "conv model with residuals \n";
    $rstor += $residuals;
    &Vtools::view($rstor)  if $viewout;


    $rstor -> sethdr($hrstor);

    
    return $rstor;
    

}

sub restorematrix {
### global: $data
    my $f = shift;
    my ($N, $M) = $data->dims;
    my $i = sequence($M)*$N; 
    $i = $i->dummy(0);  # last fix DEV DEV simon jun 2019
    my $b=  $f->range($i,$N)->transpose;
    return $b;
}

sub aux_inv_exp {
    my $x = shift;
    $x -> where($x <= $minpix) .= $minpix;
    my $I = log($x);
    return $I;
}

sub aux_exp {
    my $x = shift;
    my $I = exp($x);
    return $I;
}

sub identity {
    my $x = shift;
    my $I = $x;
    return $I;
}

sub didentity {
    my $x = shift;
    my $dI = ones($x->dims);
    return $dI;
}

sub daux_exp {
    my $x = shift;
    my $dI = exp($x);
    return $dI;
}


#sub conv {
#    my  ($im, $k, $opts) = @_;
#
#    
#
##print "DEGRADING \n";
##my $im_1 = $im_0->copy;
##my $dum = $Br->copy;
##fftconvolve($im_1,$dum);
##$im_1 = &reord($im_1, $N);
#
#
#
#}


sub f {
### global: $kernel, $data, $residual_global
    my $x = shift;
#    print "in perl f : x(0) is ",$x(0),"\n";
#    my $b = restorematrix($x);
#    &Vtools::view($b);
    if ( ($doclip) && ($entropy_switch >= 0) && ($aux eq 'identity') ) {
	$x->where($x <= $minpix) .= $minpix;
	$clips++;
    } 
    my $I = &$aux($x);
    

    my $b = restorematrix($I);
#    print "chi2 : input model \n";
#    print $b->info,"\n";
#    &Vtools::view($b);
#    print $PB->info,"\n";
#    print "PB\n";

#    print "real beam \n";

#    &Vtools::view($PB);

#    print "fudge beam \n";
#    $PB = exp(-0.5 * ($PB->rvals**2)/10**2);

    undef @residuals_global;
    my $chi2 = 0;
    foreach my $field (@datapass) {
	my $bfield = $b->copy; # copy model pdl
	if (@pointsources) {
	    foreach (@pointsources) {
		my $pixpos = &CelCoord::pos($tref,$_->{Pos}->[0],$_->{Pos}->[1]);
		print "ADDING point source at pixel $pixpos \n";
#		&Vtools::view($bfield);
		$bfield(list $pixpos) .= $_->{Flux} / $pixsolidangle_inbeams;
#		&Vtools::view($bfield);
	    }
	}

	my $PB = $field->{PB};
	my $kernel = $field->{kernel};
	my $data = $field->{data};

	$bfield *= $PB;

	my $model = convolveND $bfield, $kernel, {Boundary => $boundary};



#	&Vtools::view($bfield);
#	&Vtools::view($model);
	
	my $residual = $data - $model;
	push @residuals_global, $residual;
	
	my $sigma = $field->{sigma};
	$chi2 += sum(($model-$data)**2/$sigma**2);
#    my $norm  = sum($model);

    }
    my $mem = 0;
    if ($entropy_switch == 5) {
	$b->where($b <=  $minpix) .= $minpix;
	$mem = sum(log($b));
    } elsif ($entropy_switch == 6) {
#	print "ENTROPY 6 MODEL  \n";
#	&Vtools::view($b);
	$b->where($b <=  $minpix) .= $minpix;
	$mem = sum($b*log($b/$minpix));
    } elsif ($entropy_switch == 7) {
	$b->where($b <=  $minpix) .= $minpix;
	if ($globaldebug) {
	    print "ENTROPY 7 MODEL  \n";
	    &Vtools::view($b);
	    print "prior info ",$prior->info,"\n";
	    &Vtools::view($prior);
	}
	$mem = sum($b*log($b/$prior));
	print "max b ",max($b)," minpix $minpix \n";
	print "mem = $mem \n";
    }
    my $L = $chi2-$lambda*$mem; 
    print "perl side red chi2 = ",$chi2/$omega,"  mem = $mem  L=$L normfree = ",sqrt(sum($x**2)),"  \n";
    print "  ITER $iterdf, lambda $lambda entropy_switch $entropy_switch  \n";
#    if ( abs($chi2/$omega - 1.0) < 1.1) {
#	wfits $b, "MEM_$iterdf.fits";
#    }
    
    @bestpass = ($L, $chi2, $omega);
#    print "perl side chi2 = ",$chi2/$omega,"  omega=$omega mem = $mem \n";
    return $L;
}



sub df {
### global: $kernel, $data
    my $x = shift;


    if ( ($doclip) && ($entropy_switch >= 0) && ($aux eq 'identity') ) {
	$x->where($x <= $minpix) .= $minpix;
    } 
    my $I = &$aux($x);


    my $b = restorematrix($I);

    my $chi2 = 0;
    my $deriv;

    my $ifield = 0;
    foreach my $field (@datapass) {
	my $bfield = $b->copy; # copy model pdl
	if (@pointsources) {
	    foreach (@pointsources) {
		my $pixpos = &CelCoord::pos($tref,$_->{Pos}->[0],$_->{Pos}->[1]);
		$bfield(list $pixpos) .= $_->{Flux} / $pixsolidangle_inbeams;
	    }
	}
	my $PB = $field->{PB};
	$bfield *= $PB;
	my $kernel = $field->{kernel};
	my $data = $field->{data};
	my $sigma = $field->{sigma};
	my $model = convolveND $bfield, $kernel, {Boundary => $boundary};



	my $residual = $data - $model;
	wfits $residual, $outdir."RES$ifield\_$iterdf.fits";
	$ifield++;

	my $A =  2 * ($model-$data) / $sigma**2 ;
	
#    print "in df, viewing deriv A \n";
#    &Vtools::view($A);
#    print "in df, viewing kernel  \n";
#    &Vtools::view($kernel);
	
	my $deriv2D = convolveND $A,$kernel, {Boundary => $boundary};
	
	$deriv2D *= $PB;
#    print "deriv2D \n";
#    &Vtools::view($deriv2D); 
	
	$deriv += $deriv2D->flat;
	
    }

    my $memderiv = 0;
    print "entropy_switch $entropy_switch \n";

    if ($entropy_switch == 5) {
	$I->where($I <=  $minpix) .= $minpix;
	$memderiv = 1/$I;
    } elsif ($entropy_switch == 6) {
	$I->where($I <=  $minpix) .= $minpix;
	$memderiv = log($I/$minpix) + 1.;
    } elsif ($entropy_switch == 7) {
	$I->where($I <=  $minpix) .= $minpix;
	$memderiv = log($I/$fprior) + 1.;
    }

    $memderiv *= &$daux($x);

    my $retvect = $deriv - $lambda* $memderiv; 
    print "in df, lambda $lambda normgrad ",sqrt(sum($retvect**2))," \n";
    print "in df, input model at iter $iterdf \n";
#    &Vtools::view($deriv2D);
    
    
    my $save = 0;
    if ($save) {

	my $grad2D = restorematrix($retvect);
    	wfits $grad2D,$outdir.'grad_out.fits';
    }

#    &Vtools::view($b);
    if ($iterdf > $iter_threshold) {
	$entropy_switch = $entropy_switch0;
	$minpix /= 10. unless ($minpix <= $target_minpix);
	$lambda *= 10. unless ($lambda >= $target_lambda);
	$lambda = $target_lambda if  ($lambda > $target_lambda);
    }
    wfits $b, $outdir."MEM_$iterdf.fits";

    
    my $PBmosaic_flat = $PBmosaic->flat;
    $retvect->where($PBmosaic_flat <= 0. ) .= 0;

    $iterdf++;
    return $retvect;
}



