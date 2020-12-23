
#
# GENERATED WITH PDL::PP! Don't modify!
#
package WCS;

@EXPORT_OK  = qw( PDL::PP proj_meat );
%EXPORT_TAGS = (Func=>[@EXPORT_OK]);

use PDL::Core;
use PDL::Exporter;
use DynaLoader;



   
   @ISA    = ( 'PDL::Exporter','DynaLoader' );
   push @PDL::Core::PP, __PACKAGE__;
   bootstrap WCS ;





=head1 NAME

WCS                                                                    
                                                                                
=head1 DESCRIPTION

interface for libwcs, which can be found as part of the wcstools
package at http://tdc-www.harvard.edu/software/wcstools/

libwcs is the interface written by Doug Mink to access the projections
routines of Mark Calabretta, 

http://www.atnf.csiro.au/people/mcalabre/WCS/index.html

based on the latest FITS standard (Calabretta & Griesen astro-ph/0207413).


=head1 SYNOPSIS                                      
            
         
=head1 FUNCTIONS

=head2 wcs2pix

  Signature: (a(n,m); d(n,m); [o] x(n,m); [o] y(n,m))

input WCS in deg. 
wcs2pix can also take single scalar arguments.

=head2 pix2wcs

         Signature: (x(n,m); y(n,m); [o] a(n,m); [o] d(n,m))

output WCS in deg.
pix2wcs can also take single scalar arguments.

=for usage

Usage:


($x,$y) = &WCS::wcs2pix($filename,$ra,$dec);

($ra2,$dec2) = &WCS::pix2wcs($filename,$x,$y);


   



=head1 SEE ALSO

 .....

=head1 AUTHOR

simon casassus

simon@das.uchile.cl

=cut







=head1 FUNCTIONS



=cut





use PDL::Basic;
use strict;

sub mkhead {
    my $h = shift;
    my $hdr;
#    print "INPUT PDL HEADER \n";
#    print Dumper($h);
    foreach   my $key (keys %$h) { 
#       	  my $value = sprintf("%21s",$$h{$key}); # 21+8+2+2  80-33
#       	  my $value = sprintf("%67s",$$h{$key}); # 80 - 8 -2 -2 
       	  my $value = sprintf("%21s",$$h{$key}); # 
          unless ( ($key =~ /COMMENT/) || ($key =~ /AUTHOR/)|| ($key =~ /END/) || ($key =~ /HISTORY/) || ($key =~ /HIERARCH/) || ($key =~ /^\s*$/) ) {		 
#	    $hdr .= sprintf("%-8s",$key)."= $value "."/ ".sprintf("%46s");
	    my $nwhitespaces = (80 -8 - 2 -2 -1 - length($value));
	    $hdr .= sprintf("%-8s",$key)."= $value / ".sprintf("%$nwhitespaces\s");
 	  }
    }
    return $hdr;
}

sub proj_prep { 
    my ($ra,$dec,$x,$y,$direction,$filename0) = @_;
    my $filename = $filename0;
    my $dohead = 0;
    use Data::Dumper;
    if  (ref($filename) =~ /HASH/) {
        my $h = $filename;
	my $hdr = &mkhead($h);    	       
#	$hdr .= sprintf("%-3s","END");
	$hdr .= sprintf("%-s","END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ");
	$dohead = 1;
    	&proj_meat($ra,$dec,$x,$y,$direction,$dohead,$hdr);	
    } elsif ( (ref($filename) =~ /ARR/) && (ref($filename->[0]) =~ /HASH/ ) ) {
        my @hhs = @$filename;
	my $hdr;
        foreach my $h (@hhs) {
	       if (ref($h) =~ /HASH/) {
	         $hdr .= &mkhead($h);    	       
	       } 
	}
	$hdr .= sprintf("%-s","END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ");
	$dohead = 1;
	&proj_meat($ra,$dec,$x,$y,$direction,$dohead,$hdr);	
    } else {
	&proj_meat($ra,$dec,$x,$y,$direction,$dohead,$filename);	
    }
}


sub wcs2pix {
   my ($filename, $ra0, $dec0)  = @_;
   my ($ra, $dec, $x, $y);
   if (ref($ra0) =~ /PDL/) {
      $ra = $ra0;
      $dec = $dec0;
   } else {
      $ra = pdl [ [$ra0] ];
      $dec = pdl [ [$dec0] ];
   } 
   $x = zeroes($ra ->dims);
   $y = zeroes($ra ->dims);
   my $direction = -1;
   &proj_prep($ra,$dec,$x,$y,$direction,$filename);	
   if (ref($ra0) !~ /PDL/) {
      $x = $x->sclr;
      $y = $y->sclr;
   }
   return ($x, $y);
}	

sub pix2wcs {
   my ($filename, $x0, $y0)  = @_;
   my ($ra, $dec, $x, $y);
   if (ref($x0) =~ /PDL/) {
      $x = $x0;
      $y = $y0;
   } else {
      $x = pdl[ [$x0] ];
      $y = pdl[ [$y0] ];
   } 
   $ra = zeroes($x ->dims);
   $dec = zeroes($y ->dims);
   my $direction = 1;
   &proj_prep($ra,$dec,$x,$y,$direction,$filename);	
   if (ref($x0) !~ /PDL/) {
      $ra = $ra->sclr;
      $dec = $dec->sclr;
   }
   return ($ra, $dec);
}	




=head2 proj_meat

=for sig

  Signature: (ra(n,m);  dec(n,m); x(n,m); y(n,m); ; int direction; int dohead; char* filename; )


=for ref

info not available


=for bad

proj_meat does not process bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*proj_meat = \&WCS::proj_meat;



;



# Exit with OK status

1;

		   