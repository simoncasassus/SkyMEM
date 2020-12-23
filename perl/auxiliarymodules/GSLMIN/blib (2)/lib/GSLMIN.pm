
#
# GENERATED WITH PDL::PP! Don't modify!
#
package GSLMIN;

@EXPORT_OK  = qw( PDL::PP fr_meat );
%EXPORT_TAGS = (Func=>[@EXPORT_OK]);

use PDL::Core;
use PDL::Exporter;
use DynaLoader;



   
   @ISA    = ( 'PDL::Exporter','DynaLoader' );
   push @PDL::Core::PP, __PACKAGE__;
   bootstrap GSLMIN ;







=head1 NAME
                                                                    
GSLMIN

                                                                                
=head1 DESCRIPTION
       
error codes in /usr/include/gsl/gsl_errno.h

or use const char * gsl_strerror (const int gsl_errno);
                                                               

=head1 SYNOPSIS                                      
            
   use PDL;
   use GSLMIN;
         
=head1 FUNCTIONS

=head2 fr

fletcher-reeves

=for usage

Usage:
   

=for ref

=head1 SEE ALSO


 .....

=head1 AUTHOR


=cut







=head1 FUNCTIONS



=cut





sub fr{	
 	my ($x, , $f1, $f2, $opt) = @_;	
        my $stepsize = pdl $opt ->{InitialStep};	 
        my $linmintol = pdl $opt ->{LinMinTol};	 
	my $gradtol = pdl $opt->{GradTol};

#	print "in testf, input is $x $stepsize $linmintol, $gradtol \n";
        my $out = pdl(0);
	fr_meat($x,$out, $stepsize, $linmintol, $gradtol,$f1,$f2);
#	print "back to testf: output is $out, x is $x \n";
        return $out;
}




=head2 fr_meat

=for sig

  Signature: (double xfree(n);  double [o] out(); double stepsize();double linmintol(); double gradtol(); ; SV* funcion1;SV* funcion2)


=for ref

info not available


=for bad

fr_meat does not process bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*fr_meat = \&GSLMIN::fr_meat;



;



# Exit with OK status

1;

		   