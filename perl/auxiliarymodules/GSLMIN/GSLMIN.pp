pp_bless('GSLMIN');

pp_addhdr('
#include <math.h>

#include "FUNC.c"
//extern void conjgrad(double*,int*, double* out);
');     
pp_addpm('
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
');


pp_def('fr_meat',
        Pars => 'double xfree(n);  double [o] out(); double stepsize();double linmintol(); double gradtol(); ',
        OtherPars => 'SV* funcion1;SV* funcion2;',
        Inplace => ['out'],
        Docs => undef,
        Code =>'
int n; int i;
n = $SIZE(n);
//printf ("in fr_meat, input is \n");
//for (i=0;i<n; i++) {
//  printf (" %g ",$xfree(n=> i));
//}
printf (" \n");
ene = $SIZE(n);
ext_funname1 = $COMP(funcion1);
ext_funname2 = $COMP(funcion2);
//printf ("in fr_meat, input n is %d  \n",n);
conjgrad($P(xfree), &n, $P(out), $P(stepsize), $P(linmintol),  $P(gradtol)); 
//printf(" OUT, in fr_meat %g \n",$out());
');  


#
#pp_addpm('
#sub dummy{
#    my ($x) = @_;
#    return  trace($x);
#}
#');
#
#pp_def('trace',
#     Pars => 'a(n); [o] b();', 
#     Code => 'int i; n=SIZE($n); 
#	 for(i=0; i<n; i++) { 
#			$b  = $b + $a(n => i);
#		}
#	');
#

pp_addpm({At=>Top},<<'EOD');       



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

EOD

pp_done();
