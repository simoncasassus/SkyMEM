# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

#BEGIN { $| = 1; print "1..1\n"; }
#END {print "not ok 1\n" unless $loaded;}
#use PDL;
#use Vor;
#$loaded = 1;
#print "ok 1\n";
#

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

use blib;
use PDL;
use PDL::NiceSlice;
use GSLMIN;
use Vtools;
use strict;

print "TESTING \n";

my $x=random(2);
my $offset = sequence(nelem($x));

#print "DUMMY: ",&PDL::Dev::dummy($x),"\n";


#print "OUTPUT f: ",&GSLMIN::fr($x,\&f,\&df),"\n";

my $stepsize = 1;
my $linmintol = 1E-2;
my $gradtol = 1E-9;
my $xfree = random(2);
print "init xfree: $xfree \n";
print "OUTPUT f: ",&GSLMIN::fr($xfree,\&f,\&df,{InitialStep=> $stepsize, LinMinTol => $linmintol, GradTol => $gradtol}),"\n";

print "output xfree: $xfree \n";
print "DONE \n";

sub f {
    my $y = shift;
#    print $y->info,"\n";
#    print "offset : $offset \n";
    my $output = sum(($y-$offset)**2);
#    print "IN PERLFUNC f, input is $y \n";
#    print "IN PERLFUNC f, output is $output \n";
    return $output;
}



sub df {
    my $y = shift;
#    print $y->info,"\n";
    my $output = sum($y**2);
#    print "IN PERLFUNC df, input is $y \n";
    my $deriv = 2*($y-$offset);
#    foreach my $i (0..nelem($y)-1) {
#	$y($i) .= $deriv($i);
#    }
#    print "IN PERLFUNC df, output deriv  is $deriv \n";
    return $deriv;
}










