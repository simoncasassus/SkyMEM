use ExtUtils::MakeMaker;
use PDL::Core::Dev;

@pack = (["GSLMIN.pp",GSLMIN,GSLMIN]);

%hash = pdlpp_stdargs(@::pack);
#use Data::Dumper;
#print Dumper(\%hash);

$hash{'LIBS'} = [' -L/usr/local/lib   -lgsl -lgslcblas -lm'];

$hash{INC} .= " -I/usr/local/include/";
#$hash{INC} .= " -I/Developer/SDKs/MacOSX10.6.sdk/usr/include/sys/";
#$hash{CCFLAGS} .= " -g";

# $hash{'OPTIMIZE'}  = '-g'; # If you want to debug, uncomment this.
# $hash{INC} .= " -I/usr/local/include";  # uncomment as required
# $hash{LIBS}[0] .= " -L/usr/local/lib -lmylib "; # uncomment as required
#$hash{LIBS}[0] .= " -L/opt/SUNWspro/WS6U2/lib -lF77 -lM77 -lsunmath -lm";
#$hash{MYEXTLIB} .= " libminuit.a ";
#$hash{MYEXTLIB} .= " ccode.o ";
#$hash{PREREQ_PM} .= ('ExtUtils::F77' => 1.13);
#$hash{DEFINE} .= ExtUtils::F77->trail_ ? "" : "-DNO_TRAILING_USCORE";
#$has{DEFINE} .= ExtUtils::F77->runtimeok ? "" : "-DF77_NOK";

WriteMakefile(%hash);

# Add genpp rule
# add other makefile additions as required (see also ExtUtils::MakeMaker)
sub MY::postamble {
	pdlpp_postamble(@::pack);
}

