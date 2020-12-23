use ExtUtils::MakeMaker;
use PDL::Core::Dev;
my %hash = pdlpp_stdargs(['WCS.pp','WCS','WCS']);
$hash{'LIBS'} = " -Llibwcs  -lwcs  -lm";
#$hash{'INC'} .= " -I/sw/include/";
$hash{CCFLAGS} .= " -g";
$hash{'MYEXTLIB'} = ' wcs_call.o';
WriteMakefile(%hash);
sub MY::postamble{
    pdlpp_postamble(['WCS.pp','WCS','WCS'])
    }
