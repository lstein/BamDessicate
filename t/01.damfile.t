#-*-Perl-*-

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use File::Temp qw(tempfile);
use FindBin '$Bin';
use constant TEST_COUNT => 1;

use lib "$Bin/../lib","$Bin/../blib/lib","$Bin/../blib/arch";

BEGIN {
  use Test::More tests => TEST_COUNT;
  use_ok('Bio::DB::DamFile');

}


{
    # no tests yet
}

__END__
