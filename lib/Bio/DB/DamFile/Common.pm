package Bio::DB::DamFile::Common;
 
use strict;
use base 'Exporter';
our @EXPORT = qw(HEADER HEADER_STRUCT MAGIC FORMAT_VERSION BLOCKSIZE DEFAULT_BLOCK_CACHE_SIZE);

# header format:
#   4 bytes   -- magic number 'DAM1'
#   4 bytes   -- version number * 100
#   8 bytes   -- offset to beginning of SAM header data
#   8 bytes   -- offset to beginning of bunzip2 block data
#   8 bytes   -- offset to beginning of read name index
#   N+1 bytes -- full path to original BAM/SAM file, zero terminated string
#   512-(36+N+1) -- reserved for future expansion

use constant HEADER        => 512;
use constant FORMAT_VERSION=> 1.01;
use constant HEADER_STRUCT => 'a4LQQQZ*';
use constant MAGIC         => 'DAM1';
use constant BLOCKSIZE     => 1_048_576;   # A megabyte
use constant DEFAULT_BLOCK_CACHE_SIZE => BLOCKSIZE * 100; # hold ~100 blocks in memory
1;
