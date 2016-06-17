# BamDessicate
Reduce the size of BAM files by stripping out read and quality information, while preserving mapping data.

# Command-line usage

```shell 
 # compress a BAM file by stripping read and quality information, creating a .dam file
 dessicate.pl big_file.bam small_file.dam

 # rehydrate the .dam file by adding back the read and quality info
 # (read information can come from a SAM, BAM or FASTQ file)
 hydrate.pl small_file.dam reads.fastq.gz out.bam
```

# Description
This Perl library was created to solve the issue of maintaining
multiple aligned BAM files from the same set of reads. This happens
when a BAM file is remapped onto different genome builds or using
different alignment software/settings. Rather than have multiple
copies of the same read and quality score information, one wishes to
maintain a single BAM (or FASTQ) file with the read information, and
store the alternative alignments in separate data files.

The DAM ("dessicated BAM") format is very simple, and consists of the
standard SAM header followed by a series of bzip2-compressed chunks of
alpha-sorted SAM (text format) alignment lines from which the read and
quality score fields have been removed. This is followed by a
alphabetically-sorted index of the first read name in each block
followed by the offset to that block. An individual read line can be
found by performing a binary search on the read name index, fetching
and uncompressing the corresponding block, and then performing a
binary search on the SAM lines contained within the uncompressed
block.

The DAM file is typically reduced by a factor of 6-8 relative to the
size of the original BAM file, making it an economical alternative to
storing multiple remapped BAM files.

Here is an ASCII text representation of the DAM file:

```
 --------------------------------------------------------------------------------------
 |                       HEADER  (512 bytes)                                          |
 |DAM1                -- magic number,                                      4 bytes   |
 |SAM_header_offset   -- byte position of the start of the SAM header,      16 bytes  |
 |block_offset        -- byte position of the first compressed block,       16 bytes  |
 |index_offset        -- byte position of the start of the read name index, 16 bytes  |
 |path_name           -- full path to the source BAM/SAM file                variable |
 |                        (this is a zero-terminated string)                          |
 --------------------------------------------------------------------------------------
 |                      SAM HEADER (variable length)                                  |
 | Uncompressed text version of the SAM/BAM header                                    |
 --------------------------------------------------------------------------------------
 |  BLOCK 1 (variable length) -bzip2 compressed SAM lines, newline terminated,        |
 |                              up to 1 MB in length prior to compression.            |
 --------------------------------------------------------------------------------------
 |  BLOCK 2 (variable length)                                                         |
 --------------------------------------------------------------------------------------
 ~
 --------------------------------------------------------------------------------------
 |  BLOCK N (variable length)                                                         |
 --------------------------------------------------------------------------------------
 |  READ Index (variable length) -- 
              bzip2-compressed data. When uncompressed looks like this:               |  
 | ReadName1\0Offset
 | ReadName2\0Offset
 |...
 | ReadNameN\Offset
 |
 | ReadName is a null terminated string. Offset is the position of the compressed block
 | that begins with the ReadName SAM line. Offset is 8 bytes in length.
 --------------------------------------------------------------------------------------
```

# Installation

After cloning or downloading the repository, please run the following
commands from the shell:

```shell
> perl Build.PL
> ./Build test
> sudo ./Build install
```

# Perl Library Documentation

The Perl documentation is embedded in the file
lib/Bio/DB/DamFile.pm. It will also be available as a manual page once
you install the module on your system. Use:

```shell
perldoc Bio::DB::DamFile
```