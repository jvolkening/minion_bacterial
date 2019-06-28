#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use File::Temp;
use Getopt::Long;
use List::Util qw/sum/;

my $fn_ref;
my $min_depth = 5;
my $min_af    = 0.8;
my $threads   = 1;

GetOptions(
    'ref=s'       => \$fn_ref,
    'threads=i'   => \$threads,
    'min_depth=i' => \$min_depth,
    'min_af=i'    => \$min_af,
);

# remaining trailing arguments are read names
my @reads = @ARGV;
die "Must supply at least one reads filename\n"
    if (scalar @reads < 1);
die "Too many filenames specified\n"
    if (scalar @reads > 2);

my $ret;

# index assembly
say STDERR "indexing assembly";
$ret = system( "bwa index $fn_ref 2> /dev/null" );
die "Error during bwa index: $!\n"
    if ($ret);

# perform read mapping
say STDERR "mapping reads";
my $bam = File::Temp->new(UNLINK => 1);
$ret = system( "bwa mem -t $threads $fn_ref @reads 2> /dev/null | samtools sort -\@$threads -o $bam - 2> /dev/null" );
die "Error during bwa mem: $!\n"
    if ($ret);
$ret = system( "samtools index $bam ");
die "Error during samtools index: $!\n"
    if ($ret);

# LoFreq indel quals
say STDERR "adding indel qualities";
my $lofreq_dir = File::Temp->newdir(CLEANUP => 0);
my $dindel = "$lofreq_dir/tmp.bam";
$ret = system( "lofreq indelqual --dindel -f $fn_ref -o $dindel $bam 2> /dev/null" );
die "Error during lofreq indelqual: $!\n"
    if ($ret);
$ret = system( "samtools index $dindel ");
die "Error during samtools index dindel: $!\n"
    if ($ret);

# LoFreq variant calling
say STDERR "calling variants";
my $vcf = "$lofreq_dir/tmp.vcf";
$ret = system( "lofreq call-parallel --pp-threads $threads -f $fn_ref -o $vcf --call-indels $dindel 2> /dev/null" );
die "Error during lofreq call: $!\n"
    if ($ret);

# calculate "covered length" -- i.e. the length of the genome covered by reads
# at a minimum depth to allow variant calling. This will be used as the genome
# length in downstream calculations

say STDERR "calculating coverage";
open my $cov, '-|',
    "bedtools",
    'genomecov',
    '-d',
    '-ibam' => $bam,
;

my %covered;
my %total;
while (<$cov>) {
    chomp;
    my @fields = split "\t";
    ++$total{$fields[0]};
    ++$covered{$fields[0]}
        if ($fields[2] >= $min_depth);
}

# calculate the number of mismatched bases
my %mm;
my %snv;
open my $in, '<', $vcf;
while (my $line = <$in>) {
    
    next if ($line =~ /^#/);
    chomp $line;

    my (
        $chr,
        $pos,
        $id,
        $given,
        $alt,
        $qual,
        $filt,
        $info,
    ) = split "\t", $line;

    my %fields;
    for my $pair (split ';', $info) {
        my ($key, $val) = ($pair =~ /^([^=]+)(?:=(.+))?$/);
        die "Invalid pair $pair"
            if (! defined $key);
        $fields{$key} = $val; # $val can be undefined
    }
    next if ($fields{AF} < $min_af);
    next if ($fields{DP} < $min_depth);

    my $delta = length($alt) - length($given);

    # deletions in assembly
    if ($delta > 0) {
        # add a position to the assembly length, since this essentially
        # represents positions called as '-'
        ++$covered{$chr};
        $mm{$chr} += $delta;
    }
    # insertions in assembly
    elsif ($delta < 0) {
        $mm{$chr} += abs($delta);
    }
    else {
        $mm{$chr} += 1;
        $snv{$chr} += 1;
    }
}

my $cov_tot = sum values %covered;
my $mm_tot  = sum values %mm;
my $snv_tot = sum values %snv;

my $err_rate = $mm_tot/$cov_tot;
my $snv_rate = $snv_tot/$cov_tot;
my $indel_rate = $err_rate - $snv_rate;
my $accuracy = 1 - $err_rate;

say join "\t",
    sprintf("%0.5f", $accuracy),
    $snv_rate * 1000,
    $indel_rate * 1000,
;

#for my $chr (sort keys %covered) {
    #my $err_rate   = $mm{$chr}/$covered{$chr};
    #my $snv_rate   = $snv{$chr}/$covered{$chr};
    #my $indel_rate = $err_rate - $snv_rate;
    #my $accuracy   = 1 - $err_rate;
    #say join "\t",
        #$chr,
        #map {sprintf "%0.2f", $_*100} (
            #$accuracy,
            #$snv_rate,
            #$indel_rate,
        #),
    #;
#}
