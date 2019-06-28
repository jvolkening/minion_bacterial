#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use autodie;
use BioX::Seq::Stream;
use Cwd qw/abs_path/;
use Getopt::Long;
use PerlIO::gzip;

my $in_map;
my $in_dir;
my $in_trim;
my $in_untrim;
my $out_dir = '.';
my @durations; # in minutes

GetOptions(
    'read_map=s'  => \$in_map,
    'in_fast5=s'  => \$in_dir,
    'in_trim=s'   => \$in_trim,
    'in_untrim=s' => \$in_untrim,
    'durations=s' => \@durations,
    'out_dir=s'   => \$out_dir,
);

# allow durations to be specified separately or as comma-joined string
@durations = split ',', join(',', @durations);

# set up output directories and output fastq filehandles
my %fhs;
die "$out_dir exists and won't overwrite"
    if (-e $out_dir);
mkdir $out_dir;
for my $d (@durations) {
    my $tgt = "$out_dir/$d";
    die "$tgt exists and won't overwrite"
        if (-e $tgt);
    mkdir $tgt;
    mkdir "$tgt/fast5";
    if (defined $in_trim) {
        open my $fh1, '>:gzip', "$tgt/trimmed.fq.gz";
        open my $fh2, '>:gzip', "$tgt/untrimmed.fq.gz";
        $fhs{$d} = [$fh1,$fh2];
    }
}

my %times;

open my $fh_map, '<:gzip', $in_map;
my $head = <$fh_map>;
while (my $line = <$fh_map>) {
    chomp $line;
    my ($fn, $id, $elapsed) = split "\t", $line;
    $times{$id} = $elapsed;
    my $fn_in = "$in_dir/$fn";
    # create symlinks input FAST5
    if (-e $fn_in) {
        for my $t (@durations) {
            if ($elapsed <= $t*60) {
                my $fn_out = "$out_dir/$t/fast5/$fn";
                symlink abs_path($fn_in), $fn_out;
            }
        }
    }
    else {
        die "Failed to find input FAST5: $fn_in\n";
    }
}

# write FASTQ output if asked
if (defined $in_trim) {

    my $p1 = BioX::Seq::Stream->new($in_trim);
    while (my $seq = $p1->next_seq) {
        my $id = $seq->id;
        my $elapsed = $times{ $id }
            // die "Trimmed read $id missing in the file map\n";
        for my $t (@durations) {
            if ($elapsed <= $t*60) {
                print {$fhs{$t}->[0]} $seq->as_fastq;
            }
        }

    }
    my $p2 = BioX::Seq::Stream->new($in_untrim);
    while (my $seq = $p2->next_seq) {
        my $id = $seq->id;
        my $elapsed = $times{ $id }
            // die "Untrimmed read $id missing in the file map\n";
        for my $t (@durations) {
            if ($elapsed <= $t*60) {
                print {$fhs{$t}->[1]} $seq->as_fastq;
            }
        }

    }

}

