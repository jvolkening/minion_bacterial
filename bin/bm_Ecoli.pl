#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use File::Basename qw/basename/;
use File::Copy qw/move/;
use Getopt::Long;

my $dir_in;
my $dir_out;
my $dir_meta;
my $threads = 1;

GetOptions(
    'in=s'      => \$dir_in,
    'out=s'     => \$dir_out,
    'meta=s'    => \$dir_meta,
    'threads=s' => \$threads,
);

if (! -e $dir_out) {
    mkdir $dir_out;
}

my @durations = sort {basename($a) <=> basename($b)} glob "$dir_in/*";

# add last run with Pilon
#push @durations, "$durations[-1]_pilon";

for my $dur (@durations) {

    my $base = basename($dur);

    say "Running workflow on $base";

    mkdir "$dir_out/$base";

    $dur =~ s/_pilon$//;

    my @call = (
        'nextflow'       => 'run',
        '--outdir'       => "$dir_out/$base",
        '--threads'      => $threads,
        '--trimmed'      => "$dur/trimmed.fq.gz",
        '--fast5dir'     => "$dur/fast5",
        '--id_map'       => "$dir_meta/Ecoli.map.tsv.gz",
        '--untrimmed'    => "$dur/untrimmed.fq.gz",
        '--depth'        => '75',
        '--starts'       => "$dir_meta/ref/Ecoli_starts.fa",
        '--buscolineage' => "$dir_meta/enterobacteriales_odb9",
        '--illumina'     => 'raw/illumina/Ecoli/\*.fastq.gz',
        "minion_bacterial.nf",
        '-resume',
        '-with-trace',
    );
    if ($base =~ /_pilon$/) {
        push @call, '--run_pilon';
    }

    my $ret = system(@call);
    die "workflow failed for $base: $!"
        if ($ret);
    move 'trace.txt', "$dir_out/$base/trace.txt";
}
