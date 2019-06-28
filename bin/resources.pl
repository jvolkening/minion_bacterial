#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use autodie;
use Cwd qw/getcwd abs_path/;
use File::Basename qw/basename/;
use Getopt::Long;
use List::Util qw/any/;
use Excel::Writer::XLSX;

GetOptions(
);

my @skip_processes = qw/
    dnadiff
    busco
    lofreq
/;

my @headers = qw/
    duration
    reads
    sampled_reads
    assembly_size
    circular_contigs
    linear_contigs
    longest_contig
    longest_circular
    NG50
/;
my @rep_headers = qw/
    reference_coverage
    avg_identity
    relocations
    translocations
    inversions
    insertions
    insertion_sum
    SNPs_per_kb
    indels_per_kb
    BUSCO_complete
    BUSCO_fragmented
    BUSCO_missing
    time(min)
    max_memory(GB)
/;

my $ref;

my $h = <STDIN>;
chomp $h;
while (my $line = <STDIN>) {
    chomp $line;
    my @fields = split "\t", $line;
    my $proc = $fields[3];
    my $cpu  = $fields[9];
    $cpu =~ s/\%//;
    $cpu /= 100;
    $proc =~ s/\s.+//g;
    next if (any {$proc eq $_} @skip_processes);
    my $dur = parse_duration($fields[8]);
    my $cpu_hr = $dur * $cpu;
    my $mem = parse_mem($fields[10]);
    if ($proc =~ /^(?:nanopolish|vcf2fasta)/) {
        my $tproc = 'nanopolish';
        $ref->{$tproc}->{'cpu_hr'}->[2] += $cpu_hr;
        $ref->{$tproc}->{'max_memory(GB)'}->[2] = $mem
            if ($mem > ($ref->{$tproc}->{'max_memory(GB)'}->[2] // 0));
        if ($proc !~ /_2$/) {
            $ref->{$tproc}->{'cpu_hr'}->[1] += $cpu_hr;
            $ref->{$tproc}->{'max_memory(GB)'}->[1] = $mem
                if ($mem > ($ref->{$tproc}->{'max_memory(GB)'}->[1] // 0));
        }
    }
    else {
        my $tproc = $proc;
        if ($proc =~ /^(?:circularize|orient|preprocess|filtlong)/) {
            $tproc = 'other';
        }
            
        for my $rnd (0..2) {
            $ref->{$tproc}->{'cpu_hr'}->[$rnd] += $cpu_hr;
            $ref->{$tproc}->{'max_memory(GB)'}->[$rnd] = $mem
                if ($mem > ($ref->{$tproc}->{'max_memory(GB)'}->[$rnd] // -1));
        }
    }
}
my @procs = qw/
    unicycler
    nanopolish
    other
/;

say "# CPU usage";
say join "\t",
    "process",
    "no Nanopolish",
    "one rnd Nanopolish",
    "two rnds Nanopolish",
;
for my $proc (@procs) {
    my @vals;
    for my $rnd (0..2) {
         push @vals, sprintf '%0.1f',
            ($ref->{$proc}->{'cpu_hr'}->[$rnd] // 0)
    }
    say join "\t",
        $proc,
        @vals,
    ;
}

say "# max memory usage";
say join "\t",
    "process",
    "no Nanopolish",
    "one rnd Nanopolish",
    "two rnds Nanopolish",
;
for my $proc (@procs) {
    my @vals;
    for my $rnd (0..2) {
         push @vals, sprintf '%0.1f',
            ($ref->{$proc}->{'max_memory(GB)'}->[$rnd] // 0)
    }
    say join "\t",
        $proc,
        @vals,
    ;
}

sub parse_duration {
    
    my ($str) = @_;
    my @parts = split ' ', $str;
    my $t = 0;
    for (@parts) {
        my ($f, $l) = (/^([\d\.]+)([a-z]+)$/);
        die "bad duration string: $str"
            if (! defined $l);
        $t += $l eq 'd'  ? $f * 24
            : $l eq 'h'  ? $f
            : $l eq 'm'  ? $f/60
            : $l eq 's'  ? $f/60/60
            : $l eq 'ms' ? $f/60/60/1000
            : die "Invalid time unit $l in string $str";
    }
    return $t;

}

sub parse_mem {
    
    my ($str) = @_;
    if ($str eq '0') {
        return 0;
    }
    elsif ($str =~ /^([\d\.]+) (GB|MB)/) {
        my $mem = $1;
        $mem /= 1024 if ($2 eq 'MB');
        return $mem;
    }
    else {
        die "bad memory string: $str";
    }

}
