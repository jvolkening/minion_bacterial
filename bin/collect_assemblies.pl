#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use File::Copy qw/copy/;
use File::Path qw/make_path/;

my $dir_out = $ARGV[0]
    // die "output directory undefined";

my @bcs = qw/
    BC01
    BC02
/;

my @times = (
    15,
    30,
    60,
    120,
    240,
    480,
    960,
    1500,
);

for my $bc (@bcs) {
for my $t (@times) {

make_path("$dir_out/$bc/$t")
    or die "Failed to make output path: $@\n";

my %fns = (
    rnd0 => "$bc/results.prev/$t/orient/oriented.fasta",
    rnd1 => "$bc/results.prev/$t/nanopolish/oriented.polished.fa",
    rnd2 => "$bc/results.prev/$t/nanopolish/oriented.polished.polished.fa",
);

while (my ($rnd, $fn_in) = each %fns) {
    die "Failed to find $fn_in\n"
        if (! -e $fn_in);
    my $fn_out = "$dir_out/$bc/$t/$rnd.fa";
    say "copying $fn_in to $fn_out";
    copy( $fn_in => $fn_out )
        or die "Failed to copy file: $@\n";
}

}}
