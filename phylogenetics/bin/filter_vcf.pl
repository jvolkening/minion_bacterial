#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use BioX::Seq::Fetch;
use List::Util qw/any/;

use Getopt::Long;

my $in_vcf;
my $out_vcf;
my $min_freq    = 0.9;
my $drop_indels = 0;
my $drop_snps   = 0;

GetOptions(
    'in=s'        => \$in_vcf,
    'out=s'       => \$out_vcf,
    'min_freq=f'  => \$min_freq,
    'drop_indels' => \$drop_indels,
    'drop_snps'   => \$drop_snps,
);

open my $in, '<', $in_vcf;
open my $out, '>', $out_vcf;

LINE:
while (my $line = <$in>) {

    if ($line =~ /^\s*#/) {
        print {$out} $line;
        next LINE;
    }

    chomp $line;

    my @fields = split "\t", $line;
    my $info = $fields[7];
    my @annots = split ';', $info;
    my $is_indel = any {$_ eq 'INDEL'} @annots;
    next LINE if ($is_indel && $drop_indels);
    next LINE if ((! $is_indel) && $drop_snps);

    for my $annot (@annots) {
        if ($annot =~ /^AF=([\d\.]+)/) {
            say STDERR $1;
            next LINE if ($1 < $min_freq);
            say STDERR 'FOO';
        }
    }

    say {$out} $line;

}
