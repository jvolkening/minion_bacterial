#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use BioX::Seq::Fetch;

use Getopt::Long;
use List::Util qw/sum/;

my $in_vcf;
my $in_fa;
my $out_vcf;
my $min_dist    = 2;
my $homo_n      = 3;
my $homo_dist   = 3;
my $methyl_dist = 10;
my $tag_agtc    = 0;
my $min_af      = 0.9;

GetOptions(
    'in=s'          => \$in_vcf,
    'ref=s'         => \$in_fa,
    'out=s'         => \$out_vcf,
    'min_dist=i'    => \$min_dist,
    'homo_n=i'      => \$homo_n,
    'homo_dist=i'   => \$homo_dist,
    'methyl_dist=i' => \$methyl_dist,
    'tag_agtc'      => \$tag_agtc,
    'min_af=f'      => \$min_af,
);

my $fa = BioX::Seq::Fetch->new($in_fa);

say STDERR "Parsing $in_vcf";
open my $in, '<', $in_vcf;

my %indels;
my %snps;
my %chrs;
my $i_chr = 0;

open my $out, '>', $out_vcf;

LINE:
while (my $line = <$in>) {

    if ($line =~ /^\s*#/) {
        print {$out} $line;
        next LINE;
    }

    chomp $line;

    my @fields = split "\t", $line;
    my $chr  = $fields[0];
    my $pos  = $fields[1];
    my $info = $fields[7];
    if ($info =~ /INDEL/) {
        $indels{$chr}->{$pos} = $line;
    }
    else {
        $snps{$chr}->{$pos} = $line;
    }

    $chrs{$chr} //= $i_chr++;

}

for my $chr ( sort {$chrs{$a} <=> $chrs{$b}} keys %chrs) {

    my $chr_len = $fa->length($chr);


    for my $pos (sort {$a <=> $b} keys %{ $snps{$chr} }) {
        my $methyl;
        my $cluster;
        my $hrun;
        my $agtc;
        my $minaf;

        my $line = $snps{$chr}->{$pos};
        my @fields = split "\t", $line;

        # flag A->G and T->C SNPs
        my $ref = $fields[3];
        my $alt = $fields[4];
        my $tag = "$ref$alt";
        if ($tag_agtc && ($tag eq 'AG' || $tag eq 'TC')) {
            $agtc = "AG_TC";
        }

        # check alternative allele frequency
        my $info = $fields[7];
        if ($info =~ /DP4=([\d,]+)/) {
            my @dps = split ',', $1;
            my $ratio = sum(@dps[2..3])/sum(@dps);
            if ($ratio < $min_af) {
                $minaf = 'MINAF';
            }
        }
        else {
            die "Unexpected or missing DP4 string";
        }

        # flag SNPs near methylation sites
        my $start = $pos - $methyl_dist;
        my $end   = $pos + $methyl_dist;
        $start = 1 if ($start < 1);
        $end = $chr_len if ($end > $chr_len);
        my $sub = $fa->fetch_seq($chr, $start, $end);
        for my $motif (qw/
            GATC
            [GC]C[AT]G[GC]
        /)  {
            if ($sub =~ /$motif/) {
                $methyl = 'METHYL';
            }
        }

        # flag SNPs within a maximum distance of other variants
        for my $p ( ($pos-$min_dist)..($pos+$min_dist) ) {
            next if ($p == $pos);
            if ( defined $snps{$chr}->{$p}
             || defined $indels{$chr}->{$p} ) {
                $cluster = 'CLUSTER';
            }
        }

        # flag SNPs near homopolymer stretches
        $start = $pos - ($homo_n + $homo_dist - 1);
        $end   = $pos + ($homo_n + $homo_dist - 1);
        $start = 1 if ($start < 1);
        $end = $chr_len if ($end > $chr_len);
        $sub = $fa->fetch_seq($chr, $start, $end);
        for my $base (qw/ A T G C/) {
            if ($sub =~ /${base}{$homo_n,}/) {
                $hrun = 'HRUN';
            }
        }

        my $status = join '|', grep {defined $_} ($methyl,$cluster,$hrun,$agtc,$minaf);
        $status = 'OK' if (! $status);
        
        say {$out} $snps{$chr}->{$pos}, ";STATUS=$status";

    }

}

close $out; 
