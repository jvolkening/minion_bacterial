#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use BioX::Seq::Stream;
use File::Basename qw/basename/;
use Getopt::Long;
use List::Util qw/uniq/;

my $fn_ref;
my $use_all       = 0; # include conserved sites in output
my $min_diff      = 2; # minimum number of sequences a SNP must be present in
my $only_matching = 0; # only use SNPs that are found in at least one Illumina dataset
my $min_ctg_len   = 1000;

GetOptions(
    'ref=s'  => \$fn_ref,
    'use_all' => \$use_all,
    'min_diff=i' => \$min_diff,
    'only_matching' => \$only_matching,
    'min_ctg_len=i' => \$min_ctg_len,
);
my @fns_vcf = @ARGV;

my %snps;
my $skip;
my $keep;

say STDERR "Parsing VCFs...";

for my $vcf (@fns_vcf) {

    my ($base) = (basename($vcf) =~ /^(\w+)/);
    die "base not found for $vcf\n"
        if (! defined $base);

    $snps{$base} = {};

    open my $in, '<', $vcf;

    LINE:
    while (my $line = <$in>) {
        chomp $line;
        next if ($line =~ /^\s*#/);
        my (
            $chr,
            $pos,
            $ig1,
            $ref,
            $alt,
            $score,
            $ig2,
            $info,
        ) = split "\t", $line;
        my $state = 'OK';
        if ($info =~ /STATUS=(\w+)/) {
            if ($1 ne 'OK') {
                $skip->{$chr}->{$pos} = 1;
                next LINE;
            }
        }
        else {
            $keep->{$chr}->{$pos} = $base;
        }
            
        $snps{$base}->{$chr}->{$pos} = [$ref, $alt];
    }

}

my @samples = sort keys %snps;
my @seqs;
for (@samples) {
    push @seqs, BioX::Seq->new('', $_);
}

say STDERR "Walking genome...";

my $p = BioX::Seq::Stream->new($fn_ref);
while (my $seq = $p->next_seq) {
    next if (length($seq) < $min_ctg_len);
    my $chr = $seq->id;
    for my $i (1..length($seq)) {
        next if ($skip->{$chr}->{$i});
        #next if ($only_matching && ! $keep->{$chr}->{$i});
        my $ref = substr $seq, $i-1, 1;
        next if ($ref eq 'N');
        my @nas;
        for my $sample (@samples) {
            my $snp = $snps{$sample}->{$chr}->{$i};
            if (defined $snp) {
                die "Ref base mismatch with $sample $chr $i\n"
                    if ($snp->[0] ne $ref);
                push @nas, $snp->[1];
            }
            else {
                push @nas, $ref;
            }
        }

        # filter out SNPs that occur in less than $min_diff samples
        my $multis = 0;
        for my $c (uniq @nas) {
            if ((scalar grep {$_ eq $c} @nas) >= $min_diff) {
                ++$multis;
            }
        }
        my $n_diff = scalar uniq @nas;
           
        if ( ($use_all && $n_diff==1)
          || ($multis >= 2 && (! $only_matching || $keep->{$chr}->{$i}) ) ) {
            my $ill = $keep->{$chr}->{$i} // 'none'; 
            say STDERR "Using $chr $i ($ill)";
            for (0..$#seqs) {
                $seqs[$_]->seq .= $nas[$_];
            }
        }
    }
}
                
print $_->as_fasta for (@seqs);    

