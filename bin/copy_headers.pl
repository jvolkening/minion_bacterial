#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use BioX::Seq::Stream;

my ($fn_from, $fn_to) = @ARGV;

my $p_from = BioX::Seq::Stream->new($fn_from);
my $p_to   = BioX::Seq::Stream->new($fn_to);

while (my $seq_from = $p_from->next_seq) {
    my $seq_to = $p_to->next_seq
        or die "Not enough sequences in target file";
    $seq_to->id   = $seq_from->id;
    $seq_to->desc = $seq_from->desc;
    print $seq_to->as_fasta;
}

die "Too many sequences in target file"
    if ($p_to->next_seq);

exit;
    
