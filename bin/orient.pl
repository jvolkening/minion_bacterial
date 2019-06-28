#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use BioX::Seq::Stream;
use Getopt::Long;
use File::Temp;

my $fn_assembly;
my $fn_starts;
my $threads = 1;

GetOptions(
    'assembly=s' => \$fn_assembly,
    'starts=s'   => \$fn_starts,
    'threads=i'  => \$threads,
);


die "Missing or unreadable assembly file"
    if (! -r $fn_assembly);
die "Missing or unreadable starts file"
    if (! -r $fn_starts);
die "Bad thread count"
    if ($threads =~ /\D/);

my %seqs;
my $p = BioX::Seq::Stream->new($fn_assembly);
while (my $seq = $p->next_seq) {
    $seqs{$seq->id} = $seq;
}

my $ret;

my $meta = {};

bwa(
    $fn_assembly,
    $fn_starts,
    $meta,
);

# rev com and re-map if needed
if (grep {$_} values %{ $meta->{rc} }) {
    my $dir = File::Temp->newdir(CLEANUP => 1);
    my $fn_tmp = "$dir/tmp.fa";
    open my $tmp, '>', $fn_tmp;
    for my $id (
        grep {$meta->{rc}->{$_}} keys %{ $meta->{rc} }
    ) {
        $seqs{$id} = $seqs{$id}->rev_com();
        say STDERR "Reverse complemented $id";
        print {$tmp} $seqs{$id}->as_fasta;
    }
    close $tmp;
    bwa(
        $fn_tmp,
        $fn_starts,
        $meta,
    );

}

# perform rotation
for my $id (keys %{$meta->{rotate} }) {
    next if ($seqs{$id}->desc !~ /circular=true/);
    my $str = $seqs{$id}->seq;
    my $begin = substr $str, 0, $meta->{rotate}->{$id}, '';
    $str .= $begin;
    $seqs{$id}->seq = $str;
    say STDERR "Rotated $id by $meta->{rotate}->{$id}";
}


for my $chr (sort {
    length($seqs{$b}) <=> length($seqs{$a})
} keys %seqs) {
    print $seqs{$chr}->as_fasta;
}

sub bwa {

    my ($ass, $starts, $meta) = @_;

    # index assembly
    $ret = system( "bwa index $ass 2> /dev/null" );
    die "bwa index failed: $!" if $ret;

    open my $bwa, '-|', "bwa mem -t $threads $ass $starts 2> /dev/null";

    while (my $line = <$bwa>) {
        chomp $line;
        next if ($line =~ /^@/);
        my @fields = split "\t", $line;
        if ($fields[1] == 16) {
            $meta->{rc}->{$fields[2]} = 1;
        }
        elsif ($fields[1] == 0) {
            $meta->{rotate}->{$fields[2]} = $fields[3] - 1;
            $meta->{rc}->{$fields[2]} = 0;
        }
    }

}
