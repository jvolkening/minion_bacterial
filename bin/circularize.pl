#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use BioX::Seq::Stream;
use Cwd qw/getcwd abs_path/;
use File::Temp;

my $fn_minimus2 = abs_path( $ARGV[0] );

my $max_split_size = 200000;

my %seqs;
while (my $line = <STDIN>) {

    chomp $line;
    my @fields = split "\t", $line;

    if ($fields[0] eq 'S') {
        my $seq = BioX::Seq->new(
            $fields[2],
            $fields[1],
        );
        $seqs{$seq->id} = {
            seq => $seq,
            len => length($seq),
        }
    }

    elsif ($fields[0] eq 'L') {
        if (
            $fields[1] eq $fields[3]
         && $fields[2] eq $fields[4]
         && $fields[5] eq '0M'
        ) {
            $seqs{$fields[1]}->{is_circular} = 1;
        }

    }

}

for my $chr (sort {
    $seqs{$b}->{len} <=> $seqs{$a}->{len}
} keys %seqs) {
    if ($seqs{$chr}->{is_circular}) {
        my $seq = $seqs{$chr}->{seq};
        $seq->desc = defined $seq->desc
            ? $seq->desc . ' circular=true'
            : 'circular=true';
        print $seq->as_fasta;
    }
    else {
        print circularize(
            $seqs{$chr}->{seq}
        )->as_fasta;
    }
}

sub circularize {

    my ($seq) = @_;
    my $wd = File::Temp->newdir(CLEANUP => 1);

    my $cwd = abs_path( getcwd());
    chdir $wd;

    my $id = $seq->id;

    my $l = length($seq);
    my $split_size = $l/2 < $max_split_size
        ? int($l/2)
        : $max_split_size;

    open my $split, '>', "split.fa";
    print {$split} $seq->range(1 => $split_size)->as_fasta;
    print {$split} $seq->range($split_size+1 => $l)->as_fasta;
    close $split;

    my $ret = system( "toAmos -s split.fa -o split.afg > /dev/null" );
    die "toAmos failed: $!" if ($ret);

    # Try to run minimus2 with increasingly lower errors until it succeeds.
    # This is to deal with the issue of minimus2 crashing with some input seqs
    for my $err (0.2,0.1,0.05) {
        my $ret = system( "runAmos -C $fn_minimus2 split -D CONSERR=$err > /dev/null");
        last if (! $ret);
    }

    # run succeeded
    if (-e 'split.fasta') {

        # was circularized
        if (-s 'split.fasta') {
            my @seqs;
            my $p = BioX::Seq::Stream->new('split.fasta');
            while (my $s = $p->next_seq) {
                push @seqs, $s;
            }
            if (scalar(@seqs) == 1 && (! -s 'split.singletons')) {
                my $id = $seq->id;
                $seq = $seqs[0];
                $seq->id = $id;
                $seq->desc = defined $seq->desc
                    ? $seq->desc . ' circular=true'
                    : 'circular=true';
                say STDERR "Circularized $id";
            }
            else {
                warn "Unexpected number of output contigs or singletons for $id\n";
            }

        }
        else {
            say STDERR "No circularization for $id";
        }

    }

    else {
        warn "Minimus2 failed completely on $id, returning original sequence\n";
    }

    chdir $cwd;

    return $seq;

}

