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

use BioX::Seq::Stream;

my @bcs = qw/
    BC01
    BC02
/;

my $n_rnds = 3;
my $out_fmt = 'xlsx';
my $dir_in = 'results';
my $orig_in;

GetOptions(
    'out_fmt=s'     => \$out_fmt,
    'dir_in=s'      => \$dir_in,
    'orig_in=s'     => \$orig_in,
);

my @skip_processes = qw/
    dnadiff
    busco
    lofreq
/;

my $home = abs_path(getcwd());

binmode(STDOUT);

my ($wb, %ws, $fmt_bc, $fmt_c);

if ($out_fmt eq 'xlsx') {
    $wb = Excel::Writer::XLSX->new( \*STDOUT );
    for (@bcs) {
        $ws{$_} = $wb->add_worksheet($_);
    }
    $fmt_bc = $wb->add_format();
    $fmt_bc->set_bold();
    $fmt_bc->set_align( 'center' );
    $fmt_c = $wb->add_format();
    $fmt_c->set_align( 'center' );
}

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

for my $bc (@bcs) {

    # Pre-process reference sizes
    my $ref_len;
    my $fn_ref = "$bc/ref.fa";
    die "Missing reference file"
        if (! -e $fn_ref);
    my $p = BioX::Seq::Stream->new($fn_ref);
    while (my $seq = $p->next_seq) {
        $ref_len += length $seq;
    }

    # Pre-process original read counts
    my %orig_counts;
    for my $dir (
        sort {basename($a) <=> basename($b)}
        glob "$orig_in/$bc/*"
    ) {
        my $dur = basename($dir);
        my $str = `wc -l $dir/reads.fq`;
        #my $str = '12345 foo'; # for testing
        if ($str =~ /^(\d+)/) {
            $orig_counts{$dur} = $1/4;
        }
    }

    my @results;

    for my $dir (
        sort {
            (basename($a)  =~ /(\d+)/)[0] <=> (basename($b) =~ /(\d+)/)[0]
        || basename($a) cmp basename($b)
        }
        glob "$bc/$dir_in/*"
    ) {

        my $dur = basename($dir);

        #next if ($dur > 480);

        warn "Parsing $dur\n";


        chdir $dir;

        # account for tagged durations
        my $sub = $dur;
        $sub =~ s/\D//g;

        my $ref = {
            duration => $dur,
            reads    => $orig_counts{$sub},
        };

        parse_filtlong($ref);

        parse_assembly($ref, $ref_len);

        parse_dnadiff($ref)
            if (-e 'dnadiff/oriented.report');

        parse_lofreq($ref)
            if (-e 'lofreq/calc_accuracy_oriented.tsv');

        parse_busco($ref)
            if (-e 'busco/full_table_oriented.tsv');

        parse_traces($ref);

        push @results, $ref;

        chdir $home;

    }

    if ($out_fmt eq 'tsv') {
        say join "\t", qw/
            barcode
            time
            polishings
            variable
            value
        /;
    }

    # filter out unused columns
    my @used_headers;
    my $col = 0;
    for my $h (@headers) {
        next if (! defined $results[0]->{$h});
        push @used_headers, $h;
        if ($out_fmt eq 'xlsx') {
            $ws{$bc}->write(
                1,
                $col,
                $h,
                $fmt_bc,
            );
            $ws{$bc}->set_column(
                $col,
                $col++,
                length($h)+4,
            );
        }
    }

    my @used_rep_headers;
    my @rnd_labs = (
        "No polishing",
        "One round Nanopolish",
        "Two rounds Nanopolish",
    );
    for (0..($n_rnds-1)) {
        my $col_start = $col;
        for my $h (@rep_headers) {
            next if (! defined $results[0]->{$h});
            push @used_rep_headers, $h if ($_ == 1);
            if ($out_fmt eq 'xlsx') {
                $ws{$bc}->write(
                    1,
                    $col,
                    $h,
                    $fmt_bc,
                );
                $ws{$bc}->set_column(
                    $col,
                    $col++,
                    length($h)+4,
                );
            }
        }
        if ($out_fmt eq 'xlsx') {
            $ws{$bc}->merge_range(0,$col_start,0,$col-1,$rnd_labs[$_],$fmt_bc);
        }
    }

    # print out results table
    #say join "\t", @used_headers, ((@used_rep_headers) x $n_rnds);
    my $row = 2;
    for my $res (@results) {

        my @final;
        push @final,
            map {$res->{$_}} @used_headers;
        for my $rnd (0..($n_rnds-1)) {
            push @final,
                map {$res->{$_}->[$rnd]} @used_rep_headers;
        }
        if ($out_fmt eq 'xlsx') {
            $ws{$bc}->write_row($row++, 0, \@final, $fmt_c);
        }
        elsif ($out_fmt eq 'tsv') {
            my $dur = $res->{duration};
            delete $res->{duration};
            for my $key (keys %{ $res }) {
                for my $rnd (0..($n_rnds-1)) {
                    my $val = ref($res->{$key}) eq 'ARRAY'
                        ? $res->{$key}->[$rnd]
                        : $res->{$key};
                    say join "\t",
                        $bc,
                        $dur,
                        $rnd,
                        $key,
                        $val,
                    ;
                }
            }
        }
                    
    }

}

sub parse_duration {
    
    my ($str) = @_;
    my @parts = split ' ', $str;
    my $t = 0;
    for (@parts) {
        my ($f, $l) = (/^([\d\.]+)([a-z]+)$/);
        die "bad duration string: $str"
            if (! defined $l);
        $t += $l eq 'd'  ? $f * 60 * 24
            : $l eq 'h'  ? $f * 60
            : $l eq 'm'  ? $f
            : $l eq 's'  ? $f/60
            : $l eq 'ms' ? $f/60/1000
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

sub parse_filtlong {

    my ($ref) = @_;

    # count sampled reads
    my $p = BioX::Seq::Stream->new('filtlong/reads.filt.fq');
    while (my $seq = $p->next_seq) {
        ++$ref->{sampled_reads};
    }

}

sub parse_assembly {

    my ($ref, $ref_len) = @_;

    # parse Unicycler log
    $ref->{linear_contigs}   = 0;
    $ref->{circular_contigs} = 0;
    $ref->{assembly_size}    = 0;
    $ref->{longest_contig}   = 0;
    $ref->{longest_circular} = 0;

    # parse assembly
    #my $fn = -r 'pilon/oriented.polished.fasta'   ? 'pilon/oriented.polished.fasta'
           #: -r 'orient/oriented.fasta'           ? 'orient/oriented.fasta'
           #: -r 'circularized/circularized.fasta' ? 'circularized/circularized.fasta'
           #: die "no assembly found!";
    my $fn = -r 'nanopolish/oriented.polished.polished.fa'
        ? 'nanopolish/oriented.polished.polished.fa'
        : die "no assembly found!";
    my $fn2 = -r 'circularized/circularized.fasta'
        ? 'circularized/circularized.fasta'
        : die "no unpolished assembly found!";
    my @lens;
    my $p2 = BioX::Seq::Stream->new($fn2);
    my %circ;
    while (my $seq = $p2->next_seq) {
        if ($seq->desc && $seq->desc =~ /circular=true/) {
            $circ{$seq->id} = 1;
        }
    }
    my $p = BioX::Seq::Stream->new($fn);
    while (my $seq = $p->next_seq) {
        my $l = length $seq;
        $ref->{assembly_size} += $l;
        $ref->{longest_contig} = $l
            if ($l > $ref->{longest_contig});
        push @lens, $l;
        if ($circ{$seq->id}) {
            ++$ref->{circular_contigs};
            $ref->{longest_circular} = $l
                if ($l > $ref->{longest_circular});
        }
        else {
            ++$ref->{linear_contigs};
        }
    }

    # Calculate NG50
    my $cum  = 0;
    $ref->{NG50} = 0;
    for (sort {$b <=> $a} @lens) {
        $cum += $_;
        if ($cum >= $ref_len/2) {
            $ref->{NG50} = $_;
            last;
        }
    }

}

sub parse_dnadiff {

    my ($ref) = @_;

    my @assemblies = qw{
        dnadiff/oriented.report
        dnadiff/oriented.polished.report
        dnadiff/oriented.polished.polished.report
    };

    for my $rnd (0..2) {

        my $fn = $assemblies[$rnd];

        # parse dnadiff output
        open my $in, '<', $fn;
        while (my $line = <$in>) {
            chomp $line;
            if ($line =~ /^AlignedBases\s+(\d+)\(([\d\.]+)\%\)/) {
                $ref->{aln_len}->[$rnd] = $1;
                $ref->{reference_coverage}->[$rnd] = $2;
            }
            for (qw/
                AvgIdentity
                Relocations
                Translocations
                Inversions
                Insertions
                InsertionSum
                TotalSNPs
                TotalIndels
            /) {
                if ($line =~ /^$_\s+\S+\s+([\d\.]+)/) {
                    my $tag = $_;
                    my $val = $1;
                    $tag =~ s/(?<=[a-z])([A-Z])/_$1/g;
                    $tag = lc $tag;
                    $ref->{$tag}->[$rnd] //= $val;
                }
            }
        }

        $ref->{SNPs_per_kb}->[$rnd] = sprintf(
            "%0.2f",
            $ref->{total_snps}->[$rnd]/$ref->{aln_len}->[$rnd]*1000,
        );
        $ref->{indels_per_kb}->[$rnd] = sprintf(
            "%0.2f",
            $ref->{total_indels}->[$rnd]/$ref->{aln_len}->[$rnd]*1000,
        );

        # These count 'ends', and we want event counts
        for (qw/translocations inversions relocations/) {
            $ref->{$_}->[$rnd] /= 2;
        }

        # Remove differences in circular start from relocation count
        $ref->{relocations}->[$rnd] = $ref->{relocations}->[$rnd]-2;
        $ref->{relocations}->[$rnd] = $ref->{relocations}->[$rnd] > 0
            ? $ref->{relocations}->[$rnd]
            : 0;

    }

}

sub parse_busco {

    my ($ref) = @_;

    my @assemblies = qw{
        busco/full_table_oriented.tsv
        busco/full_table_oriented.polished.tsv
        busco/full_table_oriented.polished.polished.tsv
    };

    for my $rnd (0..2) {

        my $fn = $assemblies[$rnd];

        # parse BUSCO output
        my %counts;
        open my $in, '<', $fn;
        while (my $line = <$in>) {

            next if ($line =~ /^#/);

            chomp $line;
            my @fields = split "\t", $line;
            ++$counts{lc($fields[1])};
            ++$counts{total};
        
        }

        $ref->{BUSCO_complete}->[$rnd]   = sprintf(
            "%0.2f",
            ($counts{complete}   // 0) / $counts{total},
        );
        $ref->{BUSCO_fragmented}->[$rnd]   = sprintf(
            "%0.2f",
            ($counts{fragmented}   // 0) / $counts{total},
        );
        $ref->{BUSCO_missing}->[$rnd]   = sprintf(
            "%0.2f",
            ($counts{missing}   // 0) / $counts{total},
        );

    }

}

sub parse_traces {

    my ($ref) = @_;

    if (! -e 'trace.txt') {
        warn "trace file missing\n";
        $ref->{'time(min)'}      = 'NA';
        $ref->{'max_memory(GB)'} = 'NA';
        return;
    }

    open my $in, '<', 'trace.txt';
    my $h = <$in>;
    chomp $h;
    while (my $line = <$in>) {
        chomp $line;
        my @fields = split "\t", $line;
        my $proc = $fields[3];
        $proc =~ s/\s.+//g;
        next if (any {$proc eq $_} @skip_processes);
        my $dur = parse_duration($fields[8]);
        my $mem = parse_mem($fields[10]);
        if ($proc =~ /^(?:nanopolish|vcf2fasta)/) {
            $ref->{'time(min)'}->[2] += $dur;
            $ref->{'max_memory(GB)'}->[2] = $mem
                if ($mem > ($ref->{'max_memory(GB)'}->[2] // 0));
            if ($proc !~ /_2$/) {
                $ref->{'time(min)'}->[1] += $dur;
                $ref->{'max_memory(GB)'}->[1] = $mem
                    if ($mem > ($ref->{'max_memory(GB)'}->[1] // 0));
            }
        }
        else {
            for my $rnd (0..2) {
                $ref->{'time(min)'}->[$rnd] += $dur;
                $ref->{'max_memory(GB)'}->[$rnd] = $mem
                    if ($mem > ($ref->{'max_memory(GB)'}->[$rnd] // 0));
            }
        }
    }

    for my $rnd (0..2) {
        $ref->{'time(min)'}->[$rnd] = sprintf '%0.1f',
            $ref->{'time(min)'}->[$rnd];
        $ref->{'max_memory(GB)'}->[$rnd] = sprintf '%0.1f',
            $ref->{'max_memory(GB)'}->[$rnd];
    }

}

sub parse_lofreq {

    my ($ref) = @_;

    my @assemblies = qw{
        lofreq/calc_accuracy_oriented.tsv
        lofreq/calc_accuracy_oriented.polished.tsv
        lofreq/calc_accuracy_oriented.polished.polished.tsv
    };

    for my $rnd (0..2) {

        my $fn = $assemblies[$rnd];

        open my $in, '<', $fn;
        my $line = <$in>;
        chomp $line;
        my @fields = split "\t", $line;
        close $in;

        $ref->{avg_identity}->[$rnd]  = sprintf "%0.2f", $fields[0] * 100;
        $ref->{SNPs_per_kb}->[$rnd]   = sprintf "%0.2f", $fields[1];
        $ref->{indels_per_kb}->[$rnd] = sprintf "%0.2f", $fields[2];

    }

}
