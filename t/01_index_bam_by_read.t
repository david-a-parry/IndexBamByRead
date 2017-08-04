use strict;
use warnings;
use Test::More;# tests => 2;
use FindBin qw($RealBin $Script);
use lib "$RealBin/../";
BEGIN 
{ 
    use_ok("IndexBamByRead");
}
my $n_tests = 1;

my $bam = "$RealBin/test_data/test.bam";
(my $sort_out = $bam) =~ s/\.bam/_rid_sorted.bam/;
IndexBamByRead::sort_bam($bam);
ok(-e $sort_out, "Produced expected sorted BAM");
$n_tests++;
is(1010, get_read_count($sort_out), "Sort output has correct number of reads");
$n_tests++;
is(0, reads_are_sorted($bam), "Input reads are not sorted be ID");
$n_tests++;
ok(reads_are_sorted($sort_out), "Sort output reads are sorted be ID");
$n_tests++;
ok(IndexBamByRead::index_bam($sort_out, records_per_chunk => 100), "Index sorted BAM");
$n_tests++;
ok(my %idx = IndexBamByRead::read_index(bam => $sort_out), "Read index");
$n_tests++;
ok(find_reads($sort_out, \%idx), "Can find reads by ID");

done_testing($n_tests);

#################################################
sub find_reads{
    my $f = shift;
    my $idx = shift;
    my $bam = Bio::DB::Bam->open($f);
    my $header = $bam->header;#required to call header before reading alignments
    my @targets = qw/
        SRR043348.7107241
        SRR043348.14284382
        SRR043354.3858155
        SRR043354.12939592
        SRR043366.2208972
        SRR043372.1247100
        SRR043372.11318637
        SRR043378.2924881
        SRR043378.13722892
        SRR043386.7232182
        SRR043386.16375485
        SRR043396.6476248
        SRR043396.16790502
    /; 
    foreach my $t (@targets){
        my @matches = IndexBamByRead::get_by_id($bam, $t, $idx);
        my $m = grep {$_->qname eq $t} @matches; 
        if ($m != 2){
            return 0;
        }
    }
    return 1;
}

#################################################
sub reads_are_sorted{
    my $f = shift;
    my $bam = Bio::DB::Bam->open($f);
    my $n = 0;
    my $header = $bam->header;#required to call header before reading alignments
    my $prev_read = undef;
    while (my $align = $bam->read1){
        my $this_read = $align->qname;
        $n++;
        if (defined $prev_read){
            my $cmp = IndexBamByRead::_cmp_id($prev_read, $this_read);
            if ($cmp > 0){
                return 0;
            }
        }
        $prev_read = $this_read;
    }
    return $n;#require at least one read to have been read
}
#################################################
sub get_read_count{
    my $f = shift;
    my $bam = Bio::DB::Bam->open($f);
    my $n = 0;
    my $header = $bam->header;#required to call header before reading alignments
    while (my $align = $bam->read1){
        $n++;
    }
    return $n;
}
