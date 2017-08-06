=head1 NAME

IndexBamByRead.pm - sort, index and retrieve BAM records by read ID

=head1 VERSION

version 0.1

=head1 SYNOPSIS

 use IndexBamByRead;
 
 my $bam = 'example.bam';
 #sort by Read ID 
 IndexBamByRead::sort_bam($bam, 'rid_sorted');

 #create index (.ibbr) file
 IndexBamByRead::index_bam('rid_sorted.bam');
 
 #read index into memory
 my %index = IndexBamByRead::read_index('rid_sorted.bam.ibbr')

 #retrieve by read ID
 my $bam = Bio::DB::Bam->open('rid_sorted.bam');
 my @reads = IndexBamByRead::get_by_id($bam, 'READ_ID' \%index);

=cut

package IndexBamByRead;
use strict;
use warnings;
use IO::Uncompress::Gunzip qw/ gunzip $GunzipError /;
use IO::Compress::Gzip qw/ gzip $GzipError / ;
use Fcntl ':seek';
use File::Copy;
use File::Temp qw/ tempfile /;
use Data::Dumper;
use Carp;
use Bio::DB::Sam;
our $VERSION = 0.1;

=head1 FUNCTIONS

=over 8

=item B<sort_bam>

Sort BAM file ascibetically by read ID. First argument is filename of your 
input BAM file. 

=over 12

=item B<Optional args:>

=over 16

=item B<output>:

Prefix for output filename. Defaults to basename of input plus '_rid_sorted' 
(e.g. for input.bam it will default to input_rid_sorted.bam).

=item B<max_mem>:

Amount of memory to use. Default=500 (MB).

=back

=item B<Examples:>

 #output filename may be specified
 IndexBamByRead::sort_bam($bam, 'rid_sorted');

 # if output prefix is not specified, for input.bam it will default to 
 # input_rid_sorted.bam and will be returned
 my $sorted = IndexBamByRead::sort_bam($bam);

=back

=cut 

sub sort_bam{
    my $bam = shift;
    my %args = @_;
    my $out = $args{output};
    if (not $out){
        ($out = $bam) =~ s/\.bam/_rid_sorted/i;
    }
    my @args  = (1, $bam, $out);
    if ($args{max_mem}){
        push @args, $args{max_mem};
    }
    Bio::DB::Bam->sort_core(@args);
    return "$out.bam";
}


sub _cmp_id{
    my $r1 = shift;
    my $r2 = shift;
    return 0 if $r1 eq $r2;
    my @s1 = split(/[\.\:]/, $r1);
    my @s2 = split(/[\.\:]/, $r2);
    for (my $i = 0; $i < @s1; $i++){
        my $cmp = 0;
        if ($s1[$i] =~ /^\d+$/ and $s2[$i] =~ /^\d+$/){
            $cmp = $s1[$i] <=> $s2[$i];
        }else{
            $cmp = $s1[$i] cmp $s2[$i];
        }
        if ($cmp){
            return $cmp;
        }
    }
    return 0;
}

=item B<index_bam>

Create a bam index. The first argument must be a BAM file sorted ascibetically 
by read ID.

=over 12

=item B<Optional args:>

=over 16

=item B<index_name>:

Name for the index file (defaults to the input name + '.ibbr').

=item B<records_per_chunk>:

Number of records spanning each indexed segment. Higher values will result in a
larger index size (taking up more memory when read into memory in the get_by_id
method) but potentially faster read retrieval. Default=50,000.

 IndexBamByRead::index_bam('rid_sorted.bam');

=back
=cut

sub index_bam{
    local $Data::Dumper::Terse = 1;
    local $Data::Dumper::Useqq = 1;
    my $bam = shift;
    my %args = @_;
    if (not $args{index_name}){
        $args{index_name} = $bam . '.ibbr';
    } 
    my (undef, $tmp_index) = tempfile
    ( 
        "tmp_ibbrXXXXX", 
        UNLINK => 1,
        TMPDIR => 1 
    );
    my $gz = new IO::Compress::Gzip $tmp_index
        or croak "IO::Compress::Gzip failed to write to temporary file for ".
        "index: $GzipError\n";
    my $chunk_size = $args{records_per_chunk} || 50000;
    my $b = Bio::DB::Bam->open($bam);
    my $n = 0;
    my $header = $b->header;#required to call header before reading alignments
    my %idx = ();
    my $pos = $b->tell();
    $idx{start} = $pos;
    $idx{chunk_size} = $chunk_size; 
    my $prev = undef;
    my $prev_pos = undef;
    while (my $align = $b->read1){
        my $this_read =  $align->qname;
        if (defined $prev){
            my $cmp = IndexBamByRead::_cmp_id($prev, $this_read);
            if ($cmp > 0){
                croak("Reads are unsorted - cannot index $bam ");
            }
        }
        if ($n % $chunk_size == 0){
            $idx{pos}->{$this_read} = $pos;
            push @{$idx{ids}}, $this_read; #keep a sorted array of indexed IDs
        }
        $prev_pos = $pos;
        $pos = $b->tell();
        $n++;
        $prev = $this_read;
    }
    $idx{pos}->{$prev} = $prev_pos;
    push @{$idx{ids}}, $prev;
    print $gz Dumper \%idx;
    close $gz;
    move($tmp_index, $args{index_name}) or croak 
      "Could not create index '$args{index_name}' from temporary index file: $! ";
    return $args{index_name}; 
}

=item B<read_index>

Read an .ibbr bam index. If 'bam' argument is given, will create index if it 
does not exist.

=over 12

=item B<Args:>

=over 16

=item B<bam>:

Name of bam file. If used, the index is assumed to be the bam name + '.ibbr'.

=item B<index>:

Name of the index file. Required if not specifying 'bam' argument.

=item B<records_per_chunk>:

Number of records spanning each indexed segment. Ignored if index exists, but 
will be passed to 'index_bam' method if it does not. Default=50,000

 my %index = IndexBamByRead::read_index('rid_sorted.bam');

=back
=cut 

sub read_index{
    my %args = @_;
    if (not $args{index}){
        croak "Either 'index' or 'bam' argument is required " if not $args{bam};
        $args{index} = $args{bam} . '.ibbr';
    }
    my $chunk_size = $args{records_per_chunk} || 50000;
    if (not -e $args{index}){
        croak "$args{index} does not exist! " if not $args{bam};
        index_bam($args{bam}, index_name => $args{index}, $chunk_size);
    }
    if ($args{bam}){
        carp "\nWARNING: index $args{index} is older than $args{bam} " 
          if (-M $args{bam}) < (-M $args{index}); 
    }
    my $block_dump;
    my $z = new IO::Uncompress::Gunzip $args{index}
        or die "gunzip failed to read index $args{index}: $GunzipError\n";
    {
        local $/;
        $block_dump = <$z>;
    }
    close ($z);
    my %i = %{ eval $block_dump };
    return %i;
}
 


=item B<get_by_id>

Retrieve read (or read pair) by read ID. The first argument must be a 
Bio::DB::BAM object from an indexed file. The second argument must be the 
read ID to retrieve. The third argument must be a reference to an index hash 
as retrieved by the IndexBamByRead::read_index method.

If no matching read can be found returns an empty array.

 my $bam = Bio::DB::Bam->open('rid_sorted.bam');
 my @matches = IndexBamByRead::get_by_id($bam, 'READ_ID', \%index);

=cut

sub get_by_id{
    my ($bam, $id, $idx) = @_;
    if (_cmp_id($id, $idx->{ids}->[0]) < 0){ #not in VCF - lt first read ID
        return ();
    }
    if (_cmp_id($id, $idx->{ids}->[-1]) > 0){ #not in VCF - gt last read ID
        return ();
    }
    my @matches = (); 
    my ($before, $after) = _get_nearest_indices($id, $idx);
    if ($before == -1){ #not in bounds = shouldn't happen after initial checks?
        return ();
    }
    if ($before == $after){
        #get read at this index and look up and down for same read id
        push @matches, _search_either_side($bam, $id, $idx, $before);
    }else{
        #look between the two indices for any matching reads
        my $start = $idx->{pos}->{$idx->{ids}->[$before]};
        my $stop = $idx->{pos}->{$idx->{ids}->[$after]};
        push @matches, _get_matching($bam, $id, $start, $stop);
    }
    return @matches;
}

sub _search_either_side{
    my ($bam, $id, $idx, $i) = @_;
    my @matches = ();
    push @matches, _get_read_at_offset($bam, $idx->{$id});
    if ($i > 0){
        my $start = $idx->{pos}->{$idx->{ids}->[$i - 1]};
        my $stop = $idx->{pos}->{$idx->{ids}->[$i]};
        my @reads_before = _get_reads_in_region($bam, $start, $stop);
        for (my $j = @reads_before -1; $j >= 0; $j--){
            if ($reads_before[$j]->qname eq $id){
                unshift @matches, $reads_before[$j];
            }
        }
    }
    if ($i < @{$idx->{ids}} - 1){
        my $start = $idx->{pos}->{$idx->{ids}->[$i]};
        my $stop = $idx->{pos}->{$idx->{ids}->[$i + 1]};
        my @reads_after = _get_reads_in_region($bam, $start, $stop);
        for (my $j = 0; $j < @reads_after; $j++){
            if ($reads_after[$j]->qname eq $id){
                push @matches, $reads_after[$j];
            }
        }
    }
    return @matches;
}

sub _get_matching{
    my ($bam, $id, $start, $stop) = @_;
    my @reads = _get_reads_in_region($bam, $start, $stop);
    my $hit = _binsearch_reads(\@reads, $id);
    return () if $hit < 0; #shouldn't happen
    my @matches = ($reads[$hit]);
    for (my $i = $hit - 1; $i >= 0; $i--){
        if ($reads[$i]->qname eq $id){
            unshift @matches, $reads[$i];
        }else{
            last;
        }
    }
    for (my $i = $hit + 1; $i < @reads; $i++){
        if ($reads[$i]->qname eq $id){
            push @matches, $reads[$i];
        }else{
            last;
        }
    }
    return @matches;
}

sub _binsearch_reads{
    my $r = shift;
    my $id = shift;
    my $l = 0;
    my $u = @$r;
    while ( $l <= $u ) {
        my $i = int( ( $u + $l ) / 2 );
        my $rid = $r->[$i]->qname;
        my $cmp = _cmp_id($id, $rid);
        if ($cmp == 0){ #exact match
            return $i;
        }elsif($cmp < 0){ #id is less than $rid
            $u = $i -1;
        }else{
            $l = $i + 1;
        }
    }
    return -1;
} 

sub _get_reads_in_region{
    my ($bam, $start, $stop) = @_;
    $bam->header;
    $bam->seek($start, SEEK_SET); #this is not working - bug in Bio::DB::Bam?
    my @reads = ();
    push @reads, $bam->read1;
    my $pos = $bam->tell();
    while ($pos < $stop){
        push @reads, $bam->read1;
        $pos = $bam->tell();
    }
    return @reads;
}

sub _get_read_at_offset{
    my ($bam, $pos) = @_;
    $bam->seek($pos, SEEK_SET);
    return $bam->read1;
}

sub _get_nearest_indices{
    my ($read, $idx) = @_;
    #binary search of $idx->{ids}
    my $l = 0;
    my $u = @{$idx->{ids}} - 1;
    my $upper = $u;
    while ( $l <= $u ) {
        my $i = int( ( $u + $l ) / 2 );
        my $iid = $idx->{ids}->[$i];
        my $cmp = _cmp_id($read, $iid);
        if ($cmp == 0){ #exact match
            return ($iid, $iid);
        }elsif($cmp < 0){ #read ID is less than iid
            my $pre = $idx->{ids}->[$i - 1];
            if (_cmp_id($read, $pre) > 0){ #read ID is greater than prev in index
                return ($i-1 , $i);
            }
            $u = $i - 1;
        }else{ #read ID is greater than iid
            my $nxt = $idx->{ids}->[$i + 1];
            if (_cmp_id($read, $nxt) < 0){ #read ID is less than next in index
                return ($i, $i+1);
            }
            $l = $i + 1;
        }
    }
    return (-1, -1);
}

=back

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

MIT License

Copyright (c) 2017 David A. Parry

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=cut

1;
