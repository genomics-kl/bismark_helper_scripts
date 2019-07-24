#!/usr/bin/perl
 
use strict;
use warnings;
 
use Bio::DB::Sam;
#use Data::Dumper;

# Parses Bismark BAM and retains only the reads where at least one alignment (mate) has GpC methylation based on the 'maowar' file from NOMe_filter.

my $usage = "Usage: $0 <original bismark bam> <outfile.bam> <manOwar.txt.gz>\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $manowar_file = shift or die $usage;

die "Input file does not exist." unless (-e $infile && -e $manowar_file);
die unless $infile =~ /\.bam$/;
die unless $manowar_file =~ /any_C_context_.*manOwar\.txt\.gz$/;

my $in_bam = Bio::DB::Bam->open("$infile" , 'r');
open(my $manowar_fh, "gunzip -c $manowar_file |") or die $!;

die "Output file exists." if -e $outfile;

my $out_bam = Bio::DB::Bam->open("$outfile" , 'w');
 
#my @targets = $sam->seq_ids;

my %ids_to_keep;
my @alignments;

# save the header from the input BAM
my $header = $in_bam->header();
# Given an open BAM file, return a Bio::DB::Bam::Header object containing information about the reference sequence(s). Note that you must invoke header() at least once before calling read1().

# output the header to the output BAM.
my $status_code = $out_bam->header_write($header);
# Given a Bio::DB::Bam::Header object and a BAM file opened in write mode, write the header to the file. If the write fails the process will be terminated at the C layer. The result code is (currently) always zero.

# Read one alignment from the BAM file and return it as a Bio::DB::Bam::Alignment object. Note that you must invoke header() at least once before calling read1().
my $total_Z_calls = 0;
my $Z_read_ct = 0;
while (my $alignment = $in_bam->read1()){
    my $readname = $alignment->qname;
    die if $readname =~ /R[12]$/; # must be using original bismark BAM with original names

    my $xm_tag = $alignment->get_tag_values('XM');

    # count number of CpG calls per read
    my $read_Z_calls = ($xm_tag =~ tr/Z//);

    # add read CpG count to running count.
    $total_Z_calls += $read_Z_calls;

    # don't save the read id if it contains a methylated CpG because this script filters only for GpC meth. Do keep count for output summary stats.
    if ($read_Z_calls > 0){
#        $ids_to_keep{$readname}++;
        $Z_read_ct++;
    }

    # save all alignment
    push(@alignments, $alignment);
}

my $num_total_reads = scalar(@alignments);
my $cpg_meth_read_rate = $Z_read_ct / $num_total_reads;# each mate/alignment is treated as a separate count

print "###\n";
print "$0 summary stats:\n\n";
print "Total reads\t$num_total_reads\n";
print "Total Z calls\t$total_Z_calls\n";
print "Total Z reads\t$Z_read_ct\n";
print "%reads with CpG meth\t$cpg_meth_read_rate\n";

# remove header of manowar table
my $manowar_header = <$manowar_fh>;
my $expected_manowar_format = '^ReadID\tChr\tStart\tEnd\tmeth_CG\tunmeth_CG\tmeth_GC\tunmeth_GC$';
die 'manowar file format is incorrect.' unless $manowar_header =~ /$expected_manowar_format/;

# Look for GpC methylation according to manowar file
# assumes read names were renamed with _R1 and _R2 suffices
my $total_gpc_calls = 0;
my $gpc_read_ct = 0;
while (<$manowar_fh>){
    chomp;
    my @line = split("\t", $_);
    my $meth_gc_ct = $line[6];
    my $readname = $line[0];

    die unless $readname =~ /_R[12]$/;# this allows the bismark original bam to be filtered by read name, since that does not have the R1 R2 suffices.
    $readname =~ s/_R[12]$//;

    die unless $meth_gc_ct =~ /^\d+$/;

    # add to running total
    $total_gpc_calls += $meth_gc_ct;

    # keep read if at least one GpC call
    if ($meth_gc_ct > 0){
        $ids_to_keep{$readname}++;
        $gpc_read_ct++;
    }
}
close $manowar_fh;

my $gpc_meth_read_rate = $gpc_read_ct / $num_total_reads;
#my $num_reads_w_cg_or_gc_meth = (scalar keys %ids_to_keep) * 2;

print "Total GpC calls (counts only GCC, GCT, GCA)\t$total_gpc_calls\n";
print "Total GpC reads (counts only GCC, GCT, GCA)\t$gpc_read_ct\n";
print "%reads with GpC meth\t$gpc_meth_read_rate\n";
#print "Expected number of reads/alignments to be retained\t$num_reads_w_cg_or_gc_meth\n";

# output alignment if the id was stored
foreach my $alignment (@alignments){
    my $readname = $alignment->qname;
    if (exists $ids_to_keep{$readname}){
        my $bytes = $out_bam->write1($alignment);
    }
}



