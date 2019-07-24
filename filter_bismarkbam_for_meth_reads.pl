#!/usr/bin/perl
 
use strict;
use warnings;
 
use Bio::DB::Sam;
#use Data::Dumper;

# Parses Bismark BAM and retains only the reads where at least one alignment (mate) has methylated CpG.
# 'skip_flanking_calls' is an integer specifying the number of bases to ignore on either side of a read

my $usage = "Usage $0 <infile.bam> <outfile.bam>\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
#my $skip_left_flanking_calls = shift or die $usage;

die "Input file does not exist." unless -e $infile;

my $in_bam = Bio::DB::Bam->open("$infile" , 'r');

# die unless $infile =~ /([^\/]+)\.bam$/;
# my $outfile = "${1}.CpGfilt.bam";

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
while (my $alignment = $in_bam->read1()){
    my $readname = $alignment->qname;
    my $xm_tag = $alignment->get_tag_values('XM');

#    my $untrimmed_xm_len = length($xm_tag);
#    $xm_tag = substr ($xm_tag, $skip_flanking_calls, -$skip_flanking_calls); 

    # save the read id if it contains a methylated CpG
    if ($xm_tag =~ /Z/){
        $ids_to_keep{$readname}++;
    }

    # save all alignment
    push(@alignments, $alignment);
}

# output alignment if the id was stored
foreach my $alignment (@alignments){
    my $readname = $alignment->qname;
    if (exists $ids_to_keep{$readname}){
        my $bytes = $out_bam->write1($alignment);
    }
}



