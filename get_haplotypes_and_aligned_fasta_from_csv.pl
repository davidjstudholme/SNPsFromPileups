#!/usr/bin/perl

use strict;
use warnings ;

my $required_consensus = 0.95;

my %genome2seq;

my $csv_file = shift or die "Usage: $0 <tab-delimited .csv file>\n" ;
open(FILE, "<$csv_file") or die "Failed to open file '$csv_file' for reading\n$!\n";

### Read the heading line of the .csv file 
my $header_line = <FILE>;
chomp $header_line;
my @headings = split /\t/, $header_line;
foreach my $expected_heading ('chromosome', 'pos', 'ref', 'alt') {
    my $heading = shift @headings;
    die "Expected '$expected_heading' but got '$heading'" unless $expected_heading eq $heading;
}
my @genomes = @headings;

### Read the data from the .csv file
my $count_lines;
while (<FILE>) {
    chomp;
    my @fields = split /\t/;
    my $chromosome = shift @fields;
    my $pos = shift @fields;
    my $ref = shift @fields;
    my $alt = shift @fields;
    
    ### Are there any ambiguous nucleotides calls (Ns)?
    my $unambiguous = 1;
    foreach my $alt_freq (@fields) {
	if ($alt_freq > $required_consensus or
	    $alt_freq < (1 - $required_consensus) ) {
	    ### OK
	} else {
	    ### Ambiguity
	    $unambiguous = 0;
	}
    }

    if ($unambiguous) {
	foreach my $genome (@genomes) {
	    my $alt_freq = shift @fields;
	    if ($alt_freq < (1 - $required_consensus) ) {
		$genome2seq{$genome} .= $ref;
	    } 
	    elsif ( $alt_freq > $required_consensus) {
		$genome2seq{$genome} .= $alt;
	    } else {
		$genome2seq{$genome} .= 'N';
	    }
	}
        $count_lines++;
    }
}
close FILE;


### Write pseudosequences to FastA file
my $pseudosequence_fasta_file = "$csv_file.haplotype-as-fasta.fna";
open(PSEUDOSEQ, ">$pseudosequence_fasta_file ") or die $!;
foreach my $genome (sort keys %genome2seq) {
    my $pseudoseq = $genome2seq{$genome};
    print PSEUDOSEQ ">$genome\n$pseudoseq\n";
}
close PSEUDOSEQ;
warn "Wrote pseudosequences to file '$pseudosequence_fasta_file'\n";
    
### Write pseudosequences to Nexus file
my $pseudosequence_nexus_file = "$csv_file.haplotype.nex";
my @taxa = sort keys %genome2seq;
my $taxa_count = scalar(@taxa);
open(PSEUDOSEQ, ">$pseudosequence_nexus_file ") or die $!;
print PSEUDOSEQ "#NEXUS\n\n";

### Print taxa block
print PSEUDOSEQ "BEGIN TAXA;\n";
print PSEUDOSEQ "DIMENSIONS NTAX=$taxa_count;\n\n";
print PSEUDOSEQ "TAXLABELS\n";
foreach my $taxon (@taxa) {
    print PSEUDOSEQ "$taxon\n";
}
print PSEUDOSEQ ";\n\n";
print PSEUDOSEQ "END;\n\n";

### Print characters block
my $nchar= $count_lines;
print PSEUDOSEQ "BEGIN CHARACTERS;\n";
print PSEUDOSEQ	"DIMENSIONS NCHAR=$nchar;\n";
print PSEUDOSEQ	"FORMAT DATATYPE=DNA MISSING=N GAP=- ;\n";
print PSEUDOSEQ	"MATRIX\n";
foreach my $genome (sort keys %genome2seq) {
    my $pseudoseq = $genome2seq{$genome};
    print PSEUDOSEQ "$genome\t$pseudoseq\n";
}
print PSEUDOSEQ ";\n";
print PSEUDOSEQ "END;\n\n";

if (0) {
    ### Print traits block
    print PSEUDOSEQ "Begin Traits;\n";
    print PSEUDOSEQ "Dimensions NTraits=2;\n";
    print PSEUDOSEQ "Format labels=yes missing=? separator=Comma;\n";
    print PSEUDOSEQ "TraitLabels trait1 trait2;\n";
    print PSEUDOSEQ "Matrix\n";
    
    foreach my $genome (keys %genome2seq) {
	my ($trait1, $trait2) = (0, 0);
	print PSEUDOSEQ "$genome\t$trait1,$trait2\n";
    }
    
    print PSEUDOSEQ ";\n";
    print PSEUDOSEQ "End;\n";
}

close PSEUDOSEQ;
warn "Wrote pseudosequences to file '$pseudosequence_nexus_file'\n";

