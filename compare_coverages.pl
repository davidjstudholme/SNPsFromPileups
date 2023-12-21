#!/usr/bin/perl -w

use strict;
use warnings;

my $usage = "Usage: $0 <annotation file> <csv filenames>";

my $max_breadth = 0.00;
my $min_breadth = 1.00;

#my $max_breadth = 1.00;
#my $min_breadth = 0.00;

### Parse annotation GFF file
my $annotation_file = shift or die "$usage\n";
die "$usage\n" unless $annotation_file =~ m/\.gff$/i; 
my %id2location;
my %id2start;
open (FILE, "<$annotation_file ") or die "Failed to open file '$annotation_file'\n$!\n";
while (<FILE>) {
    unless (m/^#/) {
	chomp;
	s/\"//g;
	my @fields = split /\t/;
	my ($contig_id,
	    $source,
	    $feature,
	    $start,
	    $end,
	    $score,
	    $frame,
	    $strand,
	    $attributes,
	    ) = @fields;
	if ($contig_id =~ m/gb\|(\S+)\|/){
	    $contig_id=$1;
	}

	if ( $feature eq 'CDS' and $attributes =~ m/protein_id=([\w\d\.]+)/) {
	    my $id = $1;
	    $id2location{$id} = "$contig_id:$start..$end";
	    $id2start{$id} = "$start";
	    #warn "\$id2location{'$id'} = $contig_id:$start..$end\n";
	}
    }
}
close FILE;



### Now read the coverageBed output .csv files
my %gene2breadths;
@ARGV or die "$usage\n";
foreach my $file (@ARGV) {
    
    open(FILE, "<$file") or die "Failed to open file '$file'\n";
    while (<FILE>) {
	chomp;
	my @fields = split /\t/, $_;
	if (@fields == 13) {
	    my ($contig_id,
		$b,
		$c,
		$start,
		$stop,
		$f,
		$g,
		$h,
		$gene,
		$j,
		$k,
		$l,
		$breadth,
		)
		= @fields;

	    if ($c eq 'CDS' and $gene =~ m/protein_id=([\w\d\.]+)/) {
		my $id = $1;

		$gene2breadths{$id}{$file} = $breadth if $breadth =~ m/^[\d+\.]+$/;
		#warn "\$gene2breadths{'$id'}{$file} = $breadth\n";
	    }
	} else {
	    warn "Ignoring: $_\n";
	}
    }
    close FILE;
}

print "\"Gene\"";
print "\t";
print "\"Location\"";
print "\t";
print "\"Start\"";
foreach my $file (@ARGV) {
    my $heading = $file;

    $heading =~ s/.aln.sorted.bam.csv//;
    $heading =~ s/_filtered.subsampled.bam.csv//;
    $heading =~ s/_filtered.sorted.bam.csv//;
    
    print "\t\"$heading\"";
}
print "\n";

foreach my $gene (sort keys %gene2breadths) {
    my @breadths;
    foreach my $file (@ARGV) {
	push @breadths, $gene2breadths{$gene}{$file};
    }

    my $no_difference = 1;
    while (my $breadth1 = shift @breadths) {
	foreach my $breadth2 (@breadths) {
	    if ( ($breadth1 >= $min_breadth and $breadth2 <= $max_breadth) or
		 ($breadth2 >= $min_breadth and $breadth1 <= $max_breadth) ) {
		$no_difference = 0;
	    }
	}
    }
	    
    unless ($no_difference ) {
	
	### Look up location of gene
	my $location = $id2location{$gene};
	my $start = $id2start{$gene};
	die "Couldn't find location of '$gene'" unless defined $location;

	my @breadths;
	foreach my $file (@ARGV) {
	    push @breadths, $gene2breadths{$gene}{$file};
	}

	### clean up the gene name
	if ($gene =~ m/name=(\S+?);.*/i) {
	    
	    #warn "gene= '$1'\n";
	    
	    my $cleaned_gene = $1;
	    if ($gene =~ m/product=(.*?);.*/i) {
		$cleaned_gene .= " $1";
	    } else {
		$cleaned_gene = 'exclude';
	    }
	    $gene = $cleaned_gene;
	} else {
	    #"\$gene does not contain gene= = $gene\n";
	}
	$gene =~s/;Name=/ /gi;
	$gene =~ s/;Ontology_term.*//gi;
	$gene =~ s/ID=//gi;


	unless ($gene eq 'exclude') {
	    
	    print "\"$gene\"";
	    print "\t";
	    print "\"$location\"";
	    print "\t";
	    print "\"$start\"";
	    foreach my $breadth (@breadths) {
		print "\t$breadth";
	    }
	    print "\n";
	}
    }
}
