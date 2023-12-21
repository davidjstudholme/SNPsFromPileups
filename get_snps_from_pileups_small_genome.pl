#!/usr/bin/perl -w

use strict;
use warnings;
my $usage = "Usage: $0  <minimum depth of read-coverage> <vcf_files> <pileup files>";

my $min_depth = shift or die "$usage\n";

die $usage unless $min_depth =~ m/^\d+$/;

warn "Only considering sites where depth >= $min_depth in all pileup files\n";

### Read locations of candidate SNPs from VCF file to reduce search space
my $no_search_space = 1;
my %search_space;
foreach my $vcf_file (@ARGV) {
    if ($vcf_file =~ m/\.vcf$/) {
	$no_search_space = 0;
	open (FILE, "<$vcf_file") and
	    warn "Parsing VCF file '$vcf_file'\n" or
	    die "Failed to open file '$vcf_file'\n$!\n";
	while(<FILE>) {
	    chomp;
	    if (m/^\#/ or m/^$/) {
		### ignore comment line
	    } else {
		### Read the variant
		my @fields = split /\t/;
		my ($chrom, $pos, $variant_id, $ref, $alt, $qual, $filter, $info, $format, $na00001, $na00002, $na00003) = @fields;
		my ($id,$position,$from,$to) = ($chrom, $pos, lc($ref), uc($alt));
		
		$search_space{$chrom}{$pos}++;
		#warn "\$search_space{$chrom}{$pos}++\n";
		
	    }
	}
	close FILE;
	warn "finished reading VCF file '$vcf_file'\n";
    }
}
my $search_space_count = 0;
foreach my $chrom (keys %search_space) {
    $search_space_count += keys %{ $search_space{$chrom}};
}
warn "There are $search_space_count candidate SNP sites in the search space\n";

### Get a list of input pileup files
my %filehandles;
foreach my $infile (@ARGV) {
    if ($infile =~ m/pileup$/i) {
	$filehandles{$infile}++;
    }
}

my %chrom_to_pos_to_file_to_base_to_freq;
my %chrom_to_pos_to_file_to_depth;
my %chrom_to_pos_to_refbase;

### Print header line
print "chromosome";
print "\t";
print "pos";
print "\t";
print "ref";
print "\t";
print "alt";
foreach my $infile (sort keys %filehandles) {
	my $heading = $infile;
	$heading =~ s/\.combined.*\.pileup//;
	$heading =~ s/\.versus.*\.pileup//;
	$heading =~ s/.*\///;
	print "\t";
	print "$heading";
}
print "\n";

### open the input pileup files                                                                                                                         
foreach my $infile (keys %filehandles) {
    
    open(my $fh, "<", $infile) and
            warn "Opened file < $infile\n" or
            die "Can't open < $infile: $!";
    
    while (defined (my $readline = <$fh>)) {
	#warn "\nRead a line from $infile: $readline";
	chomp $readline;
	my @fields = split /\t/, $readline;
	my ($chromosome, $pos, $ref_base, $reads, $alignment ) = @fields; 
	
	### Infer the genotype at this site
	if ((defined $search_space{$chromosome}{$pos} or $no_search_space) and
	    $reads and 
	    defined $alignment) {
	    #warn "fields: @fields\n";
	    my %base2freq = %{ get_genotype_from_alignment_string(
				   \$chromosome,
				   \$pos, 
				   \$ref_base,
				   \$reads, 
				   \$alignment) };
	    
	    $chrom_to_pos_to_refbase{$chromosome}{$pos} = lc($ref_base);
	    $chrom_to_pos_to_file_to_depth{$chromosome}{$pos}{$infile} = $reads;
	    foreach my $base (keys %base2freq) {
		$chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{$base} = $base2freq{$base};
		#warn "\$chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{$base} = $base2freq{$base}\n";
		
	    }
	}
    }

    close $fh;
    warn "\tFinished parsing file '$infile'\n";

    
}


### Now print the current data
my %insufficient_depth;
my %not_defined;

foreach my $chromosome (keys %chrom_to_pos_to_file_to_base_to_freq) {

    foreach my $key (sort keys %insufficient_depth) {
	warn "$key => $insufficient_depth{$key} insufficient depth\n";
    }
    foreach my $key (sort keys %not_defined) {
	warn "$key => $not_defined{$key} not defined\n";
    }
    
    
    warn "For each chromosome ($chromosome)";
    foreach my $pos (keys %{ $chrom_to_pos_to_file_to_base_to_freq{$chromosome} }) {
	#warn "\tExamining $chromosome: $pos\n";
	
	### Check whether data for this position is defined from all pileup infiles and is non-zero in at least one
	my $all_defined = 1;
	
	foreach my $infile (sort keys %filehandles) {
	    if(defined $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{'a'} and
	       defined $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{'c'} and
	       defined $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{'g'} and
	       defined $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{'t'}) {
		### Data are defined for this locus
		
	    } else {
		$all_defined = 0;
		$not_defined{$infile}++;
	    }
	}
	if ($all_defined) {
	    
	    ### Check that we have sufficient depth in all pileup infiles at this position
	    my $sufficient_depth = 1;
	    foreach my $infile (keys %filehandles) {
		my $depth = $chrom_to_pos_to_file_to_depth{$chromosome}{$pos}{$infile};
		if (defined $depth) {
		    #warn "depth for $chromosome: $pos $infile is $depth\n";
		} else {
		    die warn "depth for $chromosome: $pos $infile is UNDEFINED; this should never happen\n";
		}
		
		if ($depth < $min_depth) {
		    $sufficient_depth = 0;
		    #warn "Insufficient depth ($depth x) for $chromosome: $pos in $infile\n"; 
		    $insufficient_depth{$infile}++;
		}
		
	    }
	    
	    if ($sufficient_depth) {
		
		foreach my $base ('a', 'c', 'g', 't') {
		    
		    my $all_zero = 1;
		    
		    #warn "Examining $base at $chromosome: $pos (we have sufficient depth and all defined)\n";
		    my $print_line = '';	
		    $print_line .= "$chromosome";
		    $print_line .=  "\t";
		    $print_line .=  "$pos";
		    $print_line .= "\t";
		    $print_line .= $chrom_to_pos_to_refbase{$chromosome}{$pos};
		    $print_line .= "\t";
		    $print_line .=  "$base";
		    
		    foreach my $infile (sort keys %filehandles) {
			my $freq = $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{$base};
			unless (defined $freq) {
			    die "Frequency is UNDEFINED for $chromosome: $pos $infile $base; this should never happen";
			}
			$print_line .= "\t";
			$print_line .= "$freq";
			$all_zero = 0 if $freq > 0.3;
		    }
		    
		    ### Print the data
		    if (
			$base ne $chrom_to_pos_to_refbase{$chromosome}{$pos} and
			$all_zero == 0
			) {
			print "$print_line\n";
			#warn "$chromosome: $pos  $chrom_to_pos_to_refbase{$chromosome}{$pos} -> $base\n";
			my @pos_count = keys %{ $chrom_to_pos_to_file_to_base_to_freq{$chromosome} };
			my $pos_count = @pos_count;
			@pos_count = sort {$a<=>$b} @pos_count;
			#warn "There are $pos_count positions to consider on $chromosome: $pos_count[0] ... $pos_count[-1]\n";
		    }
		    
		    
		}
		
	    }    
	}
    }
}

foreach my $key (sort keys %insufficient_depth) {
    warn "$key => $insufficient_depth{$key} insufficient depth\n";
}

foreach my $key (sort keys %not_defined) {
    warn "$key => $not_defined{$key} not defined\n";
}


exit;


sub get_genotype_from_alignment_string{
	
    my $chromosome_ref = shift or die;
    my $chromosome = $$chromosome_ref;
    
    my $pos_ref = shift or die;
    my $pos = $$pos_ref;
    
    my $ref_base_ref = shift or die;
    my $ref_base= $$ref_base_ref;
    
    my $reads_ref = shift or die;
    my $reads = $$reads_ref;

    my $alignment_ref = shift or die;
    my $alignment = $$alignment_ref; 

    #warn "Alignment: '$alignment'\n";
	
    ### Remove special symbols for read ends and starts
    $alignment =~ s/\^.//gi; 
    $alignment =~ s/\$//gi; 
	
    ### Count the dots and commas
    my @agree_hits = ($alignment =~ m/[\,\.]/g );
    #warn "Done counting dots and commas\n";
    
    ### Count the deletions
    my @deletion;
    while  ($alignment =~ m/\-(\d+)/) {
	my $n = $1;
	if ( $alignment =~ m/(\-$n[ACGTNacgtn]{$n})/) {
	    my $deletion = $1;
	    #warn "\n$_\n$deletion\n";
	    push @deletion, $deletion;
	    $alignment =~ s/$deletion//;
	}
    }
    #warn "Done counting deletions\n";
    
    ### Count the insertions
    my @insertion;
    while  ($alignment =~ m/\+(\d+)/) {
	my $n = $1;
	if ( $alignment =~ m/(\+$n[ACGTNacgtn]{$n})/) {
	    my $insertion = "\\$1";
	    #warn "\n$_\n$insertion\n";
	    push @insertion, $insertion;
	    $alignment =~ s/$insertion//;
	}
    }
    #warn "Done counting insertions\n";
    
    ### Count substitutions
    my @a = ($alignment =~ m/a/gi );
    my @c = ($alignment =~ m/c/gi );	
    my @g = ($alignment =~ m/g/gi );
    my @t = ($alignment =~ m/t/gi );
    my @n = ($alignment =~ m/n/gi );
    my @star = ($alignment =~ m/\*/gi );
    #warn "Done counting substitutions\n";
    
    ### Do a sanity check that everything adds up!
    my $depth = @agree_hits + @a + @c + @g + @t + @n + @star;
    my $length = length($alignment);
    unless ($depth == $reads) {
	warn "$_\nlength($alignment)=$length\n$depth != $reads\n";
    }
    
    ### What are the frequencies of each base?
    my %base2freq;

    if ($reads) {
	
	$base2freq{a} = @a / $reads if @a or 1;
	$base2freq{c} = @c / $reads if @c or 1;
	$base2freq{g} = @g / $reads if @g or 1;
	$base2freq{t} = @t / $reads if @t or 1;
	$base2freq{ lc($ref_base) } = @agree_hits / $reads if @agree_hits;
	#$base2freq{'*'} = @star / $reads;
	#$base2freq{'ins'} = @insertion / $reads;
	#$base2freq{'del'} = @deletion / $reads;
	#$base2freq{'ref'} = @agree_hits / $reads;
    }
    
    return \%base2freq;
}
