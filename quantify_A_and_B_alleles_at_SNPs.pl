#!/usr/bin/perl -w

use strict;
use warnings;
my $usage = "Usage: $0 <chromosome, e.g. chr01>  <minimum depth of read-coverage>  <maximum depth of read-coverage>";

my $chr = shift or die "$usage\n";

my @vcf_files = (
    './filtered_vcf/M_balbisiana_S71.versus.musa_acuminata_v2_pseudochromosome.aln.sorted.bam.filtered.vcf',
    #'./filtered/M_balbisiana_S75.versus.musa_acuminata_v2_pseudochromosome.aln.sorted.bam.filtered.vcf'
    );

my @pileup_files_B = (
    "./pileup/balbisiana_IITA.hiseq.$chr.pileup",
    "./pileup/SRR6996488.$chr.pileup",
    "./pileup/SRR6996489.$chr.pileup",
    "./pileup/2001-1027.$chr.pileup", # Eden balbisiana
    );
my @pileup_files_A = (
    "./pileup/SRR6996490.$chr.pileup",
    "./pileup/SRR6996491.$chr.pileup",
    "./pileup/SRR6996492.$chr.pileup",
    "./pileup/SRR6996493.$chr.pileup",
    
    "./pileup/2012-1161.$chr.pileup",
    "./pileup/2012-1173.$chr.pileup",
    "./pileup/2012-1154.$chr.pileup",
    "./pileup/cavendish.hiseq.$chr.pileup",
    "./pileup/2012-1156.$chr.pileup",
#    "./pileup/1998-2307.$chr.pileup", # Pisang Mas 
    );
my @pileup_files_unknowns = (

    "./pileup/SRR6996486.$chr.pileup",
    "./pileup/SRR6996487.$chr.pileup",
    "./pileup/SRR6996488.$chr.pileup",
    "./pileup/SRR6996494.$chr.pileup",

#    "./pileup/2011-0950.$chr.pileup", # Congo 2
    "./pileup/gonja_manjaya.hiseq.$chr.pileup",
#    "./pileup/paradisiaca.1999-2846.$chr.pileup", # 'paradisiaca'

#    "./pileup/1999-0524.$chr.pileup", # textilis
    "./pileup/sukali_ndiizi.hiseq.$chr.pileup", 

    
    #"./pileup/2011-0952.novaseq.$chr.pileup", # One hand planty
    #"./pileup/2012-1152.novaseq.$chr.pileup", # Safet Velchi
    "./pileup/2012-1164.$chr.pileup", # Calypso
    
    "./pileup/kayinja.hiseq.$chr.pileup",
    
    "./pileup/2012-1152.$chr.pileup", # Safet Velchi
    
    "./pileup/2011-0952.$chr.pileup", # One Hand Planty
    
    #"./pileup/1999-0158.$chr.pileup", # trogladytarum
    #"./pileup/2012-1166.$chr.pileup", # velutina
    );

my $min_depth = shift or die "$usage\n";
my $max_depth = shift or die "$usage\n";

my $absence_max_freq = 0.01;

die $usage unless $min_depth =~ m/^\d+$/ and  $max_depth =~ m/^\d+$/;

warn "Only considering sites where depth >= $min_depth and <= $max_depth in all pileup files\n";

warn "\nConsidering these genomes as 'A': @pileup_files_A\n";
warn "\nConsidering these genomes as 'B': @pileup_files_B\n";
warn "\nConsidering these genomes as unknown: @pileup_files_unknowns\n\n";


### Read locations of candidate SNPs from VCF file; to reduce search space, we will examine only those sites
my $no_search_space = 1;
my %search_space;
foreach my $vcf_file (@vcf_files) {
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
    
### open the input pileup files
my %filehandles;
my %pileup_files_lists = ('A' => \@pileup_files_A , 
			'B' => \@pileup_files_B, 
			'unknowns' => \@pileup_files_unknowns
    );
foreach my $list_name (keys %pileup_files_lists) {
    my $list_ref = $pileup_files_lists{$list_name};
    my @list = @{$list_ref};
    foreach my $infile (@list) {
	open(my $fh, "<", $infile) and
	    warn "Opened file < $infile\n" or
	    die "Can't open < $infile: $!";
	$filehandles{$infile} = $fh;
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
	$heading =~ s/\.pileup//;
	$heading =~ s/.*\///;
	$heading =~ s/\.chr\d+$//;
	print "\t";
	print "$heading\_percentage_B";
	print "\t";
	print "$heading\_depth";
	
}
print "\n";

### Keep track of progress in each infile
my %current_position;
my @chromosomes_seen_so_far;
my %chromosomes_seen_so_far;
my $current_chromosome = '?';
my $least_progressed_file;

### Keep going until we have reached the end of all the infiles
my $eof = 0;
until ($eof) {
    
    ### Figure out which file we need to read from next
    my %pos_to_infile;
    foreach my $infile (keys %filehandles) {
	my $pos = $current_position{$infile}{$current_chromosome};
	#warn "Current position of $infile on $current_chromosome is $pos\n"; 
	$pos = 0 unless defined $pos;
	$pos_to_infile{$pos} = $infile;
    }
    my @ascending_pos = sort {$a<=>$b} keys %pos_to_infile;
    $least_progressed_file = $pos_to_infile{ $ascending_pos[0]  };
    #warn "Least progressed file is $least_progressed_file\n";
    
    ### Read a line from each file
    my $infile = $least_progressed_file;
    my $fh = $filehandles{$infile};
    if (defined (my $readline = <$fh>)) {
	#warn "\nRead a line from $infile: $readline";
	chomp $readline;
	my @fields = split /\t/, $readline;
	my ($chromosome, $pos, $ref_base, $reads, $alignment ) = @fields; 
	
	### Keep track of progress in this infile
	unless (defined $chromosomes_seen_so_far{$chromosome}) {
	    $chromosomes_seen_so_far{$chromosome}++;
	    push @chromosomes_seen_so_far, $chromosome;
	    $current_chromosome = $chromosome;
	}
	$current_position{$infile}{$chromosome} = $pos;
	
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
    
    ### Now print the current data
    foreach my $chromosome (keys %chrom_to_pos_to_file_to_base_to_freq) {
	
	#warn "For each chromosome ($chromosome)";
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
		}
	    }
	    if ($all_defined) {
		
		### Release some memory by deleting from memory everything up to the line that has been printed
		foreach my $i (keys %{$chrom_to_pos_to_file_to_base_to_freq{$chromosome}} ){
		    if ($i < $pos) {
			#warn "\tWe are currently at $chromosome: $pos\tDeleting $chromosome: $i\n";
			delete $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$i};
			delete $chrom_to_pos_to_file_to_depth{$chromosome}{$i};
			delete $chrom_to_pos_to_refbase{$chromosome}{$i};
		    }
		}
              
		### Delete from memory data on previous chromosomes
		my $previous = 1;
		foreach my $seen_chromosome (@chromosomes_seen_so_far) {
		    if ($chromosome eq $seen_chromosome) {
			$previous = 0;
		    } elsif ($previous and defined $chrom_to_pos_to_file_to_base_to_freq{$seen_chromosome}) {
			delete $chrom_to_pos_to_file_to_base_to_freq{$seen_chromosome};
			delete $chrom_to_pos_to_file_to_depth{$seen_chromosome};
			delete $chrom_to_pos_to_refbase{$seen_chromosome};
			warn "Finished with $seen_chromosome\n";
		    }
		}

  
		### Check that we have sufficient depth in all pileup infiles at this position
		my $sufficient_depth = 1;
		foreach my $infile (keys %filehandles) {
		    my $depth = $chrom_to_pos_to_file_to_depth{$chromosome}{$pos}{$infile};
		    if (defined $depth) {
			#warn "depth for $chromosome: $pos $infile is $depth\n";
		    } else {
			die warn "depth for $chromosome: $pos $infile is UNDEFINED; this should never happen\n";
		    }
		    
		    if ($depth < $min_depth or
			$depth > $max_depth) {
			$sufficient_depth = 0;
			#warn "depth for $chromosome: $pos $infile is $depth\n";
		    }
		    
		}

		

		if ($sufficient_depth) {
		    
		    foreach my $base ('a', 'c', 'g', 't') {

			### Check whether this site is homozygous and matching reference in all of the A genomes
			my $homozygous_in_A = 1;
			my $files_list_ref = $pileup_files_lists{'A'};
			my @files_list = @{$files_list_ref};
			foreach my $infile( @files_list ) {
			    my $freq = $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{$base};
			    $homozygous_in_A = 0 unless ($freq == 0);   
			}
			
			### Check whether this site is homozygous and different from reference in all of the B genomes
			my $homozygous_in_B = 1;
			$files_list_ref = $pileup_files_lists{'B'};
			@files_list = @{$files_list_ref};
			foreach my $infile( @files_list ) {
			    my $freq = $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{$base};
			    $homozygous_in_B = 0 unless ($freq == 1);   
			}
				
			my $all_zero = 1;
			
			#warn "Examining $base at $chromosome: $pos (we have sufficient depth and all defined) (Loop $loop_count)\n";
			my $print_line = '';	
			$print_line .= "$chromosome";
			$print_line .=  "\t";
			$print_line .=  "$pos";
			$print_line .= "\t";
			$print_line .= $chrom_to_pos_to_refbase{$chromosome}{$pos};
			$print_line .= "\t";
			$print_line .=  "$base";
			
			foreach my $infile (sort keys %filehandles) {
			    my $depth = $chrom_to_pos_to_file_to_depth{$chromosome}{$pos}{$infile};
			    my $freq = $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}{$infile}{$base};
			    unless (defined $freq) {
				die "Frequency is UNDEFINED for $chromosome: $pos $infile $base; this should never happen";
			    }
			    $print_line .= "\t";
			    $print_line .= "$freq";
			    $print_line .= "\t";
			    $print_line .= "$depth";
			    
			    $all_zero = 0 if $freq > $absence_max_freq; 
			}
			
			### Print the data
			if (
			    $base ne $chrom_to_pos_to_refbase{$chromosome}{$pos} and
			    $all_zero == 0 and 
			    $homozygous_in_A and $homozygous_in_B
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
		
		delete $chrom_to_pos_to_file_to_base_to_freq{$chromosome}{$pos}; ### Don't print this location again
	    }
	}
	
    }

    ### Check for end of file
    foreach my $infile (keys %filehandles) {
	my $fh = $filehandles{$infile};
	if ( eof($fh) ) {
	    $eof++;
	    warn "Reached EOF for $infile\n";
	}
    }
}

    






foreach my $infile (keys %filehandles) {
    my $fh = $filehandles{$infile};
    close $fh;
    warn "\tFinished parsing file '$infile'\n";
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
