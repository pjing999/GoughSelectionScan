#!/usr/bin/perl -w

use strict;
use IO::File;
use POSIX;                          # floor, ceil
use List::Util qw(sum min max);
# use Statistics::Basic qw(:all);
# use Math::Complex;                  # sqrt

#########
# Calculate Pi, Theta, H12, Haplotype Number, etc
# Needs to set window size, frequency first
#########

# Input file is snps_gough_domestics_phased.vcf
# fileds #10 to #41 are snp calls
my $infile = IO::File -> new ("<$ARGV[0]") or die "Couldn't open file $ARGV[0]: $!";
# outfile is only for allele frequencies
my $outfile = IO::File -> new (">$ARGV[1]") or die "Couldn't open file $ARGV[1]: $!";

my $chr = $ARGV[2];
chomp($chr);

my $size1 = 14;
my $size2 = 8;
my $windownum = 0;
my $windowsize = 5000;
my $windowstart = 0;          # chr1 from 3000000
my @fields = ();
my @subfields = ();
my %freqdist = ();          # 1 to 16
my @stats = ();
my $allelecount = 0;
my $snpcount = 0;          # total number of SNPs in a window
my $snpUnique = 0;

# pi = (2 * j * (n - j)) / (n * (n - 1)), n is number samples (chromosomes), j is count for alternative allele, j = AC & n = AN in vcf file
my $pi = 0;
my $theta = 0;
my @haplotypes = ();
my %hapHash = ();
my $hapCount = 0;
my $hapFreq = 0;
my $hapHeter = 1;         # haplotype heterozygosity

my %freqdist2 = ();          # 1 to 16
my @stats2 = ();
my $allelecount2 = 0;
my $snpcount2 = 0;          # total number of SNPs in a window
my $snpUnique2 = 0;
my $snpShared = 0;
# my $recrate = 0;
my $snpTotal = 0;

# pi = (2 * j * (n - j)) / (n * (n - 1)), n is number samples (chromosomes), j is count for alternative allele, j = AC & n = AN in vcf file
my $pi2 = 0;
my $theta2 = 0;
my @haplotypes2 = ();
my %hapHash2 = ();
my $hapCount2 = 0;
my $hapFreq2 = 0;
my $hapHeter2 = 1;          # haplotype heterozygosity

my @genotypes = ();
my $snpID = 0;
my @snps = ();              # For Rsquare
my @snpsvalues1 = ();        # Used to determine singletons
my @snpsvalues2 = (); 
my @fields2 = ();
my @fields3 = ();

my $hapfreqhigh1 = 0;
my $hapfreqmedian1 = 0;
my $hapfreqhigh2 = 0;
my $hapfreqmedian2 = 0;

my $snptmp = "";
my $snpvalue1tmp = 0;
my $snpvalue2tmp = 0;
my $tmp = 0;

print $outfile "Window\tChr\tStart\tSNPGermany\tSNPGough\tSNPShared\tHapNumGermany\tHapNumGough\tHapFreqHighGermany\tHapFreqHighGough\tHapHeterGermany\tHapHeterGough\tH12Germany\tH12Gough\n";

my $h11 = 0;
my $h121 = 0;
my $h12 = 0;
my $h122 = 0;
my $id = 0;
my $p1 = 0;
my $p2 = 0;

sub median {
  sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}

while (<$infile>) {
  if (!($_ =~ /^\#/) && !($_ =~ "INDEL")) {
    $allelecount = 0;
    $allelecount2= 0;
	
    chomp();
    @fields = split(/\t/, $_);
    
    if ($fields[0] eq $chr) {
      if (length($fields[4]) == 1) {        # Only one Alt allele
	if ($_ =~ /0\|1/ || $_ =~ /0\|0/) {      # Not fixed difference, for one population version
	  for (my $i = 9; $i < 23; $i++) {
	    if ($fields[$i] =~ /^\.\/\./) {
	      $fields[$i] =~ s/\.\/\./0\|0/;
	    }
	    @subfields = split(/[\|\:\,]/, $fields[$i]);
	    $genotypes[($i - 9) * 2] = $subfields[0];
	    $genotypes[($i - 9) * 2 + 1] = $subfields[1];
	    
	    if ($fields[$i] =~ /^0\|0/) {
	      $snps[$snpID] .= 0;
	      $snpsvalues1[$snpID] += 0;
	    }
	    elsif ($fields[$i] =~ /^0\|1/ || $fields[$i] =~ /^1\|0/) {
	      $snps[$snpID] .= 1;
	      $snpsvalues1[$snpID] += 1;
	    }
	    elsif ($fields[$i] =~ /^1\|1/) {
	      $snps[$snpID] .= 2;
	      $snpsvalues1[$snpID] += 2;
	    }
	    
	    if ($subfields[0] eq "1") {
	      $allelecount++; 
	    } 
	    if ($subfields[1] eq "1") {  
	      $allelecount++;  
	    }  
	  }
	  
	  for (my $i = 23; $i < 31; $i++) {
	    if ($fields[$i] =~ /^\.\/\./) {
	      $fields[$i] =~ s/\.\/\./0\|0/;
	    }
	    @subfields = split(/[\|\:\,]/, $fields[$i]);
	    $genotypes[($i - 9) * 2] = $subfields[0];
	    $genotypes[($i - 9) * 2 + 1] = $subfields[1];
	    
	    if ($fields[$i] =~ /^0\|0/) {
	      $snps[$snpID] .= 0;
	      $snpsvalues2[$snpID] += 0;
	    }
	    elsif ($fields[$i] =~ /^0\|1/ || $fields[$i] =~ /^1\|0/) {
	      $snps[$snpID] .= 1;
	      $snpsvalues2[$snpID] += 1;
	    }
	    elsif ($fields[$i] =~ /^1\|1/) {
	      $snps[$snpID] .= 2;
	      $snpsvalues2[$snpID] += 2;
	    }
	    
	    if ($subfields[0] eq "1") {
	      $allelecount2++;
	    } 
	    if ($subfields[1] eq "1") {  
	      $allelecount2++;
	    }  
	  }
	  
	  $snpID++;
	  
	  if ($windowstart == 0 || ($fields[1] - $windowstart > $windowsize)) {
	    if ($windowstart == 0) {
	      $windowstart = $fields[1] - $fields[1] % $windowsize;
	    }
	    else {
	      if ($snpcount > 0 || $snpcount2 > 0) {  
		# $theta = $snpcount / $a1;
		
		$hapHash{$_}++ for @haplotypes;
		$hapCount = keys %hapHash;
		
		for my $key (keys %hapHash) {
		  $hapFreq += $hapHash{$key} / (2 * $size1);
		}
		
		$hapFreq = $hapFreq / $hapCount;
		
		$id = 1;
		$p1 = 0;
		$p2 = 0;
		for my $key (sort { $hapHash{$b} <=> $hapHash{$a} } keys %hapHash) {
		  if ($id == 1) {
		    $p1 = $hapHash{$key} / (2 * $size1);
		  }
		  elsif ($id == 2) {
		    $p2 = $hapHash{$key} / (2 * $size1);
		  }
		  
		  $hapHeter = $hapHeter - ($hapHash{$key} * $hapHash{$key}) / (2 * $size1 * 2 * $size1);
		  $h11 += ($hapHash{$key} * $hapHash{$key}) / (2 * $size1 * 2 * $size1);
		  $id++;
		}
		
		$h121 = $h11 + 2 * $p1 * $p2;
		$hapfreqhigh1 = max(values %hapHash) / (2 * $size1);
		$hapfreqmedian1 = median(values %hapHash) / (2 * $size1);
		
		# $theta2 = $snpcount2 / $a12;  
		$hapHash2{$_}++ for @haplotypes2;
		$hapCount2 = keys %hapHash2;
		
		for my $key (keys %hapHash2) {
		  $hapFreq2 += $hapHash2{$key} / (2 * $size2);
		}
		
		$hapFreq2 = $hapFreq2 / $hapCount2;
		
		$id = 1;
		$p1 = 0;
		$p2 = 0;
		for my $key (sort { $hapHash2{$b} <=> $hapHash2{$a} } keys %hapHash2) {
		  if ($id == 1) {
		    $p1 = $hapHash2{$key} / (2 * $size2);
		  }
		  elsif ($id == 2) {
		    $p2 = $hapHash2{$key} / (2 * $size2);
		  }
		  
		  $hapHeter2 = $hapHeter2 - ($hapHash2{$key} * $hapHash2{$key}) / (2 * $size2 * 2 * $size2);
		  $h12 += ($hapHash2{$key} * $hapHash2{$key}) / (2 * $size2 * 2 * $size2);
		  $id++;
		}
		
		$h122 = $h12 + 2 * $p1 * $p2;
		$hapfreqhigh2 = max(values %hapHash2) / (2 * $size2);
		$hapfreqmedian2 = median(values %hapHash2) / (2 * $size2);
		
		print $outfile "$windownum\t$chr\t$windowstart\t$snpUnique2\t$snpUnique\t$snpShared\t$hapCount2\t$hapCount\t$hapfreqhigh2\t$hapfreqhigh1\t$hapHeter2\t$hapHeter\t$h122\t$h121\n";
		
	      }
	    }
	    
	    if ($allelecount == 2 * $size1 || $allelecount == 0) {
	      $snpcount = 0;
	    }
	    else {
	      $snpcount = 1;
	    }
	    
	    if ($allelecount2 == 2 * $size2 || $allelecount2 == 0) {
	      $snpcount2 = 0;
	    }
	    else {
	      $snpcount2 = 1;
	    }
	    
	    if ($allelecount > 0 && $allelecount < 2 * $size1 && $allelecount2 == 0) {
	      $snpUnique = 1;
	    }
	    else {
	      $snpUnique = 0;
	    }
	    
	    if ($allelecount == 0 && $allelecount2 > 0 && $allelecount2 < 2 * $size2) {
	      $snpUnique2 = 1;
	    }
	    else {
	      $snpUnique2 = 0;
	    }
	    
	    if ($allelecount > 0 && $allelecount2 > 0) {
	      $snpShared = 1;
	    }
	    else {
	      $snpShared = 0;
	    }
	    
	    $windowstart = $windowstart + 5000;
	    $windownum++;
	    
	    @haplotypes = ();
	    # if (!($allelecount == 2 * $size1 || $allelecount == 0)) {        # Fixed differences don't matter
	    for (my $i = 9; $i < 23; $i++) {
	      @subfields = split(/[\|\:\,]/, $fields[$i]);
	      $haplotypes[($i - 9) * 2] .= $subfields[0];
	      $haplotypes[($i - 9) * 2 + 1] .= $subfields[1];
			    }
	    # }
	    %hapHash = ();
	    $hapCount = 0;
	    $hapFreq = 0;
	    $hapHeter = 1;
	    $h11 = 0;
	    $h121 = 0;
	    
	    $pi = 2 * $allelecount * (2 * $size1 - $allelecount) / (2 * $size1 * (2 * $size1 - 1));
	    $pi2 = 2 * $allelecount2 * (2 * $size2 - $allelecount2) / (2 * $size2 * (2 * $size2 - 1));
	    
	    @haplotypes2 = ();
	    # if (!($allelecount2 == 2 * $size2 || $allelecount2 == 0)) {        # Fixed differences don't matter
	    for (my $i = 23; $i < 31; $i++) {
	      @subfields = split(/[\|\:\,]/, $fields[$i]);
	      $haplotypes2[($i - 23) * 2] .= $subfields[0];
	      $haplotypes2[($i - 23) * 2 + 1] .= $subfields[1];
	    }
	    # }
	    %hapHash2 = ();
	    $hapCount2 = 0;
	    $hapFreq2 = 0;
	    $hapHeter2 = 1;
	    $h12 = 0;
	    $h122 = 0;
	    
	    $snptmp = $snps[$#snps];
	    $snpvalue1tmp = $snpsvalues1[$#snpsvalues1];
	    $snpvalue2tmp = $snpsvalues2[$#snpsvalues2];
	    @snps = ();
	    @snpsvalues1 = ();
	    @snpsvalues2 = ();
	    $snpID = 0;
	    $snps[$snpID] = $snptmp;
	    $snpsvalues1[$snpID] = $snpvalue1tmp;
	    $snpsvalues2[$snpID] = $snpvalue2tmp;
	    $snpID++;
	    $snpTotal = 1;
	    
	  }
	  
	  else {
	    if (!($allelecount == 2 * $size1 || $allelecount == 0)) {
	      $snpcount++;
	    }
	    
	    if (!($allelecount2 == 2 * $size2 || $allelecount2 == 0)) {
	      $snpcount2++;
	    }
	    
	    # if (!($allelecount == 2 * $size1 || $allelecount == 0)) {        # Fixed differences don't matter
	    for (my $i = 9; $i < 23; $i++) {
	      @subfields = split(/[\|\:\,]/, $fields[$i]);
	      $haplotypes[($i - 9) * 2] .= $subfields[0];
	      $haplotypes[($i - 9) * 2 + 1] .= $subfields[1];
	    }
	    # }
	    
	    # if (!($allelecount2 == 2 * $size2 || $allelecount2 == 0)) {        # Fixed differences don't matter
	    for (my $i = 23; $i < 31; $i++) {
	      @subfields = split(/[\|\:\,]/, $fields[$i]);
	      $haplotypes2[($i - 23) * 2] .= $subfields[0];
	      $haplotypes2[($i - 23) * 2 + 1] .= $subfields[1];
	    }
	    # }
	    
	    if ($allelecount > 0 && $allelecount < 2 * $size1 && $allelecount2 == 0) {
	      $snpUnique++;
	    }
	    
	    if ($allelecount == 0 && $allelecount2 > 0 && $allelecount2 < 2 * $size2) {
	      $snpUnique2++;
	    }
	    
	    if ($allelecount > 0 && $allelecount2 > 0) {
	      $snpShared++;
	    }
	    
	    $snpTotal++;
	    
	    $pi += 2 * $allelecount * (2 * $size1 - $allelecount) / (2 * $size1 * (2 * $size1 - 1));
	    $pi2 += 2 * $allelecount2 * (2 * $size2 - $allelecount2) / (2 * $size2 * (2 * $size2 - 1));
	    
	  }
	  
	}
      }
    }	
  }
}
