#!/usr/bin/perl -w

use strict;
use IO::File;
use POSIX;                          # floor, ceil
use List::Util qw(sum min max);
# use Math::Complex;                  # sqrt

#########
# Calculate Pi, Theta, Tajima's D, etc
# Needs to set window size, frequency first
#########

# Input file is unphased vcf file, for FWH SNPs are polarized against Spretus
# fileds #10 to #41 are snp calls
my $infile = IO::File -> new ("<$ARGV[0]") or die "Couldn't open file $ARGV[0]: $!";
# outfile is only for allele frequencies
my $outfile = IO::File -> new (">$ARGV[1]") or die "Couldn't open file $ARGV[1]: $!";

my $chr = $ARGV[2];
chomp($chr);

my $size1 = 14;
my $size2 = 8;
my $windowsize = 5000;
my $windowstart = 0;          # chr1 from 3000000
my $windownum = 0;
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
# a1 = 1/1 + 1/2 + 1/3 ... + 1/S, S is number of segregating sites
my $a1 = 0;
my $a2 = 0;
my $b1 = 0;
my $b2 = 0;
my $c1 = 0;
my $c2 = 0;
my $e1 = 0;
my $e2 = 0;
# D(tajima) = average number of polymorphism - snp number / a1 
my $tajimasd = 0;

my %freqdist2 = ();          # 1 to 16
my @stats2 = ();
my $allelecount2 = 0;
my $snpcount2 = 0;          # total number of SNPs in a window
my $snpUnique2 = 0;
my $snpShared = 0;
my $snpTotal = 0;

# pi = (2 * j * (n - j)) / (n * (n - 1)), n is number samples (chromosomes), j is count for alternative allele, j = AC & n = AN in vcf file
my $pi2 = 0;
my $theta2 = 0;
# a1 = 1/1 + 1/2 + 1/3 ... + 1/S, S is number of segregating sites
my $a12 = 0;
my $a22 = 0;
my $b12 = 0;
my $b22 = 0;
my $c12 = 0;
my $c22 = 0;
my $e12 = 0;
my $e22 = 0;
# D(tajima) = average number of polymorphism - snp number / a1 
my $tajimasd2 = 0;

my @genotypes = ();
my $pairb = 0;
my $pairw1 = 0;
my $pairw2 = 0;
my $fst = 0;
my $snpID = 0;
my @snps = ();              # For Rsquare
my @snpsvalues1 = ();        # Used to determine singletons
my @snpsvalues2 = (); 
my @fields2 = ();
my @fields3 = ();
my $xy = 0;
my $xsquare = 0;
my $ysquare = 0;
my $xy2 = 0;
my $xsquare2 = 0;
my $ysquare2 = 0;
my $pairs1 = 0;
my $pairs2 = 0;
my $rsquare = 0;
my $rsquare2 = 0;

my $snptmp = "";
my $snpvalue1tmp = 0;
my $snpvalue2tmp = 0;
my $tmp = 0;

# For Fay & Wu's H
my $thetasquare = 0;
my $thetasquare2 = 0;
my $thetal = 0;
my $thetal2 = 0;
my $fwh = 0;
my $fwh2 = 0;
my $a2plus1 = 0;
my $a22plus1 = 0;

# Sample size (Consider missing genotypes)
my $sample1 = 0;
my $sample2 = 0;

print $outfile "Window\tChr\tStart\tSNP1\tSNP2\tSNPShared\tPi1\tPi2\tTheta1\tTheta2\tTajimasD1\tTajimasD2\tFWH1\tFWH2\tRsquare1\tRsquare2\tFst\n";

# For both Tajima'a & Fu and Li's
for (my $i = 1; $i < 2 * $size1; $i++) {
    $a1 += 1 / $i;
    $a2 += 1 / ($i * $i);
}
$a2plus1 = $a2 + (1 / (2 * $size1 * 2 * $size1));

for (my $i = 1; $i < 2 * $size2; $i++) {
    $a12 += 1 / $i;
    $a22 += 1 / ($i * $i);
}
$a22plus1 = $a22 + (1 / (2 * $size2 * 2 * $size2));

while (<$infile>) {
  if (!($_ =~ /^\#/) && !($_ =~ "INDEL")) {
    $allelecount = 0;
    $allelecount2= 0;
    
    chomp();
    @fields = split(/\t/, $_);
    
    if ($fields[0] eq $chr) {
      if (length($fields[4]) == 1) {        # Only one Alt allele
	if ($_ =~ /0\/1/ || $_ =~ /0\/0/) {      # Not fixed difference, for one population version
	  for (my $i = 9; $i < 23; $i++) {
	    if ($fields[$i] =~ /^\.\/\./) {
	      $fields[$i] =~ s/\.\/\./0\/0/;
	    }
	    @subfields = split(/[\/\:\,]/, $fields[$i]);
	    $genotypes[($i - 9) * 2] = $subfields[0];
	    $genotypes[($i - 9) * 2 + 1] = $subfields[1];
	    
	    if ($fields[$i] =~ /^0\/0/) {
	      $snps[$snpID] .= 0;
	      $snpsvalues1[$snpID] += 0;
	    }
	    elsif ($fields[$i] =~ /^0\/1/ || $fields[$i] =~ /^1\/0/) {
	      $snps[$snpID] .= 1;
	      $snpsvalues1[$snpID] += 1;
	    }
	    elsif ($fields[$i] =~ /^1\/1/) {
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
	      $fields[$i] =~ s/\.\/\./0\/0/;
	    }
	    @subfields = split(/[\/\:\,]/, $fields[$i]);
	    $genotypes[($i - 9) * 2] = $subfields[0];
	    $genotypes[($i - 9) * 2 + 1] = $subfields[1];
	    
	    if ($fields[$i] =~ /^0\/0/) {
	      $snps[$snpID] .= 0;
	      $snpsvalues2[$snpID] += 0;
	    }
	    elsif ($fields[$i] =~ /^0\/1/ || $fields[$i] =~ /^1\/0/) {
	      $snps[$snpID] .= 1;
	      $snpsvalues2[$snpID] += 1;
	    }
	    elsif ($fields[$i] =~ /^1\/1/) {
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
		if ($snpcount > 1) {
		  # Tajima's
		  $b1 = (2 * $size1 + 1) / (3 * (2 * $size1 - 1)); 
		  $b2 = 2 * (2 * $size1 * 2 * $size1 + 2 * $size1 + 3) / (9 * 2 * $size1 * (2 * $size1 - 1)); 
		  $c1 = $b1 - 1 / $a1; 
		  $c2 = $b2 - (2 * $size1 + 2) / ($a1 * 2 * $size1) + $a2 / ($a1 * $a1); 
		  $e1 = $c1 / $a1; 
		  $e2 = $c2 / ($a1 * $a1 + $a2);
		  
		  $tajimasd = ($pi - $snpcount / $a1) / sqrt(abs($e1 * $snpcount + $e2 * $snpcount * ($snpcount - 1)));
		  $thetasquare = ($snpcount * ($snpcount - 1)) / ($a1 * $a1 + $a2);
                  $fwh = ($pi - $thetal) / sqrt(($theta * ($size1 - 2)) / (6 * ($size1 - 1)) + ($thetasquare * (18 * $size1 * $size1 * (3 * $size1 + 2) * $a2plus1 - (88 * $size1 * $size1 * $size1 + 9 * $size1 * $size1 - 13 * $size1 + 6))) / (9 * $size1 * ($size1 - 1) * ($size1 - 1)));

		}
		else {
		  $tajimasd = "NA";
		  $fwh = "NA";
		}
		
		if ($snpcount2 > 1) {
		  $b12 = (2 * $size2 + 1) / (3 * (2 * $size2 - 1)); 
		  $b22 = 2 * (2 * $size2 * 2 * $size2 + 2 * $size2 + 3) / (9 * 2 * $size2 * (2 * $size2 - 1)); 
		  $c12 = $b12 - 1 / $a12; 
		  $c22 = $b22 - (2 * $size2 + 2) / ($a12 * 2 * $size2) + $a22 / ($a12 * $a12); 
		  $e12 = $c12 / $a12; 
		  $e22 = $c22 / ($a12 * $a12 + $a22);
		  
		  $tajimasd2 = ($pi2 - $snpcount2 / $a12) / sqrt(abs($e12 * $snpcount2 + $e22 * $snpcount2 * ($snpcount2 - 1)));
		  $thetasquare2 = ($snpcount2 * ($snpcount2 - 1)) / ($a12 * $a12 + $a22);
		  $fwh2 = ($pi2 - $thetal2) / sqrt(($theta2 * ($size2 - 2)) / (6 * ($size2 - 1)) + ($thetasquare2 * (18 * $size2 * $size2 * (3 * $size2 + 2) * $a22plus1 - (88 * $size2 * $size2 * $size2 + 9 * $size2 * $size2 - 13 * $size2 + 6))) / (9 * $size2 * ($size2 - 1) * ($size2 - 1)));

		}
		else {
		  $tajimasd2 = "NA";
		  $fwh2 = "NA";
		}
		
		$theta = $snpcount / $a1;
		$theta2 = $snpcount2 / $a12;
		
		# Fst and snn
		$pairb = $pairb / 448;
		$pairw1 = $pairw1 / 378;
		$pairw2 = $pairw2 / 120;
		
		if ($pairb > 0) {
		  $fst = ($pairb - (7 * $pairw1 / 11  + 4 * $pairw2 / 11)) / $pairb;
		}
		else {
		  $fst = "NA";
		}
		
		# -------- Rsquare ----------
		$pairs1 = 0;
		$pairs2 = 0;
		
		for (my $i = 0; $i < $#snps - 1; $i++) {
		  @fields2 = split(//, $snps[$i]);
				    
		  for (my $j = $i + 1; $j < $#snps; $j++) {
		    @fields3 = split(//, $snps[$j]);
		    $xy = 0;
		    $xsquare = 0;
		    $ysquare = 0;
		    $xy2 = 0;
		    $xsquare2 = 0;
		    $ysquare2 = 0;
		    
		    for (my $k = 0; $k <= $#fields2; $k++) {
		      if ($k < $size1) {
			$xsquare += ($fields2[$k] - $snpsvalues1[$i] / $size1) * ($fields2[$k] - $snpsvalues1[$i] / $size1);
			$ysquare += ($fields3[$k] - $snpsvalues1[$j] / $size1) * ($fields3[$k] - $snpsvalues1[$j] / $size1);
			$xy += ($fields2[$k] - $snpsvalues1[$i] / $size1) * ($fields3[$k] - $snpsvalues1[$j] / $size1);
		      }
		      else {
			$xsquare2 += ($fields2[$k] - $snpsvalues2[$i] / $size2) * ($fields2[$k] - $snpsvalues2[$i] / $size2);
			$ysquare2 += ($fields3[$k] - $snpsvalues2[$j] / $size2) * ($fields3[$k] - $snpsvalues2[$j] / $size2);
			$xy2 += ($fields2[$k] - $snpsvalues2[$i] / $size2) * ($fields3[$k] - $snpsvalues2[$j] / $size2);
		      }
		    }
		    
		    # if ($xsquare * $ysquare != 0 && $snpsvalues1[$i] != 1 && $snpsvalues1[$i] != (2 * $size1 - 1) && $snpsvalues1[$j] != 1 && $snpsvalues1[$j] != (2 * $size1 - 1)) {
		    if ($xsquare * $ysquare != 0) {
		      $rsquare += ($xy * $xy) / ($xsquare * $ysquare);
		      $pairs1++;
		    }
		    
		    # if ($xsquare2 * $ysquare2 != 0 && $snpsvalues2[$i] != 1 && $snpsvalues2[$i] != (2 * $size2 - 1) && $snpsvalues2[$j] != 1 && $snpsvalues2[$j] != (2 * $size2 - 1)) {
		    if ($xsquare2 * $ysquare2 != 0) {
		      $rsquare2 += ($xy2 * $xy2) / ($xsquare2 * $ysquare2);
		      $pairs2++;
		    }
		  }
		}
				
		if ($pairs1 > 0) {  
		  $rsquare = $rsquare / $pairs1;
		}
		else {
		  $rsquare = "NA";
		}
		
		if ($pairs2 > 0) {  
		  $rsquare2 = $rsquare2 / $pairs2;
		}
		else {
		  $rsquare2 = "NA";
		} 
		
		print $outfile "$windownum\t$chr\t$windowstart\t$snpUnique2\t$snpUnique\t$snpShared\t$pi2\t$pi\t$theta2\t$theta\t$tajimasd2\t$tajimasd\t$fwh2\t$fwh$rsquare2\t$rsquare\t$fst\n";
		
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
	    
	    $pi = 2 * $allelecount * (2 * $size1 - $allelecount) / (2 * $size1 * (2 * $size1 - 1));
	    $pi2 = 2 * $allelecount2 * (2 * $size2 - $allelecount2) / (2 * $size2 * (2 * $size2 - 1));
	    
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
	    $rsquare = 0;
	    $rsquare2 = 0;
	    $snpTotal = 1;
	    
	    # $recrate = $fields[41];
	    $pairb = 0;
	    $pairw1 = 0;
	    $pairw2 = 0;
	    
	  }
	  
	  else {
	    if (!($allelecount == 2 * $size1 || $allelecount == 0)) {
	      $snpcount++;
	    }
	    
	    if (!($allelecount2 == 2 * $size2 || $allelecount2 == 0)) {
	      $snpcount2++;
	    }
	    
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
	  
	  for (my $i = 0; $i < (2 * ($size1 + $size2) - 1); $i++) {
	    for (my $j = $i + 1; $j < 2 * ($size1 + $size2); $j++) {
	      if ($genotypes[$i] ne $genotypes[$j]) {
		if ($i < 2 * $size1 && $j < 2 * $size1) {
		  $pairw1++;
		}
		elsif ($i >= 2 * $size1 && $j >= 2 * $size1) {
		  $pairw2++;
		}
		else {
		  $pairb++;
		}
	      }
	    }
	  }
	}
      }
    }	
  }
}
