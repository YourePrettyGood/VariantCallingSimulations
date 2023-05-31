#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

##########################################################################
# combineSampleStats.pl                                                  #
# Version 1.0 (2019/10/21)                                               #
# Description:                                                           #
# This script adds up equivalent bins from VCF stat histograms from      #
# different sample IDs of a jointly genotyped sample, and outputs the    #
# combined counts for use in generating PDFs and CDFs.                   #
##########################################################################

my $SCRIPTNAME = "combineSampleStats.pl";
my $VERSION = "1.0";

=pod

=head1 NAME

combineSampleStats.pl - Combine VCF stat histogram bins across samples

=head1 SYNOPSIS

combineSampleStats.pl [options]

 Options:
  --help,-h,-?		Display this help documentation
  --input_tsv,-i        Uncompressed concatenated VCF stat histograms
                        (Default: STDIN)
  --version,-v          Output version string

=head1 DESCRIPTION
This script adds up equivalent bins from VCF stat histograms from
different sample IDs of a jointly genotyped sample, and outputs the
combined counts for use in generating PDFs and CDFs.

=cut

#Function for numerical comparison accounting for NAs:
sub numerical_comp_wNA {
   my $a = shift @_;
   my $b = shift @_;
   if ($a =~ /^NA|\.$/ or $b =~ /^NA|\.$/) {
      if ($a eq $b) {
         return 0;
      } elsif ($a =~ /^NA|\.$/) {
         return -1;
      } else {
         return 1;
      }
   } else { #No NAs, so handle like standard numeric comparison
      return $a <=> $b;
   }
}

my $display_version = 0;
my $help = 0;
my $man = 0;
my $input_path = "STDIN";
my $dispversion = 0;
GetOptions('input_tsv|i=s' => \$input_path, 'help|h|?+' => \$help, man => \$man, 'version|v' => \$dispversion) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

my $tsvfh;
if ($input_path ne "STDIN") {
   unless(open($tsvfh, "<", $input_path)) {
      print STDERR "Error opening input stat histogram TSV file ${input_path}.\n";
      exit 2;
   }
} else {
   open($tsvfh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $tsvfh so we can seamlessly handle piping
}

#Feed the header line through, excluding the sample ID column
my $headerline = <$tsvfh>;
chomp $headerline;
my @headers = split /\t/, $headerline;
shift @headers; #Skip the first column (sample ID)
print STDOUT join("\t", @headers, "Fraction", "CumulativeFraction"), "\n";

#Hash of hashes to store the counts for each statistic bin
#First level of hashing is on a combined key representing
# the variant caller, site category (ER, FN, FP, TN, or TP),
# and statistic (e.g. DP, MQ, QUAL, etc.).
#Second level of hashing is on the specific value of the
# statistic (hopefully read in as a number).
my %counts = ();
#Hash to store total counts for each statistic, used for
# calculating fraction and cumulative fraction
#Hash key is combined key as described above.
my %total_counts = ();
while (my $line = <$tsvfh>) {
   chomp $line;
   #Decompose the input line into columns:
   my ($sampleid, $caller, $category, $statistic, $count, $statvalue) = split /\t/, $line, 6;
   #Compose the first-level hash key:
   my $mainkey = join("\t", $caller, $category, $statistic);
   #Make sure we create the second-level hash in the first-level hash
   # if it doesn't already exist:
   $counts{$mainkey} = {} unless exists($counts{$mainkey});
   $total_counts{$mainkey} = 0 unless exists($total_counts{$mainkey});
   #Fix the potential broken input case where count is blank and
   # statvalue is a left-padded space-separated pair of count and
   # statvalue:
   if ($count eq "") {
      #Trim padding:
      $statvalue =~ s/^\s+|\s+$//g;
      #Split on whitespace:
      my ($fixed_count, $fixed_statvalue) = split/\s+/, $statvalue, 2;
      $count = $fixed_count;
      $statvalue = $fixed_statvalue;
   }
   #Initialize the count if the stat value wasn't seen before:
   unless (exists($counts{$mainkey}{$statvalue})) {
      $counts{$mainkey}{$statvalue} = $count;
   } else {
      #Otherwise, add to the existing count:
      $counts{$mainkey}{$statvalue} += $count;
   }
   #Add to total counts:
   $total_counts{$mainkey} += $count unless $statvalue =~ /^NA|\.$/;
#   print STDERR join("\t", $mainkey, $count, $total_counts{$mainkey}), "\n";
}

#Now output in somewhat sorted order
#Outer order is arbitrary, we don't care if input order is preserved
# for the first level of the hash
for my $mainkey (keys %counts) {
   #Keep track of a cumulative count for each statistic to calculate
   # cumulative fraction:
   my $cumulative_count = 0;
   #However, the sort order of the statistic values is important to preserve
   for my $statvalue (sort { numerical_comp_wNA($a, $b) } keys %{$counts{$mainkey}}) {
      $cumulative_count += $counts{$mainkey}{$statvalue} unless $statvalue =~ /^NA|\.$/;
      print join("\t", $mainkey, $counts{$mainkey}{$statvalue}, $statvalue, $counts{$mainkey}{$statvalue}/$total_counts{$mainkey}, $cumulative_count/$total_counts{$mainkey}), "\n" unless $statvalue =~ /^NA|\.$/;
      print join("\t", $mainkey, $counts{$mainkey}{$statvalue}, $statvalue, "NA", "NA"), "\n" if $statvalue =~ /^NA|\.$/;
   }
}

#Close the input file if it was indeed opened:
if ($input_path ne "STDIN") {
   close($tsvfh);
}

exit 0;
