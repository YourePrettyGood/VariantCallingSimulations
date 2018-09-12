#!/usr/bin/perl
use POSIX;
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

=pod

=head1 NAME

closestIndelDistance.pl - Calculate the distance from a SNP to its closest indel

=head1 SYNOPSIS

closestIndelDistance.pl [options]

 Options:
  --help,-h,-?         Display this help documentation
  --vcf,-v             Input VCF used to find indel positions
                       (default: -, which is STDIN)
  --insnp,-i           Input INSNP (TSV) of SNPs
  --debug,-d           Output extra information to STDERR

=head1 DESCRIPTION

closestIndelDistance.pl calculates the distance from each SNP in an
INSNP TSV file (of the style used as input to seqtk mutfa) to the
nearest indel, identified from the supplied VCF.

You can then pass the output into subsetVCFstats.pl to partition
the distances by FP vs. TP SNP.

=cut


my $vcf_path = "";
my $insnp_path = "";
my $help = 0;
my $man = 0;
my $debug = 0;

GetOptions('vcf|v=s' => \$vcf_path, 'insnp|i=s' => \$insnp_path, 'debug|d' => \$debug, 'help|h|?' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -output => \*STDERR, -verbose => 2) if $man;

#Hash of arrays of positions of indels (hash keys are scaffolds):
my %indel_positions = ();

#Open the VCF to parse out indel positions:
my $input_vcf;
if ($vcf_path eq '-') {
   open $input_vcf, "<&", \*STDIN or die "Failed to duplicate STDIN file handle for input VCF due to error $!";
} else {
   open $input_vcf, "<", $vcf_path or die "Failed to open input VCF ${vcf_path} due to error $!";
}

print STDERR "Parsing indel positions from input VCF ${vcf_path}\n" if $debug;
while (my $line = <$input_vcf>) {
   chomp $line;
   next if $line =~ /^#/; #Skip header and comment lines
   my ($scaf, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = split /\t/, $line, 10;
   #Note that we only care about pos, so the number of samples doesn't matter
   my @alts = split /,/, $alt;
   my $indel = 0; #Sentinel to detect indel records
   $indel = 1 unless length($ref) == 1; #Detect deletions
   for my $altseq (@alts) {
      $indel = 1 unless length($altseq) == 1; #Detect insertions
   }
   #Initialize array for scaffold if not already present in hash:
   $indel_positions{$scaf} = [] unless exists($indel_positions{$scaf});
   push @{$indel_positions{$scaf}}, $pos if $indel; #Add indel position to list
}
close($input_vcf);

my $input_insnp;
open $input_insnp, "<", $insnp_path or die "Failed to open input INSNP file ${insnp_path} due to error $!";

print STDERR "Reading INSNP file ${insnp_path}\n" if $debug;
my $prev_scaf = "";
my $last_indel = 0;
my $indel_index = 0;
while (my $line = <$input_insnp>) {
   chomp $line;
   my ($scaf, $pos, $from, $to) = split /\t/, $line, 4;
   #Reset index and retrieve number of indels if on a new scaffold:
   $indel_index = 0 unless $scaf eq $prev_scaf;
   if (exists($indel_positions{$scaf}) and scalar(@{$indel_positions{$scaf}}) > 0) {
      $last_indel = scalar(@{$indel_positions{$scaf}}) - 1 unless $scaf eq $prev_scaf;
      #Move index forward until current indel is at or after the current SNP, or last indel in list:
      $indel_index++ while $indel_positions{$scaf}[$indel_index] < $pos and $indel_index < $last_indel;
      #Calculate distances to the two surrounding indels:
      my $dist_prior = abs($indel_positions{$scaf}[$indel_index - 1] - $pos) + 1 unless $indel_index == 0;
      my $dist_current = abs($indel_positions{$scaf}[$indel_index] - $pos) + 1;
      #Shorter of the two distances, or dist_current if there is no prior indel:
      my $min_dist = $indel_index == 0 ? $dist_current : ($dist_current <= $dist_prior ? $dist_current : $dist_prior);
      #Print the location of the SNP and the distance to nearest indel:
      print STDOUT join("\t", $scaf, $pos, $min_dist), "\n";
   } else { #No indels on this scaffold, therefore distance is NA
      print STDOUT join("\t", $scaf, $pos, "NA"), "\n";
   }
   $prev_scaf = $scaf;
}
close($input_insnp);
