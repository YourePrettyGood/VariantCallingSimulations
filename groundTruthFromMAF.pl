#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

################################################################
#                                                              #
################################################################

#First pass script to extract ground truth INSNP from pairwise MAF
# of 1:1 alignments produced by LAST
#2023/05/31: Adjusted to output each a block at the first p line
#            instead of the last p line to handle newer MAFs from
#            LAST -- Thanks to Clair Han for catching this

my $SCRIPTNAME = "groundTruthFromMAF.pl";
my $VERSION = "1.1";

=pod

=head1 NAME

groundTruthFromMAF.pl - Extract ground truth INSNPs from pairwise 1:1 MAF

=head1 SYNOPSIS

groundTruthFromMAF.pl [options] <Species 1 prefix> <Species 2 prefix>

 Options:
  --help,-h,-?          Print this help documentation
  --input_MAF,-i        Path to input MAF 1:1 alignment file (default: STDIN)
  --debug,-d            Output debugging information to STDERR
  --version,-v          Output version string

=head1 DESCRIPTION

This script generates INSNP files of the ground truth variant calls
in each species' coordinate space based on a MAF file of 1:1 pairwise
alignments produced by LAST. The scaffold IDs in the MAF file should
follow the format [species prefix].[scaffold name], where the species
prefix and scaffold name both are alphanumeric (plus underscores). As
long as the scaffold IDs follow this convention, and you provide the
species prefixes precisely as found in the MAF, the script will work
as intended.
The script will also generate BED files to indicate the extents of each
aligned segment in each appropriate coordinate space.

For example:

=begin text

a score=10 mismap=1e-05
s DyakTai18E2.2L 0 10 + 29953808 AATAACGGCT
s DyakNY73PB.2L  0 10 + 24234981 AACAACGGGT
p                                !!!!!!!!!!
p                                !!!!!!!!!!

=end text

should output in DyakTai18E2_INSNP.tsv:

=begin text

2L	3	T	C
2L	9	C	G

=end text

and in DyakNY73PB_INSNP.tsv:

=begin text

2L	3	C	T
2L	9	G	C

=end text

as well as the following BED files:

DyakTai18E2_aligned_regions.bed:

=begin text

2L	0	10

=end text

DyakNY73PB_aligned_regions.bed:

=begin text

2L	0	10

=end text

=cut

my $help = 0;
my $man = 0;
my $maf_path = "STDIN";
my $debug = 0;
my $dispversion = 0;
GetOptions('input_maf|i=s' => \$maf_path, 'debug|d+' => \$debug, 'version|v' => \$dispversion, 'help|h|?+' => \$help, man => \$man) or pod2usage(2);
pod2usage(-exitval => 1, -verbose => $help, -output => \*STDERR) if $help;
pod2usage(-exitval => 0, -verbose => 2, -output => \*STDERR) if $man;

print STDERR "${SCRIPTNAME} version ${VERSION}\n" if $dispversion;
exit 0 if $dispversion;

#Open the MAF file, or set it up to be read from STDIN:
print STDERR "Opening MAF file\n" if $debug;
my $maf_fh;
if ($maf_path ne "STDIN") {
   unless(open($maf_fh, "<", $maf_path)) {
      print STDERR "Error opening MAF file: ${maf_path}.\n";
      exit 2;
   }
} else {
   open($maf_fh, "<&", "STDIN"); #Duplicate the file handle for STDIN to $maf_fh so we can seamlessly handle piping
}

#Extract the positional arguments, which are prefixes for the species:
my @prefixes = @ARGV;
#Make a hash so we can quickly check if the species from an "s" line matches:
print STDERR "Constructing prefix hash\n" if $debug;
my %prefix_hash = map {$_ => 1} @prefixes;

if (scalar(@prefixes) != 2) {
   close($maf_fh);
   print STDERR "It appears you have provided more than 2 prefixes. We only support pairwise alignments for the moment.\n";
   exit 3;
}

#Establish file handles for the output INSNPs and BEDs for each prefix:
print STDERR "Opening output INSNP and BED handles\n" if $debug;
my %spp_fh = ();
my %spp_bed_fh = ();
for my $prefix (@prefixes) {
   unless(open($spp_fh{$prefix}, ">", "${prefix}_INSNP.tsv")) {
      print STDERR "Could not open INSNP file for prefix ${prefix} for writing.\n";
      close($maf_fh);
      exit 4;
   }
   unless(open($spp_bed_fh{$prefix}, ">", "${prefix}_aligned_regions.bed")) {
      print STDERR "Could not open BED file for prefix ${prefix} for writing.\n";
      close($maf_fh);
      exit 4;
   }
}

sub revcomp($) {
   my $input_sequence = shift @_;
   my $reverse_sequence = reverse $input_sequence; #Reverse
   $reverse_sequence =~ tr/ACGTNacgtn/TGCANtgcan/; #Complement
   return $reverse_sequence;
}

sub stranded_allele($$) {
   my $input_allele = shift @_;
   my $strand = shift @_;
   return $input_allele unless $strand eq "-";
   return revcomp($input_allele) if $strand eq "-";
}

#Hash of hashes to store the elements of a MAF alignment:
my %aligned_seqs = ();
#Hash of offsets to stay in the appropriate coordinate space:
my %offsets = ();
#Array of species in the alignment:
my @species_aligned = ();
#Number of p lines we've encountered, which should equal the number of s lines:
my $num_p_lines = 0;

#Keep track of the total number of alignments for debugging:
my $num_a_lines = 0;

print STDERR "Parsing MAF ${maf_path}\n" if $debug;
while (my $line = <$maf_fh>) {
   chomp $line;
   next if $line =~ /^#/ or $line eq ""; #Skip header and empty lines
   my @maf_elems = split /\s+/, $line; #MAF uses padded spacing
   @species_aligned = () if $maf_elems[0] eq "a"; #Reset species on new alignment
   $num_p_lines = 0 if $maf_elems[0] eq "a"; #Reset p count on new alignment
   $num_a_lines++ if $maf_elems[0] eq "a";
   print STDERR "Parsed ${num_a_lines} alignments\n" if $debug and ${num_a_lines} % 1000 == 1 and $maf_elems[0] eq "a";
   if ($maf_elems[0] eq "s") { #For aligned Sequence records, store the details
      my ($species, $scaffold) = split /\./, $maf_elems[1], 2;
      unless (exists($prefix_hash{$species})) {
         print STDERR "Alignment has species not found in provided prefixes: ${line}\n";
         close($maf_fh);
         exit 5;
      }
      push @species_aligned, $species;
      my @seq_arr = split //, uc($maf_elems[6]);
      $aligned_seqs{$species} = {'scaffold' => $scaffold,
         'start' => $maf_elems[4] eq "-" ? $maf_elems[5]-$maf_elems[2] : $maf_elems[2]+1,
         'size' => $maf_elems[3],
         'strand' => $maf_elems[4],
         'seq' => \@seq_arr};
   }
   if ($maf_elems[0] eq "p") { #Identify SNPs once we're in the Probability lines
      $num_p_lines++;
      next unless $num_p_lines == 1; #Skip unless we're on the first p line
      #Make a hash of the aligned species so we can take set differences:
      my %aligned_spp = map {$_ => 1} @species_aligned;

      #Sequence length for each "s" record should be identical, so just take
      # the length from the first species in the alignment:
      my $aln_length = scalar(@{$aligned_seqs{$species_aligned[0]}{'seq'}});

      #Initialize the offsets for each coordinate space:
      for my $species (@species_aligned) {
         $offsets{$species} = 0;
      }

      #Construct the BED line for this alignment for each species:
      for my $species (@species_aligned) {
         next unless exists($prefix_hash{$species});
         my $scaffold = $aligned_seqs{$species}{'scaffold'};
         my $BEDstart = $aligned_seqs{$species}{'strand'} eq "-" ? $aligned_seqs{$species}{'start'}-$aligned_seqs{$species}{'size'} : $aligned_seqs{$species}{'start'}-1;
         my $BEDend = $aligned_seqs{$species}{'strand'} eq "-" ? $aligned_seqs{$species}{'start'} : $aligned_seqs{$species}{'start'}+$aligned_seqs{$species}{'size'}-1;
         print {$spp_bed_fh{$species}} join("\t", $scaffold, $BEDstart, $BEDend), "\n";
      }
      for (my $i = 0; $i < $aln_length; $i++) {
         #Find SNP, assuming there are 2 spp in the alignment:
         if ($aligned_seqs{$species_aligned[0]}{'seq'}[$i] ne $aligned_seqs{$species_aligned[1]}{'seq'}[$i] and $aligned_seqs{$species_aligned[0]}{'seq'}[$i] ne '-' and $aligned_seqs{$species_aligned[1]}{'seq'}[$i] ne '-') {
            for my $species (@species_aligned) {
               #Hacky way of getting the ID of the other of the two species:
               my @remaining_species = grep {$_ ne $species} @species_aligned;
               if (scalar(@remaining_species) == 0) {
                  print STDERR "Something went wrong trying to identify the other species' (not ${species}) prefix\n" if $debug;
                  for my $output_species (keys %aligned_spp) {
                     print STDERR "Species: ${output_species}\n" if $debug;
                  }
                  close($maf_fh);
                  for my $prefix (@prefixes) {
                     close($spp_fh{$prefix});
                     close($spp_bed_fh{$prefix});
                  }
                  exit 6;
               }
               my $other_species = $remaining_species[0];

               #Print the SNP to the appropriate INSNP:
               my $insnp_position = $aligned_seqs{$species}{'strand'} eq "-" ? $aligned_seqs{$species}{'start'}-$offsets{$species} : $aligned_seqs{$species}{'start'}+$offsets{$species};
               my $species1_allele = stranded_allele($aligned_seqs{$species}{'seq'}[$i], $aligned_seqs{$species}{'strand'});
               my $species2_allele = stranded_allele($aligned_seqs{$other_species}{'seq'}[$i], $aligned_seqs{$species}{'strand'});
               print {$spp_fh{$species}} join("\t", $aligned_seqs{$species}{'scaffold'}, $insnp_position, $species1_allele, $species2_allele), "\n";
            }
         }
         for my $species (@species_aligned) {
            $offsets{$species}++ unless $aligned_seqs{$species}{'seq'}[$i] eq "-";
         }
      }
   }
}
close($maf_fh);
print STDERR "Done parsing all ${num_a_lines} alignments from MAF file\n" if $debug;

print STDERR "Closing INSNP and BED files\n" if $debug;
for my $prefix (@prefixes) {
   close($spp_fh{$prefix});
   close($spp_bed_fh{$prefix});
}

exit 0;
