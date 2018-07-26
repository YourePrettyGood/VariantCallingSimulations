/**********************************************************************************
 * diploidizeSNPlog.cpp                                                           *
 * Written by Patrick Reilly                                                      *
 * Version 1.0 written 2017/01/23                                                 *
 * Description:                                                                   *
 *                                                                                *
 * Syntax: diploidizeSNPlog [haploid 1 merged SNP log] [haploid 2 merged SNP log] *
 **********************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cctype>
#include <vector>
#include <map>
#include <sstream>
#include <array>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define VERSION "1.0"

//Usage/help:
#define USAGE "diploidizeSNPlog\nUsage:\n diploidizeSNPlog -i [FASTA .fai] -a [haploid 1 merged SNP log] -b [haploid 2 merged SNP log]\n"

using namespace std;

vector<string> splitString(string line_to_split, char delimiter) {
   vector<string> line_vector;
   string element;
   istringstream line_to_split_stream(line_to_split);
   while (getline(line_to_split_stream, element, delimiter)) {
      line_vector.push_back(element);
   }
   return line_vector;
}

long baseToLong(string &base) {
   switch(base[0]) {
      case 'A':
         return 0;
      case 'C':
         return 1;
      case 'G':
         return 2;
      case 'T':
         return 3;
      case 'a':
         return 0;
      case 'c':
         return 1;
      case 'g':
         return 2;
      case 't':
         return 3;
      default:
         return 4;
   }
}

long degenerateBases(long a, long b) {
   if (a == b) { //Homozygous
      return a;
   } else if (a > 3 || b > 3) { //Either is an N
      return 4;
   } else if ((a == 0 && b == 1) || (a == 1 && b == 0)) { //Het A/C
      return 5; //M
   } else if ((a == 0 && b == 2) || (a == 2 && b == 0)) { //Het A/G
      return 6; //R
   } else if ((a == 0 && b == 3) || (a == 3 && b == 0)) { //Het A/T
      return 7; //W
   } else if ((a == 1 && b == 2) || (a == 2 && b == 1)) { //Het C/G
      return 8; //S
   } else if ((a == 1 && b == 3) || (a == 3 && b == 1)) { //Het C/T
      return 9; //Y
   } else if ((a == 2 && b == 3) || (a == 3 && b == 2)) { //Het G/T
      return 10; //K
   } else {
      return 4;
   }
}

int main(int argc, char **argv) {
   //Numbers to bases map:
   char int2bases[] = {'A', 'C', 'G', 'T', 'N', 'M', 'R', 'W', 'S', 'Y', 'K'};
   
   //Log file paths:
   string branch1snplog_path, branch2snplog_path, fai_path;
   
   //Option for debugging:
   bool debug = 0;
   
   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"input_fai", required_argument, 0, 'i'},
      {"hap1_snp_log", required_argument, 0, 'a'},
      {"hap2_snp_log", required_argument, 0, 'b'},
      {"debug", no_argument, 0, 'd'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "i:a:b:dvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'i':
            cerr << "Using FASTA .fai index: " << optarg << endl;
            fai_path = optarg;
            break;
         case 'a':
            cerr << "Using haploid 1 merged SNP log: " << optarg << endl;
            branch1snplog_path = optarg;
            break;
         case 'b':
            cerr << "Using haploid 2 merged SNP log: " << optarg << endl;
            branch2snplog_path = optarg;
            break;
         case 'd':
            cerr << "Debugging mode enabled." << endl;
            debug = 1;
            break;
         case 'v':
            cerr << "diploidizeSNPlog version " << VERSION << endl;
            return 0;
            break;
         case 'h':
            cerr << USAGE;
            return 0;
            break;
         default:
            cerr << "Unknown option " << (unsigned char)optchar << " supplied." << endl;
            cerr << USAGE;
            return 1;
            break;
      }
   }
   
   //Ignore positional arguments
   if (optind < argc) {
      cerr << "Ignoring extra positional arguments starting at " << argv[optind++] << endl;
   }
   
   //Check that log paths are set:
   if (branch1snplog_path.empty() || branch2snplog_path.empty() || fai_path.empty()) {
      cerr << "Missing one of the input logs.  Quitting." << endl;
      return 2;
   }
   
   //Open the FASTA .fai index:
   ifstream fasta_fai;
   fasta_fai.open(fai_path);
   if (!fasta_fai) {
      cerr << "Error opening FASTA .fai index file " << fai_path << ".  Quitting." << endl;
      return 3;
   }
   
   //Read in the scaffold order from the .fai file:
   vector<string> scaffolds;
   string failine;
   while (getline(fasta_fai, failine)) {
      vector<string> line_vector;
      line_vector = splitString(failine, '\t');
      scaffolds.push_back(line_vector[0]);
   }
   fasta_fai.close();
   
   //Open the haploid 1 merged SNP log:
   ifstream branch1_snp_log;
   branch1_snp_log.open(branch1snplog_path);
   if (!branch1_snp_log) {
      cerr << "Error opening haploid 1 merged SNP log " << branch1snplog_path << ".  Quitting." << endl;
      return 5;
   }
   
   //Open the haploid 2 merged SNP log:
   ifstream branch2_snp_log;
   branch2_snp_log.open(branch2snplog_path);
   if (!branch2_snp_log) {
      cerr << "Error opening haploid 2 merged SNP log " << branch2snplog_path << ".  Quitting." << endl;
      branch1_snp_log.close();
      return 6;
   }
   
   //Read branch 1 log into map (keyed by scaffold) of vectors of 3-element arrays (pos, oldallele, newallele):
   cerr << "Reading haploid 1 merged SNP log " << branch1snplog_path << endl;
   map<string, vector<array<long, 3>>> branch1_log;
   string b1line;
   while (getline(branch1_snp_log, b1line)) {
      vector<string> line_vector;
      line_vector = splitString(b1line, '\t');
      long oldallele = baseToLong(line_vector[2]);
      long newallele = baseToLong(line_vector[3]);
      if (debug && (oldallele > 3 || newallele > 3)) {
         cerr << "Found non-ACGT base in branch 1 SNP log at " << line_vector[0] << " position " << line_vector[1] << endl;
      }
      array<long, 3> log_record;
      log_record[0] = stoul(line_vector[1]);
      log_record[1] = oldallele;
      log_record[2] = newallele;
      branch1_log[line_vector[0]].push_back(log_record);
   }
   
   branch1_snp_log.close();
   cerr << "Done reading haploid 1 merged SNP log" << endl;
   
   //Read branch 2 log into map (keyed by scaffold) of vectors of 3-element arrays (pos, oldallele, newallele):
   cerr << "Reading haploid 2 merged SNP log " << branch2snplog_path << endl;
   map<string, vector<array<long, 3>>> branch2_log;
   string b2line;
   while (getline(branch2_snp_log, b2line)) {
      vector<string> line_vector;
      line_vector = splitString(b2line, '\t');
      long oldallele = baseToLong(line_vector[2]);
      long newallele = baseToLong(line_vector[3]);
      if (debug && (oldallele > 3 || newallele > 3)) {
         cerr << "Found non-ACGT base in branch 1 SNP log at " << line_vector[0] << " position " << line_vector[1] << endl;
      }
      array<long, 3> log_record;
      log_record[0] = stoul(line_vector[1]);
      log_record[1] = oldallele;
      log_record[2] = newallele;
      branch2_log[line_vector[0]].push_back(log_record);
   }
   
   branch2_snp_log.close();
   cerr << "Done reading haploid 2 merged SNP log" << endl;
   
   //Now iterate over scaffolds, outputting diploidized SNPs at any sites where either haploid deviates from ref:
   cerr << "Diploidizing SNP logs" << endl;
   for (auto scaffold_iterator = scaffolds.begin(); scaffold_iterator != scaffolds.end(); ++scaffold_iterator) {
      if (branch1_log.count(*scaffold_iterator) == 0) {
         if (branch2_log.count(*scaffold_iterator) > 0) { //Scaffold is only represented in one of the two haploids
            //Output haploid 2/ref degenerate base:
            for (auto b2_iterator = branch2_log[*scaffold_iterator].begin(); b2_iterator != branch2_log[*scaffold_iterator].end(); ++b2_iterator) {
               cout << *scaffold_iterator << '\t' << (*b2_iterator)[0] << '\t' << int2bases[(*b2_iterator)[1]] << '\t' << int2bases[degenerateBases((*b2_iterator)[2], (*b2_iterator)[1])] << endl;
            }
         }
      } else if (branch2_log.count(*scaffold_iterator) == 0) { //Scaffold must be represented in haploid 1
         //Output haploid 1/ref degenerate base
         for (auto b1_iterator = branch1_log[*scaffold_iterator].begin(); b1_iterator != branch1_log[*scaffold_iterator].end(); ++b1_iterator) {
            cout << *scaffold_iterator << '\t' << (*b1_iterator)[0] << '\t' << int2bases[(*b1_iterator)[1]] << '\t' << int2bases[degenerateBases((*b1_iterator)[2], (*b1_iterator)[1])] << endl;
         }
      } else { //Scaffold is represented in both haploids, so diploidize the scaffold
         auto b1_iterator = branch1_log[*scaffold_iterator].begin();
         auto b2_iterator = branch2_log[*scaffold_iterator].begin();
         while (b1_iterator != branch1_log[*scaffold_iterator].end() && b2_iterator != branch2_log[*scaffold_iterator].end()) {
            if ((*b1_iterator)[0] < (*b2_iterator)[0]) {
               //Output haploid 1/ref degenerate base:
               cout << *scaffold_iterator << '\t' << (*b1_iterator)[0] << '\t' << int2bases[(*b1_iterator)[1]] << '\t' << int2bases[degenerateBases((*b1_iterator)[2], (*b1_iterator)[1])] << endl;
               ++b1_iterator;
            } else if ((*b1_iterator)[0] > (*b2_iterator)[0]) {
               //Output haploid 2/ref degenerate base:
               cout << *scaffold_iterator << '\t' << (*b2_iterator)[0] << '\t' << int2bases[(*b2_iterator)[1]] << '\t' << int2bases[degenerateBases((*b2_iterator)[2], (*b2_iterator)[1])] << endl;
               ++b2_iterator;
            } else {
               //Check that ref alleles match:
               if ((*b1_iterator)[1] != (*b2_iterator)[1]) {
                  cerr << "Old alleles for site " << (*b1_iterator)[0] << " on scaffold " << *scaffold_iterator << " do not match between haploids." << endl;
                  cerr << "Haploid 1 says " << int2bases[(*b1_iterator)[1]] << " while haploid 2 says " << int2bases[(*b2_iterator)[1]] << endl;
               }
               //Diploidize the SNP:
               cout << *scaffold_iterator << '\t' << (*b1_iterator)[0] << '\t' << int2bases[(*b1_iterator)[1]] << '\t' << int2bases[degenerateBases((*b1_iterator)[2], (*b2_iterator)[2])] << endl;
               ++b1_iterator;
               ++b2_iterator;
            }
         }
         //Output the remainder of the scaffold from whichever branch still hasn't reached its end:
         while (b1_iterator != branch1_log[*scaffold_iterator].end()) {
            cout << *scaffold_iterator << '\t' << (*b1_iterator)[0] << '\t' << int2bases[(*b1_iterator)[1]] << '\t' << int2bases[degenerateBases((*b1_iterator)[2], (*b1_iterator)[1])] << endl;
            ++b1_iterator;
         }
         while (b2_iterator != branch2_log[*scaffold_iterator].end()) {
            cout << *scaffold_iterator << '\t' << (*b2_iterator)[0] << '\t' << int2bases[(*b2_iterator)[1]] << '\t' << int2bases[degenerateBases((*b2_iterator)[2], (*b2_iterator)[1])] << endl;
            ++b2_iterator;
         }
      }
   }
   cerr << "Done diploidizing SNP logs" << endl;
   
   return 0;
}
