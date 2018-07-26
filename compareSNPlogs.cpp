/**********************************************************************************
 * compareSNPlogs.cpp                                                             *
 * Written by Patrick Reilly                                                      *
 * Version 1.0 written 2017/01/24                                                 *
 * Version 1.1 written 2017/01/24                                                 *
 * Version 1.2 written 2017/03/02                                                 *
 * Version 1.3 written 2017/09/21                                                 *
 * Version 1.4 written 2018/05/03                                                 *
 * Description:                                                                   *
 *                                                                                *
 * Syntax: compareSNPlogs -i [.fai] -e [expected SNP log] -o [in.snp file]        *
 *         -n [output false negatives to this file] -p [output false positives]   *
 *         -t [output true positives] -r [erroneous sites]                        *
 **********************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cctype>
#include <vector>
#include <map>
#include <sstream>
#include <array>
#include <unordered_set>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define VERSION "1.4"

//Usage/help:
#define USAGE "compareSNPlogs\nUsage:\n compareSNPlogs -i [FASTA .fai] -e [expected SNP log] -o [observed in.snp]\n\t-n [output false negative in.snp] -p [output false positive in.snp]\n\t-t [output true positive in.snp] -r [output erroneous call in.snp]\n\t--min_depth [minimum callable depth]\n"

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
      case 'a':
         return 0;
      case 'C':
      case 'c':
         return 1;
      case 'G':
      case 'g':
         return 2;
      case 'T':
      case 't':
         return 3;
      case 'M':
      case 'm':
         return 5;
      case 'R':
      case 'r':
         return 6;
      case 'W':
      case 'w':
         return 7;
      case 'S':
      case 's':
         return 8;
      case 'Y':
      case 'y':
         return 9;
      case 'K':
      case 'k':
         return 10;
      default:
         return 4;
   }
}

void splitBase(long base_value, vector<long> &output) {
   switch(base_value) {
      case 0:
         output.push_back(0);
         output.push_back(0);
         break;
      case 1:
         output.push_back(1);
         output.push_back(1);
         break;
      case 2:
         output.push_back(2);
         output.push_back(2);
         break;
      case 3:
         output.push_back(3);
         output.push_back(3);
         break;
      case 4:
         output.push_back(4);
         output.push_back(4);
         break;
      case 5:
         output.push_back(0);
         output.push_back(1);
         break;
      case 6:
         output.push_back(0);
         output.push_back(2);
         break;
      case 7:
         output.push_back(0);
         output.push_back(3);
         break;
      case 8:
         output.push_back(1);
         output.push_back(2);
         break;
      case 9:
         output.push_back(1);
         output.push_back(3);
         break;
      case 10:
         output.push_back(2);
         output.push_back(3);
         break;
      default:
         output.push_back(4);
         output.push_back(4);
         break;
   }
   return;
}


int main(int argc, char **argv) {
   //Numbers to bases map:
   char int2bases[] = {'A', 'C', 'G', 'T', 'N', 'M', 'R', 'W', 'S', 'Y', 'K'};
   
   //Log file paths:
   string expected_path, observed_path, fai_path;

   //False negative and positive output file paths:
   string fn_path = "", fp_path = "";
   //True positive output file path:
   string tp_path = "";
   //Erroneous call output file path:
   string error_path = "";
   
   //Minimum depth to consider a SNP callable:
   unsigned long min_depth = 0;
      
   //Option for debugging:
   bool debug = 0;
   
   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"input_fai", required_argument, 0, 'i'},
      {"expected_snps", required_argument, 0, 'e'},
      {"observed_insnp", required_argument, 0, 'o'},
      {"output_fns", required_argument, 0, 'n'},
      {"output_fps", required_argument, 0, 'p'},
      {"output_tps", required_argument, 0, 't'},
      {"output_errors", required_argument, 0, 'r'},
      {"min_depth", required_argument, 0, 'm'},
      {"debug", no_argument, 0, 'd'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "i:e:o:n:p:t:r:m:dvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'i':
            cerr << "Using FASTA .fai index: " << optarg << endl;
            fai_path = optarg;
            break;
         case 'e':
            cerr << "Using expected SNP log: " << optarg << endl;
            expected_path = optarg;
            break;
         case 'o':
            cerr << "Using observed in.snp: " << optarg << endl;
            observed_path = optarg;
            break;
         case 'n':
            cerr << "Outputting false negative sites to: " << optarg << endl;
            fn_path = optarg;
            break;
         case 'p':
            cerr << "Outputting false positive sites to: " << optarg << endl;
            fp_path = optarg;
            break;
         case 't':
            cerr << "Outputting true positive sites to: " << optarg << endl;
            tp_path = optarg;
            break;
         case 'r':
            cerr << "Outputting erroneous call sites to: " << optarg << endl;
            error_path = optarg;
            break;
         case 'm':
            cerr << "Ignoring true SNPs with raw depth less than " << optarg << endl;
            min_depth = stoul(optarg);
            break;
         case 'd':
            cerr << "Debugging mode enabled." << endl;
            debug = 1;
            break;
         case 'v':
            cerr << "compareSNPlogs version " << VERSION << endl;
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
   if (expected_path.empty() || observed_path.empty() || fai_path.empty()) {
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
   
   //Read in the scaffold order and scaffold lengths from the .fai file:
   vector<string> scaffolds;
   map<string, unsigned long> scaffold_lengths;
   unsigned long genome_size = 0;
   string failine;
   while (getline(fasta_fai, failine)) {
      vector<string> line_vector;
      line_vector = splitString(failine, '\t');
      scaffolds.push_back(line_vector[0]);
      unsigned long scaffold_length = stoul(line_vector[1]);
      scaffold_lengths[line_vector[0]] = scaffold_length;
      genome_size += scaffold_length;
   }
   fasta_fai.close();
   
   //Open the expected SNP log:
   ifstream expected;
   expected.open(expected_path);
   if (!expected) {
      cerr << "Error opening expected SNP log " << expected_path << ".  Quitting." << endl;
      return 5;
   }
   
   //Open the observed in.snp file:
   ifstream observed;
   observed.open(observed_path);
   if (!observed) {
      cerr << "Error opening observed in.snp " << observed_path << ".  Quitting." << endl;
      observed.close();
      return 6;
   }
   
   //Read expected SNP log into map (keyed by scaffold) of vectors of 3-element arrays (pos, oldallele, newallele):
   cerr << "Reading expected SNP log " << expected_path << endl;
   map<string, vector<array<long, 3>>> expected_log;
   unordered_set<string> uncallable_sites;
   string eline;
   while (getline(expected, eline)) {
      vector<string> line_vector;
      line_vector = splitString(eline, '\t');
      long oldallele = baseToLong(line_vector[2]);
      long newallele = baseToLong(line_vector[3]);
      if (debug && (oldallele > 3 || newallele > 3)) {
         cerr << "Found non-ACGT base in branch 1 SNP log at " << line_vector[0] << " position " << line_vector[1] << endl;
      }
      if (min_depth > 0) {
         if (line_vector.size() < 5) {
            cerr << "Error: Used non-zero minimum callable depth, but no depths provided in expected log." << endl;
            observed.close();
            expected.close();
            return 7;
         }
         unsigned long curdepth = stoul(line_vector[4]);
         if (curdepth < min_depth) { //Skip sites that wouldn't be callable based on the raw sequencing depth
            uncallable_sites.insert(line_vector[0] + ":" + line_vector[1]);
            continue;
         }
      }
      array<long, 3> log_record;
      log_record[0] = stol(line_vector[1]);
      log_record[1] = oldallele;
      log_record[2] = newallele;
      expected_log[line_vector[0]].push_back(log_record);
   }
   
   expected.close();
   cerr << "Done reading expected SNP log" << endl;
   
   //Read observed in.snp into map (keyed by scaffold) of vectors of 3-element arrays of strings (pos, oldallele, newallele):
   cerr << "Reading observed in.snp file " << observed_path << endl;
   map<string, vector<array<string, 3>>> observed_log;
   string oline;
   while (getline(observed, oline)) {
      vector<string> line_vector;
      line_vector = splitString(oline, '\t');
      array<string, 3> log_record;
      log_record[0] = line_vector[1];
      log_record[1] = line_vector[2];
      log_record[2] = line_vector[3];
      observed_log[line_vector[0]].push_back(log_record);
   }
   
   observed.close();
   cerr << "Done reading observed in.snp file" << endl;
   
   //If the false negative output file path was input, open that up:
   ofstream fn_file;
   if (!fn_path.empty()) {
      fn_file.open(fn_path);
      if (!fn_file) {
         cerr << "Unable to open false negative output file, so ignoring that function." << endl;
         fn_path = "";
      }
   }

   //If the false positive output file path was input, open that up:
   ofstream fp_file;
   if (!fp_path.empty()) {
      fp_file.open(fp_path);
      if (!fp_file) {
         cerr << "Unable to open false positive output file, so ignoring that function." << endl;
         fp_path = "";
      }
   }
   
   //If the true positive output file path was input, open that up:
   ofstream tp_file;
   if (!tp_path.empty()) {
      tp_file.open(tp_path);
      if (!tp_file) {
         cerr << "Unable to open true positive output file, so ignoring that function." << endl;
         tp_path = "";
      }
   }
   
   //If the erroneous call output file path was input, open that up:
   ofstream error_file;
   if (!error_path.empty()) {
      error_file.open(error_path);
      if (!error_file) {
         cerr << "Unable to open erroneous call output file, so ignoring that function." << endl;
         error_path = "";
      }
   }
      
   //Now iterate over scaffolds, counting FP and FN variant calls, ignoring masking and indels in in.snp:
   cerr << "Comparing SNP logs" << endl;
   unsigned long tps = 0, fps = 0, tns = 0, fns = 0, wrong_calls = 0;
   //Further categorize into match and mismatch types (first letter is call, second is truth, R=ref, H=het, A=alt):
   unsigned long masked_bases = 0, indel_sites = 0;
   unsigned long RR_match = 0, RH_mismatch = 0, RA_mismatch = 0;
   unsigned long HR_mismatch = 0, HH_match = 0, HH_mismatch = 0, HA_mismatch = 0;
   unsigned long AR_mismatch = 0, AH_mismatch = 0, AA_match = 0, AA_mismatch = 0;
   unsigned long NR_masked = 0, NH_masked = 0, NA_masked = 0;
   unsigned long IR_masked = 0, IH_masked = 0, IA_masked = 0;
   for (auto scaffold_iterator = scaffolds.begin(); scaffold_iterator != scaffolds.end(); ++scaffold_iterator) {
      unsigned long scaffold_tps = 0, scaffold_tns = 0, scaffold_fps = 0, scaffold_fns = 0, scaffold_wrong_calls = 0;
      if (expected_log.count(*scaffold_iterator) == 0) {
         if (observed_log.count(*scaffold_iterator) > 0) { //Scaffold is only represented in observed in.snp file
            //Count false positives:
            for (auto o_iterator = observed_log[*scaffold_iterator].begin(); o_iterator != observed_log[*scaffold_iterator].end(); ++o_iterator) {
               //Skip indels or masked bases:
               if ((*o_iterator)[1].length() > 1 || (*o_iterator)[2].length() > 1) {
                  indel_sites += 1; //Directly at indel site, guaranteed to be ref since scaffold not in expected in.snp
                  IR_masked += 1;
                  continue;
               } else if ((*o_iterator)[2] == "N") { //Not indel site, but masked (possibly in window around indel)
                  masked_bases += 1;
                  NR_masked += 1;
                  continue;
               }
               if (baseToLong((*o_iterator)[2]) > 4) { //Het call, truth is hom ref
                  HR_mismatch += 1;
               } else { //Hom alt call, truth is hom ref
                  AR_mismatch += 1;
               }
               scaffold_fps += 2; //Every non-N, non-indel record in observed in.snp not in the expected SNP log is an FP
               if (!fp_path.empty()) { //Record false positive site to log if requested
                  fp_file << *scaffold_iterator << '\t' << (*o_iterator)[0] << '\t' << (*o_iterator)[1] << '\t' << (*o_iterator)[2] << endl;
               }
            }
         }
      } else if (observed_log.count(*scaffold_iterator) == 0) { //Scaffold must be represented in expected SNP log
         //Count false negatives:
         for (auto e_iterator = expected_log[*scaffold_iterator].begin(); e_iterator != expected_log[*scaffold_iterator].end(); ++e_iterator) {
            if ((*e_iterator)[2] > 4) { //Hom ref call, truth is het
               RH_mismatch += 1;
            } else { //Hom ref call, truth is hom alt
               RA_mismatch += 1;
            }
            scaffold_fns += 2; //Every record in the expected SNP log not in the observed in.snp file is a false negative
            if (!fn_path.empty()) { //Record false negative site to log if requested
               fn_file << *scaffold_iterator << '\t' << (*e_iterator)[0] << '\t' << int2bases[(*e_iterator)[1]] << '\t' << int2bases[(*e_iterator)[2]] << endl;
            }
         }
      } else { //Scaffold is represented in both logs, so compare contents:
         auto e_iterator = expected_log[*scaffold_iterator].begin();
         auto o_iterator = observed_log[*scaffold_iterator].begin();
         while (e_iterator != expected_log[*scaffold_iterator].end() && o_iterator != observed_log[*scaffold_iterator].end()) {
            if ((*e_iterator)[0] < stol((*o_iterator)[0])) { //Observed in.snp file is missing this SNP
               if ((*e_iterator)[2] > 4) { //Hom ref call, truth is het
                  RH_mismatch += 1;
               } else { //Hom ref call, truth is hom alt
                  RA_mismatch += 1;
               }
               //Count false negative:
               scaffold_fns += 2;
               if (!fn_path.empty()) { //Record false negative site to log if requested
                  fn_file << *scaffold_iterator << '\t' << (*e_iterator)[0] << '\t' << int2bases[(*e_iterator)[1]] << '\t' << int2bases[(*e_iterator)[2]] << endl;
               }
               ++e_iterator;
            } else if ((*e_iterator)[0] > stol((*o_iterator)[0])) { //Expected SNP log does not contain this SNP
               if ((*o_iterator)[1].length() > 1 || (*o_iterator)[2].length() > 1) { //Indel getting masked
                  IR_masked += 1;
                  indel_sites += 1;
               } else if ((*o_iterator)[2] == "N") { //Masked base
                  NR_masked += 1;
                  masked_bases += 1;
               } else if (baseToLong((*o_iterator)[2]) > 4) { //Het call, truth is hom ref
                  HR_mismatch += 1;
               } else { //Hom alt call, truth is hom ref
                  AR_mismatch += 1;
               }
               if (uncallable_sites.count(*scaffold_iterator + ":" + (*o_iterator)[0]) == 0) { //This is a callable site
                  //Count false positive if non-indel and not masked:
                  if ((*o_iterator)[1].length() == 1 && (*o_iterator)[2].length() == 1 && (*o_iterator)[2] != "N") {
                     scaffold_fps += 2;
                     if (!fp_path.empty()) { //Record false positive site to log if requested
                        fp_file << *scaffold_iterator << '\t' << (*o_iterator)[0] << '\t' << (*o_iterator)[1] << '\t' << (*o_iterator)[2] << endl;
                     }
                  }
               }
               ++o_iterator;
            } else { //Both files have this record, so compare the values
               //Check that ref alleles match:
               if (to_string(int2bases[(*e_iterator)[1]]) != (*o_iterator)[1] && debug) {
                  cerr << "Ref alleles for site " << (*e_iterator)[0] << " on scaffold " << *scaffold_iterator << " do not match between SNP logs." << endl;
                  cerr << "Expected SNP log says " << int2bases[(*e_iterator)[1]] << " while observed in.snp says " << (*o_iterator)[1] << endl;
               }
               //Compare the values:
               if ((*e_iterator)[2] == baseToLong((*o_iterator)[2])) { //True positive
                  if ((*e_iterator)[2] > 4) { //Matching het call
                     HH_match += 1;
                  } else { //Matching hom alt call
                     AA_match += 1;
                  }
                  scaffold_tps += 2;
                  if (!tp_path.empty()) { //Record true positive site to log if requested
                     tp_file << *scaffold_iterator << '\t' << (*e_iterator)[0] << '\t' << int2bases[(*e_iterator)[2]] << '\t' << (*o_iterator)[2] << endl;
                  }
               } else if (baseToLong((*o_iterator)[2]) == 4) { //Masked base => false negative
                  if ((*o_iterator)[1].length() > 1 || (*o_iterator)[2].length() > 1) { //Indel masking
                     if ((*e_iterator)[2] > 4) { //Indel masked het site
                        IH_masked += 1;
                     } else { //Indel masked hom alt site
                        IA_masked += 1;
                     }
                     indel_sites += 1;
                  } else {
                     if ((*e_iterator)[2] > 4) { //Masked het site
                        NH_masked += 1;
                     } else { //Masked hom alt site
                        NA_masked += 1;
                     }
                     masked_bases += 1;
                  }
                  scaffold_fns += 2;
                  if (!fn_path.empty()) { //Record false negative site to log if requested
                     fn_file << *scaffold_iterator << '\t' << (*e_iterator)[0] << '\t' << int2bases[(*e_iterator)[1]] << '\t' << int2bases[(*e_iterator)[2]] << endl;
                  }
               } else { //Error (Does this count as FP or FN?)
                  if ((*e_iterator)[2] > 4) { //Truth is het
                     if (baseToLong((*o_iterator)[2]) > 4) { //Wrong het
                        HH_mismatch += 1;
                     } else { //Called hom alt, truth is het
                        AH_mismatch += 1;
                     }
                  } else { //Truth is hom alt
                     if (baseToLong((*o_iterator)[2]) > 4) { //Called het, truth is hom alt
                        HA_mismatch += 1;
                     } else { //Wrong hom alt
                        AA_mismatch += 1;
                     }
                  }
                  //Split the expected and observed bases into two:
                  vector<long> x, y;
                  splitBase((*e_iterator)[2], x);
                  splitBase(baseToLong((*o_iterator)[2]), y);
                  long ref_base = (*e_iterator)[1];
                  if (y[0] == x[0] || y[0] == x[1] || y[1] == x[0] || y[1] == x[1]) { //At least one match, so add a tp
                     scaffold_tps++;
                     if (y[0] == x[0] || y[0] == x[1]) { //Condition on y0 being the base that matched
                        if (y[1] == ref_base) {
                           scaffold_fns++;
                        } else {
                           scaffold_wrong_calls++;
                        }
                     } else { //Condition on y1 being the base that matched
                        if (y[1] == ref_base) {
                           scaffold_fns++;
                        } else {
                           scaffold_wrong_calls++;
                        }
                     }
                  } else if (y[0] == ref_base || y[1] == ref_base) {
                     scaffold_fns++;
                     scaffold_wrong_calls++;
                  } else {
                     scaffold_wrong_calls += 2;
                  }
                  if (!error_path.empty()) { //Record erroneous call site to log if requested
                     error_file << *scaffold_iterator << '\t' << (*e_iterator)[0] << '\t' << int2bases[(*e_iterator)[2]] << '\t' << (*o_iterator)[2] << endl;
                  }
               }
               ++e_iterator;
               ++o_iterator;
            }
         }
         //Count the remainder of the scaffold from whichever log still hasn't reached its end:
         while (e_iterator != expected_log[*scaffold_iterator].end()) {
            if ((*e_iterator)[2] > 4) { //Called hom ref, truth is het
               RH_mismatch += 1;
            } else { //Called hom ref, truth is hom alt
               RA_mismatch += 1;
            }
            //Count false negatives:
            scaffold_fns += 2;
            if (!fn_path.empty()) { //Record false negative site to log if requested
               fn_file << *scaffold_iterator << '\t' << (*e_iterator)[0] << '\t' << int2bases[(*e_iterator)[1]] << '\t' << int2bases[(*e_iterator)[2]] << endl;
            }
            ++e_iterator;
         }
         while (o_iterator != observed_log[*scaffold_iterator].end()) {
            if ((*o_iterator)[1].length() > 1 || (*o_iterator)[2].length() > 1) { //Indel masked hom ref
               IR_masked += 1;
               indel_sites += 1;
            } else if ((*o_iterator)[2] == "N") { //Masked hom ref
               NR_masked += 1;
               masked_bases += 1;
            } else if (baseToLong((*o_iterator)[2]) > 4) { //Called het, truth is hom ref
               HR_mismatch += 1;
            } else { //Called hom alt, truth is hom ref
               AR_mismatch += 1;
            }
            //Count false positive if non-indel and not masked:
            if (uncallable_sites.count(*scaffold_iterator + ":" + (*o_iterator)[0]) == 0) { //This is a callable site
               if ((*o_iterator)[1].length() == 1 && (*o_iterator)[2].length() == 1 && (*o_iterator)[2] != "N") {
                  scaffold_fps += 2;
                  if (!fp_path.empty()) { //Record false positive site to log if requested
                     fp_file << *scaffold_iterator << '\t' << (*o_iterator)[0] << '\t' << (*o_iterator)[1] << '\t' << (*o_iterator)[2] << endl;
                  }
               }
            }
            ++o_iterator;
         }
      }
      tps += scaffold_tps;
      fps += scaffold_fps;
      fns += scaffold_fns;
      wrong_calls += scaffold_wrong_calls;
      scaffold_tns = 2*scaffold_lengths[*scaffold_iterator] - scaffold_tps - scaffold_fns - scaffold_wrong_calls - scaffold_fps;
      tns += scaffold_tns;
   }
   tns -= 2*uncallable_sites.size(); //Don't count uncallable sites as anything, including TN.
   unsigned long R_mismatches = RH_mismatch + RA_mismatch;
   unsigned long H_mismatches = HR_mismatch + HH_mismatch + HA_mismatch;
   unsigned long A_mismatches = AR_mismatch + AH_mismatch + AA_mismatch;
   unsigned long mismatches = R_mismatches + H_mismatches + A_mismatches;
   RR_match = genome_size - indel_sites - masked_bases - mismatches - HH_match - AA_match;
   if (!fn_path.empty()) {
      fn_file.close();
   }
   if (!fp_path.empty()) {
      fp_file.close();
   }
   if (!tp_path.empty()) {
      tp_file.close();
   }
   cerr << "Done comparing SNP logs" << endl;
   
   cout << setprecision(15);
   cout << "True positives\t" << (double)tps/2.0 << endl;
   cout << "False positives\t" << (double)fps/2.0 << endl;
   cout << "True negatives\t" << (double)tns/2.0 << endl;
   cout << "False negatives\t" << (double)fns/2.0 << endl;
   cout << "Wrong calls\t" << (double)wrong_calls/2.0 << endl;
   cout << "FPR\t" << (double)fps/(double)(fps+tns) << endl;
   cout << "FNR\t" << (double)fns/(double)(fns+tps) << endl;
   cout << "FNR+wrong\t" << (double)(fns+wrong_calls)/(double)(fns+wrong_calls+tps) << endl;
   cout << "Wrong call rate (wrong calls out of all calls)\t" << (double)(wrong_calls)/(double)(wrong_calls+tps+fps) << endl;
   cout << "Sensitivity\t" << (double)tps/(double)(tps+fns) << endl;
   cout << "Specificity\t" << (double)tns/(double)(tns+fps) << endl;
   cout << "FDR\t" << (double)fps/(double)(tps+fps) << endl;
   cout << endl;
   cout << "Call types:" << endl;
   cout << "Masked\t" << (double)masked_bases << endl;
   cout << "Indel site\t" << (double)indel_sites << endl;
   cout << "Homozygous ref\t" << (double)(RR_match + R_mismatches) << endl;
   cout << "Heterozygous\t" << (double)(HH_match + H_mismatches) << endl;
   cout << "Homozygous alt\t" << (double)(AA_match + A_mismatches) << endl;
   cout << endl;
   cout << "Matches:" << endl;
   cout << "Homozygous ref\t" << (double)RR_match << endl;
   cout << "Heterozygous\t" << (double)HH_match << endl;
   cout << "Homozygous alt\t" << (double)AA_match << endl;
   cout << endl;
   cout << "Mismatches:" << endl;
   cout << "Het->RR\t" << (double)RH_mismatch << endl;
   cout << "Alt->RR\t" << (double)RA_mismatch << endl;
   cout << "RR->Het\t" << (double)HR_mismatch << endl;
   cout << "Het->Other Het\t" << (double)HH_mismatch << endl;
   cout << "Alt->Het\t" << (double)HA_mismatch << endl;
   cout << "RR->Alt\t" << (double)AR_mismatch << endl;
   cout << "Het->Alt\t" << (double)AH_mismatch << endl;
   cout << "Alt->Other Alt\t" << (double)AA_mismatch << endl;
   cout << endl;
   cout << "Masking:" << endl;
   cout << "RR->N\t" << (double)NR_masked << endl;
   cout << "Het->N\t" << (double)NH_masked << endl;
   cout << "Alt->N\t" << (double)NA_masked << endl;
   cout << endl;
   cout << "Indel Sites:" << endl;
   cout << "RR->Indel\t" << (double)IR_masked << endl;
   cout << "Het->Indel\t" << (double)IH_masked << endl;
   cout << "Alt->Indel\t" << (double)IA_masked << endl;
   
   return 0;
}
