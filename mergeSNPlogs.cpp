/**********************************************************************************
 * mergeSNPlogs.cpp                                                               *
 * Written by Patrick Reilly                                                      *
 * Version 1.0 written 2017/01/16                                                 *
 * Version 1.1 written 2018/09/07 Variety of bug fixes                            *
 * Description:                                                                   *
 *                                                                                *
 * Syntax: mergeSNPlogs [branch 1 indel log] [branch 1 SNP log] [branch 2 SNP log]*
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
#define VERSION "1.1"

//Usage/help:
#define USAGE "mergeSNPlogs\nUsage:\n mergeSNPlogs -i [branch 1 indel log] -b [branch 1 SNP log] -c [branch 2 SNP log]\n"

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

bool constructIndelMap(ifstream &indel_log, map<string, vector<pair<long, long>>> &indel_map) {
   bool readfail = 0;
   string logline;
   long cumulativechange = 0;
   while(getline(indel_log, logline)) {
      vector<string> line_vector;
      line_vector = splitString(logline, '\t');
      if (indel_map.count(line_vector[0]) == 0) {
         cumulativechange = 0; //Make sure to reset the change on a new scaffold
         indel_map[line_vector[0]].push_back(make_pair(0, 0)); //Every mapping starts with 0,0
      }
      long indel_size = stol(line_vector[3]);
      if (indel_size == 0) { //Skip indels of size 0, they don't affect coordinate space mapping
         continue;
      }
      long indel_change = line_vector[2] == "ins" ? indel_size : -indel_size;
      cumulativechange += indel_change;
      long ref_position = stol(line_vector[1]) + 1;
      long new_position = ref_position + cumulativechange;
      indel_map[line_vector[0]].push_back(make_pair(ref_position, new_position));
   }
   readfail = (indel_log.fail() || indel_log.bad()) && !indel_log.eof(); //Only signal failure if fail bit or bad bit are set, but EOF bit is not.
   return readfail;
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

int main(int argc, char **argv) {
   //Numbers to bases map:
   char int2bases[] = {'A', 'C', 'G', 'T', 'N', 'M', 'R', 'W', 'S', 'Y', 'K'};
   
   //Map for coordinate space change due to indels:
   map<string, vector<pair<long, long>>> indelmap;
   
   //Log file paths:
   string branch1snplog_path, branch2snplog_path, indellog_path;
   
   //Option for debugging:
   bool debug = 0;
   
   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"indel_log", required_argument, 0, 'i'},
      {"branch1_snp_log", required_argument, 0, 'b'},
      {"branch2_snp_log", required_argument, 0, 'c'},
      {"debug", no_argument, 0, 'd'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "i:b:c:dvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'i':
            cerr << "Using branch 1 indel log: " << optarg << endl;
            indellog_path = optarg;
            break;
         case 'b':
            cerr << "Using branch 1 SNP log: " << optarg << endl;
            branch1snplog_path = optarg;
            break;
         case 'c':
            cerr << "Using branch 2 SNP log: " << optarg << endl;
            branch2snplog_path = optarg;
            break;
         case 'd':
            cerr << "Debugging mode enabled." << endl;
            debug = 1;
            break;
         case 'v':
            cerr << "mergeSNPlogs version " << VERSION << endl;
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
   if (indellog_path.empty() || branch1snplog_path.empty() || branch2snplog_path.empty()) {
      cerr << "Missing one of the input logs.  Quitting." << endl;
      return 2;
   }
   
   //Open the branch 1 indel log:
   ifstream indellog;
   indellog.open(indellog_path);
   if (!indellog) {
      cerr << "Error opening branch 1 indel log " << indellog_path << ".  Quitting." << endl;
      return 3;
   }
   
   bool indelmapfail = constructIndelMap(indellog, indelmap);
   if (indelmapfail) {
      cerr << "Failed to construct coordinate-space mapping.  Quitting." << endl;
      return 4;
   }
   indellog.close();
   if (debug) {
      for (auto imprint_iterator = indelmap.begin(); imprint_iterator != indelmap.end(); ++imprint_iterator) {
         for (auto indel_iterator = imprint_iterator->second.begin(); indel_iterator != imprint_iterator->second.end(); ++indel_iterator) {
            cerr << imprint_iterator->first << '\t' << indel_iterator->first << '\t' << indel_iterator->second << endl;
         }
      }
   }
   
   //Open the SNP logs:
   ifstream branch1_snp_log, branch2_snp_log;
   branch1_snp_log.open(branch1snplog_path);
   if (!branch1_snp_log) {
      cerr << "Error opening branch 1 SNP log " << branch1snplog_path << ".  Quitting." << endl;
      return 5;
   }
   branch2_snp_log.open(branch2snplog_path);
   if (!branch2_snp_log) {
      cerr << "Error opening branch 2 SNP log " << branch2snplog_path << ".  Quitting." << endl;
      branch1_snp_log.close();
      return 6;
   }
   
   //Read branch 1 log into map (keyed by scaffold) of vectors of 3-element arrays (pos, oldallele, newallele):
   cerr << "Reading branch 1 SNP log " << branch1snplog_path << endl;
   string firstScaffold;
   map<string, vector<array<long, 3>>> branch1_log;
   string b1line;
   while (getline(branch1_snp_log, b1line)) {
      vector<string> line_vector;
      line_vector = splitString(b1line, '\t');
      if (firstScaffold.empty()) {
         firstScaffold = line_vector[0];
      }
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
   cerr << "Done reading branch 1 SNP log" << endl;
   
   //Now iterate over branch 2 log, adjusting new position back to old position using indel map, 
   //then comparing to branch 1 log to look for overlapping changes that need to be transitively reduced:
   cerr << "Reading branch 2 SNP log " << branch2snplog_path << endl;
   string b2line;
   auto indelmap_iterator = indelmap[firstScaffold].begin();
   auto left_iterator = indelmap[firstScaffold].begin();
   auto b1log_iterator = branch1_log[firstScaffold].begin();
   while (getline(branch2_snp_log, b2line)) {
      vector<string> line_vector;
      line_vector = splitString(b2line, '\t');
      if (line_vector[0] != firstScaffold) {
         indelmap_iterator = indelmap[line_vector[0]].begin();
         left_iterator = indelmap[line_vector[0]].begin();
         b1log_iterator = branch1_log[line_vector[0]].begin();
         firstScaffold = line_vector[0];
      }
      long oldallele = baseToLong(line_vector[2]);
      long newallele = baseToLong(line_vector[3]);
      if (debug && (oldallele > 3 || newallele > 3)) {
         cerr << "Found non-ACGT base in branch 2 SNP log at " << line_vector[0] << " position " << line_vector[1] << endl;
      }
      //Adjust the position back into the branch 1 source coordinate space:
      long newref_position = stol(line_vector[1]);
      while (newref_position > indelmap_iterator->second && indelmap_iterator != indelmap[line_vector[0]].end()) {
         left_iterator = indelmap_iterator;
         ++indelmap_iterator;
      }
      auto right_iterator = indelmap_iterator;
      if (right_iterator == indelmap[line_vector[0]].end()) {
         --right_iterator;
      }
      if (indelmap_iterator == indelmap[line_vector[0]].end() || (newref_position < indelmap_iterator->second && indelmap_iterator != indelmap[line_vector[0]].begin())) {
         --indelmap_iterator;
      }
      long left_ins_flank = left_iterator->second + right_iterator->first - left_iterator->first;
      long right_ins_flank = right_iterator->second;
      if (debug) {
         cerr << newref_position << '\t' << indelmap_iterator->first - indelmap_iterator->second << '\t' << left_iterator->first << '\t' << right_iterator->first << '\t' << left_iterator->second << '\t' << right_iterator->second << '\t' << left_ins_flank << '\t' << right_ins_flank << endl;
      }
      if (newref_position > left_ins_flank && newref_position < right_ins_flank) {
         //Mutation along branch 2 is within insertion on branch 1
         if (debug) {
            cerr << "Mutation along branch 2 is within insertion on branch 1 at unadjusted position " << line_vector[0] << ":" << newref_position << endl;
         }
         continue;
      }
      long adjusted_position = newref_position + indelmap_iterator->first - indelmap_iterator->second;
      //if (adjusted_position < indelmap_iterator->second && indelmap_iterator->second > 0) { // We went too far in the indel map
      //   --indelmap_iterator;
      //} else {
      //   while (adjusted_position <= indelmap_iterator->second && indelmap_iterator != indelmap[line_vector[0]].end()) {
      //      ++indelmap_iterator;
      //   }
      //   //Make sure to check the left boundary condition, drop back to last interval with start <= adjusted_position:
      //   if (adjusted_position > indelmap_iterator->second && indelmap_iterator != indelmap[line_vector[0]].begin()) {
      //      --indelmap_iterator;
      //   }
      //}
      //adjusted_position += indelmap_iterator->first - indelmap_iterator->second;
      //If this scaffold isn't present in the branch 1 SNP log, output the adjusted record:
      if (branch1_log.count(line_vector[0]) == 0) {
         cout << line_vector[0] << '\t' << adjusted_position << '\t' << int2bases[oldallele] << '\t' << int2bases[newallele] << endl;
      } else {
         while ((*b1log_iterator)[0] < adjusted_position && b1log_iterator != branch1_log[line_vector[0]].end()) { //Output branch 1-exclusive events
            cout << line_vector[0] << '\t' << (*b1log_iterator)[0] << '\t' << int2bases[(*b1log_iterator)[1]] << '\t' << int2bases[(*b1log_iterator)[2]] << endl;
            ++b1log_iterator;
         }
         if (b1log_iterator == branch1_log[line_vector[0]].end()) { //No more branch 1 records, so short-circuit outputting branch 2 records
            cout << line_vector[0] << '\t' << adjusted_position << '\t' << int2bases[oldallele] << '\t' << int2bases[newallele] << endl;
         } else if ((*b1log_iterator)[0] == adjusted_position) { //Transitively reduce this record
            if (debug && (*b1log_iterator)[2] != oldallele) { //Transitive mismatch, output an error if debug mode is on
               cerr << "Allele mismatch during transitive reduction at " << line_vector[0] << " position " << adjusted_position << endl;
               cerr << "Branch 1 says " << int2bases[(*b1log_iterator)[1]] << "->" << int2bases[(*b1log_iterator)[2]] << endl;
               cerr << "Branch 2 says " << int2bases[oldallele] << "->" << int2bases[newallele] << endl;
            }
            cout << line_vector[0] << '\t' << adjusted_position << '\t' << int2bases[(*b1log_iterator)[1]] << '\t' << int2bases[newallele] << endl;
         } else { //Only branch 2 record at this position, so output it
            cout << line_vector[0] << '\t' << adjusted_position << '\t' << int2bases[oldallele] << '\t' << int2bases[newallele] << endl;
         }
      }
   }
   
   branch2_snp_log.close();
   cerr << "Done reading branch 2 SNP log" << endl;

   
   return 0;
}
