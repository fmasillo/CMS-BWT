
#ifndef UTILS_H
#define UTILS_H

#include <string>

// struct containing command line parameters and other globals
struct Args {
  std::string filename = "";
  std::string outname = "";
  bool format = 0; // output file format: 0 full BWT | 1 run-length encoded BWT
  bool memory_saving = 0; // memory saving mode: 0 no | 1 yes
  unsigned int buffer = 2; // buffer for internal 
  long long unsigned int prefixLength = UINT64_MAX; // prefix length for input collection file
  bool check = 0; // check the correctness of the output
};

#endif

