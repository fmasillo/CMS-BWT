
#ifndef UTILS_H
#define UTILS_H

#include <string>

// struct containing command line parameters and other globals
struct Args {
  std::string filename = "";
  std::string outname = "";
  bool format = 0; // output file format: 0 full BWT | 1 run-length encoded BWT
  unsigned int buffer = 2; // buffer for internal 
  long long unsigned int prefixLength = UINT64_MAX; // prefix length for input collection file
};

#endif

