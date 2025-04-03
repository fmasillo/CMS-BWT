
#ifndef CMS_BWT_H
#define CMS_BWT_H

#include <stdlib.h>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <string.h>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include "utils.h"
#include "rmq_tree.h"
#include "match.h"

#include <libsais.h>

using data_type = uint8_t;
using filelength_type = uint64_t;

#define sequenceSeparator (char)2
#define MAXIMUM_UINT32 std::numeric_limits<uint32_t>::max()

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

static std::string outputFile;

//static data_type *_x;
static std::string _x;
static int *_SA;
static uint32_t *_ISA;
static int32_t *_LCP;
static int32_t *_PLCP;
static rmq_tree *_rmq;
static uint8_t *_BWT;
static uint32_t _n;
//static data_type *_sx;
static filelength_type _sn = 0;
static bool _isMismatchingSymbolNeeded;
static std::vector<uint32_t> headBoundaries;

static bool verbose;

static uint16_t sizeChars = 256;
static uint64_t D = 1;


void constructISA(int32_t *sa, uint32_t *isa, uint32_t n);
std::pair<int,int> adjustInterval(int lo, int hi, int offset);
std::pair<int,int> fasterContractLeft(int lo, int hi, int offset);
std::pair<int,int> contractLeft(int lo, int hi, int offset);
void computeMSFactorAt(const std::string &_sx, filelength_type i, uint32_t *pos, uint32_t *len, int32_t & leftB, int32_t & rightB, bool & isSmallerThanMaxMatch, unsigned char &mismatchingSymbol);
inline int32_t binarySearchLB(int32_t lo, int32_t hi, uint32_t offset, data_type c);
inline int32_t binarySearchRB(int32_t lo, int32_t hi, uint32_t offset, data_type c);


void initialize_reference(Args arg, std::string &refFileName, std::string &collFileName, uint64_t prefixLength);
void process_collection_small_reference(Args arg, std::string &collFileName);
void process_collection_large_reference(Args arg, std::string &collFileName);

void computeBWT(Args arg, std::string &refFileName, std::string &collFileName);

#endif