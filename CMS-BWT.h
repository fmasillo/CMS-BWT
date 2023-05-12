
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
#include <string.h>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include "utils.h"
#include "rmq_tree.h"
#include "match.h"
#include "predecessor.h"
#include "libsais/include/libsais.h"

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

static uint64_t timesLength = 0;
static uint64_t timesISA = 0;
static uint64_t timesRank = 0;

static std::vector<Match> phrases;
static std::vector<MatchSA> headsSA;
static std::vector<MatchSANN> headsSANN;


void constructLCP(std::string t, int32_t n, int32_t *sa, uint32_t *lcp, uint32_t *temp);
void constructISA(int32_t *sa, uint32_t *isa, uint32_t n);
std::pair<int,int> adjustInterval(int lo, int hi, int offset);
std::pair<int,int> fasterContractLeft(int lo, int hi, int offset);
std::pair<int,int> contractLeft(int lo, int hi, int offset);
bool compareMatchSA(const MatchSA &a, const MatchSA &b);
inline bool compareSufToHead(const MatchSA &headA, const uint64_t distanceLeft, const MatchSA &headB);
void computeLZFactorAt(const std::string &_sx, filelength_type i, uint32_t *pos, uint32_t *len, int32_t & leftB, int32_t & rightB, bool & isSmallerThanMaxMatch, unsigned char &mismatchingSymbol);
inline int32_t binarySearchLB(int32_t lo, int32_t hi, uint32_t offset, data_type c);
inline int32_t binarySearchRB(int32_t lo, int32_t hi, uint32_t offset, data_type c);


void lzInitialize(char *refFileName, char *collFileName, uint64_t prefixLength);
int lzFactorize(Args arg, char *collFileName);
int lzFactorizeMemorySaving(Args arg, char *collFileName);

void computeBWT(Args arg, char *refFileName, char *collFileName);
void computeBWTMemorySaving(Args arg, char *refFileName, char *collFileName);

#endif