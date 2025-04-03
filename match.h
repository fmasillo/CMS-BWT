#ifndef __MATCH_H
#define __MATCH_H
#include <stdint.h>

struct MatchInVectorBigger {
  MatchInVectorBigger() : start(0), pos(0), len(0), smaller(false), toNext(0) {}
  MatchInVectorBigger(uint32_t s, uint32_t p, uint32_t l, bool sm, uint32_t i)
      : start(s), pos(p), len(l), smaller(sm), toNext(0), idx(i) {}

  uint32_t start; // position in the collection
  uint32_t pos;   // position of match in _reference
  uint32_t len;   // length of match
  bool smaller;
  uint32_t toNext;
  uint32_t idx;
};

struct MatchInSet {
  MatchInSet() : len(0), smaller(0), isaNext(0) {}
  MatchInSet(uint32_t l, bool sm, uint32_t iN)
      : len(l), smaller(sm), isaNext(iN) {}

  bool operator==(const MatchInSet &other) const {
    return (len == other.len) && (isaNext == other.isaNext);
  }

  bool operator<(const MatchInSet &other) const {
    if (len != other.len) {
      return smaller * (len < other.len) + !other.smaller * (len > other.len);
    } else {
      return isaNext < other.isaNext;
    }
  }

  uint32_t len;
  bool smaller;
  uint32_t isaNext;
};

struct ItemMatchInSet {
  ItemMatchInSet() : untilNext(0), rank(0) {}
  ItemMatchInSet(uint32_t u, uint32_t i) : untilNext(u), rank(0) {
  idxs.push_back(i);
  }

  void addIdx(uint32_t i) { idxs.push_back(i); }

  void changeUntilNext(uint32_t u) { untilNext = u; }

  void changeRank(uint32_t r) { rank = r; }

  uint32_t untilNext;
  uint32_t rank;
  std::vector<uint32_t> idxs;
};


#endif