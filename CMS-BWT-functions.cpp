#include "CMS-BWT.h"

//LCP array construction method of J. Kärkkäinen, G. Manzini, and S. J. Puglisi, CPM 2009
void constructLCP(std::string t, int32_t n, int32_t *sa, uint32_t *lcp, uint32_t *temp) {
//void constructLCP(unsigned char *t, int32_t n, uint32_t *sa, uint32_t *lcp, uint32_t *temp) {
   //fprintf(stderr,"\tComputing LCP...\n");
   std::cerr << "\tComputing LCP..." << std::endl;
   int32_t *phi = (int32_t *)lcp, *plcp = (int32_t *)temp, l = 0;
   for (int32_t i = 1; i < n; ++i)
      phi[sa[i]] = sa[i-1];
   phi[sa[0]] = -1;
   for (int32_t i = 0; i < n; ++i) {
      int32_t j = phi[i];
      if (j == -1) { plcp[i] = 0; continue; }
      else {
         while (i + l < n && j + l < n && t[i + l] == t[j + l]) ++l;
         plcp[i] = l;
         l = std::max(l - 1, 0);
      }
   }
   for (int32_t i = 0; i < n; ++i)
     lcp[i] = plcp[sa[i]];
}

void constructISA(int32_t *sa, uint32_t *isa, uint32_t n){
   //fprintf(stderr,"\tComputing ISA...\n");
   std::cerr << "\tComputing ISA..." << std::endl;
   for(uint32_t i=0; i<n; i++){
      isa[sa[i]] = i;
   }
}

std::pair<int,int> adjustInterval(int lo, int hi, int offset) {
  int psv = _rmq->psv(lo,offset);
  if(psv == -1){
    psv = 0;
  }else{
    //psv++;
  }
  int nsv = _rmq->nsv(hi+1,offset);
  if(nsv == -1){
    nsv = _n-1;
  }else{
    nsv--;
  }
  return {psv,nsv};
}

std::pair<int,int> fasterContractLeft(int lo, int hi, int offset){
   // uint32_t suflo = _SA[lo];
   uint32_t sufhi = _SA[hi];
   // if(suflo == _n-1 || sufhi == _n-1){ //if true we must be at depth 1
   //    return std::make_pair(0,_n-1); //root
   // }
   uint32_t tmplo = _ISA[lo+1];
   uint32_t tmphi = _ISA[sufhi+1];
   return adjustInterval(tmplo, tmphi, offset);
}

std::pair<int,int> contractLeft(int lo, int hi, int offset){
   uint32_t suflo = _SA[lo];
   uint32_t sufhi = _SA[hi];
   if(suflo == _n-1 || sufhi == _n-1){ //if true we must be at depth 1
      return std::make_pair(0,_n-1); //root
   }
   uint32_t tmplo = _ISA[suflo+1];
   uint32_t tmphi = _ISA[sufhi+1];
   return adjustInterval(tmplo, tmphi, offset);
}

bool compareMatchSA(const MatchSA &a, const MatchSA &b){
   if((a.smaller == 0) & (b.smaller == 0)){
      return (a.next < b.next)*(a.len == b.len) + (a.len < b.len)*(a.len != b.len);
   }
   else if((a.smaller == 1) & (b.smaller == 1)){
      return (a.next < b.next)*(a.len == b.len) + (a.len > b.len)*(a.len != b.len);
   }
   return a.smaller < b.smaller;
}

inline bool compareSufToHead(const MatchSA &headA, const uint64_t distanceLeft, const MatchSA &headB){
   if(distanceLeft != headB.len){
      return headA.smaller*(distanceLeft < headB.len) +
            !headB.smaller*(distanceLeft > headB.len);
   }
   else{
      //MatchSA headA2 = headsSA[headA.start];
      //return (headA2.pos != headB.pos)*(headA2.pos < headB.pos) + (headA2.pos == headB.pos)*(headA2.start < headB.start);
      return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
   }
}

void computeLZFactorAt(const std::string &_sx, filelength_type i, uint32_t *pos, uint32_t *len, int32_t & leftB, int32_t & rightB, bool & isSmallerThanMaxMatch, unsigned char &mismatchingSymbol){//, uint32_t *backup) {
   uint32_t offset = *len;
   filelength_type j = i + offset;

   //int64_t nlb = 0, nrb = _n - 1;
   int32_t nlb = leftB, nrb = rightB, maxMatch;
   unsigned int match = _SA[nlb];
   while (j < _sn) { //scans the string from j onwards until a maximal prefix is found between _x and _sx
      if (nlb == nrb) {
         //std::cout << "Search finished with nlb == nrb " << nlb << " " << nrb << "\n";
         if (_x[_SA[nlb] + offset] != _sx[j]) {
               //std::cout << "mismatch " << _x[_SA[nlb] + offset] << " and " << _sx[j] << "\n";
               //fprintf(stderr,"Breaking from 1\n");
               isSmallerThanMaxMatch = (_x[_SA[nlb] + offset] > _sx[j]);
               mismatchingSymbol = _sx[j];
               break;
         }
         leftB = nlb;
         rightB = nrb;
         maxMatch = nlb;
      } else { //refining the bucket in which the match is found, from left and then from right
         nlb = binarySearchLB(nlb, nrb, offset, _sx[j]);
         // nlb = std::lower_bound(rangeSA+nlb, rangeSA+nrb+1, _sx[j], [offset](const uint32_t &a, const unsigned char &b){
         //    return _x[_SA[a]+offset] < b;
         // }) - rangeSA;
         if(nlb < 0) {
         //if(_x[_SA[nlb] + offset] != _sx[j] | nlb == nrb+1){
            //no match, the game is up
            maxMatch = -nlb-1;
            //maxMatch = nlb;
            isSmallerThanMaxMatch = true;
            mismatchingSymbol = _sx[j];
            if(maxMatch == nrb+1){
               maxMatch--;
               isSmallerThanMaxMatch = false;
            }
            match = _SA[maxMatch];
            break;
         }
         nrb = binarySearchRB(nlb, nrb, offset, _sx[j]);
         // nrb = std::upper_bound(rangeSA+nlb, rangeSA+nrb+1, _sx[j], [offset](const unsigned char &a, const uint32_t &b){
         //    return a < _x[_SA[b]+offset];
         // }) - rangeSA - 1;
         leftB = nlb;
         rightB = nrb;
      }
      //std::cerr << "After if nlb: " << nlb << "\n";
      match = _SA[nlb];
      j++;
      offset++;
   }
   *pos = match;
   *len = offset;
   //std::cout << "Match " << match << "\n";
   //std::cout << "Len " << offset << "\n";
}

//Returns the leftmost occurrence of the element if it is present or (if not present)
//then returns -(x+1) where x is the index at which the key would be inserted into the
//array: i.e., the index of the first element greater than the key, or hi+1 if all elements
//in the array are less than the specified key.
inline int32_t binarySearchLB(int32_t lo, int32_t hi, uint32_t offset, data_type c) {
   int32_t low = lo, high = hi;
   while (low <= high) {
      int32_t mid = (low + high) >> 1;
      data_type midVal = _x[_SA[mid] + offset];
      if (midVal < c){
         low = mid + 1;
         __builtin_prefetch(&_x[_SA[(low + high) >> 1] + offset], 0, 0);
      }
      else if (midVal > c){
         high = mid - 1;
         __builtin_prefetch(&_x[_SA[(low + high) >> 1] + offset], 0, 0);
      }
      else { //midVal == c
         if (mid == lo)
            return mid; // leftmost occ of key found
         data_type midValLeft = _x[_SA[mid - 1] + offset];
         if (midValLeft == midVal) {
            high = mid - 1; //discard mid and the ones to the right of mid
            __builtin_prefetch(&_x[_SA[(low + high) >> 1] + offset], 0, 0);
         } else { //midValLeft must be less than midVal == c
            return mid; //leftmost occ of key found
         }
      }
   }
   return -(low + 1);  // key not found.
}


inline int32_t binarySearchRB(int32_t lo, int32_t hi, uint32_t offset, data_type c) {
   int32_t low = lo, high = hi;
   while (low <= high) {
      int32_t mid = (low + high) >> 1;
      data_type midVal = _x[_SA[mid] + offset];
      if (midVal < c){
         low = mid + 1;
         __builtin_prefetch(&_x[_SA[(low + high) >> 1] + offset], 0, 0);
      }
      else if (midVal > c){
         high = mid - 1;
         __builtin_prefetch(&_x[_SA[(low + high) >> 1] + offset], 0, 0);
      }
      else { //midVal == c
         if (mid == hi)
            return mid; // rightmost occ of key found
         data_type midValRight = _x[_SA[mid + 1] + offset];
         if (midValRight == midVal) {
            low = mid + 1; //discard mid and the ones to the left of mid
            __builtin_prefetch(&_x[_SA[(low + high) >> 1] + offset], 0, 0);
         } else { //midValRight must be greater than midVal == c
            return mid; //rightmost occ of key found
         }
      }
   }
   return -(low + 1);  // key not found.
}


void lzInitialize(char *refFileName, char *collFileName, uint64_t prefixLength) {
//void lzInitialize(data_type *ax, unsigned int an, bool isMismatchingSymbolNeeded, char *refFileName, char *collFileName) {
   auto t1 = std::chrono::high_resolution_clock::now();
   errno = 0;
   FILE *infileRef = fopen(refFileName, "r");
   if (!infileRef) {
      //fprintf(stderr, "Error opening file of base sequence %s, errno=%d\n", refFileName,errno);
      std::cerr << "Error opening file of base sequence " << refFileName << ", errno=" << errno << '\n';
      exit(1);
   }
   //fprintf(stderr, "About to read ref\n");
   std::cerr << "About to read ref\n";

   unsigned int n = 0;
   //data_type *x = new data_type[n + 1];
   fseek(infileRef, 0, SEEK_END);
   n = ftell(infileRef) / sizeof(data_type);
   std::cerr << "n = " << n << '\n';
   fseek(infileRef, 0, SEEK_SET);
   if(n){
      char *firstChar = new char[1];
      int o = fread(firstChar, sizeof(char), 1, infileRef);
      //std::cerr << "firstChar: " << firstChar << '\n';
      if(firstChar[0] == '>'){
         _x.reserve(n);
         //std::cerr << "Yes FASTA\n";
         std::ifstream streamInfileRef(refFileName);
         std::string line, content;
         while(std::getline(streamInfileRef, line).good()){
            if( line.empty() || line[0] == '>' ){
               _x += content;
               _x += '$';
               //std::cerr << _x << "size: " << _x.size() << "\n";
               std::string().swap(content);
            }
            else if (!line.empty()) {
               content += line;
            }
         }
         if(content.size()) _x += content;
         //std::cerr << _x << "size: " << _x.size() << "\n";
         std::string().swap(content);
         _x.resize(_x.size());
         streamInfileRef.close();
      }
      else{
         //std::cerr << "No FASTA\n";
         fseek(infileRef, 0, SEEK_SET);
         _x.resize(n);
         if (n != fread(&_x[0], sizeof(data_type), n, infileRef)) {
            //fprintf(stderr, "Error reading %u bytes from file %s\n", n, refFileName);
            std::cerr << "Error reading " << n << " bytes from file " << refFileName << '\n';
            exit(1);
         }
      }
      //x[n] = 0;
      delete[] firstChar;
   }
   else{
      std::cerr << "Reference file is empty!\n";
      exit(1);
   }
   fclose(infileRef);
   //fprintf(stderr, "Reference (size = %lu):\n\t", _x.size());
   std::cerr << "Reference (size = " << _x.size() << "):\n\t";

   auto t01 = std::chrono::high_resolution_clock::now();
   //_x = std::string(reinterpret_cast<char const *>(x), n);
   if((_x[_x.size()-1] == '\n') | (_x[_x.size()-1] == '\r') | (_x[_x.size()-1] == 0)){
      _x.erase(_x.size()-1);
   }
   if(_x[_x.size()-1] == '$'){
      _x.erase(_x.size()-1);
   }
   
   FILE *infile = fopen(collFileName, "r");
   if(!infile){
      //fprintf(stderr, "Error opening file of sequence (%s)\n", collFileName);
      std::cerr << "Error opening file of sequence (" << collFileName << ")\n";
      exit(1);
   }
   filelength_type sn = 0;
   fseek(infile, 0L, SEEK_END);
   std::cerr << "ftello(infile): " << ftello(infile) << '\n';
   sn = std::min(ftello(infile) / sizeof(data_type), prefixLength);
   std::cerr << "sn: " << sn << '\n';
   //std::string refAug(_x);
   fclose(infile);
   _sn = sn;
   //augmenting the reference with Ns and other characters
   if(_x.find('N') == std::string::npos) _x.append(_x.size()/10, 'N');

   for(uint16_t i = 3; i < sizeChars; i++){
      //if(verbose) std::cerr << (char)i << " ref: " << maxRunsReference[i] << ", coll: " << maxRunsCollection[i] << "\n"; 
      _x.append(1, (data_type)i);     
   }

   _x += (char)1;
   _x += (char)0;

   
   //std::cerr << refAug << "\n";
   auto t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Augmenting reference done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   int32_t *sa = new int32_t[_x.size()];
   libsais(reinterpret_cast<unsigned char*>(const_cast<char*>(_x.c_str())), sa, _x.size(), 0, NULL);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing SA done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   _n = _x.size();

   t01 = std::chrono::high_resolution_clock::now();
   _SA = sa;
   _ISA = new uint32_t[_n];
   _PLCP = new int32_t[_n];
   _LCP = new int32_t[_n + 1];
   _LCP[_n] = -1;
   t01 = std::chrono::high_resolution_clock::now();
   constructISA(_SA,_ISA,_n);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing ISA done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   libsais_plcp(reinterpret_cast<unsigned char*>(const_cast<char*>(_x.c_str())), _SA, _PLCP, _x.size());
   libsais_lcp(_PLCP, _SA, _LCP, _x.size());
   //constructLCP(_x,_n,_SA,_LCP,_PLCP);
   for(uint32_t i = 0; i < _n; i++){
      _PLCP[i] = std::max(_LCP[_ISA[i]],_LCP[_ISA[i]+1]);
   }
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing LCP and PLCP done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   //fprintf(stderr,"\tComputing RMQ...\n");
   std::cerr << "Computing RMQ...\n";
   _rmq = new rmq_tree((int *)_LCP, _n, 7);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing RMQ done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   //fprintf(stderr,"\tComputing BWT...\n");
   std::cerr << "Computing BWT...\n";
   _BWT = new uint8_t[_n];
   for(uint32_t i = 0; i < _n; i++){
      _BWT[i] = _SA[i]-1 >= 0 ? _x[_SA[i]-1] : char(0);
   }
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing BWT done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";


   //_isMismatchingSymbolNeeded = isMismatchingSymbolNeeded;
   std::cerr << "Finished pre-processing\n";
   // rangeSA = new int32_t[_x.size()];
   // for(uint32_t i = 0; i < _x.size(); i++){
   //    rangeSA[i] = i;
   // }
   auto t2 = std::chrono::high_resolution_clock::now();
   uint64_t preprocTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Preprocessing done in " << preprocTime << " ms\n";
   //sprintf(collFileName, "%s.gsa", collFileName);
   //outputFile = std::string(collFileName);
   //return _sx;
}


int lzFactorize(Args arg, char *collFileName) {
//int lzFactorize(char *fileToParse, int seqno, char* outputfilename, const bool v) {
   verbose = 0;
   //omp_set_num_threads(nThreads);
   std::cerr << "File is in memory\n";
   //std::cerr << "_x: " <<  _x << '\n';
   //std::cerr << "_sx: " <<  _sx << '\n';
   std::cerr << "_x.len: " <<  _n << '\n';
   std::cerr << "_sx.len: " <<  _sn << '\n';
   auto t1 = std::chrono::high_resolution_clock::now();

   unsigned int numfactors = 0;

   unsigned int inc = 100000000;
   uint64_t mark = inc;

   std::cerr << "About to start main parsing loop...\n";

   
   uint32_t ndoc = 1;
   int64_t *bucketsForExpandedBWT = new int64_t[_n]();
   bucketsForExpandedBWT[_n - 1] = 0;
   uint64_t *bucketISAPosHeads = new uint64_t[_n]();

   headBoundaries.push_back(0);
   uint64_t i = 0;

   uint64_t maxValue = 0;
   phrases.reserve(_sn / sizeof(Match));
   std::vector<uint8_t> listOfChars;
   listOfChars.reserve(_sn / sizeof(Match));
   //listOfChars.push_back('%');
   
   std::ifstream streamInfile(collFileName,std::ios::in);
   std::string line, content;
   content.reserve(_n*2);
   uint64_t charactersRead = 0;
   //std::string _sx_just_to_check;
   //_sx_just_to_check.reserve(_sn);
   while(std::getline(streamInfile, line).good()){
      //std::cerr << "ndoc: " << ndoc << "\n";
      if( line.empty() || line[0] == '>' ){
         //std::cerr << "content: " << content << "\n";
         content += sequenceSeparator;
         charactersRead++;
         //_sx_just_to_check += content;
         int64_t i = 0;
         int32_t leftB = 0;
         int32_t rightB = _n-1;
         bool isSmallerThanMatch;
         unsigned char mismatchingSymbol;
         uint64_t prevPos = -2;
         uint32_t pos = _n - 1, len = 0;
         uint32_t iCurrentDoc = 0;
         D++;
         //std::cerr << "content.size(): " << content.size() << "\n";
         //std::cerr << "((int64_t)content.size())-1: " << ((int64_t)content.size())-1 << "\n";
         //std::cerr << "i < ((int64_t)content.size())-1: " << (i < ((int64_t)content.size())-1) << "\n";
         
         // if(content.size() >= 5){
         //    //std::cerr << "content[0]: " << content[0] << " content[1]: " << content[1] << " content[2]: " << content[2] << " content[3]: " << content[3] << "\n";
         //    std::pair<int32_t, int32_t> pair = lookupTableSA[convertNucleotideQuadrupletHash(content[0],content[1],content[2],content[3])];
         //    //std::cerr << "pair.first: " << pair.first << " pair.second: " << pair.second << "\n";
         //    if(pair.first != -1){
         //       leftB = pair.first;
         //       rightB = pair.second;
         //       len += 4;
         //    }
         // }
         while(i < ((int64_t)content.size())-1){
            // if(len < 4){
            //    //std::cerr << "len is less than 4!\n";
            //    std::pair<int32_t, int32_t> pair = lookupTableSA[convertNucleotideQuadrupletHash(content[i],content[i+1],content[i+2],content[i+3])];
            //    if(pair.first != -1){
            //       leftB = pair.first;
            //       rightB = pair.second;
            //       len = 4;
            //    }
            // }
            computeLZFactorAt(content, i, &pos, &len, leftB, rightB, isSmallerThanMatch, mismatchingSymbol);//, &backup);
            if(pos != prevPos+1){
               phrases.push_back(Match(iCurrentDoc, pos, len, isSmallerThanMatch, mismatchingSymbol));
               if(i == 0) listOfChars.push_back(sequenceSeparator);
               else listOfChars.push_back(content[i-1]);
               if(bucketsForExpandedBWT[pos] > 0) bucketsForExpandedBWT[pos] = -bucketsForExpandedBWT[pos]-1;
               else bucketsForExpandedBWT[pos]--;
            }
            else{
               if(bucketsForExpandedBWT[pos] >= 0) bucketsForExpandedBWT[pos]++;
               else bucketsForExpandedBWT[pos]--;
            }
            iCurrentDoc++;
            len--;

            if(leftB == rightB){
               while(len > _PLCP[pos+1]){
                  i++;
                  iCurrentDoc++;
                  len--;
                  pos++;
                  if(bucketsForExpandedBWT[pos] >= 0) bucketsForExpandedBWT[pos]++;
                  else bucketsForExpandedBWT[pos]--;
                  //heurYes++;
               }
               leftB = rightB = _ISA[pos];
               std::pair<int,int> interval = adjustInterval(_ISA[pos+1], _ISA[pos+1], len);//fasterContractLeft(pos,rightB,len);
               //fasterAdjustInterval++;
               leftB = interval.first;
               rightB = interval.second;
            }
            else{
               std::pair<int,int> interval = contractLeft(leftB, rightB, len);//fasterContractLeft(pos,rightB,len);
               //slowerAdjustInterval++;
               leftB = interval.first;
               rightB = interval.second;
            }
            
            //heurNo++;
            i++;
            prevPos = pos;
         }
         len = 0;
         pos = _n - 1;
         phrases.push_back(Match(iCurrentDoc, pos, len, 0, sequenceSeparator));
         bucketsForExpandedBWT[pos]--;
         headBoundaries.push_back(phrases.size());
         if(maxValue < iCurrentDoc) maxValue = iCurrentDoc;
         iCurrentDoc = 0;
         ndoc++;
         if(i == 0) {
            listOfChars.push_back(sequenceSeparator);
         }
         else{
            listOfChars.push_back(content[content.size()-2]);
         }
         //std::cerr << "Finished parsing document " << ndoc << "\n";
         content.erase();
         //std::string().swap(content);
      }
      else if (!line.empty()) {
         charactersRead += line.size();
         if(charactersRead >= _sn - 1){ //if string is filled up to sn (useful for prefixLength)
            //std::cerr << "read too much\n";
            //std::cerr << "charactersRead: " << charactersRead << "\n";
            //std::cerr << "line.size(): " << line.size() << "\n";
            //std::cerr << "_sn - charactersRead + line.size() - 1 = " << line.size() - (charactersRead - _sn) - 1 << "\n";
            //std::cerr << "content.size(): " << content.size() << "\n";
            content += line.substr(0, line.size() - (charactersRead - _sn) - 1);
            //std::cerr << "content.size(): " << content.size() << "\n";
            //std::cerr << _sx_just_to_check.size() << "\n";
            break;
         }
         else{
            content += line;
         }
      }
   }
   streamInfile.close();
   if(content.size() != 0){
      content += sequenceSeparator;
      charactersRead++;
      if(charactersRead < _sn - 1){
         _sn = charactersRead;
      }
      //_sx_just_to_check += content;
      //std::cerr << "Finishing up with " << content.size() << " characters\n";
      //std::cerr << "content: " << content << "\n";
      D++;
      int64_t i = 0;
      int32_t leftB = 0;
      int32_t rightB = _n-1;
      bool isSmallerThanMatch;
      unsigned char mismatchingSymbol;
      uint64_t prevPos = -2;
      uint32_t pos = _n - 1, len = 0;
      uint32_t iCurrentDoc = 0;
      while(i < content.size()){
         if(content[i] == sequenceSeparator){
            leftB = 0;
            rightB = _n-1;
            len = 0;
            pos = _n - 1;
            phrases.push_back(Match(iCurrentDoc, pos, len, 0, sequenceSeparator));
            bucketsForExpandedBWT[pos]--;
            headBoundaries.push_back(phrases.size());
            if(maxValue < iCurrentDoc) maxValue = iCurrentDoc;
            iCurrentDoc = 0;
            ndoc++;
            if(i == 0) {
               listOfChars.push_back(sequenceSeparator);
            }
            else{
               listOfChars.push_back(content[i - 1]);
            }
         }
         else{
            computeLZFactorAt(content, i, &pos, &len, leftB, rightB, isSmallerThanMatch, mismatchingSymbol);//, &backup);
            if(pos != prevPos+1){
               phrases.push_back(Match(iCurrentDoc, pos, len, isSmallerThanMatch, mismatchingSymbol));
               if(i == 0) listOfChars.push_back(sequenceSeparator);
               else listOfChars.push_back(content[i-1]);
               if(bucketsForExpandedBWT[pos] > 0) bucketsForExpandedBWT[pos] = -bucketsForExpandedBWT[pos]-1;
               else bucketsForExpandedBWT[pos]--;
            }
            else{
               if(bucketsForExpandedBWT[pos] >= 0) bucketsForExpandedBWT[pos]++;
               else bucketsForExpandedBWT[pos]--;
            }
            iCurrentDoc++;
            len--;

            if(leftB == rightB){
               while(len > _PLCP[pos+1]){
                  i++;
                  iCurrentDoc++;
                  len--;
                  pos++;
                  if(bucketsForExpandedBWT[pos] >= 0) bucketsForExpandedBWT[pos]++;
                  else bucketsForExpandedBWT[pos]--;
                  //heurYes++;
               }
               leftB = rightB = _ISA[pos];
               std::pair<int,int> interval = adjustInterval(_ISA[pos+1], _ISA[pos+1], len);//fasterContractLeft(pos,rightB,len);
               //slowerAdjustInterval++;
               leftB = interval.first;
               rightB = interval.second;
            }
            else{
               std::pair<int,int> interval = contractLeft(leftB, rightB, len);//fasterContractLeft(pos,rightB,len);
               //slowerAdjustInterval++;
               leftB = interval.first;
               rightB = interval.second;
            }
         }
         //heurNo++;
         i++;
         prevPos = pos;
      }
      std::string().swap(content);
   }
   std::cerr << "Finished parsing file of length: " << _sn << "\n";
   //std::cerr << "_sx_just_to_check = " << _sx_just_to_check << "\n";
   //std::cerr << "_sx_just_to_check.size(): " << _sx_just_to_check.size() << "\n";
   // std::cerr << "fasterAdjustInterval: " << fasterAdjustInterval << "\n";
   // std::cerr << "slowerAdjustInterval: " << slowerAdjustInterval << "\n";
   // std::cerr << "heurYes: " << heurYes << "\n";
   delete[] _LCP;
   delete[] _PLCP;
   delete _rmq;
   phrases.shrink_to_fit();
   listOfChars.shrink_to_fit();
   //std::string().swap(_x);
   std::cerr << "headBoundaries.size(): " << headBoundaries.size() << ", ndoc: " << ndoc << "\n";
   //for(uint64_t i = 0; i < headBoundaries.size(); i++) std::cerr << "i: " << i << " - " << headBoundaries[i] << "\n";
   //for(uint64_t i = 0; i < docBoundaries.size(); i++) std::cerr << "i: " << i << " - " << docBoundaries[i] << "\n";
   std::cerr << "docBoundaries.size(): " << D << "\n";
   std::cerr << "phrases.size() = " << phrases.size() << "\n";
   if(verbose) for(uint64_t i = 0; i < phrases.size(); i++) std::cerr << "i:" << i  << " - " << phrases[i].start << "," << phrases[i].pos << "," << phrases[i].len << "," << phrases[i].smaller << ", " << phrases[i].next << '\n';
   predecessor2 *pHeads = new predecessor2(phrases, headBoundaries, ndoc, maxValue);
   uint64_t *bucketPosHeads = new uint64_t[_n]();
   for(uint64_t i = 0; i < phrases.size(); i++){
      bucketPosHeads[phrases[i].pos]++;
      bucketISAPosHeads[_ISA[phrases[i].pos]]++;
   }
   // docBoundaries.push_back(_sn);
   if(verbose) std::cerr << "Printing docBoundaries" << "\n";
   //if(verbose) for(size_t i = 0; i < docBoundaries.size(); i++){ std::cerr << docBoundaries[i] << ", letter: " << _sx[docBoundaries[i]] << "\n";}
   auto t2 = std::chrono::high_resolution_clock::now();
   uint64_t lzFactorizeTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to compute matching statistics: " << lzFactorizeTime << " milliseconds\n";

   t1 = std::chrono::high_resolution_clock::now();
   //Sort suffixes corrisponding to heads
   uint64_t *prefSumBucketPosHeads = new uint64_t[_n]();
   uint64_t *prefSumBucketPosHeadsCopy = new uint64_t[_n+1];
   uint32_t *prefSumBucketISAPosHeads = new uint32_t[_n];
   uint32_t *prefSumBucketISAPosHeadsCopy = new uint32_t[_n+1];
   prefSumBucketISAPosHeads[0] = 0;
   prefSumBucketISAPosHeadsCopy[0] = 0;
   for(size_t i = 1; i < _n; i++){
      prefSumBucketPosHeads[i] = prefSumBucketPosHeads[i-1] + bucketPosHeads[i-1];
      //std::cerr << "prefSumBucketPosHeads[" << i << "] = " << prefSumBucketPosHeads[i] << "\n";
      prefSumBucketPosHeadsCopy[i] = prefSumBucketPosHeads[i];
      prefSumBucketISAPosHeads[i] = prefSumBucketISAPosHeads[i-1] + bucketISAPosHeads[i-1];
      prefSumBucketISAPosHeadsCopy[i] = prefSumBucketISAPosHeads[i];
   }
   delete[] bucketPosHeads;
   delete[] bucketISAPosHeads;
   prefSumBucketPosHeadsCopy[_n] = phrases.size();
   prefSumBucketISAPosHeadsCopy[_n] = phrases.size();

   headsSA.resize(phrases.size());
   //headsSA.reserve(phrases.size());
   i = 0, ndoc = 0;
   auto t01 = std::chrono::high_resolution_clock::now();
   for(std::vector<Match>::iterator it = phrases.begin(); it < phrases.end(); it++){
      if(it->start == 0){
         ndoc++;
      }
      //headsSA[prefSumBucketISAPosHeads[_ISA[it->pos]]++] = MatchSA(i++, it->pos, it->len, !it->smaller, _sx[it->start + docBoundaries[ndoc-1] + it->len]);
      headsSA[prefSumBucketISAPosHeads[_ISA[it->pos]]++] = MatchSA(i++, it->pos, it->len, !it->smaller, it->next);
   }
   auto t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Bucketing heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   if(verbose) std::cerr << "Outputting headsSA bucketed (size=" << headsSA.size() << ")\n";
   std::vector<MatchSA>::iterator begHeads = headsSA.begin();
   t01 = std::chrono::high_resolution_clock::now();

   for(size_t i = 1; i < _n + 1; i++){
      std::sort(begHeads + prefSumBucketISAPosHeadsCopy[i-1], begHeads + prefSumBucketISAPosHeadsCopy[i], compareMatchSA);
   }
   delete[] prefSumBucketISAPosHeads;
   //delete[] prefSumBucketISAPosHeadsCopy;

   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Head sort took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to re-rank heads\n";
   int32_t *stringHeads = new int32_t[headsSA.size()+1];
   uint64_t rank = 1;
   int64_t prevLen = -1, prevSmall = -1, prevNext = -1;
   uint64_t prevPos = -1;
   i = 0;

   uint32_t *uniqueHeadsPos = new uint32_t[_n]();
   std::vector<uint32_t> uniqueHeadsPosVec;
   std::vector<uint32_t> offSet;
   std::vector<uint32_t> uniqueHeads;
   uniqueHeadsPosVec.reserve(headsSA.size());
   //uniqueHeadsPosVec.push_back(headsSA[0].pos);
   offSet.reserve(headsSA.size());
   //offSet.push_back(0);
   uniqueHeads.reserve(headsSA.size());
   //uniqueHeads.push_back(headsSA[0].len);
   //uniqueHeadsPos[headsSA[0].pos] = 1;
   uint32_t equalHeads = 0;
   for(std::vector<MatchSA>::iterator it = headsSA.begin(); it < headsSA.end(); it++){
      if(i >= headBoundaries.size()-1){
         if(it->pos != prevPos){
            equalHeads = 0;
         }
         if(it->pos != prevPos || it->len != prevLen || it->smaller != prevSmall){
            offSet.push_back(equalHeads);
            uniqueHeadsPos[it->pos]++;
            uniqueHeadsPosVec.push_back(it->pos);
            uniqueHeads.push_back(it->smaller ? MAXIMUM_UINT32 - it->len : it->len);
            //std::cerr << "Unique head: " << it->pos << " " << it->len << " " << it->smaller << " " << it->next << ". Storing: " << uniqueHeads.back() <<"\n";
            //std::cerr << "OffSet[" << offSet.end()-offSet.begin()-1 << "] = " << offSet.back() << "\n";
         }
         equalHeads++;
         //std::cerr << "it->start: " << it->start << " phrases[it->start].start: " << phrases[it->start].start << "\n";
         //std::cerr << "currentMaxNoOverlap: " << currentMaxNoOverlap << "\n";
         if(it->pos != prevPos || it->len != prevLen || it->smaller != prevSmall || it->next != prevNext){
            rank++;
            prevPos = it->pos;
            prevLen = it->len;
            prevSmall = it->smaller;
            prevNext = it->next;
         }
      }
      stringHeads[i++] = rank;
   }

   std::cerr << "offSet.back() = " << offSet.back() << " equalHeads = " << equalHeads << "\n";
   if(offSet.back() != equalHeads-1){
      offSet.push_back(equalHeads-1);
      //uniqueHeadsPosVec.push_back(headsSA.back().pos);
      //uniqueHeadsPos[headsSA.back().pos]++;
      //uniqueHeads.push_back(headsSA.back().smaller ? MAXIMUM_UINT32 - headsSA.back().len : headsSA.back().len);
      //std::cerr << "Unique head: " << headsSA.back().pos << " " << headsSA.back().len << " " << headsSA.back().smaller << " " << headsSA.back().next << ". Storing: " << uniqueHeads.back() <<"\n";
      //std::cerr << "OffSet[" << offSet.end()-offSet.begin()-1 << "] = " << offSet.back() << "\n";
   }
   uniqueHeadsPosVec.shrink_to_fit();
   offSet.shrink_to_fit();
   uniqueHeads.shrink_to_fit();
   std::cerr << "Unique heads: " << uniqueHeads.size() << "\n";
   std::cerr << "Unique heads positions: " << uniqueHeadsPosVec.size() << "\n";
   std::cerr << "OffSet size: " << offSet.size() << "\n";
   std::cerr << "Re-ranked heads\n";
   std::cerr << "Maximum rank: " << rank << "\n";

   uint32_t *uniqueHeadsPosPrefSum = new uint32_t[_n];
   uniqueHeadsPosPrefSum[0] = 0;
   for(size_t i = 1; i < _n; i++){
      uniqueHeadsPosPrefSum[i] = uniqueHeadsPosPrefSum[i-1] + uniqueHeadsPos[i-1];
   }

   //permute uniqueHeads and offSet according to uniqueHeadsPosPrefSum using uniqueHeadsPosVec
   std::vector<uint32_t> uniqueHeadsPerm;
   std::vector<uint32_t> offSetPerm;
   uniqueHeadsPerm.resize(uniqueHeads.size());
   offSetPerm.resize(offSet.size());
   for(size_t i = 0; i < uniqueHeadsPosVec.size(); i++){
      uniqueHeadsPerm[uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]] = uniqueHeads[i];
      offSetPerm[uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]] = offSet[i];
      //uniqueHeadsPerm.push_back(uniqueHeads[uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]]);
      //offSetPerm.push_back(offSet[uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]]);
      uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]++;
   }
   uniqueHeads = uniqueHeadsPerm;
   offSet = offSetPerm;
   // //print uniqueHeads and offSet
   // std::cerr << "Unique heads:\n";
   // for(size_t i = 0; i < uniqueHeads.size(); i++){
   //    std::cerr << uniqueHeads[i] << " ";
   // }
   // std::cerr << "\nOffSet:\n";
   // for(size_t i = 0; i < offSet.size()-1; i++){
   //    std::cerr << offSet[i] << " ";
   // }
   // std::cerr << "\n";
   std::vector<uint32_t>().swap(uniqueHeadsPerm);
   std::vector<uint32_t>().swap(offSetPerm);
   std::vector<uint32_t>().swap(uniqueHeadsPosVec);


   uniqueHeadsPosPrefSum[0] = 0;
   for(size_t i = 1; i < _n; i++){
      uniqueHeadsPosPrefSum[i] = uniqueHeadsPosPrefSum[i-1] + uniqueHeadsPos[i-1];
   }
   delete[] uniqueHeadsPos;

   if(verbose) for(uint64_t i = 0; i < phrases.size(); i++){
      std::cerr << headsSA[i].pos << "," << headsSA[i].len << " --> " << stringHeads[i] << "\n";
   }
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Reranking took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to invert stringHeads\n";
   int64_t hLen = headsSA.size();

   int32_t *stringHeads2 = new int32_t[headsSA.size()+1];
   for(uint32_t i = 0; i < headsSA.size(); i++){
      stringHeads2[headsSA[i].start] = stringHeads[i];
   }
   for(uint32_t i = 1; i < headsSA.size()+1; i++){
      stringHeads[i] = stringHeads2[i];
   }
   delete[] stringHeads2;
   std::vector<MatchSA>().swap(headsSA);
   stringHeads[hLen] = 0;
   // //print stringHeads
   // for(uint64_t i = 0; i < phrases.size()+1; i++){
   //    std::cerr << stringHeads[i] << "\n";
   // }
   // exit(0);
   std::cerr << "Inverted stringHeads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Inverting stringHeads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   uint32_t bufferLibSais = 60000;
   int32_t *indicesSA = new int32_t[hLen+1+bufferLibSais]();
   // for(uint64_t i = 0; i < phrases.size(); i++){
   //    headsSufArr[i] = phrases[i].pos << 31 + phrases[i].len;
   // }
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute suffix array of MS-heads\n";
   libsais_int(stringHeads, indicesSA, hLen+1, rank+1, bufferLibSais);
   //sacak_int(stringHeads, indicesSA, hLen+1, rank+1);
   std::cerr << "Computed suffix array of MS-heads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing MS-heads SA took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   if(verbose) for(uint64_t i = 0; i < phrases.size()+1; i++){
      std::cerr << stringHeads[i] << ", " << indicesSA[i] << "\n";
   }
   delete[] stringHeads;
   int32_t *indicesSA2 = new int32_t[hLen+1];
   uint8_t *BWTheads = new uint8_t[hLen+1];
   for(uint32_t i = 0; i < hLen+1; i++){
      //std::cerr << "indicesSA[" << i << "] = " << indicesSA[i] << "\n";
      indicesSA2[indicesSA[i]] = i;
      BWTheads[i] = listOfChars[indicesSA[i+1]];
      //std::cerr << "BWTheads[" << i << "] = " << BWTheads[i] << "\n";
   }
   //exit(0);
   std::vector<uint8_t>().swap(listOfChars);
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to permute headsSA\n";

   //for(uint32_t i = 1; i < hLen+1; i++){
   //   std::cerr << phrases[indicesSA[i]].pos << '\n';
   //}

   const std::vector<Match>::iterator begPhrases = phrases.begin();
   //auto firstHead = std::begin(phrases);
   std::cerr << "New for loop\n";
   uint64_t currentSequence = 1;
   uint32_t *backupPos = new uint32_t[hLen+1];
   for(size_t i = 0; i < hLen; i++){
      //std::cerr << "i = " << i << "\n";
      //std::cerr << "phrases[i].start = " << phrases[i].start << " phrases[i].len = " << phrases[i].len << " \n";
      if(phrases[i].len == 0){
         backupPos[i] = phrases[i].pos;
         currentSequence++;
         continue;
      }
      uint32_t nextStart = phrases[i].start + phrases[i].len;
      const Match headNextStart = Match(nextStart, 0, currentSequence);
      std::vector<Match>::iterator nextHead = pHeads->predQuery(headNextStart, begPhrases);
      //std::cerr << "nextHead - begPhrases " << nextHead - begPhrases << "\n";
      //std::cerr << "indicesSA2[nextHead - begPhrases] - 1 = " << indicesSA2[nextHead - begPhrases] - 1 << "\n";
      //std::cerr << "_ISA[nextHead->pos + (nextStart - nextHead->start)] = " << _ISA[nextHead->pos + (nextStart - nextHead->start)] << "\n";
      backupPos[i] = phrases[i].pos;
      phrases[i].pos = _ISA[nextHead->pos + (nextStart - nextHead->start)];
      indicesSA2[i] = indicesSA2[nextHead - begPhrases] - 1;
   }
   std::cerr << "New for loop end\n";

   //print backupPos
   //for(uint32_t i = 0; i < hLen+1; i++){
   //  std::cerr << "backupPos[" << i << "] = " << backupPos[i] << "\n";
   //}
   headsSANN.resize(hLen);
   // std::ofstream outHeads("headsSA.txt", std::ios::out | std::ios::binary);
   // outHeads.seekp(hLen*sizeof(MatchSA));
   for(int64_t i = 1; i < hLen+1; i++) {
      //headsSA[i-1] = MatchSA(indicesSA2[indicesSA[i]], phrases[indicesSA[i]].pos, phrases[indicesSA[i]].len, phrases[indicesSA[i]].smaller, listOfChars[indicesSA[i]]);
      //std::cerr << "backupPos[indicesSA[" << i << "]] = " << backupPos[indicesSA[i]] << "\n";
      //std::cerr << "prefSumBucketPosHeads[backupPos[indicesSA[" << i << "]]] = " << prefSumBucketPosHeads[backupPos[indicesSA[i]]] << "\n";
      //auto tmp = MatchSA(indicesSA2[indicesSA[i]], phrases[indicesSA[i]].pos, phrases[indicesSA[i]].len, phrases[indicesSA[i]].smaller, (char)0);
      headsSANN[prefSumBucketPosHeads[backupPos[indicesSA[i]]]++] = MatchSANN(indicesSA2[indicesSA[i]], phrases[indicesSA[i]].pos, phrases[indicesSA[i]].len, phrases[indicesSA[i]].smaller); //tmp;
      phrases[indicesSA[i]].pos = backupPos[indicesSA[i]];
      backupPos[indicesSA[i]] = prefSumBucketPosHeads[backupPos[indicesSA[i]]] - 1;
      // outHeads.seekp((prefSumBucketPosHeads[backupPos[indicesSA[i]]] - 1)*sizeof(MatchSA));
      // outHeads.write((char*)&tmp, sizeof(MatchSA));
      //BWTheads[i-1] = listOfChars[indicesSA[i]];
      //std::cerr << headsSA[i-1].start << "," << headsSA[i-1].pos << "," << headsSA[i-1].len << "," << headsSA[i-1].smaller << "," << headsSA[i-1].next <<"\n";
   }
   //std::vector<uint8_t>().swap(listOfChars);
   delete[] prefSumBucketPosHeads;

   // save headsSA to file
   std::ofstream outHeads("headsSA.txt", std::ios::out | std::ios::binary);
   //outHeads.sync_with_stdio(false);
   for(uint64_t i = 1; i < headsSANN.size(); i++){
      outHeads.write((char*)&headsSANN[backupPos[i]], sizeof(MatchSANN));
      //out << headsSA[i].start << "," << headsSA[i].pos << "," << headsSA[i].len << "," << headsSA[i].smaller << "," << headsSA[i].next <<"\n";
   }
   outHeads.close();
   delete[] backupPos;

   std::cerr << "headBoundaries.size(): " << headBoundaries.size() << "\n";
   std::vector<uint32_t>().swap(headBoundaries);
   std::cerr << "Permuted headsSA\n";
   //delete[] backupPos;
   delete[] indicesSA;
   delete[] indicesSA2;
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Permuting MS-heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t headSortTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to sort heads: " << headSortTime << " milliseconds\n";

   delete pHeads;
   // Now:
   // headsSA[i] = (nextHeadRank, nextHeadISA[pos], length, smaller, mmchar);

   uint64_t *counterSmallerThanHead = new uint64_t[headsSANN.size()]();

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute counterSmallerThanHead\n";

   //unsigned char *BWT_collection = new unsigned char[_sn]();
   //BWT_collection[0] = sequenceSeparator;
   ndoc = 1;
   uint64_t idxCurrentHead = 1;
   i = 1; //ATTENTION: skip first char if first char is %
   MatchSANN currentHead;
   uint64_t currentHeadEnd = 0;
   uint64_t iCurrentDoc = 0;

   //compute number of buckets that are critical
   uint32_t nCriticalBuckets = 0;
   for(uint32_t i = 0; i < _n; i++){
      if(bucketsForExpandedBWT[i] < 0){
         nCriticalBuckets++;
      }
   }
   bool useBuffer = arg.buffer ? true : false;
   std::cerr << "ratio of critical buckets = " << (double)nCriticalBuckets/(double)_n << "\n";
   //if((double)nCriticalBuckets/(double)_n > 0.3){
   //    useBuffer = true;
   //}
   //useBuffer = true;
   bool smallBuffer = _n < 65536;
   //std::cerr << "smallBuffer = " << smallBuffer << "\n";
   uint32_t suffixesInBuffer = 0;
   uint64_t capOfBufferSuffixes;
   double RAMavailable = arg.buffer*1000000000;
   if(smallBuffer){
      //capOfBufferSuffixes = RAMavailable/(sizeof(std::pair<uint32_t, MatchSANN>));//100000000;
      capOfBufferSuffixes = RAMavailable/(sizeof(MatchSANN));
   }
   else{
      //capOfBufferSuffixes = RAMavailable/(sizeof(std::pair<std::pair<uint32_t, uint32_t>, MatchSANN>));//100000000;//_n/65536 * 100;
      capOfBufferSuffixes = RAMavailable/(sizeof(std::pair<uint32_t, MatchSANN>));
   }
   std::cerr << "capOfBufferSuffixes = " << capOfBufferSuffixes << "\n";

   //compute some statistics to help with buffering the bufferSuffixes. Get a rough estimate of the number of suffixes per bucketsForExpandedBWT
   // uint64_t nSuffixesInCriticalBuckets = 0;
   // uint64_t *nSuffixesPerCriticalBucket = new uint64_t[_n]();
   // for(uint32_t i = 0; i < _n; i++){
   //    if(bucketsForExpandedBWT[i] < 0){
   //       nSuffixesInCriticalBuckets += -bucketsForExpandedBWT[i];
   //    }
   // }
   // std::cerr << "nSuffixesInCriticalBuckets = " << nSuffixesInCriticalBuckets << "\n";
   // for(uint32_t i = 0; i < _n; i++){
   //    if(bucketsForExpandedBWT[i] < 0){
   //       nSuffixesPerCriticalBucket[i] = (uint32_t)ceil((double)(-bucketsForExpandedBWT[i])/(double)nSuffixesInCriticalBuckets * 100);
   //       std::cerr << "bucketsForExpandedBWT[" << i << "] = " << bucketsForExpandedBWT[i] << "\n";
   //       std::cerr << "nSuffixesPerCriticalBucket[" << i << "] = " << nSuffixesPerCriticalBucket[i] << "\n";
   //       std::cerr << "ratio: " << (double)(-bucketsForExpandedBWT[i])/(double)nSuffixesInCriticalBuckets << "\n";
   //       std::cerr << "ceil: " << ceil((double)(-bucketsForExpandedBWT[i])/(double)nSuffixesInCriticalBuckets * 100) << "\n";
   //       //std::cerr << "nSuffixesPerCriticalBucket[" << i << "] = " << nSuffixesPerCriticalBucket[i] << "\n";
   //    }
   // }
   // std::cerr << "nSuffixesPerCriticalBucket computed\n";

   const std::vector<uint32_t>::iterator startUniqueHeads = uniqueHeads.begin();
   const std::vector<MatchSANN>::iterator startHeadsSA = headsSANN.begin();
   std::ifstream inHeads("headsSA.txt", std::ios::binary | std::ios::in);
   //inHeads.sync_with_stdio(false);
   char bufferRead[sizeof(MatchSANN)];
   //useBuffer = false;
   if(useBuffer){
      std::cerr << "Using buffer\n";
      
      if(smallBuffer){
         //std::vector<std::vector<std::pair<uint32_t, MatchSANN>>> bufferSuffixes;
         std::vector<std::vector<MatchSANN>> bufferSuffixes;
         bufferSuffixes.resize(_n);
         for(uint32_t i = 0; i < _n; i++){
            bufferSuffixes[i].reserve(capOfBufferSuffixes/_n + 1);
         }
         std::cerr << "Reserved space for bufferSuffixes\n";
         for(uint64_t i = 1; i < phrases.size(); i++){
            if(phrases[i].len == 0){
               //BWT_collection[ndoc] = BWTheads[ndoc];
               idxCurrentHead++;
               inHeads.read(bufferRead, sizeof(MatchSANN));
               ndoc++;
               continue;
            }
            //currentHead = headsSA[backupPos[idxCurrentHead++]];
            inHeads.read(bufferRead, sizeof(MatchSANN));
            currentHead = *reinterpret_cast<MatchSANN*>(bufferRead);
            uint32_t posSuffix = phrases[i].pos + 1;
            uint32_t lenSuffix = phrases[i].len - 1;
            for(uint64_t iHead = 1; iHead < phrases[i+1].start - phrases[i].start; iHead++){
               if(bucketsForExpandedBWT[posSuffix] < 0) [[likely]]{
                  //bufferSuffixes[posSuffix].push_back({lenSuffix, currentHead});
                  bufferSuffixes[posSuffix].push_back(MatchSANN(currentHead.start, currentHead.pos, lenSuffix, currentHead.smaller));
                  suffixesInBuffer++;
                  if(suffixesInBuffer == capOfBufferSuffixes){ //process bufferSuffixes
                     //in each bucket count how many suffixes are smaller than the heads
                     for(uint32_t pos = 0; pos < _n - 1; pos++){
                        for(uint32_t counter = 0; counter < bufferSuffixes[pos].size(); counter++){
                           //check if the len of the suffix is dfferent than the len of the uniqueHead
                           // //print bufferSuffixes[pos][counter]
                           // std::cerr << "bufferSuffixes[" << pos << "][" << counter << "] = " << bufferSuffixes[pos][counter].start << " " << bufferSuffixes[pos][counter].pos << " " << bufferSuffixes[pos][counter].len << " " << bufferSuffixes[pos][counter].smaller << "\n";
                           uint32_t counterHeadLen = bufferSuffixes[pos][counter].smaller ? bufferSuffixes[pos][counter].len : MAXIMUM_UINT32 - bufferSuffixes[pos][counter].len;
                           // std::cerr << "counterHeadLen = " << counterHeadLen << "\n";
                           // //print uniqueHeads[uniqueHeadsPosPrefSum[pos]] to uniqueHeads[uniqueHeadsPosPrefSum[pos+1]-1]
                           // for(uint32_t i = uniqueHeadsPosPrefSum[pos]; i < uniqueHeadsPosPrefSum[pos+1]; i++){
                           //    std::cerr << "uniqueHeads[" << i << "] = " << uniqueHeads[i] << "\n";
                           // }
                           //timesRank++;
                           auto apprxPointer = std::upper_bound(startUniqueHeads + uniqueHeadsPosPrefSum[pos], startUniqueHeads + uniqueHeadsPosPrefSum[pos+1], counterHeadLen) - startUniqueHeads;
                           if(apprxPointer == uniqueHeadsPosPrefSum[pos]){
                              counterSmallerThanHead[prefSumBucketPosHeadsCopy[pos]]++;
                              //timesLength++;
                           }
                           else if(apprxPointer == uniqueHeadsPosPrefSum[pos+1] & counterHeadLen != uniqueHeads[apprxPointer - 1]){
                              //timesLength++;
                              continue;
                           }
                           else if(apprxPointer != uniqueHeadsPosPrefSum[pos+1] & counterHeadLen != uniqueHeads[apprxPointer - 1]){ //maybe with -1
                              counterSmallerThanHead[prefSumBucketPosHeadsCopy[pos] + offSet[apprxPointer]]++;
                              //timesLength++;
                           }
                           else{
                              //timesISA++;
                              auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[pos] + offSet[apprxPointer-1], 
                                                               startHeadsSA + (prefSumBucketPosHeadsCopy[pos]+offSet[apprxPointer])*(offSet[apprxPointer]>0) + 
                                                               (prefSumBucketPosHeadsCopy[pos+1])*(offSet[apprxPointer]==0), 
                                                               bufferSuffixes[pos][counter],
                                                   [](const MatchSANN &headA, const MatchSANN &headB){
                                                   //timesRank++;
                                                   // if(headA.len != headB.len){
                                                   //    //timesLength++;
                                                   //    return headA.smaller*(headA.len < headB.len) +
                                                   //          !headB.smaller*(headA.len > headB.len);
                                                   // }
                                                   // else{
                                                      return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                                   // }
                                                   });
                              //get index of pointer
                              uint64_t toIncrement = pointer - startHeadsSA;
                              //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
                              if(toIncrement != prefSumBucketPosHeadsCopy[pos+1]) counterSmallerThanHead[toIncrement]++;
                           }
                        }
                        bufferSuffixes[pos].clear();
                     }
                     suffixesInBuffer = 0;
                  }
               }
               posSuffix++;
               lenSuffix--;
            }
         }
         std::cerr << "suffixesInBuffer left = " << suffixesInBuffer << "\n";
         if(suffixesInBuffer != 0){
            //in each bucket count how many suffixes are smaller than the heads
            for(uint32_t pos = 0; pos < _n; pos++){
               for(uint32_t counter = 0; counter < bufferSuffixes[pos].size(); counter++){
                  //check if the len of the suffix is dfferent than the len of the uniqueHead
                  // //print bufferSuffixes[pos][counter]
                  // std::cerr << "bufferSuffixes[" << pos << "][" << counter << "] = " << bufferSuffixes[pos][counter].start << " " << bufferSuffixes[pos][counter].pos << " " << bufferSuffixes[pos][counter].len << " " << bufferSuffixes[pos][counter].smaller << "\n";
                  uint32_t counterHeadLen = bufferSuffixes[pos][counter].smaller ? bufferSuffixes[pos][counter].len : MAXIMUM_UINT32 - bufferSuffixes[pos][counter].len;
                  // std::cerr << "counterHeadLen = " << counterHeadLen << "\n";
                  // //print uniqueHeads[uniqueHeadsPosPrefSum[pos]] to uniqueHeads[uniqueHeadsPosPrefSum[pos+1]-1]
                  // for(uint32_t i = uniqueHeadsPosPrefSum[pos]; i < uniqueHeadsPosPrefSum[pos+1]; i++){
                  //    std::cerr << "uniqueHeads[" << i << "] = " << uniqueHeads[i] << "\n";
                  // }           
                  //timesRank++;       
                  auto apprxPointer = std::upper_bound(startUniqueHeads + uniqueHeadsPosPrefSum[pos], startUniqueHeads + uniqueHeadsPosPrefSum[pos+1], counterHeadLen) - startUniqueHeads;
                  if(apprxPointer == uniqueHeadsPosPrefSum[pos]){
                     // std::cerr << "apprxPointer == uniqueHeadsPosPrefSum[pos]\n";
                     counterSmallerThanHead[prefSumBucketPosHeadsCopy[pos]]++;
                     //timesLength++;
                  }
                  else if(apprxPointer == uniqueHeadsPosPrefSum[pos+1] & counterHeadLen != uniqueHeads[apprxPointer - 1]){
                     // std::cerr << "apprxPointer == uniqueHeadsPosPrefSum[pos+1] & counterHeadLen != uniqueHeads[apprxPointer - 1]\n";
                     //timesLength++;
                     continue;
                  }
                  else if(apprxPointer != uniqueHeadsPosPrefSum[pos+1] & counterHeadLen != uniqueHeads[apprxPointer - 1]){ //maybe with -1
                     // std::cerr << "counterHeadLen != uniqueHeads[apprxPointer - 1]\n";
                     //timesLength++;
                     counterSmallerThanHead[prefSumBucketPosHeadsCopy[pos] + offSet[apprxPointer]]++;
                  }
                  else{
                     // std::cerr << "else\n";
                     //timesISA++;
                     auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[pos] + offSet[apprxPointer-1], 
                                                      startHeadsSA + (prefSumBucketPosHeadsCopy[pos]+offSet[apprxPointer])*(offSet[apprxPointer]>0) + 
                                                      (prefSumBucketPosHeadsCopy[pos+1])*(offSet[apprxPointer]==0), 
                                                      bufferSuffixes[pos][counter],
                                       [](const MatchSANN &headA, const MatchSANN &headB){
                                          //timesRank++;
                                          // if(headA.len != headB.len){
                                          //    //timesLength++;
                                          //    return headA.smaller*(headA.len < headB.len) +
                                          //          !headB.smaller*(headA.len > headB.len);
                                          // }
                                          // else{
                                             return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                          // }
                                       });
                     //get index of pointer
                     uint64_t toIncrement = pointer - startHeadsSA;
                     //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
                     if(toIncrement != prefSumBucketPosHeadsCopy[pos+1]) counterSmallerThanHead[toIncrement]++;
                  }
               }
            }
         }
         //std::vector<std::vector<std::pair<uint32_t, MatchSANN>>>().swap(bufferSuffixes);
         std::vector<std::vector<MatchSANN>>().swap(bufferSuffixes);
         //std::cerr << "timesLength = " << timesLength << "\n";
         //std::cerr << "timesRank = " << timesRank << "\n";
         //std::cerr << "timesLength/timesRank = " << (double)timesLength/timesRank << "\n";         
      }
      else{
         //buffer for big references (it has fixed size of 65536 and it is indexes by the first 16 bits of the suffix). It stores the head and the correct position of the suffix
         //std::vector<std::vector<std::pair<std::pair<uint32_t, uint32_t>, MatchSANN>>> bufferSuffixes;
         std::vector<std::vector<std::pair<uint32_t,MatchSANN>>> bufferSuffixes;
         bufferSuffixes.resize(65536);
         for(uint32_t i = 0; i < 65536; i++){
            bufferSuffixes[i].reserve(capOfBufferSuffixes/65536 + 1);
         }
         uint32_t max = _n;
         //compute shift to get the first 16 bits of the suffix
         uint32_t shift = 0;
         for(uint8_t i = 0; i < 48; i++){
            if(max < nOfDigitsBig[i]){
               shift = i;
               break;
            }
         } 
         for(uint64_t i = 1; i < phrases.size(); i++){
            if(phrases[i].len == 0){
               //BWT_collection[ndoc] = BWTheads[ndoc];
               idxCurrentHead++;
               inHeads.read(bufferRead, sizeof(MatchSANN));
               ndoc++;
               continue;
            }
            //currentHead = headsSA[backupPos[idxCurrentHead++]];
            inHeads.read(bufferRead, sizeof(MatchSANN));
            currentHead = *reinterpret_cast<MatchSANN*>(bufferRead);
            //std::cerr << "currentHead.start = " << currentHead.start << " currentHead.pos = " << currentHead.pos << " currentHead.len = " << currentHead.len << " currentHead.smaller = " << currentHead.smaller << "\n";
            //std::cerr << "phrases[i].pos = " << phrases[i].pos << " phrases[i].len = " << phrases[i].len << '\n';
            //std::cerr << "phrases[i+1].start - phrases[i].start = " << phrases[i+1].start - phrases[i].start << '\n';
            uint32_t posSuffix = phrases[i].pos + 1;
            uint32_t lenSuffix = phrases[i].len - 1;
            for(uint64_t iHead = 1; iHead < phrases[i+1].start - phrases[i].start; iHead++){
               if(bucketsForExpandedBWT[posSuffix] < 0){
                  //bufferSuffixes[posSuffix >> shift].push_back({{lenSuffix, posSuffix}, currentHead});
                  bufferSuffixes[posSuffix >> shift].push_back({posSuffix, MatchSANN(currentHead.start, currentHead.pos, lenSuffix, currentHead.smaller)});
                  suffixesInBuffer++;
                  if(suffixesInBuffer == capOfBufferSuffixes){ //process bufferSuffixes
                     //sort by correct position in each 16-bit bucket
                     for(uint32_t area = 0; area < 65536; area++){
                        std::sort(bufferSuffixes[area].begin(), bufferSuffixes[area].end(), 
                           [](const std::pair<uint32_t, MatchSANN> &a, const std::pair<uint32_t, MatchSANN> &b){
                              return a.first < b.first;
                           }
                        );
                     }
                     //in each bucket count how many suffixes are smaller than the heads
                     for(uint32_t pos = 0; pos < 65536; pos++){
                        for(uint32_t counter = 0; counter < bufferSuffixes[pos].size(); counter++){
                           //uint32_t len = bufferSuffixes[pos][counter].first.first;
                           //uint32_t realPos = bufferSuffixes[pos][counter].first.second;
                           uint32_t realPos = bufferSuffixes[pos][counter].first;
                           uint32_t counterHeadLen = bufferSuffixes[pos][counter].second.smaller ? bufferSuffixes[pos][counter].second.len : MAXIMUM_UINT32 - bufferSuffixes[pos][counter].second.len;    
                           auto apprxPointer = std::upper_bound(startUniqueHeads + uniqueHeadsPosPrefSum[realPos], startUniqueHeads + uniqueHeadsPosPrefSum[realPos+1], counterHeadLen) - startUniqueHeads;
                           if(apprxPointer == uniqueHeadsPosPrefSum[realPos]){
                              counterSmallerThanHead[prefSumBucketPosHeadsCopy[realPos]]++;
                           }
                           else if((apprxPointer == uniqueHeadsPosPrefSum[realPos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){
                              continue;
                           }
                           else if((apprxPointer != uniqueHeadsPosPrefSum[realPos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){ //maybe with -1
                              counterSmallerThanHead[prefSumBucketPosHeadsCopy[realPos] + offSet[apprxPointer]]++;
                           }
                           else{
                              auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[realPos] + offSet[apprxPointer-1], 
                                                               startHeadsSA + (prefSumBucketPosHeadsCopy[realPos]+offSet[apprxPointer])*(offSet[apprxPointer]>0) + 
                                                               (prefSumBucketPosHeadsCopy[realPos+1])*(offSet[apprxPointer]==0), 
                                                               bufferSuffixes[pos][counter].second,
                                                   [](const MatchSANN &headA, const MatchSANN &headB){
                                                      return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                                   });
                              //get index of pointer
                              uint64_t toIncrement = pointer - startHeadsSA;
                              //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
                              if(toIncrement != prefSumBucketPosHeadsCopy[realPos+1]) counterSmallerThanHead[toIncrement]++;
                           }
                        }
                        bufferSuffixes[pos].clear();
                     }
                     suffixesInBuffer = 0;
                  }
                  //__builtin_prefetch(&headsSA[prefSumBucketPosHeadsCopy[posSuffix+1]], 0, 3);
               }
               posSuffix++;
               lenSuffix--;
            }
         }
         std::cerr << "suffixesInBuffer left = " << suffixesInBuffer << "\n";
         if(suffixesInBuffer != 0){
            for(uint32_t area = 0; area < 65536; area++){
               std::sort(bufferSuffixes[area].begin(), bufferSuffixes[area].end(), [](const std::pair<uint32_t, MatchSANN> &a, const std::pair<uint32_t, MatchSANN> &b){
                  return a.first < b.first;
               });
            }
            for(uint32_t pos = 0; pos < 65536; pos++){
               for(uint32_t counter = 0; counter < bufferSuffixes[pos].size(); counter++){
                  //uint32_t len = bufferSuffixes[pos][counter].first.first;
                  //uint32_t realPos = bufferSuffixes[pos][counter].first.second;
                  uint32_t realPos = bufferSuffixes[pos][counter].first;
                  uint32_t counterHeadLen = bufferSuffixes[pos][counter].second.smaller ? bufferSuffixes[pos][counter].second.len : MAXIMUM_UINT32 - bufferSuffixes[pos][counter].second.len;    
                  auto apprxPointer = std::upper_bound(startUniqueHeads + uniqueHeadsPosPrefSum[realPos], startUniqueHeads + uniqueHeadsPosPrefSum[realPos+1], counterHeadLen) - startUniqueHeads;
                  if(apprxPointer == uniqueHeadsPosPrefSum[realPos]){
                     counterSmallerThanHead[prefSumBucketPosHeadsCopy[realPos]]++;
                  }
                  else if((apprxPointer == uniqueHeadsPosPrefSum[realPos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){
                     continue;
                  }
                  else if((apprxPointer != uniqueHeadsPosPrefSum[realPos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){ //maybe with -1
                     counterSmallerThanHead[prefSumBucketPosHeadsCopy[realPos] + offSet[apprxPointer]]++;
                  }
                  else{
                     auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[realPos] + offSet[apprxPointer-1], 
                                                      startHeadsSA + (prefSumBucketPosHeadsCopy[realPos]+offSet[apprxPointer])*(offSet[apprxPointer]>0) + 
                                                      (prefSumBucketPosHeadsCopy[realPos+1])*(offSet[apprxPointer]==0), 
                                                      bufferSuffixes[pos][counter].second,
                                       [](const MatchSANN &headA, const MatchSANN &headB){
                                             return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                       });
                     //get index of pointer
                     uint64_t toIncrement = pointer - startHeadsSA;
                     //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
                     if(toIncrement != prefSumBucketPosHeadsCopy[realPos+1]) counterSmallerThanHead[toIncrement]++;
                  }
               }
            }
         }
         //std::vector<std::vector<std::pair<std::pair<uint32_t, uint32_t>, MatchSANN>>>().swap(bufferSuffixes);
         std::vector<std::vector<std::pair<uint32_t, MatchSANN>>>().swap(bufferSuffixes);
      }
   }
   else{
      //just count how many suffixes are smaller than the heads because the suffixes are few
      for(uint64_t i = 1; i < phrases.size(); i++){
         if(phrases[i].len == 0){
            //BWT_collection[ndoc] = BWTheads[ndoc];
            idxCurrentHead++;
            inHeads.read(bufferRead, sizeof(MatchSANN));
            ndoc++;
            continue;
         }
         
         //currentHead = headsSA[backupPos[idxCurrentHead++]];
         inHeads.read(bufferRead, sizeof(MatchSANN));
         currentHead = *reinterpret_cast<MatchSANN*>(bufferRead);
         //std::cerr << "currentHead.start = " << currentHead.start << " currentHead.pos = " << currentHead.pos << " currentHead.len = " << currentHead.len << " currentHead.smaller = " << currentHead.smaller << "\n";
         //std::cerr << "phrases[i].pos = " << phrases[i].pos << " phrases[i].len = " << phrases[i].len << '\n';
         //std::cerr << "phrases[i+1].start - phrases[i].start = " << phrases[i+1].start - phrases[i].start << '\n';
         uint32_t posSuffix = phrases[i].pos + 1;
         uint32_t lenSuffix = phrases[i].len - 1;
         for(uint64_t iHead = 1; iHead < phrases[i+1].start - phrases[i].start; iHead++){
            // std::cerr << "posSuffix = " << posSuffix << " lenSuffix = " << lenSuffix << '\n';
            //__builtin_prefetch(&headsSA[prefSumBucketHeads_SA[posSuffix].first], 0, 3);
            if(bucketsForExpandedBWT[posSuffix] < 0){
               auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[posSuffix], startHeadsSA + prefSumBucketPosHeadsCopy[posSuffix+1], currentHead,
                                 [lenSuffix](const MatchSANN &headA, const MatchSANN &headB){
                                       if(lenSuffix != headB.len){
                                          return headA.smaller*(lenSuffix < headB.len) +
                                                !headB.smaller*(lenSuffix > headB.len);
                                       }
                                       else{
                                          return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                       }
                                 });
               //get index of pointer
               uint64_t toIncrement = pointer - startHeadsSA;
               //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
               if(toIncrement != prefSumBucketPosHeadsCopy[posSuffix+1]) counterSmallerThanHead[toIncrement]++;
               //__builtin_prefetch(&headsSA[prefSumBucketPosHeadsCopy[posSuffix+1]], 0, 3);
            }
            posSuffix++;
            lenSuffix--;
         }
      }
   }

   inHeads.close();
   std::remove("headsSA.txt");
   std::vector<Match>().swap(phrases);
   std::vector<MatchSANN>().swap(headsSANN);
   std::vector<uint32_t>().swap(uniqueHeads);
   delete[] uniqueHeadsPosPrefSum;
   std::vector<uint32_t>().swap(offSet);
   delete[] prefSumBucketPosHeadsCopy;

   int64_t *bucketsForExpandedBWT_copy = new int64_t[_n]();
   for(uint32_t i = 0; i < _n; i++){
      bucketsForExpandedBWT_copy[i] = bucketsForExpandedBWT[i];
   }
   //invert data in bucketsForExpandedBWT w.r.t. _ISA
   for(uint32_t i = 0; i < _n; i++){
      bucketsForExpandedBWT[_ISA[i]] = bucketsForExpandedBWT_copy[i];
   }
   std::cerr << "Inverted data in bucketsForExpandedBWT\n";
   //print bucketsForExpandedBWT
   //for(uint64_t i = 0; i < _n; i++){
   //  std::cerr << "_SA["<<i<<"]: " << _SA[i] << ", " << bucketsForExpandedBWT[i] << "\n";
   //}
   delete[] bucketsForExpandedBWT_copy;

   std::cerr << "Computed counterSmallerThanHead\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing counterSmallerThanHead took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t1 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute final BWT\n";
   uint64_t pointerNextChar = D - 1;
   //print bucketsForExpandedBWT
   //std::cerr << "bucketsForExpandedBWT: \n";
   if(verbose) for(uint32_t i = 0; i < _n; i++){
      std::cerr << "i: " << i << " - " << (int64_t)bucketsForExpandedBWT[i] << "\n";
   }

   //print prefSumBucketISAPosHeadsCopy
   //std::cerr << "prefSumBucketISAPosHeadsCopy: ";
   if(verbose) for(uint32_t i = 0; i < _n; i++){
      std::cerr << "i: " << i << " - " << prefSumBucketISAPosHeadsCopy[i] << "\n";
   }

   //permute counterSmallerThanHead w.r.t. _SA
   // print counterSmallerThanHead
   if(verbose) for(uint64_t i = 0; i < hLen; i++){
      std::cerr << "counterSmallerThanHead[" << i << "]: " << counterSmallerThanHead[i] << "\n";
   }

   uint64_t *counterSmallerThanHead_copy = new uint64_t[hLen]();
   for(uint64_t i = 0; i < hLen; i++){
      counterSmallerThanHead_copy[i] = counterSmallerThanHead[i];
   }
   uint64_t a = 0;
   std::vector<std::pair<uint64_t, uint64_t>> prefSumBucketHeads_SA(_n);
   for(size_t i = 0; i < _n; i++){
      prefSumBucketHeads_SA[_SA[i]] = std::make_pair(prefSumBucketISAPosHeadsCopy[i], prefSumBucketISAPosHeadsCopy[i+1]);
   }
   for(uint64_t i = 0; i < _n; i++){
      //std::cerr << "prefSumBucketHeads_SA[" << i << "].first: " << prefSumBucketHeads_SA[i].first << "\n";
      //std::cerr << "prefSumBucketHeads_SA[" << i << "].second: " << prefSumBucketHeads_SA[i].second << "\n";
      for(uint64_t counter = prefSumBucketHeads_SA[i].first; counter < prefSumBucketHeads_SA[i].second; counter++){
         //std::cerr << "a: " << a << "\n";
         counterSmallerThanHead[counter] = counterSmallerThanHead_copy[a++];
      }
   }
   delete[] counterSmallerThanHead_copy;
   std::vector<std::pair<uint64_t, uint64_t>>().swap(prefSumBucketHeads_SA);

   //print counterSmallerThanHead
   if(verbose) for(uint64_t i = 0; i < hLen; i++){
      std::cerr << "counterSmallerThanHead[" << i << "]: " << counterSmallerThanHead[i] << "\n";
   }

   bool RLEWrite = arg.format;
   //uint32_t pointerHeads = D-1;
   if(!RLEWrite){
      std::ofstream streamOutfile(arg.outname+".bwt", std::ios::out | std::ios::binary);
      //streamOutfile.sync_with_stdio(false);
      std::vector<uint8_t> bufferWrite;
      bufferWrite.reserve(1024*1024);
      //FILE *f = fopen("BWT_collection", "wb");

      //std::cerr << "pointerNextChar: " << pointerNextChar << "\n";
      //write BWT_collection to file
      streamOutfile.write((char*)&BWTheads[0], sizeof(uint8_t)*(D-1));
      //for(uint32_t d = 0; d < D - 1; d++){
         //BWT_collection[d] = BWTheads[d];
      //   streamOutfile.write((char*)&BWTheads[d], sizeof(uint8_t));
         //fwrite(&BWTheads[d], sizeof(uint8_t), 1, f);
      //}
      //write BWT_collection to file
      for(uint32_t i = 1; i < _x.size(); i++){
         //std::cerr << "i: " << i << "\n";
         //uint8_t characterEqual = _SA[i] ? _x[_SA[i]-1] : '$';
         uint8_t characterEqual = _BWT[i];
         //std::cerr << "characterEqual: " << characterEqual << "\n";
         if(bucketsForExpandedBWT[i] >= 0){
            //std::cerr << "Safe bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            for(uint64_t counter = 0; counter < bucketsForExpandedBWT[i]; counter++){
               //BWT_collection[pointerNextChar++] = characterEqual;
               //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
               bufferWrite.push_back(characterEqual);
               if(bufferWrite.size() == bufferWrite.capacity()){
                  streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
                  bufferWrite.clear();
               }
               //streamOutfile << characterEqual;
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
            }
         }
         else{
            //std::cerr << "Unsafe bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            for(uint64_t counter = prefSumBucketISAPosHeadsCopy[i]; counter < prefSumBucketISAPosHeadsCopy[i+1]; counter++){
               //std::cerr << "counter: " << counter << "\n";
               //std::cerr << "counterSmallerThanHead[counter]: " << counterSmallerThanHead[counter] << "\n";
               for(uint64_t counter2 = 0; counter2 < counterSmallerThanHead[counter]; counter2++){
                  //std::cerr << "counter2: " << counter2 << "\n";
                  //BWT_collection[pointerNextChar++] = characterEqual;
                  //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
                  bufferWrite.push_back(characterEqual);
                  if(bufferWrite.size() == bufferWrite.capacity()){
                     streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
                     bufferWrite.clear();
                  }
                  //streamOutfile << characterEqual;
                  //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                  bucketsForExpandedBWT[i]++;
               }
               //std::cerr << "headsSA[counter].next: " << headsSA[counter].next << "\n";
               //BWT_collection[pointerNextChar++] = BWTheads[counter];//headsSA[counter].next;
               //streamOutfile.write((char*)&BWTheads[counter], sizeof(uint8_t));
               bufferWrite.push_back(BWTheads[counter]);
               if(bufferWrite.size() == bufferWrite.capacity()){
                  streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
                  bufferWrite.clear();
               }
               //streamOutfile.write((char*)&BWTheads[pointerHeads++], sizeof(uint8_t));
               //streamOutfile << BWTheads[counter];
               //fwrite(&BWTheads[counter], sizeof(uint8_t), 1, f);
               bucketsForExpandedBWT[i]++;
            }
            //std::cerr << "Remaining suff in bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            for(int64_t counter = bucketsForExpandedBWT[i]; counter < 0; counter++){
            //BWT_collection[pointerNextChar++] = characterEqual;
            //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
            bufferWrite.push_back(characterEqual);
            if(bufferWrite.size() == bufferWrite.capacity()){
               streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
               bufferWrite.clear();
            }
            //streamOutfile << characterEqual;
            //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
            }
         }
      }
      if(bufferWrite.size() > 0){
         streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
         bufferWrite.clear();
      }
      streamOutfile.close();
   }
   else{
      //write only when character is different
      uint8_t prevChar = (char)0;
      uint32_t runLength = 0;

      std::ofstream streamOutfile(arg.outname+".rl_bwt", std::ios::out | std::ios::binary);
      //streamOutfile.sync_with_stdio(false);
      //FILE *f = fopen("BWT_collection", "wb");

      //std::cerr << "pointerNextChar: " << pointerNextChar << "\n";
      //write BWT_collection to file
      //streamOutfile.write((char*)&BWTheads[0], sizeof(uint8_t)*(D-1));
      for(uint64_t d = 0; d < D-1; d++){
         std::cerr << "d: " << d << " ";
         std::cerr << "BWTheads[d]: " << BWTheads[d] << "\n";
         if(prevChar != BWTheads[d]){
            if(runLength > 0){
               streamOutfile.write((char*)&runLength, sizeof(uint32_t));
               streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
               //fwrite(&runLength, sizeof(uint32_t), 1, f);
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
               runLength = 1;
               prevChar = BWTheads[d];
            }
            else{
               runLength = 1;
               prevChar = BWTheads[d];
            }
         }
         else{
            runLength++;
         }
      }
      //write BWT_collection to file
      for(uint32_t i = 1; i < _x.size(); i++){
         //std::cerr << "i: " << i << "\n";
         //uint8_t characterEqual = _SA[i] ? _x[_SA[i]-1] : '$';
         uint8_t characterEqual = _BWT[i];
         //std::cerr << "characterEqual: " << characterEqual << "\n";
         if(bucketsForExpandedBWT[i] > 0){
            //std::cerr << "Safe bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            //for(uint64_t counter = 0; counter < bucketsForExpandedBWT[i]; counter++){
            //BWT_collection[pointerNextChar++] = characterEqual;
            if(prevChar != characterEqual){
               streamOutfile.write((char*)&runLength, sizeof(uint32_t));
               streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
               //fwrite(&runLength, sizeof(uint32_t), 1, f);
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
               runLength = bucketsForExpandedBWT[i];
               prevChar = characterEqual;
            }
            else{
               runLength+=bucketsForExpandedBWT[i];
            }
            //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
            //}
         }
         else if(bucketsForExpandedBWT[i] < 0){
            //std::cerr << "Unsafe bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            for(uint64_t counter = prefSumBucketISAPosHeadsCopy[i]; counter < prefSumBucketISAPosHeadsCopy[i+1]; counter++){
               //std::cerr << "counter: " << counter << "\n";
               //std::cerr << "counterSmallerThanHead[counter]: " << counterSmallerThanHead[counter] << "\n";
               //for(uint64_t counter2 = 0; counter2 < counterSmallerThanHead[counter]; counter2++){
                  //std::cerr << "counter2: " << counter2 << "\n";
                  //BWT_collection[pointerNextChar++] = characterEqual;
               if(counterSmallerThanHead[counter]){
                  if(prevChar != characterEqual){
                     streamOutfile.write((char*)&runLength, sizeof(uint32_t));
                     streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
                     //fwrite(&runLength, sizeof(uint32_t), 1, f);
                     //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                     runLength = counterSmallerThanHead[counter];
                     prevChar = characterEqual;
                  }
                  else{
                     runLength+=counterSmallerThanHead[counter];
                  }
                  //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
                  //streamOutfile << characterEqual;
                  //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                  bucketsForExpandedBWT[i]+=counterSmallerThanHead[counter];
               }
               //}
               //std::cerr << "headsSA[counter].next: " << headsSA[counter].next << "\n";
               //BWT_collection[pointerNextChar++] = BWTheads[counter];//headsSA[counter].next;
               if(BWTheads[counter] != prevChar){
                  streamOutfile.write((char*)&runLength, sizeof(uint32_t));
                  streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
                  //fwrite(&runLength, sizeof(uint32_t), 1, f);
                  //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                  runLength = 1;
                  prevChar = BWTheads[counter];
               }
               else{
                  runLength++;
               }
               //streamOutfile.write((char*)&BWTheads[counter], sizeof(uint8_t));
               //streamOutfile.write((char*)&BWTheads[pointerHeads++], sizeof(uint8_t));
               //streamOutfile << BWTheads[counter];
               //fwrite(&BWTheads[counter], sizeof(uint8_t), 1, f);
               bucketsForExpandedBWT[i]++;
            }
            //std::cerr << "Remaining suff in bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            //for(int64_t counter = bucketsForExpandedBWT[i]; counter < 0; counter++){
            //BWT_collection[pointerNextChar++] = characterEqual;
            if(bucketsForExpandedBWT[i]!=0){
               if(prevChar != characterEqual){
                  streamOutfile.write((char*)&runLength, sizeof(uint32_t));
                  streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
                  //fwrite(&runLength, sizeof(uint32_t), 1, f);
                  //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                  runLength = -bucketsForExpandedBWT[i];
                  prevChar = characterEqual;
               }
               else{
                  runLength+=-bucketsForExpandedBWT[i];
               }
               //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
               //streamOutfile << characterEqual;
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
            }
         }
      }
      streamOutfile.write((char*)&runLength, sizeof(uint32_t));
      streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
      streamOutfile.close();
   }

   delete[] bucketsForExpandedBWT;
   delete[] counterSmallerThanHead;
   
   //fclose(f);
   std::cerr << "Computed final BWT\n";
   t2 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing final BWT took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms\n";
   delete[] BWTheads;
   //print BWT_collection
   //std::cerr << "BWT_collection: ";
   //for(uint32_t i = 0; i < _sx.size(); i++){
   //  std::cerr << BWT_collection[i];
   //}
   //std::cerr << "\n";

   delete[] _SA;
   delete[] _BWT;
   delete[] _ISA;
   delete[] prefSumBucketISAPosHeadsCopy;
   std::string().swap(_x);
   //std::string().swap(_sx);
   //std::vector<Match>().swap(phrases);

   if(arg.check){
      std::string _sx_just_to_check;
      _sx_just_to_check.reserve(_sn);

      std::ifstream streamInfile(collFileName,std::ios::in);
      std::string line, content;
      content.reserve(_n*2);
      uint64_t charactersRead = 0;
      //std::string _sx_just_to_check;
      //_sx_just_to_check.reserve(_sn);
      while(std::getline(streamInfile, line).good()){
         //std::cerr << "ndoc: " << ndoc << "\n";
         if( line.empty() || line[0] == '>' ){
            //std::cerr << "content: " << content << "\n";
            content += sequenceSeparator;
            charactersRead++;
            _sx_just_to_check += content;
            //std::cerr << "Finished parsing document " << ndoc << "\n";
            content.erase();
            //std::string().swap(content);
         }
         else if (!line.empty()) {
            charactersRead += line.size();
            if(charactersRead >= _sn - 1){ //if string is filled up to sn (useful for prefixLength)
               //std::cerr << "read too much\n";
               //std::cerr << "charactersRead: " << charactersRead << "\n";
               //std::cerr << "line.size(): " << line.size() << "\n";
               //std::cerr << "_sn - charactersRead + line.size() - 1 = " << line.size() - (charactersRead - _sn) - 1 << "\n";
               //std::cerr << "content.size(): " << content.size() << "\n";
               content += line.substr(0, line.size() - (charactersRead - _sn) - 1);
               //std::cerr << "content.size(): " << content.size() << "\n";
               //std::cerr << _sx_just_to_check.size() << "\n";
               break;
            }
            else{
               content += line;
            }
         }
      }
      streamInfile.close();
      if(content.size() != 0){
         content += sequenceSeparator;
         charactersRead++;
         if(charactersRead < _sn - 1){
            _sn = charactersRead;
         }
         _sx_just_to_check += content;
      }

      uint8_t *BWT_from_SA = new uint8_t[_sn];
      int32_t *arrayForBWT = new int32_t[_sn];
      t01 = std::chrono::high_resolution_clock::now();
      libsais_bwt(reinterpret_cast<unsigned char*>(const_cast<char*>(_sx_just_to_check.c_str())), BWT_from_SA, arrayForBWT, _sn, 0, NULL);
      t02 = std::chrono::high_resolution_clock::now();
      std::cerr << "Computing BWT_from_SA took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

      //read BWT_collection from file
      std::ifstream streamInfile2(arg.outname+".bwt",std::ios::in);
      streamInfile2.seekg(0, streamInfile2.end);
      uint64_t length = streamInfile2.tellg();
      streamInfile2.seekg(0, streamInfile2.beg);
      uint8_t *BWT_collection = new uint8_t[length];
      streamInfile2.read((char*)BWT_collection, length);
      streamInfile2.close();
      //print BWT_collection
      // std::cerr << "BWT_collection: ";
      // for(uint32_t i = 0; i < length; i++){
      //    std::cerr << BWT_collection[i];
      // }
      // std::cerr << "\n";
      std::cerr << "Read BWT_collection from file\n";
      //check this BWTCollection with BWT_from_SA
      std::cerr << "Checking BWT_collection with BWT_from_SA\n";
      for(uint32_t i = D-1; i < _sn; i++){
         if(BWT_collection[i] != BWT_from_SA[i]){
            std::cerr << "BWT_collection[" << i << "] = " << BWT_collection[i] << " != " << BWT_from_SA[i] << " = BWT_from_SA[" << i << "]\n";
            //std::cerr << "BWT_collection[" << i << "] = " << (int)BWT_collection[i] << " != " << (int)BWT_from_SA[i] << " = BWT_from_SA[" << i << "]\n";

            //break;
         }
         //else{
         //   std::cerr << "BWT_collection[" << i << "] = " << BWT_collection[i] << " == " << BWT_from_SA[i] << " = BWT_from_SA[" << i << "]\n";
         //}
      }
   }

   return numfactors;
}


int lzFactorizeMemorySaving(Args arg, char *collFileName) {
//int lzFactorize(char *fileToParse, int seqno, char* outputfilename, const bool v) {
   verbose = 0;
   //omp_set_num_threads(nThreads);
   std::cerr << "File is in memory\n";
   //std::cerr << "_x: " <<  _x << '\n';
   //std::cerr << "_sx: " <<  _sx << '\n';
   std::cerr << "_x.len: " <<  _n << '\n';
   std::cerr << "_sx.len: " <<  _sn << '\n';
   auto t1 = std::chrono::high_resolution_clock::now();

   unsigned int numfactors = 0;

   unsigned int inc = 100000000;
   uint64_t mark = inc;

   std::cerr << "About to start main parsing loop...\n";

   
   uint32_t ndoc = 1;
   int64_t *bucketsForExpandedBWT = new int64_t[_n]();
   bucketsForExpandedBWT[_n - 1] = 0;
   uint64_t *bucketISAPosHeads = new uint64_t[_n]();

   headBoundaries.push_back(0);
   uint64_t i = 0;

   uint64_t maxValue = 0;
   std::vector<Match> phrases;
   phrases.reserve(_sn / sizeof(Match));
   //uint64_t sumLen = 0, maxLen = 0;
   uint64_t heurYes = 0, heurNo = 0;
   std::vector<uint8_t> listOfChars;
   listOfChars.reserve(_sn / sizeof(Match));
   //listOfChars.push_back('%');
   uint64_t fasterAdjustInterval = 0, slowerAdjustInterval = 0;
   uint32_t backup = 0;

   
   std::ifstream streamInfile(collFileName,std::ios::in);
   //streamInfile.sync_with_stdio(false);
   std::string line, content;
   content.reserve(_n*2);
   uint64_t charactersRead = 0;
   //std::string _sx_just_to_check;
   //_sx_just_to_check.reserve(_sn);
   std::ofstream outPhrases("phrases.txt", std::ios::out | std::ios::binary);
   //outPhrases.sync_with_stdio(false);
   while(std::getline(streamInfile, line).good()){
      //std::cerr << "ndoc: " << ndoc << "\n";
      if( line.empty() || line[0] == '>' ){
         //std::cerr << "content: " << content << "\n";
         content += sequenceSeparator;
         charactersRead++;
         //_sx_just_to_check += content;
         int64_t i = 0;
         int32_t leftB = 0;
         int32_t rightB = _n-1;
         bool isSmallerThanMatch;
         unsigned char mismatchingSymbol;
         uint64_t prevPos = -2;
         uint32_t pos = _n - 1, len = 0;
         uint32_t iCurrentDoc = 0;
         D++;
         //std::cerr << "content.size(): " << content.size() << "\n";
         //std::cerr << "((int64_t)content.size())-1: " << ((int64_t)content.size())-1 << "\n";
         //std::cerr << "i < ((int64_t)content.size())-1: " << (i < ((int64_t)content.size())-1) << "\n";
         
         // if(content.size() >= 5){
         //    //std::cerr << "content[0]: " << content[0] << " content[1]: " << content[1] << " content[2]: " << content[2] << " content[3]: " << content[3] << "\n";
         //    std::pair<int32_t, int32_t> pair = lookupTableSA[convertNucleotideQuadrupletHash(content[0],content[1],content[2],content[3])];
         //    //std::cerr << "pair.first: " << pair.first << " pair.second: " << pair.second << "\n";
         //    if(pair.first != -1){
         //       leftB = pair.first;
         //       rightB = pair.second;
         //       len += 4;
         //    }
         // }
         while(i < ((int64_t)content.size())-1){
            // if(len < 4){
            //    //std::cerr << "WOAH! len is less than 4!\n";
            //    std::pair<int32_t, int32_t> pair = lookupTableSA[convertNucleotideQuadrupletHash(content[i],content[i+1],content[i+2],content[i+3])];
            //    if(pair.first != -1){
            //       leftB = pair.first;
            //       rightB = pair.second;
            //       len = 4;
            //    }
            // }
            computeLZFactorAt(content, i, &pos, &len, leftB, rightB, isSmallerThanMatch, mismatchingSymbol);//, &backup);
            if(pos != prevPos+1){
               phrases.push_back(Match(iCurrentDoc, pos, len, isSmallerThanMatch, mismatchingSymbol));
               outPhrases.write((char*)&phrases.back(), sizeof(Match));
               if(i == 0) listOfChars.push_back(sequenceSeparator);
               else listOfChars.push_back(content[i-1]);
               if(bucketsForExpandedBWT[pos] > 0) bucketsForExpandedBWT[pos] = -bucketsForExpandedBWT[pos]-1;
               else bucketsForExpandedBWT[pos]--;
            }
            else{
               if(bucketsForExpandedBWT[pos] >= 0) bucketsForExpandedBWT[pos]++;
               else bucketsForExpandedBWT[pos]--;
            }
            iCurrentDoc++;
            len--;

            if(leftB == rightB){
               while(len > _PLCP[pos+1]){
                  i++;
                  iCurrentDoc++;
                  len--;
                  pos++;
                  if(bucketsForExpandedBWT[pos] >= 0) bucketsForExpandedBWT[pos]++;
                  else bucketsForExpandedBWT[pos]--;
                  //heurYes++;
               }
               leftB = rightB = _ISA[pos];
               std::pair<int,int> interval = adjustInterval(_ISA[pos+1], _ISA[pos+1], len);//fasterContractLeft(pos,rightB,len);
               //fasterAdjustInterval++;
               leftB = interval.first;
               rightB = interval.second;
            }
            else{
               std::pair<int,int> interval = contractLeft(leftB, rightB, len);//fasterContractLeft(pos,rightB,len);
               //slowerAdjustInterval++;
               leftB = interval.first;
               rightB = interval.second;
            }
               
            //heurNo++;
            i++;
            prevPos = pos;
         }
         len = 0;
         pos = _n - 1;
         phrases.push_back(Match(iCurrentDoc, pos, len, 0, sequenceSeparator));
         outPhrases.write((char*)&phrases.back(), sizeof(Match));
         bucketsForExpandedBWT[pos]--;
         headBoundaries.push_back(phrases.size());
         if(maxValue < iCurrentDoc) maxValue = iCurrentDoc;
         iCurrentDoc = 0;
         ndoc++;
         if(i == 0) {
            listOfChars.push_back(sequenceSeparator);
         }
         else{
            listOfChars.push_back(content[content.size()-2]);
         }
         //std::cerr << "Finished parsing document " << ndoc << "\n";
         content.erase();
         //std::string().swap(content);
      }
      else if (!line.empty()) {
         charactersRead += line.size();
         if(charactersRead >= _sn - 1){ //if string is filled up to sn (useful for prefixLength)
            //std::cerr << "read too much\n";
            //std::cerr << "charactersRead: " << charactersRead << "\n";
            //std::cerr << "line.size(): " << line.size() << "\n";
            //std::cerr << "_sn - charactersRead + line.size() - 1 = " << line.size() - (charactersRead - _sn) - 1 << "\n";
            //std::cerr << "content.size(): " << content.size() << "\n";
            content += line.substr(0, line.size() - (charactersRead - _sn) - 1);
            //std::cerr << "content.size(): " << content.size() << "\n";
            //std::cerr << _sx_just_to_check.size() << "\n";
            break;
         }
         else{
            content += line;
         }
      }
   }
   streamInfile.close();
   if(content.size() != 0){
      content += sequenceSeparator;
      charactersRead++;
      if(charactersRead < _sn - 1){
         _sn = charactersRead;
      }
      //_sx_just_to_check += content;
      //std::cerr << "Finishing up with " << content.size() << " characters\n";
      //std::cerr << "content: " << content << "\n";
      D++;
      int64_t i = 0;
      int32_t leftB = 0;
      int32_t rightB = _n-1;
      bool isSmallerThanMatch;
      unsigned char mismatchingSymbol;
      uint64_t prevPos = -2;
      uint32_t pos = _n - 1, len = 0;
      uint32_t iCurrentDoc = 0;
      while(i < content.size()){
         if(content[i] == sequenceSeparator){
            leftB = 0;
            rightB = _n-1;
            len = 0;
            pos = _n - 1;
            phrases.push_back(Match(iCurrentDoc, pos, len, 0, sequenceSeparator));
            outPhrases.write((char*)&phrases.back(), sizeof(Match));
            bucketsForExpandedBWT[pos]--;
            headBoundaries.push_back(phrases.size());
            if(maxValue < iCurrentDoc) maxValue = iCurrentDoc;
            iCurrentDoc = 0;
            ndoc++;
            if(i == 0) {
               listOfChars.push_back(sequenceSeparator);
            }
            else{
               listOfChars.push_back(content[i - 1]);
            }
         }
         else{
            computeLZFactorAt(content, i, &pos, &len, leftB, rightB, isSmallerThanMatch, mismatchingSymbol);//, &backup);
            if(pos != prevPos+1){
               phrases.push_back(Match(iCurrentDoc, pos, len, isSmallerThanMatch, mismatchingSymbol));
               outPhrases.write((char*)&phrases.back(), sizeof(Match));
               if(i == 0) listOfChars.push_back(sequenceSeparator);
               else listOfChars.push_back(content[i-1]);
               if(bucketsForExpandedBWT[pos] > 0) bucketsForExpandedBWT[pos] = -bucketsForExpandedBWT[pos]-1;
               else bucketsForExpandedBWT[pos]--;
            }
            else{
               if(bucketsForExpandedBWT[pos] >= 0) bucketsForExpandedBWT[pos]++;
               else bucketsForExpandedBWT[pos]--;
            }
            iCurrentDoc++;
            len--;

            if(leftB == rightB){
               while(len > _PLCP[pos+1]){
                  i++;
                  iCurrentDoc++;
                  len--;
                  pos++;
                  if(bucketsForExpandedBWT[pos] >= 0) bucketsForExpandedBWT[pos]++;
                  else bucketsForExpandedBWT[pos]--;
                  //heurYes++;
               }
               leftB = rightB = _ISA[pos];
               std::pair<int,int> interval = adjustInterval(_ISA[pos+1], _ISA[pos+1], len);//fasterContractLeft(pos,rightB,len);
               //slowerAdjustInterval++;
               leftB = interval.first;
               rightB = interval.second;
            }
            else{
               std::pair<int,int> interval = contractLeft(leftB, rightB, len);//fasterContractLeft(pos,rightB,len);
               //slowerAdjustInterval++;
               leftB = interval.first;
               rightB = interval.second;
            }
         }
         //heurNo++;
         i++;
         prevPos = pos;
      }
      std::string().swap(content);
   }
   outPhrases.close();
   std::cerr << "Finished parsing file of length: " << _sn << "\n";
   //std::cerr << "_sx_just_to_check = " << _sx_just_to_check << "\n";
   //std::cerr << "_sx_just_to_check.size(): " << _sx_just_to_check.size() << "\n";
   // std::cerr << "fasterAdjustInterval: " << fasterAdjustInterval << "\n";
   // std::cerr << "slowerAdjustInterval: " << slowerAdjustInterval << "\n";
   // std::cerr << "heurYes: " << heurYes << "\n";
   delete[] _LCP;
   delete[] _PLCP;
   delete _rmq;
   phrases.shrink_to_fit();
   listOfChars.shrink_to_fit();
   uint64_t sizePhrases = phrases.size();
   //std::string().swap(_x);
   std::cerr << "headBoundaries.size(): " << headBoundaries.size() << ", ndoc: " << ndoc << "\n";
   //for(uint64_t i = 0; i < headBoundaries.size(); i++) std::cerr << "i: " << i << " - " << headBoundaries[i] << "\n";
   //for(uint64_t i = 0; i < docBoundaries.size(); i++) std::cerr << "i: " << i << " - " << docBoundaries[i] << "\n";
   std::cerr << "docBoundaries.size(): " << D << "\n";
   std::cerr << "phrases.size() = " << sizePhrases << "\n";
   if(verbose) for(uint64_t i = 0; i < sizePhrases; i++) std::cerr << "i:" << i  << " - " << phrases[i].start << "," << phrases[i].pos << "," << phrases[i].len << "," << phrases[i].smaller << ", " << phrases[i].next << '\n';
   
   uint64_t *bucketPosHeads = new uint64_t[_n]();
   for(uint64_t i = 0; i < sizePhrases; i++){
      bucketPosHeads[phrases[i].pos]++;
      bucketISAPosHeads[_ISA[phrases[i].pos]]++;
   }
   std::vector<Match>().swap(phrases);
   // docBoundaries.push_back(_sn);
   if(verbose) std::cerr << "Printing docBoundaries" << "\n";
   //if(verbose) for(size_t i = 0; i < docBoundaries.size(); i++){ std::cerr << docBoundaries[i] << ", letter: " << _sx[docBoundaries[i]] << "\n";}
   auto t2 = std::chrono::high_resolution_clock::now();
   uint64_t lzFactorizeTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to compute matching statistics: " << lzFactorizeTime << " milliseconds\n";

   t1 = std::chrono::high_resolution_clock::now();
   //Sort suffixes corrisponding to heads
   uint64_t *prefSumBucketPosHeads = new uint64_t[_n]();
   uint64_t *prefSumBucketPosHeadsCopy = new uint64_t[_n+1];
   uint32_t *prefSumBucketISAPosHeads = new uint32_t[_n];
   uint32_t *prefSumBucketISAPosHeadsCopy = new uint32_t[_n+1];
   prefSumBucketISAPosHeads[0] = 0;
   prefSumBucketISAPosHeadsCopy[0] = 0;
   for(size_t i = 1; i < _n; i++){
      prefSumBucketPosHeads[i] = prefSumBucketPosHeads[i-1] + bucketPosHeads[i-1];
      //std::cerr << "prefSumBucketPosHeads[" << i << "] = " << prefSumBucketPosHeads[i] << "\n";
      prefSumBucketPosHeadsCopy[i] = prefSumBucketPosHeads[i];
      prefSumBucketISAPosHeads[i] = prefSumBucketISAPosHeads[i-1] + bucketISAPosHeads[i-1];
      prefSumBucketISAPosHeadsCopy[i] = prefSumBucketISAPosHeads[i];
   }
   delete[] bucketPosHeads;
   delete[] bucketISAPosHeads;
   prefSumBucketPosHeadsCopy[_n] = sizePhrases;
   prefSumBucketISAPosHeadsCopy[_n] = sizePhrases;

   std::vector<MatchSA> headsSA;
   headsSA.resize(sizePhrases);
   //headsSA.reserve(sizePhrases);
   //i = 0, ndoc = 0;
   auto t01 = std::chrono::high_resolution_clock::now();
   //for(std::vector<Match>::iterator it = phrases.begin(); it < phrases.end(); it++){
   std::ifstream inPhrases("phrases.txt", std::ios::binary);
   //inPhrases.sync_with_stdio(false);
   Match itHead;
   char bufferPhrase[sizeof(Match)];
   for(uint64_t i = 0; i < sizePhrases; i++){
      inPhrases.read(bufferPhrase, sizeof(Match));
      itHead = *reinterpret_cast<Match*>(bufferPhrase);
      // if(itHead.start != phrases[i].start || itHead.pos != phrases[i].pos || itHead.len != phrases[i].len || itHead.smaller != phrases[i].smaller || itHead.next != phrases[i].next){
      //    std::cerr << "i: " << i << " - " << itHead.start << "," << itHead.pos << "," << itHead.len << "," << itHead.smaller << ", " << itHead.next << '\n';
      //    std::cerr << "i: " << i << " - " << phrases[i].start << "," << phrases[i].pos << "," << phrases[i].len << "," << phrases[i].smaller << ", " << phrases[i].next << '\n';
      //    exit(1);
      // }
      // std::cerr << "i: " << i << " - " << itHead.start << "," << itHead.pos << "," << itHead.len << "," << itHead.smaller << ", " << itHead.next << '\n';
      // if(it->start == 0){
      //    ndoc++;
      // }
      //headsSA[prefSumBucketISAPosHeads[_ISA[it->pos]]++] = MatchSA(i++, it->pos, it->len, !it->smaller, it->next);
      headsSA[prefSumBucketISAPosHeads[_ISA[itHead.pos]]++] = MatchSA(i, itHead.pos, itHead.len, !itHead.smaller, itHead.next);
   }
   auto t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Bucketing heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   if(verbose) std::cerr << "Outputting headsSA bucketed (size=" << headsSA.size() << ")\n";
   std::vector<MatchSA>::iterator begHeads = headsSA.begin();
   t01 = std::chrono::high_resolution_clock::now();

   for(size_t i = 1; i < _n + 1; i++){
      std::sort(begHeads + prefSumBucketISAPosHeadsCopy[i-1], begHeads + prefSumBucketISAPosHeadsCopy[i], compareMatchSA);
   }
   delete[] prefSumBucketISAPosHeads;
   //delete[] prefSumBucketISAPosHeadsCopy;

   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Head sort took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to re-rank heads\n";
   int32_t *stringHeads = new int32_t[headsSA.size()+1];
   uint64_t rank = 1;
   int64_t prevLen = -1, prevSmall = -1, prevNext = -1;
   uint64_t prevPos = -1;
   i = 0;

   uint32_t *uniqueHeadsPos = new uint32_t[_n]();
   std::vector<uint32_t> uniqueHeadsPosVec;
   std::vector<uint32_t> offSet;
   std::vector<uint32_t> uniqueHeads;
   uniqueHeadsPosVec.reserve(headsSA.size());
   //uniqueHeadsPosVec.push_back(headsSA[0].pos);
   offSet.reserve(headsSA.size());
   //offSet.push_back(0);
   uniqueHeads.reserve(headsSA.size());
   //uniqueHeads.push_back(headsSA[0].len);
   //uniqueHeadsPos[headsSA[0].pos] = 1;
   uint32_t equalHeads = 0;
   for(std::vector<MatchSA>::iterator it = headsSA.begin(); it < headsSA.end(); it++){
      if(i >= headBoundaries.size()-1){
         if(it->pos != prevPos){
            equalHeads = 0;
         }
         if(it->pos != prevPos || it->len != prevLen || it->smaller != prevSmall){
            offSet.push_back(equalHeads);
            uniqueHeadsPos[it->pos]++;
            uniqueHeadsPosVec.push_back(it->pos);
            uniqueHeads.push_back(it->smaller ? MAXIMUM_UINT32 - it->len : it->len);
            //std::cerr << "Unique head: " << it->pos << " " << it->len << " " << it->smaller << " " << it->next << ". Storing: " << uniqueHeads.back() <<"\n";
            //std::cerr << "OffSet[" << offSet.end()-offSet.begin()-1 << "] = " << offSet.back() << "\n";
         }
         equalHeads++;
         if(it->pos != prevPos || it->len != prevLen || it->smaller != prevSmall || it->next != prevNext){
            rank++;
            prevPos = it->pos;
            prevLen = it->len;
            prevSmall = it->smaller;
            prevNext = it->next;
         }
      }
      stringHeads[i++] = rank;
   }
   std::cerr << "offSet.back() = " << offSet.back() << " equalHeads = " << equalHeads << "\n";
   if(offSet.back() != equalHeads-1){
      offSet.push_back(equalHeads-1);
      //uniqueHeadsPosVec.push_back(headsSA.back().pos);
      //uniqueHeadsPos[headsSA.back().pos]++;
      //uniqueHeads.push_back(headsSA.back().smaller ? MAXIMUM_UINT32 - headsSA.back().len : headsSA.back().len);
      //std::cerr << "Unique head: " << headsSA.back().pos << " " << headsSA.back().len << " " << headsSA.back().smaller << " " << headsSA.back().next << ". Storing: " << uniqueHeads.back() <<"\n";
      //std::cerr << "OffSet[" << offSet.end()-offSet.begin()-1 << "] = " << offSet.back() << "\n";
   }
   uniqueHeadsPosVec.shrink_to_fit();
   offSet.shrink_to_fit();
   uniqueHeads.shrink_to_fit();
   std::cerr << "Unique heads: " << uniqueHeads.size() << "\n";
   std::cerr << "Unique heads positions: " << uniqueHeadsPosVec.size() << "\n";
   std::cerr << "OffSet size: " << offSet.size() << "\n";
   std::cerr << "Re-ranked heads\n";
   std::cerr << "Maximum rank: " << rank << "\n";

   uint32_t *uniqueHeadsPosPrefSum = new uint32_t[_n];
   uniqueHeadsPosPrefSum[0] = 0;
   for(size_t i = 1; i < _n; i++){
      uniqueHeadsPosPrefSum[i] = uniqueHeadsPosPrefSum[i-1] + uniqueHeadsPos[i-1];
   }

   //permute uniqueHeads and offSet according to uniqueHeadsPosPrefSum using uniqueHeadsPosVec
   std::vector<uint32_t> uniqueHeadsPerm;
   std::vector<uint32_t> offSetPerm;
   uniqueHeadsPerm.resize(uniqueHeads.size());
   offSetPerm.resize(offSet.size());
   for(size_t i = 0; i < uniqueHeadsPosVec.size(); i++){
      uniqueHeadsPerm[uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]] = uniqueHeads[i];
      offSetPerm[uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]] = offSet[i];
      //uniqueHeadsPerm.push_back(uniqueHeads[uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]]);
      //offSetPerm.push_back(offSet[uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]]);
      uniqueHeadsPosPrefSum[uniqueHeadsPosVec[i]]++;
   }
   uniqueHeads = uniqueHeadsPerm;
   offSet = offSetPerm;
   // //print uniqueHeads and offSet
   // std::cerr << "Unique heads:\n";
   // for(size_t i = 0; i < uniqueHeads.size(); i++){
   //    std::cerr << uniqueHeads[i] << " ";
   // }
   // std::cerr << "\nOffSet:\n";
   // for(size_t i = 0; i < offSet.size()-1; i++){
   //    std::cerr << offSet[i] << " ";
   // }
   // std::cerr << "\n";
   std::vector<uint32_t>().swap(uniqueHeadsPerm);
   std::vector<uint32_t>().swap(offSetPerm);
   std::vector<uint32_t>().swap(uniqueHeadsPosVec);


   uniqueHeadsPosPrefSum[0] = 0;
   for(size_t i = 1; i < _n; i++){
      uniqueHeadsPosPrefSum[i] = uniqueHeadsPosPrefSum[i-1] + uniqueHeadsPos[i-1];
   }
   delete[] uniqueHeadsPos;

   if(verbose) for(uint64_t i = 0; i < phrases.size(); i++){
      std::cerr << headsSA[i].pos << "," << headsSA[i].len << " --> " << stringHeads[i] << "\n";
   }
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Reranking took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to invert stringHeads\n";
   int64_t hLen = headsSA.size();

   int32_t *stringHeads2 = new int32_t[headsSA.size()+1];
   for(uint32_t i = 0; i < headsSA.size(); i++){
      stringHeads2[headsSA[i].start] = stringHeads[i];
   }
   for(uint32_t i = 1; i < headsSA.size()+1; i++){
      stringHeads[i] = stringHeads2[i];
   }
   delete[] stringHeads2;
   std::vector<MatchSA>().swap(headsSA);
   stringHeads[hLen] = 0;
   // //print stringHeads
   // for(uint64_t i = 0; i < phrases.size()+1; i++){
   //    std::cerr << stringHeads[i] << "\n";
   // }
   // exit(0);
   std::cerr << "Inverted stringHeads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Inverting stringHeads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   uint32_t bufferLibSais = 60000;
   int32_t *indicesSA = new int32_t[hLen+1+bufferLibSais]();
   // for(uint64_t i = 0; i < phrases.size(); i++){
   //    headsSufArr[i] = phrases[i].pos << 31 + phrases[i].len;
   // }
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute suffix array of MS-heads\n";
   libsais_int(stringHeads, indicesSA, hLen+1, rank+1, bufferLibSais);
   //sacak_int(stringHeads, indicesSA, hLen+1, rank+1);
   std::cerr << "Computed suffix array of MS-heads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing MS-heads SA took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   if(verbose) for(uint64_t i = 0; i < phrases.size()+1; i++){
      std::cerr << stringHeads[i] << ", " << indicesSA[i] << "\n";
   }
   delete[] stringHeads;
   int32_t *indicesSA2 = new int32_t[hLen+1];
   uint8_t *BWTheads = new uint8_t[hLen+1];
   for(uint32_t i = 0; i < hLen+1; i++){
      //std::cerr << "indicesSA[" << i << "] = " << indicesSA[i] << "\n";
      indicesSA2[indicesSA[i]] = i;
      BWTheads[i] = listOfChars[indicesSA[i+1]];
      //std::cerr << "BWTheads[" << i << "] = " << BWTheads[i] << "\n";
   }
   //exit(0);
   std::vector<uint8_t>().swap(listOfChars);
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to permute headsSA\n";

   inPhrases.seekg(0, std::ios::beg);
   phrases.resize(hLen);
   for(uint32_t i = 0; i < hLen; i++){
      inPhrases.read((char*)&phrases[i], sizeof(Match));
   }

   predecessor2 *pHeads = new predecessor2(phrases, headBoundaries, ndoc, maxValue);
   const std::vector<Match>::iterator begPhrases = phrases.begin();
   //auto firstHead = std::begin(phrases);
   std::cerr << "New for loop\n";
   uint64_t currentSequence = 1;
   uint32_t *backupPos = new uint32_t[hLen+1];
   for(size_t i = 0; i < hLen; i++){
      //std::cerr << "i = " << i << "\n";
      //std::cerr << "phrases[i].start = " << phrases[i].start << " phrases[i].len = " << phrases[i].len << " \n";
      if(phrases[i].len == 0){
         backupPos[i] = phrases[i].pos;
         currentSequence++;
         continue;
      }
      uint32_t nextStart = phrases[i].start + phrases[i].len;
      const Match headNextStart = Match(nextStart, 0, currentSequence);
      std::vector<Match>::iterator nextHead = pHeads->predQuery(headNextStart, begPhrases);
      //std::cerr << "nextHead - begPhrases " << nextHead - begPhrases << "\n";
      //std::cerr << "indicesSA2[nextHead - begPhrases] - 1 = " << indicesSA2[nextHead - begPhrases] - 1 << "\n";
      //std::cerr << "_ISA[nextHead->pos + (nextStart - nextHead->start)] = " << _ISA[nextHead->pos + (nextStart - nextHead->start)] << "\n";
      backupPos[i] = phrases[i].pos;
      phrases[i].pos = _ISA[nextHead->pos + (nextStart - nextHead->start)];
      indicesSA2[i] = indicesSA2[nextHead - begPhrases] - 1;
   }

   std::ofstream outPermutedPhrases("permutedPhrases.txt", std::ios::out | std::ios::binary);
   std::vector<Match> bufferPermutedPhrases;
   bufferPermutedPhrases.reserve(1024*1024/sizeof(Match));
   for(uint32_t i = 1; i < hLen+1; i++){
      bufferPermutedPhrases.push_back(phrases[indicesSA[i]]);
      if(bufferPermutedPhrases.size() == bufferPermutedPhrases.capacity()){
         outPermutedPhrases.write((char*)&bufferPermutedPhrases[0], bufferPermutedPhrases.size()*sizeof(Match));
         bufferPermutedPhrases.clear();
      }
      //outPermutedPhrases.write((char*)&phrases[indicesSA[i]], sizeof(Match));
   }
   if(bufferPermutedPhrases.size() > 0){
      outPermutedPhrases.write((char*)&bufferPermutedPhrases[0], bufferPermutedPhrases.size()*sizeof(Match));
      bufferPermutedPhrases.clear();
   }
   outPermutedPhrases.close();
   std::vector<Match>().swap(phrases);
   std::vector<Match>().swap(bufferPermutedPhrases);
   std::cerr << "New for loop end\n";

   //print backupPos
   //for(uint32_t i = 0; i < hLen+1; i++){
   //  std::cerr << "backupPos[" << i << "] = " << backupPos[i] << "\n";
   //}
   std::vector<MatchSANN> headsSANN;
   headsSANN.resize(hLen);
   // std::ofstream outHeads("headsSA.txt", std::ios::out | std::ios::binary);
   // outHeads.seekp(hLen*sizeof(MatchSA));
   std::ifstream inPermutedPhrases("permutedPhrases.txt", std::ios::in | std::ios::binary);
   for(int64_t i = 1; i < hLen+1; i++) {
      //headsSA[i-1] = MatchSA(indicesSA2[indicesSA[i]], phrases[indicesSA[i]].pos, phrases[indicesSA[i]].len, phrases[indicesSA[i]].smaller, listOfChars[indicesSA[i]]);
      //std::cerr << "backupPos[indicesSA[" << i << "]] = " << backupPos[indicesSA[i]] << "\n";
      //std::cerr << "prefSumBucketPosHeads[backupPos[indicesSA[" << i << "]]] = " << prefSumBucketPosHeads[backupPos[indicesSA[i]]] << "\n";
      //auto tmp = MatchSA(indicesSA2[indicesSA[i]], phrases[indicesSA[i]].pos, phrases[indicesSA[i]].len, phrases[indicesSA[i]].smaller, (char)0);
      //headsSANN[prefSumBucketPosHeads[backupPos[indicesSA[i]]]++] = MatchSANN(indicesSA2[indicesSA[i]], phrases[indicesSA[i]].pos, phrases[indicesSA[i]].len, phrases[indicesSA[i]].smaller); //tmp;
      inPermutedPhrases.read(bufferPhrase, sizeof(Match));
      itHead = *reinterpret_cast<Match*>(bufferPhrase);
      headsSANN[prefSumBucketPosHeads[backupPos[indicesSA[i]]]++] = MatchSANN(indicesSA2[indicesSA[i]], itHead.pos, itHead.len, itHead.smaller); //tmp;
      //phrases[indicesSA[i]].pos = backupPos[indicesSA[i]];
      backupPos[indicesSA[i]] = prefSumBucketPosHeads[backupPos[indicesSA[i]]] - 1;
      // outHeads.seekp((prefSumBucketPosHeads[backupPos[indicesSA[i]]] - 1)*sizeof(MatchSA));
      // outHeads.write((char*)&tmp, sizeof(MatchSA));
      //BWTheads[i-1] = listOfChars[indicesSA[i]];
      //std::cerr << headsSA[i-1].start << "," << headsSA[i-1].pos << "," << headsSA[i-1].len << "," << headsSA[i-1].smaller << "," << headsSA[i-1].next <<"\n";
   }
   //std::vector<uint8_t>().swap(listOfChars);
   inPermutedPhrases.close();
   delete[] prefSumBucketPosHeads;

   // save headsSA to file
   std::ofstream outHeads("headsSA.txt", std::ios::out | std::ios::binary);
   //outHeads.sync_with_stdio(false);
   std::vector<MatchSANN> bufferHeadsSA;
   bufferHeadsSA.reserve(1024*1024/sizeof(MatchSANN));
   for(uint64_t i = 1; i < headsSANN.size(); i++){
      bufferHeadsSA.push_back(headsSANN[backupPos[i]]);
      if(bufferHeadsSA.size() == bufferHeadsSA.capacity()){
         outHeads.write((char*)&bufferHeadsSA[0], bufferHeadsSA.size()*sizeof(MatchSANN));
         bufferHeadsSA.clear();
      }
      //outHeads.write((char*)&headsSANN[backupPos[i]], sizeof(MatchSANN));
      //out << headsSA[i].start << "," << headsSA[i].pos << "," << headsSA[i].len << "," << headsSA[i].smaller << "," << headsSA[i].next <<"\n";
   }
   if(bufferHeadsSA.size() > 0){
      outHeads.write((char*)&bufferHeadsSA[0], bufferHeadsSA.size()*sizeof(MatchSANN));
      bufferHeadsSA.clear();
   }
   outHeads.close();
   delete[] backupPos;

   std::cerr << "headBoundaries.size(): " << headBoundaries.size() << "\n";
   std::vector<uint32_t>().swap(headBoundaries);
   std::cerr << "Permuted headsSA\n";
   //delete[] backupPos;
   delete[] indicesSA;
   delete[] indicesSA2;
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Permuting MS-heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t headSortTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to sort heads: " << headSortTime << " milliseconds\n";

   delete pHeads;
   // Now:
   // headsSA[i] = (nextHeadRank, nextHeadISA[pos], length, smaller, mmchar);

   uint64_t *counterSmallerThanHead = new uint64_t[headsSANN.size()]();

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute counterSmallerThanHead\n";

   //unsigned char *BWT_collection = new unsigned char[_sn]();
   //BWT_collection[0] = sequenceSeparator;
   ndoc = 1;
   uint64_t idxCurrentHead = 1;
   i = 1; //ATTENTION: skip first char if first char is %
   MatchSANN currentHead;
   uint64_t currentHeadEnd = 0;

   //compute number of buckets that are critical
   uint32_t nCriticalBuckets = 0;
   for(uint32_t i = 0; i < _n; i++){
      if(bucketsForExpandedBWT[i] < 0){
         nCriticalBuckets++;
      }
   }
   bool useBuffer = arg.buffer ? true : false;
   std::cerr << "ratio of critical buckets = " << (double)nCriticalBuckets/(double)_n << "\n";
   // if((double)nCriticalBuckets/(double)_n > 0.3){
   //    useBuffer = true;
   // }
   bool smallBuffer = _n < 65536;
   //std::cerr << "smallBuffer = " << smallBuffer << "\n";
   uint32_t suffixesInBuffer = 0;
   uint64_t capOfBufferSuffixes;
   double RAMavailable = arg.buffer * 1000000000;
   if(smallBuffer){
      //capOfBufferSuffixes = RAMavailable/(sizeof(std::pair<uint32_t, MatchSANN>));//100000000;
      capOfBufferSuffixes = RAMavailable/(sizeof(MatchSANN));
   }
   else{
      //capOfBufferSuffixes = RAMavailable/(sizeof(std::pair<std::pair<uint32_t, uint32_t>, MatchSANN>));//100000000;//_n/65536 * 100;
      capOfBufferSuffixes = RAMavailable/(sizeof(std::pair<uint32_t, MatchSANN>));
   }
   std::cerr << "capOfBufferSuffixes = " << capOfBufferSuffixes << "\n";

   //compute some statistics to help with buffering the bufferSuffixes. Get a rough estimate of the number of suffixes per bucketsForExpandedBWT
   // uint64_t nSuffixesInCriticalBuckets = 0;
   // uint64_t *nSuffixesPerCriticalBucket = new uint64_t[_n]();
   // for(uint32_t i = 0; i < _n; i++){
   //    if(bucketsForExpandedBWT[i] < 0){
   //       nSuffixesInCriticalBuckets += -bucketsForExpandedBWT[i];
   //    }
   // }
   // std::cerr << "nSuffixesInCriticalBuckets = " << nSuffixesInCriticalBuckets << "\n";
   // for(uint32_t i = 0; i < _n; i++){
   //    if(bucketsForExpandedBWT[i] < 0){
   //       nSuffixesPerCriticalBucket[i] = (uint32_t)ceil((double)(-bucketsForExpandedBWT[i])/(double)nSuffixesInCriticalBuckets * 100);
   //       std::cerr << "bucketsForExpandedBWT[" << i << "] = " << bucketsForExpandedBWT[i] << "\n";
   //       std::cerr << "nSuffixesPerCriticalBucket[" << i << "] = " << nSuffixesPerCriticalBucket[i] << "\n";
   //       std::cerr << "ratio: " << (double)(-bucketsForExpandedBWT[i])/(double)nSuffixesInCriticalBuckets << "\n";
   //       std::cerr << "ceil: " << ceil((double)(-bucketsForExpandedBWT[i])/(double)nSuffixesInCriticalBuckets * 100) << "\n";
   //       //std::cerr << "nSuffixesPerCriticalBucket[" << i << "] = " << nSuffixesPerCriticalBucket[i] << "\n";
   //    }
   // }
   // std::cerr << "nSuffixesPerCriticalBucket computed\n";

   const std::vector<uint32_t>::iterator startUniqueHeads = uniqueHeads.begin();
   const std::vector<MatchSANN>::iterator startHeadsSA = headsSANN.begin();
   std::ifstream inHeads("headsSA.txt", std::ios::binary | std::ios::in);
   //inHeads.sync_with_stdio(false);
   char bufferRead[sizeof(MatchSANN)];
   //useBuffer = false;
   Match firstHead, secondHead;
   inPhrases.seekg(0, std::ios::beg);
   if(useBuffer){
      std::cerr << "Using buffer\n";
      
      if(smallBuffer){
         //std::vector<std::vector<std::pair<uint32_t, MatchSANN>>> bufferSuffixes;
         std::vector<std::vector<MatchSANN>> bufferSuffixes;
         bufferSuffixes.resize(_n);
         for(uint32_t i = 0; i < _n; i++){
            bufferSuffixes[i].reserve(capOfBufferSuffixes/_n + 1);
         }
         std::cerr << "Reserved space for bufferSuffixes\n";
         inPhrases.read(bufferPhrase, sizeof(Match));
         inPhrases.read(bufferPhrase, sizeof(Match));
         secondHead = *reinterpret_cast<Match*>(bufferPhrase);
         for(uint64_t i = 1; i < sizePhrases; i++){
            firstHead = secondHead;
            inPhrases.read(bufferPhrase, sizeof(Match));
            secondHead = *reinterpret_cast<Match*>(bufferPhrase);
            if(firstHead.len == 0){
               //BWT_collection[ndoc] = BWTheads[ndoc];
               idxCurrentHead++;
               inHeads.read(bufferRead, sizeof(MatchSANN));
               ndoc++;
               continue;
            }
            //currentHead = headsSA[backupPos[idxCurrentHead++]];
            inHeads.read(bufferRead, sizeof(MatchSANN));
            currentHead = *reinterpret_cast<MatchSANN*>(bufferRead);
            uint32_t posSuffix = firstHead.pos + 1;
            uint32_t lenSuffix = firstHead.len - 1;
            for(uint64_t iHead = 1; iHead < secondHead.start - firstHead.start; iHead++){
               if(bucketsForExpandedBWT[posSuffix] < 0) [[likely]]{
                  //bufferSuffixes[posSuffix].push_back({lenSuffix, currentHead});
                  bufferSuffixes[posSuffix].push_back(MatchSANN(currentHead.start, currentHead.pos, lenSuffix, currentHead.smaller));
                  suffixesInBuffer++;
                  if(suffixesInBuffer == capOfBufferSuffixes){ //process bufferSuffixes
                     //in each bucket count how many suffixes are smaller than the heads
                     for(uint32_t pos = 0; pos < _n - 1; pos++){
                        for(uint32_t counter = 0; counter < bufferSuffixes[pos].size(); counter++){
                           //check if the len of the suffix is dfferent than the len of the uniqueHead
                           // //print bufferSuffixes[pos][counter]
                           // std::cerr << "bufferSuffixes[" << pos << "][" << counter << "] = " << bufferSuffixes[pos][counter].start << " " << bufferSuffixes[pos][counter].pos << " " << bufferSuffixes[pos][counter].len << " " << bufferSuffixes[pos][counter].smaller << "\n";
                           uint32_t counterHeadLen = bufferSuffixes[pos][counter].smaller ? bufferSuffixes[pos][counter].len : MAXIMUM_UINT32 - bufferSuffixes[pos][counter].len;
                           // std::cerr << "counterHeadLen = " << counterHeadLen << "\n";
                           // //print uniqueHeads[uniqueHeadsPosPrefSum[pos]] to uniqueHeads[uniqueHeadsPosPrefSum[pos+1]-1]
                           // for(uint32_t i = uniqueHeadsPosPrefSum[pos]; i < uniqueHeadsPosPrefSum[pos+1]; i++){
                           //    std::cerr << "uniqueHeads[" << i << "] = " << uniqueHeads[i] << "\n";
                           // }
                           //timesRank++;
                           auto apprxPointer = std::upper_bound(startUniqueHeads + uniqueHeadsPosPrefSum[pos], startUniqueHeads + uniqueHeadsPosPrefSum[pos+1], counterHeadLen) - startUniqueHeads;
                           if(apprxPointer == uniqueHeadsPosPrefSum[pos]){
                              counterSmallerThanHead[prefSumBucketPosHeadsCopy[pos]]++;
                              //timesLength++;
                           }
                           else if((apprxPointer == uniqueHeadsPosPrefSum[pos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){
                              //timesLength++;
                              continue;
                           }
                           else if((apprxPointer != uniqueHeadsPosPrefSum[pos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){ //maybe with -1
                              counterSmallerThanHead[prefSumBucketPosHeadsCopy[pos] + offSet[apprxPointer]]++;
                              //timesLength++;
                           }
                           else{
                              //timesISA++;
                              auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[pos] + offSet[apprxPointer-1], 
                                                               startHeadsSA + (prefSumBucketPosHeadsCopy[pos]+offSet[apprxPointer])*(offSet[apprxPointer]>0) + 
                                                               (prefSumBucketPosHeadsCopy[pos+1])*(offSet[apprxPointer]==0), 
                                                               bufferSuffixes[pos][counter],
                                                [](const MatchSANN &headA, const MatchSANN &headB){
                                                   //timesRank++;
                                                   // if(headA.len != headB.len){
                                                   //    //timesLength++;
                                                   //    return headA.smaller*(headA.len < headB.len) +
                                                   //          !headB.smaller*(headA.len > headB.len);
                                                   // }
                                                   // else{
                                                      return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                                   // }
                                                });
                              //get index of pointer
                              uint64_t toIncrement = pointer - startHeadsSA;
                              //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
                              if(toIncrement != prefSumBucketPosHeadsCopy[pos+1]) counterSmallerThanHead[toIncrement]++;
                           }
                        }
                        bufferSuffixes[pos].clear();
                     }
                     suffixesInBuffer = 0;
                  }
               }
               posSuffix++;
               lenSuffix--;
            }
         }
         std::cerr << "suffixesInBuffer left = " << suffixesInBuffer << "\n";
         if(suffixesInBuffer != 0){
            //in each bucket count how many suffixes are smaller than the heads
            for(uint32_t pos = 0; pos < _n; pos++){
               for(uint32_t counter = 0; counter < bufferSuffixes[pos].size(); counter++){
                  //check if the len of the suffix is dfferent than the len of the uniqueHead
                  // //print bufferSuffixes[pos][counter]
                  // std::cerr << "bufferSuffixes[" << pos << "][" << counter << "] = " << bufferSuffixes[pos][counter].start << " " << bufferSuffixes[pos][counter].pos << " " << bufferSuffixes[pos][counter].len << " " << bufferSuffixes[pos][counter].smaller << "\n";
                  uint32_t counterHeadLen = bufferSuffixes[pos][counter].smaller ? bufferSuffixes[pos][counter].len : MAXIMUM_UINT32 - bufferSuffixes[pos][counter].len;
                  // std::cerr << "counterHeadLen = " << counterHeadLen << "\n";
                  // //print uniqueHeads[uniqueHeadsPosPrefSum[pos]] to uniqueHeads[uniqueHeadsPosPrefSum[pos+1]-1]
                  // for(uint32_t i = uniqueHeadsPosPrefSum[pos]; i < uniqueHeadsPosPrefSum[pos+1]; i++){
                  //    std::cerr << "uniqueHeads[" << i << "] = " << uniqueHeads[i] << "\n";
                  // }           
                  //timesRank++;       
                  auto apprxPointer = std::upper_bound(startUniqueHeads + uniqueHeadsPosPrefSum[pos], startUniqueHeads + uniqueHeadsPosPrefSum[pos+1], counterHeadLen) - startUniqueHeads;
                  if(apprxPointer == uniqueHeadsPosPrefSum[pos]){
                     // std::cerr << "apprxPointer == uniqueHeadsPosPrefSum[pos]\n";
                     counterSmallerThanHead[prefSumBucketPosHeadsCopy[pos]]++;
                     //timesLength++;
                  }
                  else if((apprxPointer == uniqueHeadsPosPrefSum[pos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){
                     // std::cerr << "apprxPointer == uniqueHeadsPosPrefSum[pos+1] & counterHeadLen != uniqueHeads[apprxPointer - 1]\n";
                     //timesLength++;
                     continue;
                  }
                  else if((apprxPointer != uniqueHeadsPosPrefSum[pos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){ //maybe with -1
                     // std::cerr << "counterHeadLen != uniqueHeads[apprxPointer - 1]\n";
                     //timesLength++;
                     counterSmallerThanHead[prefSumBucketPosHeadsCopy[pos] + offSet[apprxPointer]]++;
                  }
                  else{
                     // std::cerr << "else\n";
                     //timesISA++;
                     auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[pos] + offSet[apprxPointer-1], 
                                                      startHeadsSA + (prefSumBucketPosHeadsCopy[pos]+offSet[apprxPointer])*(offSet[apprxPointer]>0) + 
                                                      (prefSumBucketPosHeadsCopy[pos+1])*(offSet[apprxPointer]==0), 
                                                      bufferSuffixes[pos][counter],
                                       [](const MatchSANN &headA, const MatchSANN &headB){
                                          //timesRank++;
                                          // if(headA.len != headB.len){
                                          //    //timesLength++;
                                          //    return headA.smaller*(headA.len < headB.len) +
                                          //          !headB.smaller*(headA.len > headB.len);
                                          // }
                                          // else{
                                             return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                          // }
                                       });
                     //get index of pointer
                     uint64_t toIncrement = pointer - startHeadsSA;
                     //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
                     if(toIncrement != prefSumBucketPosHeadsCopy[pos+1]) counterSmallerThanHead[toIncrement]++;
                  }
               }
            }
         }
         //std::vector<std::vector<std::pair<uint32_t, MatchSANN>>>().swap(bufferSuffixes);
         std::vector<std::vector<MatchSANN>>().swap(bufferSuffixes);
         //std::cerr << "timesLength = " << timesLength << "\n";
         //std::cerr << "timesRank = " << timesRank << "\n";
         //std::cerr << "timesLength/timesRank = " << (double)timesLength/timesRank << "\n";         
      }
      else{
         //buffer for big references (it has fixed size of 65536 and it is indexes by the first 16 bits of the suffix). It stores the head and the correct position of the suffix
         //std::vector<std::vector<std::pair<std::pair<uint32_t, uint32_t>, MatchSANN>>> bufferSuffixes;
         std::vector<std::vector<std::pair<uint32_t,MatchSANN>>> bufferSuffixes;
         bufferSuffixes.resize(65536);
         for(uint32_t i = 0; i < 65536; i++){
            bufferSuffixes[i].reserve(capOfBufferSuffixes/65536 + 1);
         }
         uint32_t max = _n;
         //compute shift to get the first 16 bits of the suffix
         uint32_t shift = 0;
         for(uint8_t i = 0; i < 48; i++){
            if(max < nOfDigitsBig[i]){
               shift = i;
               break;
            }
         } 
         inPhrases.read(bufferPhrase, sizeof(Match));
         inPhrases.read(bufferPhrase, sizeof(Match));
         secondHead = *reinterpret_cast<Match*>(bufferPhrase);
         for(uint64_t i = 1; i < sizePhrases; i++){
            firstHead = secondHead;
            inPhrases.read(bufferPhrase, sizeof(Match));
            secondHead = *reinterpret_cast<Match*>(bufferPhrase);
            if(firstHead.len == 0){
               //BWT_collection[ndoc] = BWTheads[ndoc];
               idxCurrentHead++;
               inHeads.read(bufferRead, sizeof(MatchSANN));
               ndoc++;
               continue;
            }
            //currentHead = headsSA[backupPos[idxCurrentHead++]];
            inHeads.read(bufferRead, sizeof(MatchSANN));
            currentHead = *reinterpret_cast<MatchSANN*>(bufferRead);
            //std::cerr << "currentHead.start = " << currentHead.start << " currentHead.pos = " << currentHead.pos << " currentHead.len = " << currentHead.len << " currentHead.smaller = " << currentHead.smaller << "\n";
            //std::cerr << "phrases[i].pos = " << phrases[i].pos << " phrases[i].len = " << phrases[i].len << '\n';
            //std::cerr << "phrases[i+1].start - phrases[i].start = " << phrases[i+1].start - phrases[i].start << '\n';
            uint32_t posSuffix = firstHead.pos + 1;
            uint32_t lenSuffix = firstHead.len - 1;
            for(uint64_t iHead = 1; iHead < secondHead.start - firstHead.start; iHead++){
               if(bucketsForExpandedBWT[posSuffix] < 0){
                  //bufferSuffixes[posSuffix >> shift].push_back({{lenSuffix, posSuffix}, currentHead});
                  bufferSuffixes[posSuffix >> shift].push_back({posSuffix, MatchSANN(currentHead.start, currentHead.pos, lenSuffix, currentHead.smaller)});
                  suffixesInBuffer++;
                  if(suffixesInBuffer == capOfBufferSuffixes){ //process bufferSuffixes
                     //sort by correct position in each 16-bit bucket
                     for(uint32_t area = 0; area < 65536; area++){
                        std::sort(bufferSuffixes[area].begin(), bufferSuffixes[area].end(), [](const std::pair<uint32_t, MatchSANN> &a, const std::pair<uint32_t, MatchSANN> &b){
                           return a.first < b.first;
                        });
                     }
                     //in each bucket count how many suffixes are smaller than the heads
                     for(uint32_t pos = 0; pos < 65536; pos++){
                        for(uint32_t counter = 0; counter < bufferSuffixes[pos].size(); counter++){
                           //uint32_t len = bufferSuffixes[pos][counter].first.first;
                           //uint32_t realPos = bufferSuffixes[pos][counter].first.second;
                           uint32_t realPos = bufferSuffixes[pos][counter].first;
                           uint32_t counterHeadLen = bufferSuffixes[pos][counter].second.smaller ? bufferSuffixes[pos][counter].second.len : MAXIMUM_UINT32 - bufferSuffixes[pos][counter].second.len;    
                           auto apprxPointer = std::upper_bound(startUniqueHeads + uniqueHeadsPosPrefSum[realPos], startUniqueHeads + uniqueHeadsPosPrefSum[realPos+1], counterHeadLen) - startUniqueHeads;
                           if(apprxPointer == uniqueHeadsPosPrefSum[realPos]){
                              counterSmallerThanHead[prefSumBucketPosHeadsCopy[realPos]]++;
                           }
                           else if((apprxPointer == uniqueHeadsPosPrefSum[realPos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){
                              continue;
                           }
                           else if((apprxPointer != uniqueHeadsPosPrefSum[realPos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){ //maybe with -1
                              counterSmallerThanHead[prefSumBucketPosHeadsCopy[realPos] + offSet[apprxPointer]]++;
                           }
                           else{
                              auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[realPos] + offSet[apprxPointer-1], 
                                                               startHeadsSA + (prefSumBucketPosHeadsCopy[realPos]+offSet[apprxPointer])*(offSet[apprxPointer]>0) + 
                                                               (prefSumBucketPosHeadsCopy[realPos+1])*(offSet[apprxPointer]==0), 
                                                               bufferSuffixes[pos][counter].second,
                                                [](const MatchSANN &headA, const MatchSANN &headB){
                                                      return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                                });
                              //get index of pointer
                              uint64_t toIncrement = pointer - startHeadsSA;
                              //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
                              if(toIncrement != prefSumBucketPosHeadsCopy[realPos+1]) counterSmallerThanHead[toIncrement]++;
                           }
                        }
                        bufferSuffixes[pos].clear();
                     }
                     suffixesInBuffer = 0;
                  }
                  //__builtin_prefetch(&headsSA[prefSumBucketPosHeadsCopy[posSuffix+1]], 0, 3);
               }
               posSuffix++;
               lenSuffix--;
            }
         }
         std::cerr << "suffixesInBuffer left = " << suffixesInBuffer << "\n";
         if(suffixesInBuffer != 0){
            for(uint32_t area = 0; area < 65536; area++){
               std::sort(bufferSuffixes[area].begin(), bufferSuffixes[area].end(), [](const std::pair<uint32_t, MatchSANN> &a, const std::pair<uint32_t, MatchSANN> &b){
                  return a.first < b.first;
               });
            }
            for(uint32_t pos = 0; pos < 65536; pos++){
               for(uint32_t counter = 0; counter < bufferSuffixes[pos].size(); counter++){
                  //uint32_t len = bufferSuffixes[pos][counter].first.first;
                  //uint32_t realPos = bufferSuffixes[pos][counter].first.second;
                  uint32_t realPos = bufferSuffixes[pos][counter].first;
                  uint32_t counterHeadLen = bufferSuffixes[pos][counter].second.smaller ? bufferSuffixes[pos][counter].second.len : MAXIMUM_UINT32 - bufferSuffixes[pos][counter].second.len;    
                  auto apprxPointer = std::upper_bound(startUniqueHeads + uniqueHeadsPosPrefSum[realPos], startUniqueHeads + uniqueHeadsPosPrefSum[realPos+1], counterHeadLen) - startUniqueHeads;
                  if(apprxPointer == uniqueHeadsPosPrefSum[realPos]){
                     counterSmallerThanHead[prefSumBucketPosHeadsCopy[realPos]]++;
                  }
                  else if((apprxPointer == uniqueHeadsPosPrefSum[realPos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){
                     continue;
                  }
                  else if((apprxPointer != uniqueHeadsPosPrefSum[realPos+1]) & (counterHeadLen != uniqueHeads[apprxPointer - 1])){ //maybe with -1
                     counterSmallerThanHead[prefSumBucketPosHeadsCopy[realPos] + offSet[apprxPointer]]++;
                  }
                  else{
                     auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[realPos] + offSet[apprxPointer-1], 
                                                      startHeadsSA + (prefSumBucketPosHeadsCopy[realPos]+offSet[apprxPointer])*(offSet[apprxPointer]>0) + 
                                                      (prefSumBucketPosHeadsCopy[realPos+1])*(offSet[apprxPointer]==0), 
                                                      bufferSuffixes[pos][counter].second,
                                       [](const MatchSANN &headA, const MatchSANN &headB){
                                             return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                       });
                     //get index of pointer
                     uint64_t toIncrement = pointer - startHeadsSA;
                     //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
                     if(toIncrement != prefSumBucketPosHeadsCopy[realPos+1]) counterSmallerThanHead[toIncrement]++;
                  }
               }
            }
         }
         //std::vector<std::vector<std::pair<std::pair<uint32_t, uint32_t>, MatchSANN>>>().swap(bufferSuffixes);
         std::vector<std::vector<std::pair<uint32_t, MatchSANN>>>().swap(bufferSuffixes);
      }
   }
   else{
      //just count how many suffixes are smaller than the heads because the suffixes are few
      inPhrases.read(bufferPhrase, sizeof(Match));
      inPhrases.read(bufferPhrase, sizeof(Match));
      secondHead = *reinterpret_cast<Match*>(bufferPhrase);
      for(uint64_t i = 1; i < sizePhrases; i++){
         firstHead = secondHead;
         inPhrases.read(bufferPhrase, sizeof(Match));
         secondHead = *reinterpret_cast<Match*>(bufferPhrase);
         if(firstHead.len == 0){
            //BWT_collection[ndoc] = BWTheads[ndoc];
            idxCurrentHead++;
            inHeads.read(bufferRead, sizeof(MatchSANN));
            ndoc++;
            continue;
         }
         
         //currentHead = headsSA[backupPos[idxCurrentHead++]];
         inHeads.read(bufferRead, sizeof(MatchSANN));
         currentHead = *reinterpret_cast<MatchSANN*>(bufferRead);
         //std::cerr << "currentHead.start = " << currentHead.start << " currentHead.pos = " << currentHead.pos << " currentHead.len = " << currentHead.len << " currentHead.smaller = " << currentHead.smaller << "\n";
         //std::cerr << "phrases[i].pos = " << phrases[i].pos << " phrases[i].len = " << phrases[i].len << '\n';
         //std::cerr << "phrases[i+1].start - phrases[i].start = " << phrases[i+1].start - phrases[i].start << '\n';
         uint32_t posSuffix = firstHead.pos + 1;
         uint32_t lenSuffix = firstHead.len - 1;
         // std::cerr << "secondHead.start - firstHead.start = " << secondHead.start - firstHead.start << '\n';
         // std::cerr << "secondHead.len " << secondHead.len << '\n';
         // std::cerr << "posSuffix = " << posSuffix << " lenSuffix = " << lenSuffix << '\n';
         for(uint64_t iHead = 1; iHead < secondHead.start - firstHead.start; iHead++){
            // std::cerr << "posSuffix = " << posSuffix << " lenSuffix = " << lenSuffix << '\n';
            //__builtin_prefetch(&headsSA[prefSumBucketHeads_SA[posSuffix].first], 0, 3);
            if(bucketsForExpandedBWT[posSuffix] < 0){
               auto pointer = std::upper_bound(startHeadsSA + prefSumBucketPosHeadsCopy[posSuffix], startHeadsSA + prefSumBucketPosHeadsCopy[posSuffix+1], currentHead,
                                 [lenSuffix](const MatchSANN &headA, const MatchSANN &headB){
                                    if(lenSuffix != headB.len){
                                       return headA.smaller*(lenSuffix < headB.len) +
                                             !headB.smaller*(lenSuffix > headB.len);
                                    }
                                    else{
                                       return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
                                    }
                                 });
               //get index of pointer
               uint64_t toIncrement = pointer - startHeadsSA;
               //if the pointer is not the last element + 1 of the bucket, increment counterSmallerThanHead
               if(toIncrement != prefSumBucketPosHeadsCopy[posSuffix+1]) counterSmallerThanHead[toIncrement]++;
               //__builtin_prefetch(&headsSA[prefSumBucketPosHeadsCopy[posSuffix+1]], 0, 3);
            }
            posSuffix++;
            lenSuffix--;
            //std::cerr << "posSuffix = " << posSuffix << " lenSuffix = " << lenSuffix << '\n';
         }
      }
   }
   inHeads.close();
   inPhrases.close();
   std::vector<Match>().swap(phrases);
   std::vector<MatchSANN>().swap(headsSANN);
   std::vector<uint32_t>().swap(uniqueHeads);
   delete[] uniqueHeadsPosPrefSum;
   std::vector<uint32_t>().swap(offSet);
   delete[] prefSumBucketPosHeadsCopy;

   int64_t *bucketsForExpandedBWT_copy = new int64_t[_n]();
   for(uint32_t i = 0; i < _n; i++){
      bucketsForExpandedBWT_copy[i] = bucketsForExpandedBWT[i];
   }
   //invert data in bucketsForExpandedBWT w.r.t. _ISA
   for(uint32_t i = 0; i < _n; i++){
      bucketsForExpandedBWT[_ISA[i]] = bucketsForExpandedBWT_copy[i];
   }
   std::cerr << "Inverted data in bucketsForExpandedBWT\n";
   //print bucketsForExpandedBWT
   //for(uint64_t i = 0; i < _n; i++){
   //  std::cerr << "_SA["<<i<<"]: " << _SA[i] << ", " << bucketsForExpandedBWT[i] << "\n";
   //}
   delete[] bucketsForExpandedBWT_copy;

   std::cerr << "Computed counterSmallerThanHead\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing counterSmallerThanHead took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t1 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute final BWT\n";
   uint64_t pointerNextChar = D - 1;
   //print bucketsForExpandedBWT
   //std::cerr << "bucketsForExpandedBWT: \n";
   if(verbose) for(uint32_t i = 0; i < _n; i++){
      std::cerr << "i: " << i << " - " << (int64_t)bucketsForExpandedBWT[i] << "\n";
   }

   //print prefSumBucketISAPosHeadsCopy
   //std::cerr << "prefSumBucketISAPosHeadsCopy: ";
   if(verbose) for(uint32_t i = 0; i < _n; i++){
      std::cerr << "i: " << i << " - " << prefSumBucketISAPosHeadsCopy[i] << "\n";
   }

   //permute counterSmallerThanHead w.r.t. _SA
   // print counterSmallerThanHead
   if(verbose) for(uint64_t i = 0; i < hLen; i++){
     std::cerr << "counterSmallerThanHead[" << i << "]: " << counterSmallerThanHead[i] << "\n";
   }

   uint64_t *counterSmallerThanHead_copy = new uint64_t[hLen]();
   for(uint64_t i = 0; i < hLen; i++){
      counterSmallerThanHead_copy[i] = counterSmallerThanHead[i];
   }
   uint64_t a = 0;
   std::vector<std::pair<uint64_t, uint64_t>> prefSumBucketHeads_SA(_n);
   for(size_t i = 0; i < _n; i++){
      prefSumBucketHeads_SA[_SA[i]] = std::make_pair(prefSumBucketISAPosHeadsCopy[i], prefSumBucketISAPosHeadsCopy[i+1]);
   }
   for(uint64_t i = 0; i < _n; i++){
      //std::cerr << "prefSumBucketHeads_SA[" << i << "].first: " << prefSumBucketHeads_SA[i].first << "\n";
      //std::cerr << "prefSumBucketHeads_SA[" << i << "].second: " << prefSumBucketHeads_SA[i].second << "\n";
      for(uint64_t counter = prefSumBucketHeads_SA[i].first; counter < prefSumBucketHeads_SA[i].second; counter++){
         //std::cerr << "a: " << a << "\n";
         counterSmallerThanHead[counter] = counterSmallerThanHead_copy[a++];
      }
   }
   delete[] counterSmallerThanHead_copy;
   std::vector<std::pair<uint64_t, uint64_t>>().swap(prefSumBucketHeads_SA);

   //print counterSmallerThanHead
   if(verbose) for(uint64_t i = 0; i < hLen; i++){
     std::cerr << "counterSmallerThanHead[" << i << "]: " << counterSmallerThanHead[i] << "\n";
   }

   bool RLEWrite = arg.format;
   //uint32_t pointerHeads = D-1;
   if(!RLEWrite){
      std::ofstream streamOutfile(arg.outname+".bwt", std::ios::out | std::ios::binary);
      //streamOutfile.sync_with_stdio(false);
      std::vector<uint8_t> bufferWrite;
      bufferWrite.reserve(1024*1024);
      //FILE *f = fopen("BWT_collection", "wb");

      //std::cerr << "pointerNextChar: " << pointerNextChar << "\n";
      //write BWT_collection to file
      streamOutfile.write((char*)&BWTheads[0], sizeof(uint8_t)*(D-1));
      //for(uint32_t d = 0; d < D - 1; d++){
         //BWT_collection[d] = BWTheads[d];
      //   streamOutfile.write((char*)&BWTheads[d], sizeof(uint8_t));
         //fwrite(&BWTheads[d], sizeof(uint8_t), 1, f);
      //}
      //write BWT_collection to file
      for(uint32_t i = 1; i < _x.size(); i++){
         //std::cerr << "i: " << i << "\n";
         //uint8_t characterEqual = _SA[i] ? _x[_SA[i]-1] : '$';
         uint8_t characterEqual = _BWT[i];
         //std::cerr << "characterEqual: " << characterEqual << "\n";
         if(bucketsForExpandedBWT[i] >= 0){
            //std::cerr << "Safe bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            for(uint64_t counter = 0; counter < bucketsForExpandedBWT[i]; counter++){
               //BWT_collection[pointerNextChar++] = characterEqual;
               //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
               bufferWrite.push_back(characterEqual);
               if(bufferWrite.size() == bufferWrite.capacity()){
                  streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
                  bufferWrite.clear();
               }
               //streamOutfile << characterEqual;
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
            }
         }
         else{
            //std::cerr << "Unsafe bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            for(uint64_t counter = prefSumBucketISAPosHeadsCopy[i]; counter < prefSumBucketISAPosHeadsCopy[i+1]; counter++){
               //std::cerr << "counter: " << counter << "\n";
               //std::cerr << "counterSmallerThanHead[counter]: " << counterSmallerThanHead[counter] << "\n";
               for(uint64_t counter2 = 0; counter2 < counterSmallerThanHead[counter]; counter2++){
                  //std::cerr << "counter2: " << counter2 << "\n";
                  //BWT_collection[pointerNextChar++] = characterEqual;
                  //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
                  bufferWrite.push_back(characterEqual);
                  if(bufferWrite.size() == bufferWrite.capacity()){
                     streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
                     bufferWrite.clear();
                  }
                  //streamOutfile << characterEqual;
                  //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                  bucketsForExpandedBWT[i]++;
               }
               //std::cerr << "headsSA[counter].next: " << headsSA[counter].next << "\n";
               //BWT_collection[pointerNextChar++] = BWTheads[counter];//headsSA[counter].next;
               //streamOutfile.write((char*)&BWTheads[counter], sizeof(uint8_t));
               bufferWrite.push_back(BWTheads[counter]);
               if(bufferWrite.size() == bufferWrite.capacity()){
                  streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
                  bufferWrite.clear();
               }
               //streamOutfile.write((char*)&BWTheads[pointerHeads++], sizeof(uint8_t));
               //streamOutfile << BWTheads[counter];
               //fwrite(&BWTheads[counter], sizeof(uint8_t), 1, f);
               bucketsForExpandedBWT[i]++;
            }
            //std::cerr << "Remaining suff in bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            for(int64_t counter = bucketsForExpandedBWT[i]; counter < 0; counter++){
               //BWT_collection[pointerNextChar++] = characterEqual;
               //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
               bufferWrite.push_back(characterEqual);
               if(bufferWrite.size() == bufferWrite.capacity()){
                  streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
                  bufferWrite.clear();
               }
               //streamOutfile << characterEqual;
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
            }
         }
      }
      if(bufferWrite.size() > 0){
         streamOutfile.write((char*)&bufferWrite[0], sizeof(uint8_t)*bufferWrite.size());
         bufferWrite.clear();
      }
      streamOutfile.close();
   }
   else{
      //write only when character is different
      uint8_t prevChar = (char)0;
      uint32_t runLength = 0;

      std::ofstream streamOutfile(arg.outname+".rl_bwt", std::ios::out | std::ios::binary);
      //streamOutfile.sync_with_stdio(false);
      //FILE *f = fopen("BWT_collection", "wb");

      //std::cerr << "pointerNextChar: " << pointerNextChar << "\n";
      //write BWT_collection to file
      //streamOutfile.write((char*)&BWTheads[0], sizeof(uint8_t)*(D-1));
      for(uint64_t d = 0; d < D-1; d++){
         if(prevChar != BWTheads[d]){
            if(runLength > 0){
               streamOutfile.write((char*)&runLength, sizeof(uint32_t));
               streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
               //fwrite(&runLength, sizeof(uint32_t), 1, f);
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
               runLength = 1;
               prevChar = BWTheads[d];
            }
            else{
               runLength = 1;
               prevChar = BWTheads[d];
            }
         }
         else{
            runLength++;
         }
      }
      //write BWT_collection to file
      for(uint32_t i = 1; i < _x.size(); i++){
         //std::cerr << "i: " << i << "\n";
         //uint8_t characterEqual = _SA[i] ? _x[_SA[i]-1] : '$';
         uint8_t characterEqual = _BWT[i];
         //std::cerr << "characterEqual: " << characterEqual << "\n";
         if(bucketsForExpandedBWT[i] > 0){
            //std::cerr << "Safe bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            //for(uint64_t counter = 0; counter < bucketsForExpandedBWT[i]; counter++){
            //BWT_collection[pointerNextChar++] = characterEqual;
            if(prevChar != characterEqual){
               streamOutfile.write((char*)&runLength, sizeof(uint32_t));
               streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
               //fwrite(&runLength, sizeof(uint32_t), 1, f);
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
               runLength = bucketsForExpandedBWT[i];
               prevChar = characterEqual;
            }
            else{
               runLength+=bucketsForExpandedBWT[i];
            }
            //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
            //}
         }
         else if(bucketsForExpandedBWT[i] < 0){
            //std::cerr << "Unsafe bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            for(uint64_t counter = prefSumBucketISAPosHeadsCopy[i]; counter < prefSumBucketISAPosHeadsCopy[i+1]; counter++){
               //std::cerr << "counter: " << counter << "\n";
               //std::cerr << "counterSmallerThanHead[counter]: " << counterSmallerThanHead[counter] << "\n";
               //for(uint64_t counter2 = 0; counter2 < counterSmallerThanHead[counter]; counter2++){
                  //std::cerr << "counter2: " << counter2 << "\n";
                  //BWT_collection[pointerNextChar++] = characterEqual;
               if(counterSmallerThanHead[counter]){
                  if(prevChar != characterEqual){
                     streamOutfile.write((char*)&runLength, sizeof(uint32_t));
                     streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
                     //fwrite(&runLength, sizeof(uint32_t), 1, f);
                     //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                     runLength = counterSmallerThanHead[counter];
                     prevChar = characterEqual;
                  }
                  else{
                     runLength+=counterSmallerThanHead[counter];
                  }
                  //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
                  //streamOutfile << characterEqual;
                  //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                  bucketsForExpandedBWT[i]+=counterSmallerThanHead[counter];
               }
               //}
               //std::cerr << "headsSA[counter].next: " << headsSA[counter].next << "\n";
               //BWT_collection[pointerNextChar++] = BWTheads[counter];//headsSA[counter].next;
               if(BWTheads[counter] != prevChar){
                  streamOutfile.write((char*)&runLength, sizeof(uint32_t));
                  streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
                  //fwrite(&runLength, sizeof(uint32_t), 1, f);
                  //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                  runLength = 1;
                  prevChar = BWTheads[counter];
               }
               else{
                  runLength++;
               }
               //streamOutfile.write((char*)&BWTheads[counter], sizeof(uint8_t));
               //streamOutfile.write((char*)&BWTheads[pointerHeads++], sizeof(uint8_t));
               //streamOutfile << BWTheads[counter];
               //fwrite(&BWTheads[counter], sizeof(uint8_t), 1, f);
               bucketsForExpandedBWT[i]++;
            }
            //std::cerr << "Remaining suff in bucket[" << i << "]: " << bucketsForExpandedBWT[i] << "\n";
            //for(int64_t counter = bucketsForExpandedBWT[i]; counter < 0; counter++){
            //BWT_collection[pointerNextChar++] = characterEqual;
            if(bucketsForExpandedBWT[i]!=0){
               if(prevChar != characterEqual){
                  streamOutfile.write((char*)&runLength, sizeof(uint32_t));
                  streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
                  //fwrite(&runLength, sizeof(uint32_t), 1, f);
                  //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
                  runLength = -bucketsForExpandedBWT[i];
                  prevChar = characterEqual;
               }
               else{
                  runLength+=-bucketsForExpandedBWT[i];
               }
               //streamOutfile.write((char*)&characterEqual, sizeof(uint8_t));
               //streamOutfile << characterEqual;
               //fwrite(&characterEqual, sizeof(uint8_t), 1, f);
            }
         }
      }
      streamOutfile.write((char*)&runLength, sizeof(uint32_t));
      streamOutfile.write((char*)&prevChar, sizeof(uint8_t));
      streamOutfile.close();
   }
   //fclose(f);
   std::cerr << "Computed final BWT\n";
   t2 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing final BWT took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms\n";
   delete[] bucketsForExpandedBWT;
   delete[] counterSmallerThanHead;
   delete[] BWTheads;
   //print BWT_collection
   //std::cerr << "BWT_collection: ";
   //for(uint32_t i = 0; i < _sx.size(); i++){
   //  std::cerr << BWT_collection[i];
   //}
   //std::cerr << "\n";

   delete[] _SA;
   delete[] _BWT;
   delete[] _ISA;
   delete[] prefSumBucketISAPosHeadsCopy;
   std::string().swap(_x);
   //std::string().swap(_sx);
   //std::vector<Match>().swap(phrases);
   std::remove("phrases.txt");
   std::remove("permutedPhrases.txt");
   std::remove("headsSA.txt");
   
   if(arg.check){
      std::string _sx_just_to_check;
      _sx_just_to_check.reserve(_sn);

      std::ifstream streamInfile(collFileName,std::ios::in);
      std::string line, content;
      content.reserve(_n*2);
      uint64_t charactersRead = 0;
      //std::string _sx_just_to_check;
      //_sx_just_to_check.reserve(_sn);
      while(std::getline(streamInfile, line).good()){
         //std::cerr << "ndoc: " << ndoc << "\n";
         if( line.empty() || line[0] == '>' ){
            //std::cerr << "content: " << content << "\n";
            content += sequenceSeparator;
            charactersRead++;
            _sx_just_to_check += content;
            //std::cerr << "Finished parsing document " << ndoc << "\n";
            content.erase();
            //std::string().swap(content);
         }
         else if (!line.empty()) {
            charactersRead += line.size();
            if(charactersRead >= _sn - 1){ //if string is filled up to sn (useful for prefixLength)
               //std::cerr << "read too much\n";
               //std::cerr << "charactersRead: " << charactersRead << "\n";
               //std::cerr << "line.size(): " << line.size() << "\n";
               //std::cerr << "_sn - charactersRead + line.size() - 1 = " << line.size() - (charactersRead - _sn) - 1 << "\n";
               //std::cerr << "content.size(): " << content.size() << "\n";
               content += line.substr(0, line.size() - (charactersRead - _sn) - 1);
               //std::cerr << "content.size(): " << content.size() << "\n";
               //std::cerr << _sx_just_to_check.size() << "\n";
               break;
            }
            else{
               content += line;
            }
         }
      }
      streamInfile.close();
      if(content.size() != 0){
         content += sequenceSeparator;
         charactersRead++;
         if(charactersRead < _sn - 1){
            _sn = charactersRead;
         }
         _sx_just_to_check += content;
      }

      uint8_t *BWT_from_SA = new uint8_t[_sn];
      int32_t *arrayForBWT = new int32_t[_sn];
      t01 = std::chrono::high_resolution_clock::now();
      libsais_bwt(reinterpret_cast<unsigned char*>(const_cast<char*>(_sx_just_to_check.c_str())), BWT_from_SA, arrayForBWT, _sn, 0, NULL);
      t02 = std::chrono::high_resolution_clock::now();
      std::cerr << "Computing BWT_from_SA took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

      //read BWT_collection from file
      std::ifstream streamInfile2(arg.outname+".bwt",std::ios::in);
      streamInfile2.seekg(0, streamInfile2.end);
      uint64_t length = streamInfile2.tellg();
      streamInfile2.seekg(0, streamInfile2.beg);
      uint8_t *BWT_collection = new uint8_t[length];
      streamInfile2.read((char*)BWT_collection, length);
      streamInfile2.close();
      //print BWT_collection
      // std::cerr << "BWT_collection: ";
      // for(uint32_t i = 0; i < length; i++){
      //    std::cerr << BWT_collection[i];
      // }
      // std::cerr << "\n";
      std::cerr << "Read BWT_collection from file\n";
      //check this BWTCollection with BWT_from_SA
      std::cerr << "Checking BWT_collection with BWT_from_SA\n";
      for(uint32_t i = D-1; i < _sn; i++){
         if(BWT_collection[i] != BWT_from_SA[i]){
            std::cerr << "BWT_collection[" << i << "] = " << BWT_collection[i] << " != " << BWT_from_SA[i] << " = BWT_from_SA[" << i << "]\n";
            //std::cerr << "BWT_collection[" << i << "] = " << (int)BWT_collection[i] << " != " << (int)BWT_from_SA[i] << " = BWT_from_SA[" << i << "]\n";

            //break;
         }
         //else{
         //   std::cerr << "BWT_collection[" << i << "] = " << BWT_collection[i] << " == " << BWT_from_SA[i] << " = BWT_from_SA[" << i << "]\n";
         //}
      }
   }

   return numfactors;
}


void computeBWT(Args arg, char *refFileName, char *collFileName){
   lzInitialize(refFileName, collFileName, arg.prefixLength);
   lzFactorize(arg, collFileName);
   //return MSGSA;
}


void computeBWTMemorySaving(Args arg, char *refFileName, char *collFileName){
   lzInitialize(refFileName, collFileName, arg.prefixLength);
   lzFactorizeMemorySaving(arg, collFileName);
   //return MSGSA;
}
