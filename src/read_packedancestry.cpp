/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <fstream>
//#include <bigsnpr/bed-acc.h>

using namespace Rcpp;
using namespace std;

/******************************************************************************/


/* There are four individuals, coded as groups of 2 bits, in each byte */
/* we ave use the mask below to extract each group */
/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */

// [[Rcpp::export]]
bool read_packedancestry(const char * filename,
              Environment BM,
              const RawMatrix& tab) {
  
  XPtr<FBM_RW> xpBM = BM["address_rw"];
  unsigned char* ptr = static_cast<unsigned char*>(xpBM->matrix());
  const unsigned char* code_ptr;
  int n = xpBM->nrow(); // number of individuals
  int m = xpBM->ncol(); // number of loci
  
  int length = n / 4; // number of complete packed bytes
  int extra = n - 4 * length; //extra indivs in the last byte
  int j, k;
  
  ifstream myFile(filename, ios::in | ios::binary);
  
  // size checks
  myFile.seekg(0, std::ifstream::end);
  // file size in bytes
  long long len = (long long)myFile.tellg();
  // bytes per snp (there are additional empty bytes at the end of lines)
  long long bytespersnp = len/(m+1);

  // skip the first row (which is metadata on n_ind and n_snps)
  myFile.seekg((0+1)*bytespersnp, std::ifstream::beg);
  
  // buffer to read in a snp
  unsigned char* tmp = new unsigned char[bytespersnp];

  for (j = 0; j < m; j++) {
    // read raw genotypes
    myFile.read((char*)tmp, sizeof(char) * bytespersnp);
    
    for (k = 0; k < length; k++) {
      code_ptr = &tab(0, tmp[k]);
      ptr = std::copy(code_ptr, code_ptr + 4, ptr);

      // debug
      // note that we start reading from the left (i.e. higher bits in the byte)
      // char tmpi = tmp[k];
      // Rcout<< "\t" << ((tmpi & MASK3) >> 6); // shift the bytes to the right
      // Rcout<< "\t" << ((tmpi & MASK2) >> 4);
      // Rcout<< "\t" << ((tmpi & MASK1) >> 2);
      // Rcout<< "\t" << (tmpi & MASK0);
    }
    if (extra){
      code_ptr = &tab(0, tmp[k]);
      ptr = std::copy(code_ptr, code_ptr + extra, ptr);
      
      // debug
      // char tmpi = tmp[k];
      // Rcout<< "\t" << ((tmpi & MASK3) >> 6);
      // Rcout<< "\t" << ((tmpi & MASK2) >> 4);
      // Rcout<< "\t" << ((tmpi & MASK1) >> 2);
      // Rcout<< "\t" << (tmpi & MASK0);
    }
  }

  char c;
  bool is_eof = !(myFile.get(c));
  myFile.close();
  delete[] tmp;
  return is_eof;
}
