/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <fstream>
//#include <bigsnpr/bed-acc.h>

using namespace Rcpp;
using namespace std;

/******************************************************************************/

#define PACK_DENSITY 4

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
  int lengthExtra = length + (extra > 0); // total number of bites with information
  int j, k;
  
  unsigned char *buffer = new unsigned char[std::max(3, lengthExtra) + 1];
  
  ifstream myFile(filename, ios::in | ios::binary);
  
  // size checks
  myFile.seekg(0, std::ifstream::end);
  // file size in bytes
  long long len = (long long)myFile.tellg();
  // byets per snp (there are additional empty bytes at the end of lines)
  long long bytespersnp = len/(m+1);

  // skip the first row (which is metadata on n_ind and n_snps)
  myFile.seekg((0+1)*bytespersnp, std::ifstream::beg);
  
  // buffer to read in
  unsigned char* tmp = new unsigned char[bytespersnp];
  char tmpi;
  
  for (j = 0; j < m; j++) {
    // read from bedfile
    //myFile.read((char*)buffer, lengthExtra);
    Rcout<<"locus "<< j << endl;
    // read raw genotypes
    myFile.read((char*)tmp, sizeof(char) * bytespersnp);
    
    for (k = 0; k < length; k++) {
      code_ptr = &tab(0, tmp[k]);
      ptr = std::copy(code_ptr, code_ptr + 4, ptr);
      // Rcout<< "locus "<< j << "\t" << k << "\t" << buffer[k] << endl;
      // Rcout<< "locus "<< j << "\t" << k << "\t" << (int)buffer[k] << endl;
      
      tmpi = tmp[k];
      Rcout<< "\t" << ((tmpi & MASK3) >> 6);
      Rcout<< "\t" << ((tmpi & MASK2) >> 4);
      Rcout<< "\t" << ((tmpi & MASK1) >> 2);
      Rcout<< "\t" << (tmpi & MASK0);
    }
    if (extra){
      code_ptr = &tab(0, tmp[k]);
      ptr = std::copy(code_ptr, code_ptr + extra, ptr);
      
      tmpi = tmp[k];
      Rcout<< "\t" << ((tmpi & MASK3) >> 6);
      Rcout<< "\t" << ((tmpi & MASK2) >> 4);
      Rcout<< "\t" << ((tmpi & MASK1) >> 2);
      Rcout<< "\t" << (tmpi & MASK0);
    }
    Rcout<<endl;
  }
   
    // for (k = 0; k < length; k++) {
    //   code_ptr = &tab(0, buffer[k]);
    //   ptr = std::copy(code_ptr, code_ptr + 4, ptr);
    // }
    // if (extra) {
    //   code_ptr = &tab(0, buffer[k]);
    //   ptr = std::copy(code_ptr, code_ptr + extra, ptr);
    // }

  
  char c;
//  bool is_eof = !(myFile.get(c));
  myFile.close();
  delete[] buffer;
  delete[] tmp;
  return true;
}
