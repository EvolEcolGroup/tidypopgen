/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>
#include <fstream>
//#include <bigsnpr/bed-acc.h>

using namespace Rcpp;
using namespace std;

/******************************************************************************/

#define PACK_DENSITY 4

// [[Rcpp::export]]
bool read_packedancestry(const char * filename,
              Environment BM,
              const RawMatrix& tab) {
  
  XPtr<FBM_RW> xpBM = BM["address_rw"];
  unsigned char* ptr = static_cast<unsigned char*>(xpBM->matrix());
  const unsigned char* code_ptr;
  int n = xpBM->nrow(); // number of individuals
  int m = xpBM->ncol(); // number of loci
  
  int length = n / 4;
  int extra = n - 4 * length;
  int lengthExtra = length + (extra > 0);
  int j, k;
  
  unsigned char *buffer = new unsigned char[std::max(3, lengthExtra) + 1];
  
  ifstream myFile(filename, ios::in | ios::binary);
  
  // size checks
  myFile.seekg(0, std::ifstream::end);
  // file size in bytes
  long long len = (long long)myFile.tellg();
  Rcout<<"len is "<<len << endl;
  long long bytespersnp = len/(m+1);
  Rcout<<" bytespersnp is "<<bytespersnp<<endl;
  
  

  for (j = 0; j < m; j++) {
    // read from bedfile
    myFile.read((char*)buffer, lengthExtra);
    
    for (k = 0; k < length; k++) {
      code_ptr = &tab(0, buffer[k]);
      ptr = std::copy(code_ptr, code_ptr + 4, ptr);
    }
    if (extra) {
      code_ptr = &tab(0, buffer[k]);
      ptr = std::copy(code_ptr, code_ptr + extra, ptr);
    }
  }
  
  char c;
  bool is_eof = !(myFile.get(c));
  myFile.close();
  delete[] buffer;
  
  return is_eof;
}
