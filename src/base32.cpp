#include <Rcpp.h>
using namespace Rcpp;

static const char alphabet62[] =
  "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

// ------------------------------------------------
// base-62 encode uint64
// ------------------------------------------------
inline std::string encode62_uint64(uint64_t x) {
  if (x == 0) return "0";
  
  char buf[12]; // enough for uint64 in base-62
  int i = 11;
  buf[i] = '\0';
  
  while (x > 0) {
    buf[--i] = alphabet62[x % 62];
    x /= 62;
  }
  return std::string(buf + i);
}

// ------------------------------------------------
// base-62 decode uint64
// ------------------------------------------------
inline uint64_t decode62_uint64(const std::string& s) {
  uint64_t x = 0;
  for (char c : s) {
    uint64_t v;
    if (c >= '0' && c <= '9')       v = c - '0';
    else if (c >= 'A' && c <= 'Z')  v = c - 'A' + 10;
    else if (c >= 'a' && c <= 'z')  v = c - 'a' + 36;
    else stop("Invalid base-62 character");
    
    x = x * 62 + v;
  }
  return x;
}

// ------------------------------------------------
// scalar encode (a, b)
// ------------------------------------------------
// [[Rcpp::export]]
std::string encode_pair_cpp(uint64_t a, uint64_t b, int width = 2) {
  if (width < 1)
    stop("width must be >= 1");
  
  uint64_t max_a = 1;
  for (int i = 0; i < width; ++i)
    max_a *= 62;
  max_a -= 1;
  
  if (a > max_a)
    stop("a exceeds maximum representable value for this width");
  
  std::string A = encode62_uint64(a);
  if ((int)A.size() < width)
    A.insert(0, width - A.size(), '0');
  
  return encode62_uint64(b) + A;
}

// ------------------------------------------------
// scalar decode → numeric vector c(a, b)
// ------------------------------------------------
// [[Rcpp::export]]
NumericVector decode_pair_cpp(const std::string& s, int width = 2) {
  if (width < 1)
    stop("width must be >= 1");
  
  int n = s.size();
  if (n <= width)
    stop("string too short for given width");
  
  uint64_t b = decode62_uint64(s.substr(0, n - width));
  uint64_t a = decode62_uint64(s.substr(n - width, width));
  
  NumericVector out(2);
  out[0] = static_cast<double>(a);
  out[1] = static_cast<double>(b);
  return out;
}

// ------------------------------------------------
// vectorised encoder
// ------------------------------------------------
// [[Rcpp::export]]
CharacterVector encode_pair_vec_cpp(NumericVector a,
                                    NumericVector b,
                                    int width = 2) {
  int n = a.size();
  if (b.size() != n)
    stop("a and b must have the same length");
  
  CharacterVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = encode_pair_cpp(
      static_cast<uint64_t>(a[i]),
      static_cast<uint64_t>(b[i]),
      width
    );
  }
  return out;
}

// ------------------------------------------------
// vectorised decoder
// ------------------------------------------------
// [[Rcpp::export]]
NumericMatrix decode_pair_vec_cpp(CharacterVector s,
                                  int width = 2) {
  int n = s.size();
  NumericMatrix out(n, 2);
  
  for (int i = 0; i < n; ++i) {
    NumericVector tmp = decode_pair_cpp(
      Rcpp::as<std::string>(s[i]),
      width
    );
    out(i, 0) = tmp[0]; // a
    out(i, 1) = tmp[1]; // b
  }
  
  colnames(out) = CharacterVector::create("chr", "pos");
  return out;
}
