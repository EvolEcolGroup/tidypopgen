/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// Copy a subset of the raw bytes of a FBM.code256 into another FBM.code256
//
// Column j of `BM_to` is filled with column colInd[j] of `BM_from`, restricted
// (and ordered) to the rows rowInd, and written from row `row_offset` onwards
// (so that two FBMs can be stacked by rows, e.g. by rbind()). The raw bytes are
// recoded with `copy_lut` (a lookup table of 256 values; usually 0:255, so
// that bytes are copied as they are, and the copy works for any code256),
// except for the loci flagged in `swap`, whose bytes are recoded with
// `swap_lut` to swap the reference and alternate alleles.
//
// @param BM_from the FBM.code256 to copy from
// @param BM_to the FBM.code256 to copy to (with length(colInd) columns, and at
//   least row_offset + length(rowInd) rows)
// @param rowInd the (1-based) rows of BM_from to copy
// @param colInd the (1-based) columns of BM_from to copy
// @param col_offset the (0-based) column of BM_to where the copy starts
// @param row_offset the (0-based) row of BM_to where the copy starts
// @param swap a logical vector (one value per element of colInd) flagging
//   the loci where alleles are swapped, or an empty vector (no swapping)
// @param swap_lut a vector of 256 values giving the new raw value of each raw
//   value at swapped loci
// @param copy_lut a vector of 256 values giving the new raw value of each raw
//   value at the other loci (0:255 to copy the bytes as they are)
// [[Rcpp::export]]
void gt_copy_fbm(Environment BM_from,
                 Environment BM_to,
                 const IntegerVector& rowInd,
                 const IntegerVector& colInd,
                 const int col_offset,
                 const int row_offset,
                 const LogicalVector& swap,
                 const IntegerVector& swap_lut,
                 const IntegerVector& copy_lut) {
  XPtr<FBM> xp_from = BM_from["address"];
  XPtr<FBM_RW> xp_to = BM_to["address_rw"];
  const unsigned char* from_ptr =
    static_cast<const unsigned char*>(xp_from->matrix());
  unsigned char* to_ptr = static_cast<unsigned char*>(xp_to->matrix());
  size_t n_from = xp_from->nrow(); // number of rows of the source FBM
  size_t m_from = xp_from->ncol(); // number of columns of the source FBM
  size_t n_to = xp_to->nrow(); // number of rows of the destination FBM
  size_t m_to = xp_to->ncol(); // number of columns of the destination FBM

  size_t n = rowInd.size(); // number of rows to copy
  size_t m = colInd.size(); // number of columns to copy
  bool any_swap = swap.size() > 0;

  // check that all indices are within the FBMs
  if ((row_offset < 0) || ((size_t) row_offset + n > n_to) ||
      (col_offset < 0) || ((size_t) col_offset + m > m_to)) {
    stop("the destination FBM is too small for the copy");
  }
  std::vector<size_t> rows(n);
  for (size_t i = 0; i < n; i++) {
    if ((rowInd[i] < 1) || ((size_t) rowInd[i] > n_from)) {
      stop("'rowInd' is out of the bounds of the source FBM");
    }
    rows[i] = rowInd[i] - 1;
  }
  for (size_t j = 0; j < m; j++) {
    if ((colInd[j] < 1) || ((size_t) colInd[j] > m_from)) {
      stop("'colInd' is out of the bounds of the source FBM");
    }
  }
  if (any_swap && ((size_t) swap.size() != m)) {
    stop("'swap' should have one value per column");
  }
  if ((swap_lut.size() != 256) || (copy_lut.size() != 256)) {
    stop("'swap_lut' and 'copy_lut' should have 256 elements");
  }
  for (int b = 0; b < 256; b++) {
    if ((swap_lut[b] < 0) || (swap_lut[b] > 255) ||
        (copy_lut[b] < 0) || (copy_lut[b] > 255)) {
      stop("'swap_lut' and 'copy_lut' should only contain values from 0 to 255");
    }
  }
  // lookup tables of raw values: at loci that are not swapped (usually
  // unchanged), and with swapped alleles
  unsigned char same[256], swapped[256];
  for (int b = 0; b < 256; b++) {
    same[b] = (unsigned char) copy_lut[b];
    swapped[b] = (unsigned char) swap_lut[b];
  }

  for (size_t j = 0; j < m; j++) {
    const unsigned char* col_from = from_ptr + n_from * (colInd[j] - 1);
    unsigned char* col_to = to_ptr + n_to * (col_offset + j) + row_offset;
    const unsigned char* lut =
      (any_swap && (swap[j] == TRUE)) ? swapped : same;
    for (size_t i = 0; i < n; i++) {
      col_to[i] = lut[col_from[rows[i]]];
    }
  }
}

/******************************************************************************/
