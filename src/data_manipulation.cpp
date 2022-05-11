
#include <queue>
#include <Rcpp.h>
#include <RcppCommon.h>

using namespace Rcpp;
using namespace std;

template <int RTYPE>
Matrix<RTYPE> UpdateMat(
    Matrix<RTYPE> raw_mat,
    const Matrix<RTYPE> new_mat,
    const IntegerVector row_idx,
    const IntegerVector col_idx
) {
  int ncol = col_idx.length();
  int nrow = row_idx.length();
  for (int c = 0; c < ncol; ++c) {
    for (int r = 0; r < nrow; ++r) {
      raw_mat(row_idx[r], col_idx(c)) = new_mat(r, c);
    }
  }
  return Matrix<RTYPE>(raw_mat);
}

// [[Rcpp::export]]
NumericMatrix UpdateNumMat(
    NumericMatrix & raw_mat,
    NumericMatrix & new_mat,
    IntegerVector & row_idx,
    IntegerVector & col_idx
) {
  return UpdateMat(raw_mat, new_mat, row_idx, col_idx);
}

// [[Rcpp::export]]
IntegerMatrix UpdateIntMat(
    IntegerMatrix & raw_mat,
    IntegerMatrix & new_mat,
    IntegerVector & row_idx,
    IntegerVector & col_idx
) {
  return UpdateMat(raw_mat, new_mat, row_idx, col_idx);
}

// [[Rcpp::export]]
IntegerMatrix GetRawIndex(
    IntegerMatrix & idx_mat,
    IntegerVector & raw_idx
) {
  int ncol = idx_mat.ncol();
  int nrow = idx_mat.nrow();
  for (int c = 0; c < ncol; ++c) {
    for (int r = 0; r < nrow; ++r) {
      idx_mat(r, c) = raw_idx[idx_mat(r, c) - 1];
    }
  }
  return idx_mat;
}

// [[Rcpp::export]]
List GetSparseDist(
  IntegerMatrix & knn_index,
  NumericMatrix & knn_dist,
  int n_obs,
  int n_neighbors
) {
  IntegerVector rows(n_obs * n_neighbors);
  IntegerVector cols(n_obs * n_neighbors);
  NumericVector vals(n_obs * n_neighbors);
  
  int nrow = knn_index.nrow();
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < n_neighbors; ++j) {
      if (knn_index(i, j) == 0) continue;
      if (knn_index(i, j) == i + 1) {
        vals[i * n_neighbors + j] = 0;
      } else {
        vals[i * n_neighbors + j] = knn_dist(i, j);
      }
      cols[i * n_neighbors + j] = i;
      rows[i * n_neighbors + j] = knn_index(i, j) - 1;
    }
  }
  List res = List::create(_["i"] = rows, _["j"] = cols, _["x"] = vals);
  return res;
}

// [[Rcpp::export]]
void Trimming(
  NumericVector & x,
  IntegerVector & row_idx,
  IntegerVector & p,
  int trim
) {
  int ncol = p.length() - 1;
  
  // Find threshold for each column
  NumericVector thres(ncol);
  double val;
  for (size_t col = 0; col < ncol; ++col) {
    priority_queue< double, vector<double>, greater<double> > q;
    for (size_t i = p[col]; i < p[col + 1]; ++i) {
      if (q.size() < trim) {
        q.push(x[i]);
      } else {
        val = x[i];
        if (q.top() < val) {
          q.pop();
          q.push(val);
        }
      }
    }
    thres[col] = q.top();
  }
  
  // Trim
  int row;
  double val_col, val_row;
  for (size_t col = 0; col < ncol; ++col) {
    val_col = thres[col];
    for (size_t i = p[col]; i < p[col + 1]; ++i) {
      row = row_idx[i];
      val_row = thres[row];
      if (x[i] < val_col || x[i] < val_row) {
        x[i] = 0;
      }
    }
  }
}
