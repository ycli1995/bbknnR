
#include <queue>
#include <Rcpp.h>
#include <RcppCommon.h>
#include <RcppEigen.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd correct_idx_cpp(
    Eigen::MatrixXd nnidx, 
    std::vector<double> &raw_idx
) {
  for (int i = 0; i < nnidx.cols(); ++i) {
    for (int j = 0; j < nnidx.rows(); ++j) {
      nnidx(j, i) = raw_idx[nnidx(j, i) - 1];
    }
  }
  return nnidx;
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> trim_graph(
    Eigen::SparseMatrix<double> cnt, 
    int trim
) {
  // Find threshold for each column
  std::vector<double> thres(cnt.cols(), 0.0);
  for (size_t col = 0; col < cnt.outerSize(); ++col) {
    if (cnt.col(col).nonZeros() <= trim) {
      continue;
    }
    double val;
    std::priority_queue<double, std::vector<double>, std::greater<double>> q;
    for (Eigen::SparseMatrix<double>::InnerIterator it(cnt, col); it; ++it){
      if (q.size() < trim) {
        q.push(it.value());
      } else {
        val= it.value();
        if (q.top() < val) {
          q.pop();
          q.push(val);
        }
      }
    }
    thres[col] = q.top();
  }
  // Trim
  for (size_t col = 0; col < cnt.outerSize(); ++col) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(cnt, col); it; ++it){
      if (it.value() < thres[col] || it.value() < thres[it.index()]) {
        it.valueRef() = 0;
      }
    }
  }
  cnt.prune(0.0);
  return cnt;
}
