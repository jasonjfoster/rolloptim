// todo (roll >= 1.2.1): use namespace roll except for min_dimnames

#ifndef ROLLOPTIM_ATTR_H
#define ROLLOPTIM_ATTR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

namespace rolloptim {

// xts attributes: copy dimnames, index, .indexCLASS, .indexTZ, tclass, tzone, class
template <typename T, typename S>
inline void xts_attr(T& target, const S& source) {

  target.attr("dimnames") = source.attr("dimnames");
  target.attr("index") = source.attr("index");
  target.attr(".indexCLASS") = source.attr(".indexCLASS");
  target.attr(".indexTZ") = source.attr(".indexTZ");
  target.attr("tclass") = source.attr("tclass");
  target.attr("tzone") = source.attr("tzone");
  target.attr("class") = source.attr("class");

}

// min dimnames: conditionally set dimnames with R_NilValue rows and source column names
template <typename T, typename S>
inline void min_dimnames(T& target, const S& source) {

  List dimnames = source.attr("dimnames");

  if (dimnames.size() > 1) {
    target.attr("dimnames") = List::create(R_NilValue, dimnames[1]);
  }

}

}

#endif