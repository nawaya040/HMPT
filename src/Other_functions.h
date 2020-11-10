#ifndef OTHER_FUNCTIONS
#define OTHER_FUNCTIONS

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "helpers.h"

using namespace Rcpp;
using namespace arma;
using namespace std;



//Make the transition matrix.
arma::mat CreateXi_OTHER(const int& I, const List& model_parameters_list, const int& k){
  arma::mat xi = zeros(I,I);
  return xi;
}

void pre_compute_OTHER(const int& G, const int& I, const int& NL, const List& model_parameters_list){
    
}

//Compute the marginal likelihood
double log_ML_compute_OTHER(const int& G, const int& I, const arma::ivec& n_l_vec,
                      const arma::ivec& n_r_vec,const double& L_input, const int& NL, const int& V, const int& k){
  return 0.0;
  
}


//Compute the posterior mean
double PostMean_OTHER(const int& I, const ivec& n_l_vec, const ivec& n_r_vec, const double& L_input, const int& V, const int& g, const int& k){
  return 0;
}

#endif