#ifndef MRS_FUNCTIONS
#define MRS_FUNCTIONS

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "helpers.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


arma::mat CreateXi_MRS(const int& I, const List& model_parameters_list, const int& k){
  
  double beta = model_parameters_list["beta"];
  double gamma = model_parameters_list["gamma"];
  double eta_prune = model_parameters_list["eta_prune"];

  arma::mat xi = zeros(I,I);
  
  xi(0,0) = (1 - eta_prune) * gamma * pow(beta, (double)-k);
  xi(0,1) = (1 - eta_prune) * (1 - gamma * pow(beta, (double)-k));
  xi(0,2) = eta_prune;
  
  xi(1,0) = (1 - eta_prune) * gamma * pow(2.0, (double)-k);
  xi(1,1) = (1 - eta_prune) * (1 - gamma * pow(2.0, (double)-k));
  xi(1,2) = eta_prune;
  
  xi(2,2) = 1;
  
  return xi;
}

double precision_theta_MRS_g;

void pre_compute_MRS(const int& G, const int& I, const int& NL, const List& model_parameters_list){
  
  precision_theta_MRS_g = model_parameters_list["precision_theta"];
    
}

double log_ML_compute_MRS(const int& G, const int& I, const arma::ivec& n_l_vec,
                      const arma::ivec& n_r_vec,const double& L_input, const int& NL, const int& V, const int& k){
  
  double out = 0.0;
  int n_l, n_r;
  
  double alpha_l = precision_theta_MRS_g * L_input;
  double alpha_r = precision_theta_MRS_g * (1-L_input);
  
  if(V == 1){
    
    for(int g=0;g<G;g++){
      n_l = n_l_vec(g);
      n_r = n_r_vec(g);
      
      out += lgamma(alpha_l+n_l) + lgamma(alpha_r+n_r) - lgamma(alpha_l+n_l+alpha_r+n_r) -
                (lgamma(alpha_l) + lgamma(alpha_r) - lgamma(alpha_l + alpha_r));
    }
    
  }else{
    
    n_l = sum(n_l_vec);
    n_r = sum(n_r_vec);
    
    out = lgamma(alpha_l+n_l) + lgamma(alpha_r+n_r) - lgamma(alpha_l+n_l+alpha_r+n_r) -
               (lgamma(alpha_l) + lgamma(alpha_r) - lgamma(alpha_l + alpha_r));
  }
  
  return out;
  
}


double PostMean_MRS(const int& I, const ivec& n_l_vec, const ivec& n_r_vec, const double& L_input, const int& V, const int& g, const int& k){
  
  //In the case of MRS, we do nothing here.
  return 0;
}

#endif