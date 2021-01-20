#ifndef APT_FUNCTIONS
#define APT_FUNCTIONS

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "helpers.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

double log_ML_compute_core_APT(const int& G, const int& I, const arma::ivec& n_l_vec,const arma::ivec& n_r_vec,const double& L_input, const int& V);

arma::mat CreateXi_APT(const int& I, const List& model_parameters_list, const int& k){
  
  double beta = model_parameters_list["beta"];

  arma::mat xi = zeros(I,I);
  for(int i=0;i<I;i++){
    for(int j=i;j<I;j++){
      xi(i,j) = exp(beta * (double)(i-j));
    }
    xi.row(i) = xi.row(i) /sum(xi.row(i));
  }
  return xi;
}

//# grid points for nu (precision)
int n_grid_nu_APT_g;

//Matrix to store the precision at the grid points.
//Notice that this matrix is (# grid points) by (I-1).
arma::mat precision_mat_APT_g;

//Cube to store results of pre-computing the marginal likelihoods.
arma::cube ML_pre_cube_APT_g;
//Size of pre-computing
int size_pre_APT_g = 100;

void pre_compute_APT(const int& G, const int& I, const int& NL, const List& model_parameters_list){
  
  double L = model_parameters_list["L"];
  double U = model_parameters_list["U"];
  n_grid_nu_APT_g = model_parameters_list["n_grid_nu"];
  
  //Make a matrix to store all possible precision
  precision_mat_APT_g = zeros(n_grid_nu_APT_g, I-1);
  
  arma::vec grid_log_10 = linspace(L,U,(I-1)+1);
  
  for(int i=0;i<I-1;i++){
    double dleft = grid_log_10(i) + 1.0 / (double) (2 * n_grid_nu_APT_g);
    double dright = grid_log_10(i+1) - 1.0 / (double) (2 * n_grid_nu_APT_g);
    
    arma::vec log_10_i = linspace(dleft, dright, n_grid_nu_APT_g);
    arma::vec prec_i(n_grid_nu_APT_g);
    
    for(int j=0;j<n_grid_nu_APT_g;j++){
      prec_i(j) = pow(10.0, log_10_i(j));
    }
    
    precision_mat_APT_g.col(i) = prec_i;
  }
  
  //Make a cobe that stores results of the pre-computing.
  int n_location_options = NL-1;
  ML_pre_cube_APT_g.set_size(size_pre_APT_g, size_pre_APT_g, n_location_options*I);
  ivec n_l_vec(1);
  ivec n_r_vec(1);
  int index_slice;
  double L_input;
  
  for(int V=1;V<I+1;V++){
    for(int l=0;l<n_location_options;l++){
      index_slice = (V-1)*n_location_options + l;
      L_input = (double) (l+1) / (double) NL;
      
      for(int i=0;i<size_pre_APT_g;i++){
        
        n_l_vec(0) = i;
        
        for(int j=0;j<size_pre_APT_g;j++){
          
          n_r_vec(0) = j;
          ML_pre_cube_APT_g(i, j, index_slice) = log_ML_compute_core_APT(G, I, n_l_vec, n_r_vec, L_input, V);
        
        }
      }
      
    }
  }
}

double log_ML_compute_APT(const int& G, const int& I, const arma::ivec& n_l_vec,
                      const arma::ivec& n_r_vec,const double& L_input, const int& NL, const int& V, const int& k){
  
  double out = 0.0;
  
  int n_l = n_l_vec(0);
  int n_r = n_r_vec(0);
  
  //If the sample size is small, we can use the result stored in the cube.
  if((n_l < size_pre_APT_g) && (n_r < size_pre_APT_g)){
    int L_index = round(L_input * NL) - 1;
    int n_location_options = NL-1;
    int index_slice = (V-1)*n_location_options + L_index;
    out = ML_pre_cube_APT_g(n_l, n_r, index_slice);
  }else{
    out = log_ML_compute_core_APT(G,I,n_l_vec,n_r_vec,L_input,V);
  }
  
  return out;
}


double PostMean_APT(const int& I, const ivec& n_l_vec, const ivec& n_r_vec, const double& L_input, const int& V, const int& g, const int& k){
  
  double post_mean;
  int n_l = n_l_vec(g);
  int n_r = n_r_vec(g);
  
  if(V < I)
  {
    vec log_ML_temp(n_grid_nu_APT_g);
    vec post_means_temp(n_grid_nu_APT_g);
    
    for(int j=0;j<n_grid_nu_APT_g;j++){
      double alpha_l = precision_mat_APT_g(j,V-1) * L_input;
      double alpha_r = precision_mat_APT_g(j,V-1) * (1-L_input);
      
      log_ML_temp(j) = lgamma(alpha_l+n_l) + lgamma(alpha_r+n_r) - lgamma(alpha_l+n_l+alpha_r+n_r) -
        (lgamma(alpha_l) + lgamma(alpha_r) - lgamma(alpha_l + alpha_r));
      
      post_means_temp(j) = (alpha_l+n_l) / (alpha_l+n_l + alpha_r+n_r);
    }

    vec post_nu = Normalize_log(log_ML_temp);
    post_mean = sum(post_means_temp % post_nu);
    
  }else{
    post_mean = L_input;
  }
  
  return post_mean;
}


double log_ML_compute_core_APT(const int& G, const int& I, const arma::ivec& n_l_vec,const arma::ivec& n_r_vec,const double& L_input, const int& V){
  
  double out = 0.0;
  double alpha_l, alpha_r;
  
  if(V < I)
  {
    vec log_ML_temp(n_grid_nu_APT_g);
    
    for(int g=0;g<G;g++){
      
      int n_l = n_l_vec(g);
      int n_r = n_r_vec(g);
      
      for(int j=0;j<n_grid_nu_APT_g;j++){
        
        alpha_l = precision_mat_APT_g(j,V-1) * L_input;
        alpha_r = precision_mat_APT_g(j,V-1) * (1-L_input);
        
        log_ML_temp(j) = lgamma(alpha_l+n_l) + lgamma(alpha_r+n_r) - lgamma(alpha_l+n_l+alpha_r+n_r) -
          (lgamma(alpha_l) + lgamma(alpha_r) - lgamma(alpha_l + alpha_r));
        
      }
      
      out += ColSum_log(log(1.0 / (double) n_grid_nu_APT_g) + log_ML_temp);
    }
    
  }else{
    out	= sum(n_l_vec) * log(L_input) + sum(n_r_vec) * log(1-L_input);
  }
  return out;
  
}


#endif