smc.apt = function(data, #Data matrix
                   grid.points, #Grid points used to compute the predictive density
                   max.resol = 10, #Maximum resolution
                   n.particle = 1000, #Number of particles
                   n.grid.L = 31, #Number of grid points for the L's prior
                   eta.R = 0.1, #Hyperparameter for the L's prior
                   min.obs = 5, #Partitioning is terminated if # observations is less than min.obs
                   lognu.lb = -1, lognu.ub = 4, n.grid.nu = 5, beta = 0.1, #Hyperparameters for APT (See Ma(2017))
                   n.states = 5 #Number of states
){
  if(ncol(data) == ncol(grid.points)){
    d = ncol(data)
    n = nrow(data)
    n_grid = nrow(grid.points)
  }else{
    print("Error: The dimension of the data set and the grid points do not match")
    return(0)
  }

  model_parameters_list = list(beta,
                               lognu.lb,
                               lognu.ub,
                               n.grid.nu)
  names(model_parameters_list) = c("beta","L","U","n_grid_nu")

  #SMC
  out = SMCforPT(data,
                 1,
                 rep(1, n),
                 grid.points,
                 rep(1,n_grid),
                 eta.R,
                 n.states,
                 model_parameters_list,
                 n.particle,
                 max.resol,
                 n.grid.L+1,
                 0.5,
                 0.1,
                 min.obs,
                 1,
                 1
  )

  out_APT = list("MAPtree.left" = out$tree_left,
                 "MAPtree.right" = out$tree_right,
                 "pred.density" = out$pred_density,
                 "posterior.state.nodewise" = out$posterior_states,
                 "children.IDs" = out$children_IDs+1)

  return(out_APT)
}

