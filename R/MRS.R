smc.mrs = function(data, #Data matrix
                         groups, #Vector of groups numbers
                         max.resol = 10, #Maximum resolution
                         n.particle = 1000, #Number of particles
                         n.grid.L = 31, #Number of grid points for the L's prior
                         eta.R = 0.1, #Hyperparameter for the L's prior
                         min.obs = 5, #Partitioning is terminated if # observations is less than min.obs
                         beta = 1.0, gamma = 0.3, eta.prune = 0.3, #Hyperparameters for MRS (See Soriano and Ma(2017))
                         prec.theta = 1, #Precision parameter for theta
                         rescale = 0 #If=1, all the variables are automatically rescaled
){
  if(nrow(data) == length(groups)){
    d = ncol(data)
    n = nrow(data)
  }else{
    print("Error: The size of the data set and the group vector do not match")
    return(0)
  }

  #Rescale the data if necessary
  range.store = matrix(NA, nrow=d, ncol=2)
  for(j in 1:d){
    range_j = range(data[,j])

    if((range_j[1] < 0) | (1 < range_j[2]) | (rescale == 1)){
      min0 = range_j[1]
      max0 = range_j[2]

      data[,j] = (data[,j] - min0) / (max0 - min0) * 0.999 + 0.0005

      range.store[j,] = c(min0-0.0005, max0+0.0005)
    }else{
      range.store[j,] = c(0, 1)
    }
  }

  model_parameters_list = list(beta,
                               gamma,
                               eta.prune,
                               prec.theta
                              )

  names(model_parameters_list) = c("beta", "gamma", "eta_prune","precision_theta")

  #SMC
  out = SMCforPT(data,
                 2,
                 groups,
                 matrix(0.5,nrow=1,ncol=d),
                 c(1),
                 eta.R,
                 3,
                 model_parameters_list,
                 n.particle,
                 max.resol,
                 n.grid.L+1,
                 0.5,
                 0.1,
                 min.obs,
                 0,
                 2
  )

  tree_left = out$tree_left
  tree_right = out$tree_right

  n_nodes = ncol(tree_left)
  for(i in 1:n_nodes){
    tree_left[,i]  = range.store[,1] + tree_left[,i] * (range.store[,2] - range.store[,1])
    tree_right[,i] = range.store[,1] + tree_right[,i] * (range.store[,2] - range.store[,1])
  }

  #Compute the posterior probability of the null hypothesis
  children_IDs= out$children_IDs
  xi_post = out$posterior_xi
  psi_tilde_store = numeric(n_nodes)

  for(i in n_nodes:2){
    #Check if the node is a leaf or not
    Is_leaf = (children_IDs[1,i] == -1)

    #Leaf
    if(Is_leaf == TRUE){
      psi_tilde_store[i] = xi_post[2,2,i] + xi_post[2,3,i]
    }

    #Non-leaf
    if(Is_leaf == FALSE){
      child_ID_l = children_IDs[1,i] + 1 #Don't forget to add 1!
      child_ID_r = children_IDs[2,i] + 1

      psi_tilde_child_l = psi_tilde_store[child_ID_l]
      psi_tilde_child_r = psi_tilde_store[child_ID_r]

      psi_tilde_store[i] = xi_post[2,2,i] * psi_tilde_child_l * psi_tilde_child_r + xi_post[2,3,i]
    }
  }

  psi_tilde_store[1] = xi_post[2,2,1] * psi_tilde_store[2] * psi_tilde_store[3] + xi_post[2,3,1]

  #Compute the posterior probability of the two groups being equal on each node
  post.null.nodewise = out$posterior_states[2,] + out$posterior_states[3,]

  out_MRS = list("MAPtree.left" = tree_left,
                 "MAPtree.right" = tree_right,
                 "post.null.global" = out$post_null,
                 "post.null.nodewise" = post.null.nodewise,
                 "eff.nodewise" = out$eff_MAP,
                 "children.IDs" = out$children_IDs+1)

  return(out_MRS)
}
