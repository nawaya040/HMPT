#Generate the data
set.seed(0)
n=400
data = matrix(NA, nrow=n, ncol=2)
groups = rep(c(1,2), each=n/2)

w = c(1/3, 1/3, 1/3)
cum_w = cumsum(w)

means1 = c(-2.5,1,2)
means2 = c(1,-2,2.5)

for(i in 1:n){

  group.i = groups[i]
  u = runif(1)

  if(u < cum_w[1]){

    data[i,1] = rnorm(1, means1[1], 0.5)
    data[i,2] = rnorm(1, means2[1], 0.7)
  }else if(u < cum_w[2]){
    data[i,1] = rnorm(1, means1[2], 0.5)
    data[i,2] = rnorm(1, means2[2], 0.7)
  }else{
    if(group.i == 2){
      data[i,1] = rnorm(1, means1[3]-0.5, 0.5)
      data[i,2] = rnorm(1, means2[3]-0.5, 0.7)
    }else{
      data[i,1] = rnorm(1, means1[3], 0.5)
      data[i,2] = rnorm(1, means2[3], 0.7)
    }
  }
}

#SMC
out =  smc.mrs(data, groups, max.resol = 8, n.particle = 1000, n.grid.L = 31, eta.R = 0.1)

#Visualize a node with significant difference
data1 = data[groups==1,]
data2 = data[groups==2,]

node.ID = which.min(out$post.null.nodewise)
#node.ID = which.max(out$eff.nodewise)

left  = out$MAPtree.l[,node.ID]
right = out$MAPtree.r[,node.ID]

x_st = c(left[1],left[1],right[1],right[1])
x_end = c(left[1],right[1],right[1],left[1])
y_st = c(left[2],right[2],right[2],left[2])
y_end = c(right[2],right[2],left[2],left[2])

plot(data1[,1], data1[,2], xlab = "x1", ylab = "x2", pch=21, col="blue", main = paste("P(H0|Data)=",round(out$post.null.global, digits=3),sep=""))
points(data2[,1], data2[,2], pch=24, col="red")
segments(x_st, y_st, x_end, y_end, lwd = 2)
