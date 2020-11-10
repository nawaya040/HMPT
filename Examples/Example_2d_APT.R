#Generate the data
set.seed(6)
n = 1000
data = matrix(NA, nrow=n, ncol=2)

w = c(1/3, 1/3, 1/3)
b1.x = c(0.1, 0.45)
b1.y = c(0.35, 0.9)
b2.x = c(0.2, 0.8)
b2.y = c(0.45, 0.5)
b3.x = c(0.7, 0.9)
b3.y = c(0.05, 0.6)

for(i in 1:n){
  u = runif(1)
  if(u < w[1]){
    x.temp = b1.x[1] + runif(1) * (b1.x[2] - b1.x[1])
    y.temp = b1.y[1] + runif(1) * (b1.y[2] - b1.y[1])
  }else if(u < w[1]+w[2]){
    x.temp = b2.x[1] + runif(1) * (b2.x[2] - b2.x[1])
    y.temp = b2.y[1] + runif(1) * (b2.y[2] - b2.y[1])
  }else{
    x.temp = b3.x[1] + runif(1) * (b3.x[2] - b3.x[1])
    y.temp = b3.y[1] + runif(1) * (b3.y[2] - b3.y[1])
  }

  data[i,1] = x.temp
  data[i,2] = y.temp
}

plot(data[,1], data[,2], xlab = "x1", ylab = "x2")

#Make the grid points
x.grid = seq(0.005,0.995,by=0.01)
y.grid = seq(0.005,0.995,by=0.01)
grid.points = as.matrix(expand.grid(x.grid,y.grid))

#SMC
out_APT = smc.apt(data, grid.points, max.resol = 8, n.particle = 1000, n.grid.L = 31, eta.R = 0.01, n.states = 5)

#Plot the estimated density
library(ggplot2)
library(viridis)
densities.df = data.frame(x1=grid.points[,1],x2=grid.points[,2],density = out_APT$pred.density)
p = ggplot() + geom_tile(data = densities.df, aes(x=x1, y=x2, fill=density)) + scale_fill_viridis(discrete=FALSE)
print(p)

#Plot the MAP tree
tree_left = out_APT$MAPtree.left
tree_right = out_APT$MAPtree.right

n_all_node = ncol(tree_left)
N = (n_all_node-1)/2

x_all = numeric(N)
y_all = numeric(N)
x_end_all = numeric(N)
y_end_all = numeric(N)

for(j in 1:N){
  start = tree_left[,2*j+1]
  end   = tree_right[,2*j]

  x_all[j] = start[1]
  y_all[j] = start[2]

  x_end_all[j] = end[1]
  y_end_all[j] = end[2]
}

x_st = c(0,0,1,1)
x_end = c(0,1,1,0)
y_st = c(0,1,1,0)
y_end = c(1,1,0,0)

print(p + geom_segment(aes(x = x_all, y = y_all, xend = x_end_all, yend = y_end_all)) +
  geom_segment(aes(x = x_st, y = y_st, xend = x_end, yend = y_end)) + xlab("x1") + ylab("x2"))
