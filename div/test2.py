import fastperm
import numpy as np
import rpy2.robjects as robjects



data = robjects.r("""
n = 100

set.seed(10002)

ts = seq(-1,1, length.out = n)

U = runif(n)

I =  (U <= (1/3 * pnorm((ts+2)/0.7) + 2/3 * pnorm((ts-1)/0.7)))

xleft = ts[I]
xright = ts[!I]

nleft = length(xleft)
nright = length(xright)

a = c(rep(-Inf, nleft), xright)
b = c(xleft, rep(Inf, nright))


T = 100

alpha_param = 100

X = matrix(NA, nrow = T, ncol = n)
# generate samples from prior
for (t in 1:T) {
  X[t,1] = rnorm(1)
  Us = runif(n-1)
  for (i in 2:n) {
    if(Us[i-1] <= alpha_param / (alpha_param +i-1)){
      X[t,i] = rnorm(1)
    }else{
      if(i==2){
        X[t,i] = X[t,i-1]
      }else{
        X[t,i] = sample(x = X[t,1:(i-1)],1)
      }
      if(X[t,i]==1){
        print("her")
      }
    }
    
  }
  X[t,] = sort(X[t,])
}
X = t(X)

""")


#ts = np.array(robjects.r["ts"])
X = np.array(robjects.r["X"],dtype=np.float64,order='F')
a = np.array(robjects.r["a"],dtype=np.float64)
b = np.array(robjects.r["b"],dtype=np.float64)
T = int(robjects.r["T"][0])
n = int(robjects.r["n"][0])

# print(type(T))
# print(type(n))
# print(type(X))
# print(type(a))
# print(type(b))


# print(X[1,:])

res = fastperm.get_log_permanents(X,a,b,n,T,0)
print(res)
