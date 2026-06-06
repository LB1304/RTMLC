Draw.RTMLC <- function(H, Enh, J, Jsel, nc1, nc2, 
                       sup1resp, sup2resp, sup1time, sup2time, 
                       dif, speed0, speed1, si2, 
                       coef1resp, coef2resp, coef1time, coef2time) {
  
  nv = pmax(rpois(H,Enh),1)  # vector of cluster sizes
  nt = sum(nv)
  k1 = length(sup1resp)
  k2 = length(sup2resp)
  u = sample(1:k1,H,rep=TRUE)  # cluster-specific latent classes
  v = sample(1:k2,nt,rep=TRUE) # individual-specific latent classes
  data = matrix(NA,nt,2+nc1+nc2+2*J)
  
  ind = 0
  for(h in 1:H){
    xh = rnorm(nc1)
    for(i in 1:nv[h]){
      ind = ind+1
      zhi = rnorm(nc2)
      yhi = thi = rep(NA,J)
      dhi = sample(1:J,Jsel)
      tmp = pnorm(sup1resp[u[h]]+sup2resp[v[ind]]+c(xh%*%coef1resp+zhi%*%coef2resp)-dif[dhi])
      yhi[dhi] = 1*(runif(Jsel)<tmp)
      tmp = speed0[dhi]-speed1[dhi]*(sup1time[u[h]]+sup2time[v[ind]]+c(xh%*%coef1time+zhi%*%coef2time))
      thi[dhi] = exp(rnorm(Jsel,tmp,sqrt(si2)))
      data[ind,] = c(h,i,xh,zhi,yhi,thi)
    }
  }
  data = as.data.frame(data)
  names(data) = c("cluster", "individual", paste0("X", 1:nc1), paste0("Z", 1:nc2), 
                  paste0("resp", 1:J), paste0("time", 1:J))
  
  return(list(Data = data, U = u, V = v))
}

