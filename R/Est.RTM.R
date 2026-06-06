Est.RTM <- function(data, nc1, nc2, J, control = NULL) {

# Estimate IRT_Time model without classes

#---- fix input ----
  rel.tol = ifelse(is.null(control$rel.tol), 1e-8, control$rel.tol)
  out.se = ifelse(is.null(control$out.se), FALSE, control$out.se)

#---- preliminaries ----
  H = max(data$cluster)
  nv = tapply(data$individual,data$cluster,max)
  nt = nrow(data)
  if(nc1==0 & nc2==0) XZ = NULL else XZ = as.matrix(data[,2+(1:(nc1+nc2))])
  YY = as.matrix(data[,2+nc1+nc2+(1:J)])
  LT = log(as.matrix(data[,2+nc1+nc2+J+(1:J)]))
  if(nc1==0 & nc2==0) XZv = NULL else XZv = XZ%x%matrix(1,J,1)
  YYv = c(t(YY))
  LTv = c(t(LT))
  indobs = which(!is.na(YYv))
  if(nc1==0 & nc2==0) XZvobs = NULL else XZvobs = XZv[indobs,]
  YYvobs = YYv[indobs]
  LTvobs = LTv[indobs]
  DD.resp = cbind(XZvobs,-(matrix(1,nt,1)%x%diag(J))[indobs,])
  INT0obs.time = (matrix(1,nt,1)%x%diag(J))[indobs,]
  np1 = nc1+nc2+J
  np2 = nc1+nc2
  np3 = J+(nc1+nc2>0)*(J-1)
  ntobs = length(indobs)

#---- starting values ----
  coef.resp = rep(0,np1)
  speed0.time = rep(0,J)
  speed1.time = rep(1,J-1)
  lin.resp = c(DD.resp%*%coef.resp)
  pred.resp = pnorm(lin.resp)

# update parameters of the response model
  dnorm.resp = dnorm(lin.resp)
  p1mp.pres = pred.resp*(1-pred.resp)
  sc.resp = c(t(DD.resp)%*%((YYvobs-pred.resp)*dnorm.resp/p1mp.pres))
  Fi.resp = t(DD.resp)%*%((dnorm.resp^2/p1mp.pres)*DD.resp)
  if(rcond(Fi.resp)>10^-15){
    coef.resp = coef.resp+solve(Fi.resp)%*%sc.resp
  }else{
    coef.resp = coef.resp+ginv(Fi.resp)%*%sc.resp
    print("Singular Fi.resp matrix in updating coef.resp:")
    print(rcond(Fi.resp))
    print(coef.resp)
  }
  lin.resp = c(DD.resp%*%coef.resp)
  pred.resp = pmin(pmax(pnorm(lin.resp),10^-10),1-10^-10)

# update parameters of the time model
  if(nc1==0 & nc2==0){
    DD2.time = INT0obs.time
    offset.time = rep(0,ntobs)
    coef1.time = NULL
    if(rcond(t(DD2.time)%*%DD2.time)>10^-15){
      coef2.time = c(solve(t(DD2.time)%*%DD2.time)%*%t(DD2.time)%*%(LTvobs-offset.time))
    }else{
      coef2.time = c(ginv(t(DD2.time)%*%DD2.time)%*%t(DD2.time)%*%(LTvobs-offset.time))
      print("Singular t(DD2.time)%*%DD2.time matrix in updating coef2.time:")
      print(rcond(t(DD2.time)%*%DD2.time))
      print(coef2.time)
    }
    speed0.time = coef2.time
    speed1.time = NULL
    pred.time = c(DD2.time%*%coef2.time)
  }else{
    DD.time = rep(-c(1,speed1.time),nt)[indobs]*XZvobs
    if(rcond(t(DD.time)%*%DD.time)>10^-15){
      coef1.time = c(solve(t(DD.time)%*%DD.time)%*%t(DD.time)%*%(LTvobs-rep(speed0.time,nt)[indobs]))
    }else{
      coef1.time = c(ginv(t(DD.time)%*%DD.time)%*%t(DD.time)%*%(LTvobs-rep(speed0.time,nt)[indobs]))
      print("Singular t(DD.time)%*%DD.time matrix in updating coef1.time:")
      print(rcond(t(DD.time)%*%DD.time))
      print(coef1.time)
    }
    INT1obs.time = -c(XZvobs%*%coef1.time)*INT0obs.time
    offset.time = INT1obs.time[,1]
    DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
    if(rcond(t(DD2.time)%*%DD2.time)>10^-15){
      coef2.time = c(solve(t(DD2.time)%*%DD2.time)%*%t(DD2.time)%*%(LTvobs-offset.time))
    }else{
      coef2.time = c(ginv(t(DD2.time)%*%DD2.time)%*%t(DD2.time)%*%(LTvobs-offset.time))
      print("Singular t(DD2.time)%*%DD2.time matrix in updating coef2.time")
      print(rcond(t(DD2.time)%*%DD2.time))
      print(coef2.time)
    }
    speed0.time = coef2.time[1:J]
    speed1.time = coef2.time[J+(1:(J-1))]
    pred.time = c(DD2.time%*%coef2.time+offset.time)
  }
  res.time = LTvobs-pred.time
  si2.time = mean(res.time^2)

# compute log-likelihood
  lk = sum(dbinom(YYvobs,1,pred.resp,log=TRUE)+dnorm(LTvobs,pred.time,sqrt(si2.time),log=TRUE))
  t0 = proc.time()[3]
  cat("\n")
  cat("* Maximum likelihood estimation *\n")
  cat("------------|-------------|-------------|-------------|\n")
  cat("     it     |      lk     |    lk-lko   | lk-lko rel. |\n")
  cat("------------|-------------|-------------|-------------|\n")
  cat(sprintf("%11g", c(0, lk)), "\n", sep = " | ")

#---- iterate until convergence ----
  lko = lk
  it = 0
  nsec = 2.5
  while(((lk-lko)/abs(lko)>rel.tol & it<1000) | it==0){
    lko = lk
    it = it+1

# update parameters of the response model
    dnorm.resp = dnorm(lin.resp)
    p1mp.pres = pred.resp*(1-pred.resp)
    sc.resp = c(t(DD.resp)%*%((YYvobs-pred.resp)*dnorm.resp/p1mp.pres))
    Fi.resp = t(DD.resp)%*%((dnorm.resp^2/p1mp.pres)*DD.resp)
    if(rcond(Fi.resp)>10^-15){
      coef.resp = coef.resp+solve(Fi.resp)%*%sc.resp
    }else{
      coef.resp = coef.resp+ginv(Fi.resp)%*%sc.resp
      print("Singular Fi.resp matrix in updating coef.resp:")
      print(rcond(Fi.resp))
      print(coef.resp)
    }
    lin.resp = c(DD.resp%*%coef.resp)
    pred.resp = pmin(pmax(pnorm(lin.resp),10^-10),1-10^-10)

# update parameters of the time model
    if(nc1==0 & nc2==0){
      DD2.time = INT0obs.time
      offset.time = rep(0,ntobs)
      coef1.time = NULL
      if(rcond(t(DD2.time)%*%DD2.time)>10^-15){
        coef2.time = c(solve(t(DD2.time)%*%DD2.time)%*%t(DD2.time)%*%(LTvobs-offset.time))
      }else{
        coef2.time = c(ginv(t(DD2.time)%*%DD2.time)%*%t(DD2.time)%*%(LTvobs-offset.time))
        print("Singular matrix t(DD2.time)%*%DD2.time in updating coef2.time:")
        print(rcond(t(DD2.time)%*%DD2.time))
        print(coef2.time)
      }      
      speed0.time = coef2.time
      speed1.time = NULL
      pred.time = c(DD2.time%*%coef2.time)
    }else{
      DD.time = rep(-c(1,speed1.time),nt)[indobs]*XZvobs
      if(rcond(t(DD.time)%*%DD.time)>10^-15){
        coef1.time = c(solve(t(DD.time)%*%DD.time)%*%t(DD.time)%*%(LTvobs-rep(speed0.time,nt)[indobs]))
      }else{
        coef1.time = c(ginv(t(DD.time)%*%DD.time)%*%t(DD.time)%*%(LTvobs-rep(speed0.time,nt)[indobs]))
        print("Singular matrix DD.time in updating coef1.time:")
        print(rcond(t(DD.time)%*%DD.time))
        print(coef1.time)
      }
      INT1obs.time = -c(XZvobs%*%coef1.time)*INT0obs.time
      offset.time = INT1obs.time[,1]
      DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
      if(rcond(t(DD2.time)%*%DD2.time)>10^-15){
        coef2.time = c(solve(t(DD2.time)%*%DD2.time)%*%t(DD2.time)%*%(LTvobs-offset.time))
      }else{
        coef2.time = c(ginv(t(DD2.time)%*%DD2.time)%*%t(DD2.time)%*%(LTvobs-offset.time))
        print("Singular matrix t(DD2.time)%*%DD2.time in updating coef2.time:")
        print(rcond(t(DD2.time)%*%DD2.time))
        print(coef2.time)
      }
      speed0.time = coef2.time[1:J]
      speed1.time = coef2.time[J+(1:(J-1))]
      pred.time = c(DD2.time%*%coef2.time+offset.time)
    }
    res.time = LTvobs-pred.time
    si2.time = mean(res.time^2)

# compute log-likelihood
    lk = sum(dbinom(YYvobs,1,pred.resp,log=TRUE)+dnorm(LTvobs,pred.time,sqrt(si2.time),log=TRUE))
    flag = TRUE
    if(proc.time()[3]-t0>nsec){
      cat(sprintf("%11g", c(it,lk,lk-lko,(lk-lko)/abs(lko))), "\n", sep = " | ")
      nsec = nsec+2.5
      flag = FALSE
    }
  }
  if(flag) cat(sprintf("%11g", c(it,lk,lk-lko,(lk-lko)/abs(lko))), "\n", sep = " | ")
  cat("------------|-------------|-------------|-------------|\n")

#---- compute standard errors ----
  np = np1+np2+np3+1
  if(out.se){
    cat("\n")
    cat("* Compute standard errors *\n")
    th = c(coef.resp,coef1.time,coef2.time,si2.time)
    out = lk_RTM(th,data,nc1,nc2,J)
    scn = rep(0,np)
    Dn = matrix(0,np,np)
    t0 = proc.time()[3]
    nsec = 2.5
    for(j in 1:np){
      th1 = th; th1[j] = th1[j]+10^-6
      out1 = lk_RTM(th1,data,nc1,nc2,J)
      scn[j] = 10^6*(out1$lk-out$lk)
      Dn[,j] = 10^6*(out1$sc-out$sc)
      flag = TRUE
      if(proc.time()[3]-t0>nsec){
        cat("(",j,"/",np,")\n", sep = "")
        nsec = nsec+2.5
        flag = FALSE
      }
    }
    if(flag) cat("(",j,"/",np,")\n", sep = "")
    # print(c(out$lk,lk,out$lk/lk-1)) # to check log-likelihood
    # print(cbind(out$sc,scn,out$sc/scn-1)) # to check derivatives
    Dn = (Dn+t(Dn))/2
    if(rcond(Dn)>10^-15){
      se = sqrt(diag(solve(-Dn)))
    }else{
      se = try(sqrt(diag(ginv(-Dn))))
      if(inherits(se,"try-error")) se = rep(NA,length(th))
      print("Singular matrix Dn in computing se:")
      print(rcond(Dn))
      print(se)
    }
    secoef.resp = se[1:np1]
    secoef1.time = se[np1+(1:np2)]
    secoef2.time = se[np1+np2+(1:np3)]
    sespeed0.time = secoef2.time[1:J]
    sespeed1.time = secoef2.time[J+(1:(J-1))]
    sesi2.time = se[np1+np2+np3+1]
    cat("")
  }

#---- final output ----
  AIC = -2*lk+2*np
  BIC = -2*lk+log(nt)*np
  if(nc1==0) coefX.resp = NULL else coefX.resp = coef.resp[1:nc1]
  if(nc2==0) coefZ.resp = NULL else coefZ.resp = coef.resp[nc1+(1:nc2)]
  diff.resp = coef.resp[nc1+nc2+(1:J)]
  if(nc1==0) coefX.time = NULL else coefX.time = coef1.time[1:nc1]
  if(nc2==0) coefZ.time = NULL else coefZ.time = coef1.time[nc1+(1:nc2)]
  speed1.time = c(0,speed1.time)
  if(out.se){
    if(nc1==0) secoefX.resp = NULL else secoefX.resp = secoef.resp[1:nc1]
    if(nc2==0) secoefZ.resp = NULL else secoefZ.resp = secoef.resp[nc1+(1:nc2)]
    sediff.resp = secoef.resp[nc1+nc2+(1:J)]
    if(nc1==0) secoefX.time = NULL else secoefX.time = secoef1.time[1:nc1]
    if(nc2==0) secoefZ.time = NULL else secoefZ.time = secoef1.time[nc1+(1:nc2)]
    sespeed1.time = c(0,sespeed1.time)
  }
  out = list(lk=lk,coefX.resp=coefX.resp,coefZ.resp=coefZ.resp,diff.resp=diff.resp,
             coefX.time=coefX.time,speed0.time=speed0.time,speed1.time=speed1.time,
             coefZ.time=coefZ.time,si2.time=si2.time,nit=it,np=np,
             AIC=AIC,BIC=BIC)
  if(out.se) out = c(out,list(secoefX.resp=secoefX.resp,secoefZ.resp=secoefZ.resp,
                     sediff.resp=sediff.resp,secoefX.time=secoefX.time,
                     sespeed0.time=sespeed0.time,sespeed1.time=sespeed1.time,
                     secoefZ.time=secoefZ.time,sesi2.time=sesi2.time,Dn=Dn))
  return(out)

}