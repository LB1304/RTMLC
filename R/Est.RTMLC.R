Est.RTMLC <- function(data, nc1, nc2, J, k1, k2, param.init, control = NULL) {

# Estimate IRT_Time model without classes

#---- fix input ----
  rel.tol = ifelse(is.null(control$rel.tol), 1e-8, control$rel.tol)
  out.se = ifelse(is.null(control$out.se), FALSE, control$out.se)
  init.rand = ifelse(is.null(control$init.rand), FALSE, control$init.rand)
  coefX.resp = param.init$coefX.resp
  coefZ.resp = param.init$coefZ.resp
  diff.resp = param.init$diff.resp
  speed0.time = param.init$speed0.time
  speed1.time = param.init$speed1.time
  coefX.time = param.init$coefX.time
  coefZ.time= param.init$coefZ.time
  si2.time = param.init$si2.time
  if(nc1==0 & nc2==0 & (k1>1 | k2>1)) speed1.time = rep(1,J-1) else speed1.time = speed1.time[-1]

#---- preliminaries ----
  cluster = as.vector(data$cluster)
  H = max(cluster)
  nv = tapply(data$individual,data$cluster,max)
  nt = nrow(data)
  if(nc1==0 & nc2==0) XZ = NULL else XZ = as.matrix(data[,2+(1:(nc1+nc2))])
  YY = as.matrix(data[,2+nc1+nc2+(1:J)])
  LT = log(as.matrix(data[,2+nc1+nc2+J+(1:J)]))
  if(nc1==0 & nc2==0) XZv = NULL else XZv = XZ%x%matrix(1,J,1)
  YYv = c(t(YY))
  LTv = c(t(LT))
  indobs = which(!is.na(YYv))
  if(nc1==0 & nc2==0) XZvobs = NULL else  XZvobs = XZv[indobs,]
  YYvobs = YYv[indobs]
  LTvobs = LTv[indobs]
  DD.resp = cbind(XZvobs,-(matrix(1,nt,1)%x%diag(J))[indobs,])
  INT0obs.time = (matrix(1,nt,1)%x%diag(J))[indobs,]
  ntobs = length(indobs)
  np1 = nc1+nc2+J+k1+k2-2
  np2 = nc1+nc2+k1+k2-2
  np3 = J+(nc1+nc2+k1+k2-2>0)*(J-1)

#---- starting values ----
  if(k1==1){
    sup1.resp = sup1.time = NULL
  }else{
    if(init.rand){
      tmp = rnorm(k1,0,2)
      sup1.resp = tmp[-1]
      diff.resp = diff.resp+mean(tmp)
    }else{
      sup1.resp = seq(0,4,length.out=k1)[-1]
      diff.resp = diff.resp+2
    }
    sup1.time = rep(0,k1-1)
  }
  if(k2==1){
    sup2.resp = sup2.time = NULL
  }else{
    if(init.rand){
      tmp = rnorm(k2,0,2)
      sup2.resp = tmp[-1]
      diff.resp = diff.resp+mean(tmp)
    }else{
      sup2.resp = seq(0,4,length.out=k2)[-1]
      diff.resp = diff.resp+2
    }
    sup2.time = rep(0,k2-1)
  }
  if(init.rand){
    piv1 = runif(k1)
    piv1 = piv1/sum(piv1)
    piv2 = runif(k2)
    piv2 = piv2/sum(piv2)
  }else{
    piv1 = rep(1/k1,k1)
    piv2 = rep(1/k2,k2)
  }
  coef.resp = c(coefX.resp,coefZ.resp,diff.resp)
  coef1.time = c(coefX.time,coefZ.time)
  coef2.time = c(speed0.time,speed1.time)
  lin.resp = c(DD.resp%*%coef.resp)

# compute log-likelihood
  lQcc = array(0,c(nt,k1,k2))
  for(u in 1:k1){
    for(v in 1:k2){
      tmp = 0
      if(u>1) tmp = sup1.resp[u-1]
      if(v>1) tmp = tmp+sup2.resp[v-1]
      pred.resp = pmin(pmax(pnorm(tmp+lin.resp),10^-10),1-10^-10)
      tmp = 0
      if(u>1) tmp = sup1.time[u-1]
      if(v>1) tmp = tmp+sup2.time[v-1]
      if(nc1==0 & nc2==0 & k1==1 & k2==1){
        offset.time = rep(0,ntobs)
        DD2.time = INT0obs.time
      }else if(nc1==0 & nc2==0){
        INT1obs.time = -tmp*INT0obs.time
        offset.time = INT1obs.time[,1]
        DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
      }else{
        INT1obs.time = -c(tmp+XZvobs%*%coef1.time)*INT0obs.time
        offset.time = INT1obs.time[,1]
        DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
      }
      pred.time = c(DD2.time%*%coef2.time+offset.time)
      tmp = rep(NA,nt*J)
      tmp[indobs] = dbinom(YYvobs,1,pred.resp,log=TRUE)+dnorm(LTvobs,pred.time,sqrt(si2.time),log=TRUE)
      Tmp = matrix(tmp,J,nt)
      lQcc[,u,v] = colSums(Tmp,na.rm=TRUE)
    }
  }
  Tmp = apply(lQcc,1:2,max)
  Qcm1 = matrix(0,nt,k1)
  for(v in 1:k2) Qcm1 = Qcm1+exp(lQcc[,,v]-Tmp)*piv2[v]
  lQcm = log(Qcm1)+Tmp
  lPc = matrix(0,H,k1)
  for(u in 1:k1) lPc[,u] = tapply(lQcm[,u],cluster,sum)
  tmp = apply(lPc,1,max)
  lPm = log(c(exp(lPc-tmp)%*%piv1))+tmp
  lk = sum(lPm)
  t0 = proc.time()[3]
  cat("\n")
  cat("* Maximum likelihood estimation *\n")
  cat("------------|-------------|-------------|-------------|-------------|-------------|\n")
  cat("     it     |      k1     |      k2     |      lk     |    lk-lko   | lk-lko rel. |\n")
  cat("------------|-------------|-------------|-------------|-------------|-------------|\n")
  cat(sprintf("%11g", c(0,k1,k2,lk)), "\n", sep = " | ")
  
#---- iterate until convergence ----
  lko = lk
  it = 0
  nsec = 2.5
  while((abs(lk-lko)/abs(lko)>rel.tol & it<1000) | it==0){
    lko = lk
    it = it+1
    
# E-step
    Pp = exp(lPc-lPm)%*%diag(piv1)
    Pp = pmax(Pp,10^-10) # correction
    Pp = Pp/rowSums(Pp)
    Qcp = array(0,c(nt,k1,k2))
    for(u in 1:k1) Qcp[,u,] = exp(lQcc[,u,]-lQcm[,u])%*%diag(piv2)
    Qpp = array(0,c(nt,k1,k2))
    for(h in 1:H){
      indh = which(cluster==h)
      for(u in 1:k1) Qpp[indh,u,] = Pp[h,u]*Qcp[indh,u,] 
    }
    Qpp = pmax(Qpp,10^-10) # correction
    Qpp = Qpp/array(apply(Qpp,1,sum),c(nt,k1,k2))
    Qmp = apply(Qpp,c(1,3),sum)

# M-step
# update cluster weights
    if(k1>1){
      piv1 = colSums(Pp)/H
      piv1 = piv1/sum(piv1)
    }
    if(k2>1){
      piv2 = colSums(Qmp)/nt
      piv2 = piv2/sum(piv2)
    }

# update parameters for response
    sc.resp = rep(0,np1)
    Fi.resp = matrix(0,np1,np1)
    for(u in 1:k1){
      for(v in 1:k2){
        tmp = 0
        if(u>1) tmp = sup1.resp[u-1]
        if(v>1) tmp = tmp+sup2.resp[v-1]
        pred.resp = pmin(pmax(pnorm(tmp+lin.resp),10^-10),1-10^-10)
        dnorm.resp = dnorm(tmp+lin.resp)
        p1mp.pres = pred.resp*(1-pred.resp)
        Tmp = matrix(0,ntobs,k1+k2-2)
        if(u>1) Tmp[,u-1] = 1
        if(v>1) Tmp[,k1+v-2] = 1
        w = rep(Qpp[,u,v],each=J)[indobs]
        Tmp = cbind(Tmp,DD.resp)
        sc.resp = sc.resp+c(t(Tmp)%*%(w*(YYvobs-pred.resp)*dnorm.resp/p1mp.pres))
        Fi.resp = Fi.resp+t(Tmp)%*%((w*dnorm.resp^2/p1mp.pres)*Tmp)
      }
    }
    if(rcond(Fi.resp)>10^-15){
      tmp = solve(Fi.resp)%*%sc.resp
    }else{
      tmp = ginv(Fi.resp)%*%sc.resp
      print("Singular Fi.resp matrix in updating coef.resp:")
      print(rcond(Fi.resp))
      print(coef.resp)
    }
    dtmp = max(abs(tmp))
    if(dtmp>0.1) tmp = tmp/dtmp*0.1
    coef.resp = c(sup1.resp,sup2.resp,coef.resp)+tmp
    if(k1>1) sup1.resp = coef.resp[1:(k1-1)]
    if(k2>1) sup2.resp = coef.resp[k1-1+(1:(k2-1))]
    if(k1>1 | k2>1) coef.resp = coef.resp[-(1:(k1+k2-2))]
    lin.resp = c(DD.resp%*%coef.resp)

# update parameters of the time model
    if(nc1>0 | nc2>0 | k1>1 | k2>1){
      num = rep(0,np2)
      DEN = matrix(0,np2,np2)
      for(u in 1:k1){
        for(v in 1:k2){
          Tmp = matrix(0,ntobs,k1+k2-2)
          if(u>1) Tmp[,u-1] = 1
          if(v>1) Tmp[,k1+v-2] = 1
          Tmp = cbind(Tmp,XZvobs)
          DD.time = rep(-c(1,speed1.time),nt)[indobs]*Tmp
          w = rep(Qpp[,u,v],each=J)[indobs]
          num = num+c(t(DD.time)%*%(w*(LTvobs-rep(speed0.time,nt)[indobs])))
          DEN = DEN+t(DD.time)%*%(w*DD.time)
        }
      }
      if(rcond(DEN)>10^-15){
        coef1.time = solve(DEN)%*%num
      }else{
        coef1.time = ginv(DEN)%*%num
        print("Singular DEN matrix in updating coef1.time")
        print(rcond(DEN))
        print(coef1.time)
      }
    }
    if(k1>1) sup1.time = coef1.time[1:(k1-1)]
    if(k2>1) sup2.time = coef1.time[k1-1+(1:(k2-1))]
    if(k1>1 | k2>1) coef1.time = coef1.time[-(1:(k1+k2-2))]
    num = rep(0,np3)
    DEN = matrix(0,np3,np3)
    for(u in 1:k1){
      for(v in 1:k2){
        tmp = 0
        if(u>1) tmp = sup1.time[u-1]
        if(v>1) tmp = tmp+sup2.time[v-1]
        if(nc1==0 & nc2==0 & k1==1 & k2==1){
          offset.time = rep(0,ntobs)
          DD2.time = INT0obs.time
        }else if(nc1==0 & nc2==0){
          INT1obs.time = -tmp*INT0obs.time
          offset.time = INT1obs.time[,1]
          DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
        }else{
          INT1obs.time = -c(tmp+XZvobs%*%coef1.time)*INT0obs.time
          offset.time = INT1obs.time[,1]
          DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
        }
        w = rep(Qpp[,u,v],each=J)[indobs]
        num = num+c(t(DD2.time)%*%(w*(LTvobs-offset.time)))
        DEN = DEN+t(DD2.time)%*%(w*DD2.time)
      }
    }
    if(rcond(DEN)>10^-15){
      coef2.time = solve(DEN)%*%num
    }else{
      coef2.time = ginv(DEN)%*%num
      print("Singular DEN matrix in updating coef2.time:")
      print(rcond(DEN))
      print(coef2.time)
    }
    speed0.time = coef2.time[1:J]
    speed1.time = coef2.time[J+(1:(J-1))]
    si2.time = 0
    for(u in 1:k1){
      for(v in 1:k2){
        tmp = 0
        if(u>1) tmp = sup1.time[u-1]
        if(v>1) tmp = tmp+sup2.time[v-1]
        if(nc1==0 & nc2==0 & k1==1 & k2==1){
          offset.time = rep(0,ntobs)
          DD2.time = INT0obs.time
        }else if(nc1==0 & nc2==0){
          INT1obs.time = -tmp*INT0obs.time
          offset.time = INT1obs.time[,1]
          DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
        }else{
          INT1obs.time = -c(tmp+XZvobs%*%coef1.time)*INT0obs.time
          offset.time = INT1obs.time[,1]
          DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
        }
        w = rep(Qpp[,u,v],each=J)[indobs]
        pred.time = c(DD2.time%*%coef2.time+offset.time)
        res.time = LTvobs-pred.time
        si2.time = si2.time+sum(w*res.time^2)
      }
    }
    si2.time = si2.time/ntobs
    
# compute log-likelihood
    lQcc = array(0,c(nt,k1,k2))
    for(u in 1:k1){
      for(v in 1:k2){
        tmp = 0
        if(u>1) tmp = sup1.resp[u-1]
        if(v>1) tmp = tmp+sup2.resp[v-1]
        pred.resp = pmin(pmax(pnorm(tmp+lin.resp),10^-10),1-10^-10)
        tmp = 0
        if(u>1) tmp = sup1.time[u-1]
        if(v>1) tmp = tmp+sup2.time[v-1]
        if(nc1==0 & nc2==0 & k1==1 & k2==1){
          offset.time = rep(0,ntobs)
          DD2.time = INT0obs.time
        }else if(nc1==0 & nc2==0){
          INT1obs.time = -tmp*INT0obs.time
          offset.time = INT1obs.time[,1]
          DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
        }else{
          INT1obs.time = -c(tmp+XZvobs%*%coef1.time)*INT0obs.time
          offset.time = INT1obs.time[,1]
          DD2.time = cbind(INT0obs.time,INT1obs.time[,-1])
        }
        pred.time = c(DD2.time%*%coef2.time+offset.time)
        tmp = rep(NA,nt*J)
        tmp[indobs] = dbinom(YYvobs,1,pred.resp,log=TRUE)+dnorm(LTvobs,pred.time,sqrt(si2.time),log=TRUE)
        Tmp = matrix(tmp,J,nt)
        lQcc[,u,v] = colSums(Tmp,na.rm=TRUE)
      }
    }
    Tmp = apply(lQcc,1:2,max)
    Qcm1 = matrix(0,nt,k1)
    for(v in 1:k2) Qcm1 = Qcm1+exp(lQcc[,,v]-Tmp)*piv2[v]
    lQcm = log(Qcm1)+Tmp
    lPc = matrix(0,H,k1)
    for(u in 1:k1) lPc[,u] = tapply(lQcm[,u],cluster,sum)
    tmp = apply(lPc,1,max)
    lPm = log(c(exp(lPc-tmp)%*%piv1))+tmp
    lk = sum(lPm)
    flag = TRUE
    if(proc.time()[3]-t0>nsec){
      cat(sprintf("%11g", c(it,k1,k2,lk,lk-lko,(lk-lko)/abs(lko))), "\n", sep = " | ")
      nsec = nsec+2.5
      flag = FALSE
    }
  }
  cat(sprintf("%11g", c(it,k1,k2,lk,lk-lko,(lk-lko)/abs(lko))), "\n", sep = " | ")
  cat("------------|-------------|-------------|-------------|-------------|-------------|\n")

#---- compute standard errors ----
  np = np1+np2+np3+1+k1+k2-2
  if(out.se){
    th = c(sup1.resp,sup2.resp,coef.resp,sup1.time,sup2.time,coef1.time,coef2.time,
           si2.time)
    if(k1>1){
      L1 = cbind(-1,diag(k1-1))
      th = c(th,L1%*%log(piv1))
    }  
    if(k2>1){
      L2 = cbind(-1,diag(k2-1))
      th = c(th,L2%*%log(piv2))
    }
    out = lk_RTMLC(th,data,nc1,nc2,J,k1,k2)
    scn = rep(0,np)
    Dn = matrix(0,np,np)
    cat("\n")
    cat("* Compute standard errors *\n")
    for(j in 1:np){
      th1 = th; th1[j] = th1[j]+10^-6
      out1 = lk_RTMLC(th1,data,nc1,nc2,J,k1,k2)
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
    # print(c(out$lk,lk,out$lk/lk-1))
    # print(cbind(out$sc,scn,out$sc/scn-1)) # to check derivatives
    Dn = (Dn+t(Dn))/2
    if(rcond(Dn)>10^-15){
      se = sqrt(diag(solve(-Dn)))
    }else{
      se = try(sqrt(diag(ginv(-Dn))))
      if(inherits(se,"try-error")) se = rep(NA,length(th))
      print("Singular Fi.resp matrix in computing se:")
      print(rcond(Dn))
      print(se)
    }
    secoef.resp = se[1:np1]
    if(k1==1) sesup1.resp = NULL else sesup1.resp = secoef.resp[1:(k1-1)]
    if(k2==1) sesup2.resp = NULL else sesup2.resp = secoef.resp[(k1-1)+(1:(k2-1))]
    secoef.resp = secoef.resp[-(1:(k1+k2-2))]
    secoef1.time = se[np1+(1:np2)]
    if(k1==1) sesup1.time = NULL else sesup1.time = secoef1.time[1:(k1-1)]
    if(k2==1) sesup2.time = NULL else sesup2.time = secoef1.time[(k1-1)+(1:(k2-1))]
    secoef1.time = secoef1.time[-(1:(k1+k2-2))]
    secoef2.time = se[np1+np2+(1:np3)]
    sespeed0.time = secoef2.time[1:J]
    sespeed1.time = secoef2.time[J+(1:(J-1))]
    sesi2.time = se[np1+np2+np3+1]
    if(k1==1) selpiv1 = NULL else selpiv1 = se[np1+np2+np3+1+(1:(k1-1))]
    if(k2==1) selpiv2 = NULL else selpiv2 = se[np1+np2+np3+1+k1-1+(1:(k2-1))]
  }

#---- final output ----
  AIC = -2*lk+2*np
  BIC = -2*lk+log(nt)*np
  if(nc1==0) coefX.resp = NULL else coefX.resp = coef.resp[1:nc1]
  if(nc2==0) coefZ.resp = NULL else coefZ.resp = coef.resp[nc1+(1:nc2)]
  diff.resp = coef.resp[nc1+nc2+(1:J)]
  if(nc1==0) coefX.time = NULL else coefX.time = coef1.time[1:nc1]
  if(nc2==0) coefZ.time = NULL else coefZ.time = coef1.time[nc1+(1:nc2)]
  speed1.time = c(1,speed1.time)
  sup1.resp = c(0,sup1.resp)
  sup2.resp = c(0,sup2.resp)
  sup1.time = c(0,sup1.time)
  sup2.time = c(0,sup2.time)
  if(out.se){
    if(nc1==0) secoefX.resp = NULL else secoefX.resp = secoef.resp[1:nc1]
    if(nc2==0) secoefZ.resp = NULL else secoefZ.resp = secoef.resp[nc1+(1:nc2)]
    sediff.resp = secoef.resp[nc1+nc2+(1:J)]
    if(nc1==0) secoefX.time = NULL else secoefX.time = secoef1.time[1:nc1]
    if(nc2==0) secoefZ.time = NULL else secoefZ.time = secoef1.time[nc1+(1:nc2)]
    sespeed1.time = c(0,sespeed1.time)
    sesup1.resp = c(0,sesup1.resp)
    sesup2.resp = c(0,sesup2.resp)
    sesup1.time = c(0,sesup1.time)
    sesup2.time = c(0,sesup2.time)
  }
  out = list(lk=lk,coefX.resp=coefX.resp,coefZ.resp=coefZ.resp,diff.resp=diff.resp,
             coefX.time=coefX.time,coefZ.time=coefZ.time,speed0.time=speed0.time,
             speed1.time=speed1.time,si2.time=si2.time,piv1=piv1,piv2=piv2,
             sup1.resp=sup1.resp,sup2.resp=sup2.resp,sup1.time=sup1.time,
             sup2.time=sup2.time,nit=it,AIC=AIC,BIC=BIC)
  if(out.se) out = c(out,list(secoefX.resp=secoefX.resp,secoefZ.resp=secoefZ.resp,
                              sediff.resp=sediff.resp,secoefX.time=secoefX.time,
                              secoefZ.time=secoefZ.time,sespeed0.time=sespeed0.time,
                              sespeed1.time=sespeed1.time,sesi2.time=sesi2.time,
                              selpiv1=selpiv1,selpiv2=selpiv2,sesup1.resp=sesup1.resp,
                              sesup2.resp=sesup2.resp,sesup1.time=sesup1.time,
                              sesup2.time=sesup2.time,Dn=Dn))
  return(out)

}