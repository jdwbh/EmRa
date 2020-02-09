Dif<-function(x, d = 1, lag = 1,...){ # Delta函数, Dif(xt)=xt-xt_1,有别与R中的diff函数，保持与x长度相同
  c(rep(NA, d*lag),diff(x,lag=lag,differences=d,...))
}
Lag<-function(x, d = 1){ # Lag函数有别与R中的lag函数，可保持与x长度相同
    if (d != as.integer(d) || d < 0) 
        stop("d must be a non-negative integer")
    if (d == 0) return(x)
    else return(c(rep(NA, d), x[1:(length(x)-d)]))
}

plot.Lag<-function(y,d=1,plot=FALSE){ # plot(yt_1,yt)，可计算相关系数绘制回归线 
  yt=y
  yt_1=Lag(yt,d)
  #r=cor.test(yt_1,yt,use='complete.obs')
  plot(yt_1,yt);#,ylab=paste(names(y),"t"),xlab=paste(names(y),"t-1")
  if(plot){
    fm=lm(yt~yt_1-1)
    title(paste('ruo = ',round(fm$coef[[1]],4)))
    #deparse(substitute(y))
    abline(fm);
    summary(fm)$coef
  }
  #c(r=r[[4]],t=r[[1]],p=r[[3]])
}

D.W<-function(res){  # 计算D.W.值
  sum(diff(res)^2)/sum(res^2) 
}

DW.test<-function(formula, order.by = NULL, alternative = c("greater", 
    "two.sided", "less"), iterations = 15, exact = NULL, tol = 1e-10, 
    data = list()) 
{
    dname <- paste(deparse(substitute(formula)))
    alternative <- match.arg(alternative)
    if (!inherits(formula, "formula")) {
        if (!is.null(w <- weights(formula))) {
            if (!isTRUE(all.equal(as.vector(w), rep(1L, length(w))))) 
                stop("weighted regressions are not supported")
        }
        X <- if (is.matrix(formula$x)) 
            formula$x
        else model.matrix(terms(formula), model.frame(formula))
        y <- if (is.vector(formula$y)) 
            formula$y
        else model.response(model.frame(formula))
    }
    else {
        mf <- model.frame(formula, data = data)
        y <- model.response(mf)
        X <- model.matrix(formula, data = data)
    }
    if (!is.null(order.by)) {
        if (inherits(order.by, "formula")) {
            z <- model.matrix(order.by, data = data)
            z <- as.vector(z[, ncol(z)])
        }
        else {
            z <- order.by
        }
        X <- as.matrix(X[order(z), ])
        y <- y[order(z)]
    }
    n <- nrow(X)
    if (is.null(exact)) 
        exact <- (n < 100)
    k <- ncol(X)
    res <- lm.fit(X, y)$residuals
    dw <- sum(diff(res)^2)/sum(res^2)
    Q1 <- chol2inv(qr.R(qr(X)))
    if (n < 3) {
        warning("not enough observations for computing a p value, set to 1")
        pval <- 1
    }
    else {
        if (exact) {
            A <- diag(c(1, rep(2, n - 2), 1))
            A[abs(row(A) - col(A)) == 1] <- -1
            MA <- diag(rep(1, n)) - X %*% Q1 %*% t(X)
            MA <- MA %*% A
            ev <- eigen(MA)$values[1:(n - k)]
            if (any(Im(ev) > tol)) 
                warning("imaginary parts of eigenvalues discarded")
            ev <- Re(ev)
            ev <- ev[ev > tol]
            pdw <- function(dw) .Fortran("pan", as.double(c(dw, 
                ev)), as.integer(length(ev)), as.double(0), as.integer(iterations), 
                x = double(1), PACKAGE = "lmtest")$x
            pval <- switch(alternative, two.sided = (2 * min(pdw(dw), 
                1 - pdw(dw))), less = (1 - pdw(dw)), greater = pdw(dw))
            if (is.na(pval) || ((pval > 1) | (pval < 0))) {
                warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
                exact <- FALSE
            }
        }
        if (!exact) {
            if (n < max(5, k)) {
                warning("not enough observations for computing an approximate p value, set to 1")
                pval <- 1
            }
            else {
                AX <- matrix(as.vector(filter(X, c(-1, 2, -1))), 
                  ncol = k)
                AX[1, ] <- X[1, ] - X[2, ]
                AX[n, ] <- X[n, ] - X[(n - 1), ]
                XAXQ <- t(X) %*% AX %*% Q1
                P <- 2 * (n - 1) - sum(diag(XAXQ))
                Q <- 2 * (3 * n - 4) - 2 * sum(diag(crossprod(AX) %*% 
                  Q1)) + sum(diag(XAXQ %*% XAXQ))
                dmean <- P/(n - k)
                dvar <- 2/((n - k) * (n - k + 2)) * (Q - P * 
                  dmean)
                pval <- switch(alternative, two.sided = (2 * 
                  pnorm(abs(dw - dmean), sd = sqrt(dvar), lower.tail = FALSE)), 
                  less = pnorm(dw, mean = dmean, sd = sqrt(dvar), 
                    lower.tail = FALSE), greater = pnorm(dw, 
                    mean = dmean, sd = sqrt(dvar)))
            }
        }
    }
    alternative <- switch(alternative, two.sided = "true autocorrelation is not 0", 
        less = "true autocorrelation is less than 0", greater = "true autocorrelation is greater than 0")
    names(dw) <- "DW"
    RVAL <- list(statistic = dw, method = "Durbin-Watson test", 
        alternative = alternative, p.value = pval, data.name = dname)
    class(RVAL) <- "htest"
    return(RVAL)
}


White.test<-function(LM){  # White 异方差性检验
  n=nrow(LM$model)
  R2=summary(LM)$r.sq
  W=n*R2;W
  P=1-pchisq(W,2)
  c(W=W,P=P)
}
lm.test<-function(LM,dw=FALSE){ 
  sfm=summary(LM); #print(sfm)
  fc=sfm$coef;#cat('R2=',sfm$r.sq,'\n')
  cat("\n")
  print(round(rbind(b=fc[,1],t=fc[,3],p=fc[,4]),5))
  cat("\n R^2= ", sfm$r.sq,"\n")
  if(dw)
    cat(" D.W. = ", D.W(resid(LM)),"\n")
  cat("\n")
}

lm.plot<-function(LM,x=1:nrow(LM$model)){  # 拟合值和学生化残差
  par(mfrow=c(2,1),cex=0.75,mar=c(1,4,2,1)+0.1)
    matplot(x,cbind(LM$model[1],fitted(LM)),type='l',lty=c(3,1),lwd=2,ylab=names(LM$model[1]));
    title(paste(" ... Actual ", " --- Fitted "),cex.main=1)
    #plot(x,LM$model[1]);lines(x,fitted(LM));
    #plot(x,rstudent(LM),axes = FALSE,ylim=c(-4,4),type='b',ylab='rstudent'); abline(h=0)
    plot(x,rstudent(LM),axes = FALSE,frame.plot = TRUE,ylim=c(-4,4),type='o',ylab='rstudent'); 
    abline(h=c(-3,0,3),lty=3,col='blue');axis(2, at = seq(-4,4,by=1))
  par(mfrow=c(1,1))
}

lm.diag<-function(LM){
  n=nrow(LM$model);p=ncol(LM$model)
  X=as.matrix(cbind(rep(1,n),LM$model[,-1]));
  hi=diag(X%*%solve(t(X)%*%X)%*%t(X));
  ho=ifelse(hi>2*p/n,'*','');dh
  ri=rstudent(LM)
  ro=ifelse(abs(ri)>3,'*','')
  Di=1/(p+1)*ri^2*hi/(1-hi) 
  Do=ifelse(Di>1,'*','')
  data.frame(ri,ro,hi,ho,Di,Do)
}

acf.pacf<-function(y,plot=FALSE,...){ # 自相关与偏自相关及其检验与绘图
  AC=acf(y,plot=FALSE,...)$acf[-1]
  PAC=pacf(y,plot=FALSE,...)$acf
  k=length(AC)
  Q=rep(0,k); P=rep(0,k)
  for(i in 1:k){
     Qstat=Box.test(y,lag=i,type="Ljung-Box");
     Q[i]=Qstat$stat
     P[i]=Qstat$p.value
  }
  if(plot){
    #par(mfrow=c(2,1),cex=0.75,mar=c(4,4,2,1)+0.1)
    par(mfrow=c(2,1),cex=0.75,mar=c(0,4,2,1)+0.2)
      n=length(y)
      plot(AC,type='h',xlab='lag',...); #ylim=c(-1,1),
      abline(h=c(-1.96/sqrt(n),0,1.96/sqrt(n)),lty=c(3,1,3),col='blue');
      plot(PAC,type='h',xlab='lag',...); #ylim=c(-1,1),
      abline(h=c(-1.96/sqrt(n),0,1.96/sqrt(n)),lty=c(3,1,3),col='blue');
    par(mfrow=c(1,1))
  }
  data.frame(AC=round(AC[1:k],4),PAC=round(PAC[1:k],4),Qstat=round(Q,4),Prob=round(P,4)) 
}
acf.test<-function(y,...){ # 自相关系数及其检验与绘图
  AC=acf(y,...)$acf #AC=acf(y,...)$acf[-1]
  k=length(AC)
  Q=rep(0,k); P=rep(0,k)
  for(i in 1:k){
     Qstat=Box.test(y,lag=i,type="Ljung-Box");
     Q[i]=Qstat$stat
     P[i]=Qstat$p.value
  }
  
  data.frame(AC=round(AC,4),Qstat=round(Q,4),Prob=round(P,4)) 
}

pacf.test<-function(y,...){ # 偏自相关系数及其检验与绘图
  PAC=pacf(y,...)$acf #AC=acf(y,...)$acf[-1]
  k=length(PAC)
  Q=rep(0,k); P=rep(0,k)
  for(i in 1:k){
     Qstat=Box.test(y,lag=i,type="Ljung-Box");
     Q[i]=Qstat$stat
     P[i]=Qstat$p.value
  }
  data.frame(PAC=round(PAC,4),Qstat=round(Q,4),Prob=round(P,4)) 
}

arima.test <- function (yt,...,summary=FALSE){ # arima系数的假设检验 
     object=arima(yt,...)
     coef <- coef(object) 
     if (length(coef) > 0) { 
         mask <- object$mask 
         sdev <- sqrt(diag(vcov(object))) 
         t.rat <- rep(NA, length(mask)) 
         t.rat[mask] <- coef[mask]/sdev 
         pt <- 2 * pnorm(-abs(t.rat)) 
         setmp <- rep(NA, length(mask)) 
         setmp[mask] <- sdev 
         #sum <- signif(rbind(coef, setmp, t.rat, pt),4) 
         print(object$call);cat("\n") 
         sum <- round(rbind(coef, setmp, t.rat, pt),4) 
         dimnames(sum) <- list(c("coef", "s.e.", "t-stat", "p-value"), names(coef)) 
         if(summary)
           print(t(sum))
         else print(sum)
         e=resid(object)
         D.W=sum(diff(e)^2)/sum(e^2);D.W
         cat("\nsigma^2 estimated as ",object$sigma2," :  log likelihood = ", object$loglik, "aic = ", object$aic, "\n") 
         cat("Durbin-Watson stat  ", D.W, "\n") 
         #return(sum) 
     } else return(NA) 
} 
acf2=function(series,max.lag=NULL){
  num=length(series)
  if (is.null(max.lag)) max.lag=ceiling(10+sqrt(num))
  if (max.lag > (num-1)) stop("Number of lags exceeds number of observations")
  ACF=acf(series, max.lag, plot=FALSE)$acf[-1]
  PACF=pacf(series, max.lag, plot=FALSE)$acf
  LAG=1:max.lag/frequency(series)
  minA=min(ACF)
  minP=min(PACF)
  U=2/sqrt(num)
  L=-U
  minu=min(minA,minP,L)-.01
  old.par <- par(no.readonly = TRUE)
  par(mfrow=c(2,1), mar = c(3,3,2,0.8),oma = c(1,1.2,1,1), mgp = c(1.5,0.6,0))
  plot(LAG, ACF, type="h",ylim=c(minu,1), 
    main=paste("Series: ",deparse(substitute(series))))
    abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,4,4))
  plot(LAG, PACF, type="h",ylim=c(minu,1))
    abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,4,4))
  on.exit(par(old.par))  
  ACF<-round(ACF,2); PACF<-round(PACF,2)    
  return(cbind(ACF, PACF)) 
}
