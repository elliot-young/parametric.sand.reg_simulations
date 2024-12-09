# Reproducibility code for simulations in Section 4.1 of "Sandwich Regression: Accurate and robust inference in generalised linear multilevel and longitudinal models"

# Linear model simulation
LAMBDA=0 # Vary to 0 (homoscedastic) or 3 (heteroscedastic) accordinly
betas_unw_ALL = betas_gee_ALL = betas_mle_ALL = betas_sl_ALL = betas_sl_naive_ALL = list()
Is = c(5,10,seq(20,200,by=20))
MSE_unw = MSE_gee = MSE_mle = MSE_sl = MSE_sl_naive = rep(0,length(Is))
n=4
arc = 0
MMM=10000
for (I in Is) {
  arc=arc+1
  betas_unw = betas_gee = betas_mle = betas_sl = betas_sl_naive = numeric()
  SEEDj = 0
  Sigma0root = matrix(0.9,n,n) + diag(1-0.9,n)
  s = function(x) 1+LAMBDA*exp(-2*x^2)
  while (length(betas_gee)<MMM) {
    SEEDj = SEEDj + 1
    {
      x = rnorm(n*I)
      y = rep(0,n*I)
      for (i in 1:I) {
        x_i = x[(1+(i-1)*n):(i*n)]
        Sigma0 = diag(s(x_i)) %*% Sigma0root %*% diag(s(x_i))
        epsilon_i = expm::sqrtm(Sigma0) %*% as.vector(rnorm(n))
        y[(1+(i-1)*n):(i*n)] = x_i + epsilon_i
      }
      rdf = data.frame(x=x,y=y,id=rep(1:I,each=n))
    }

    betas_gee=append(betas_gee,   summary(geepack::geeglm(y~x-1, data=rdf, id=id, corstr="exchangeable"))$coefficients[1,1]   )
    betas_mle=append(betas_mle,   summary(nlme::gls(y~x-1, data=rdf, correlation=corCompSymm(form=~1|id)))$coefficients[1]   )
    betas_unw=append(betas_unw,   summary(nlme::gls(y~x-1, data=rdf, correlation=corCompSymm(form=~1|id, value=0, fixed=TRUE)))$coefficients[1]   )

    find_params = function(params) clubSandwich::vcovCR(  nlme::gls(y~x-1, data=rdf, correlation=corCompSymm(form=~1|id, value=c(params),fixed=TRUE)), cluster=rdf$id, type="CR3"  )[1,1]
    params_here = optim(c(0), find_params, method="Brent", lower=-0.24, upper=0.99)
    betas_sl=append(betas_sl,   summary(nlme::gls(y~x, data=rdf, correlation=corCompSymm(form=~1|id, value=params_here$par, fixed=TRUE)))$coefficients[2]   )

    find_params_naive = function(params) {
      beta_unweighted = summary(nlme::gls(y~x-1, data=rdf, correlation=corCompSymm(form=~1|id, value=0, fixed=TRUE)))$coefficients
      resids = rdf$y - beta_unweighted * rdf$x
      XWX = XWSWX = 0
      for (i in 1:I) {
        x_i = rdf$x[(1+(i-1)*n):(i*n)]
        W_i = solve( matrix(params,n,n) + diag(1-params,n) )
        r_i = resids[(1+(i-1)*n):(i*n)]
        S_i = r_i %*% t(r_i)
        XWX = XWX + t(x_i) %*% W_i %*% t(t(x_i))
        XWSWX = XWSWX + t(x_i) %*% W_i %*% S_i %*% W_i %*% t(t(x_i))
      }
      XWSWX / (XWX)^2
    }
    params_here_naive = optim(c(0), find_params_naive, method="Brent", lower=-0.24, upper=0.99)
    betas_sl_naive=append(betas_sl_naive,   summary(nlme::gls(y~x, data=rdf, correlation=corCompSymm(form=~1|id, value=params_here_naive$par, fixed=TRUE)))$coefficients[2]   )

  }
  #
  betas_unw_ALL[[arc]] = betas_unw
  betas_gee_ALL[[arc]] = betas_gee
  betas_mle_ALL[[arc]] = betas_mle
  betas_sl_ALL[[arc]] = betas_sl
  betas_sl_naive_ALL[[arc]] = betas_sl_naive
  
  MSE_unw[arc] = mean((betas_unw-1)^2)
  MSE_gee[arc] = mean((betas_gee-1)^2)
  MSE_mle[arc] = mean((betas_mle-1)^2)
  MSE_sl[arc] = mean((betas_sl-1)^2)
  MSE_sl_naive[arc] = mean((betas_sl_naive-1)^2)
  
}


# Binomial simulation (well-specified covariance)
corstructure="equicorr"
betas_unw_ALL = betas_gee_ALL = betas_mle_ALL = betas_sl_ALL = betas_sl_naive_ALL = list()
Is = c(5,10,seq(20,200,by=20));
MSE_unw = MSE_gee = MSE_mle = MSE_sl = MSE_sl_naive = rep(0,length(Is))
n=20
arc = 0
MMM=10000
for (I in Is) {
  arc=arc+1
  betas_unw = betas_gee = betas_mle = betas_sl = betas_sl_naive = numeric()
  SEEDj = 0
  sqrtCov1 = expm::sqrtm(matrix(0.5,n,n)+diag(0.5,n))
  sqrtZcov = expm::sqrtm(matrix(0.6,n,n)+diag(0.4,n))
  while (length(betas_gee)<MMM) {
    SEEDj = SEEDj + 1
    {
      x = rep(0,n*I)
      for (i in 1:I) {
        x[(1+(i-1)*n):(i*n)] = sqrtCov1 %*% as.vector(rnorm(n))
      }
      y = rep(0,n*I)
      for (i in 1:I) {
        x_i = x[(1+(i-1)*n):(i*n)]
        z_i = sqrtZcov %*% as.vector(rnorm(n))
        U_i = pnorm(z_i)
        Y_i = as.numeric( U_i >= 1-exp(x_i)/(1+exp(x_i)) )
        y[(1+(i-1)*n):(i*n)] = Y_i
      }
      rdf = data.frame(x=x,y=y,id=rep(1:I,each=n))
    }

    betas_gee=append(betas_gee,   summary(geepack::geeglm(y~x-1, data=rdf, family=binomial("logit"), id=id, corstr="exchangeable"))$coefficients[1,1]   )

    {
      cor.fixed = matrix(0,n,n) + diag(1-0,n)
      zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)
      unw.mod <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
      betas_unw = append(betas_unw, summary(unw.mod)$coefficients[1,1] )
    }

    {
      gam.mle.finder <- function(gam) {

        if (corstructure == "equicorr") cor.fixed = matrix(gam,n,n) + diag(1-gam,n)
        if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
        zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)

        mod0 <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
        resids <- residuals(mod0, "pearson")
        MLE.obj <- 0
        phi <- 0
        if (corstructure == "equicorr") R_i <- matrix(gam,n,n) + diag(1-gam,n)
        if (corstructure == "ar1") R_i <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
        R_i_inv <- solve(R_i)
        for (i in 1:I) {
          r_i <- resids[rdf$id==i]
          phi <- phi + t(r_i)%*%R_i_inv%*%r_i
        }
        phi <- phi/(dim(rdf)[1])
        for (i in 1:I) {
          r_i <- resids[rdf$id==i]
          MLE.obj <- MLE.obj + (1/phi) * t(r_i)%*%R_i_inv%*%t(t(r_i)) + log(det(R_i))
        }
        MLE.obj.solve.for.zero = MLE.obj
        MLE.obj.solve.for.zero^2
      }
      if (corstructure == "equicorr")  gam.mle = optim(0, gam.mle.finder, method="Brent", lower=-0.001, upper=0.99)$par
      if (corstructure == "ar1")  gam.mle = optim(0, gam.mle.finder, method="Brent", lower=-0.99, upper=0.99)$par
      if (corstructure == "equicorr") cor.fixed = matrix(gam.mle,n,n) + diag(1-gam.mle,n)
      if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=gam.mle,lag.max=n-1))
      zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)
      mle.mod <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
      betas_mle = append(betas_mle, summary(mle.mod)$coefficients[1,1] )
    }

    gam.sl.finder <- function(gam) {
      X = model.matrix( ~ x-1, data = rdf)
      if (corstructure == "equicorr") cor.fixed <- matrix(gam,n,n) + diag(1-gam,n)
      if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=c(gam),lag.max=n-1))
      zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)

      mod0 <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
      beta_hat_gam <- summary(mod0)$coefficients[1,1]
      resids <- residuals(mod0, "pearson")
      SL.obj <- 0
      g_inv = function(x) log(x/(1-x))
      g_dif = function(mu) 1/(mu*(1-mu))
      v = function(mu) mu*(1-mu)

      X_is = A_is = A_is.sqrt = DWD_is = M_is = S_is = B_is = SBS_is = list()
      if (corstructure == "equicorr") {
        R_i <- matrix(gam,n,n) + diag(1-gam,n)
        R_i_inv <- solve(R_i)
        R_i_deriv <- matrix(1,n,n) - diag(1,n)
      }
      if (corstructure == "ar1") {
        R_i <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
        R_i_inv <- solve(R_i)
        R_i_deriv <- matrix(0,n,n)
        for (j in 1:n) {
          for (k in 1:n) {
            if (j!=k) R_i_deriv[j,k] = abs(j-k)*(gam)^(abs(j-k)-1)
          }
        }
      }

      for (i in 1:I) {
        X_i = X[ rdf$id==i , ]
        X_is[[i]] = X_i
        A_is[[i]] = diag(v(mod0$fitted.values[rdf$id==i]))
        A_is.sqrt[[i]] = diag(sqrt(v(mod0$fitted.values[rdf$id==i])))
        DWD_is[[i]] <- t(X_is[[i]]) %*% A_is.sqrt[[i]] %*% R_i_inv %*% A_is.sqrt[[i]] %*% (X_is[[i]])
      }
      DWD <- Reduce("+", DWD_is)
      DWD_inv <- solve(DWD)
      for (i in 1:I) {
        S_is[[i]]  <-  solve(DWD - DWD_is[[i]])
        resids_i <- residuals(mod0, "pearson")
        B_is[[i]]  <-  t(X_is[[i]])%*%A_is.sqrt[[i]]%*%R_i_inv%*%((resids_i[rdf$id==i])%*%t(resids_i[rdf$id==i]))%*%R_i_inv%*%A_is.sqrt[[i]]%*%X_is[[i]]
        SBS_is[[i]] <- S_is[[i]] %*% B_is[[i]] %*% S_is[[i]]
      }
      SL <- Reduce("+", SBS_is)
      SL[1,1]
    }
    if (corstructure == "equicorr")  gam.sl = optim(0, gam.sl.finder, method="Brent", lower=-0.01, upper=0.99)$par
    if (corstructure == "ar1")  gam.sl = optim(0, gam.sl.finder, method="Brent", lower=-0.99, upper=0.99)$par
    if (corstructure == "equicorr") cor.fixed = matrix(gam.sl,n,n) + diag(1-gam.mle,n)
    if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=gam.sl,lag.max=n-1))
    zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)
    sl.mod <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
    betas_sl = append(betas_sl, summary(sl.mod)$coefficients[1,1] )


    find_params_naive <- function(gam) {
      X = model.matrix( ~ x-1, data = rdf)

      beta_unweighted <- summary(unw.mod)$coefficients[1,1]
      resids <- residuals(unw.mod, "pearson")

      v = function(mu) mu*(1-mu)

      X_is = A_is = A_is.sqrt = DWD_is = M_is = S_is = B_is = SBS_is = list()
      if (corstructure == "equicorr") {
        R_i <- matrix(gam,n,n) + diag(1-gam,n)
        R_i_inv <- solve(R_i)
        R_i_deriv <- matrix(1,n,n) - diag(1,n)
      }
      if (corstructure == "ar1") {
        R_i <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
        R_i_inv <- solve(R_i)
        R_i_deriv <- matrix(0,n,n)
        for (j in 1:n) {
          for (k in 1:n) {
            if (j!=k) R_i_deriv[j,k] = abs(j-k)*(gam)^(abs(j-k)-1)
          }
        }
      }

      for (i in 1:I) {
        X_i = X[ rdf$id==i , ]
        X_is[[i]] = X_i
        A_is[[i]] = diag(v(unw.mod$fitted.values[rdf$id==i]))
        A_is.sqrt[[i]] = diag(sqrt(v(unw.mod$fitted.values[rdf$id==i])))
        DWD_is[[i]] <- t(X_is[[i]]) %*% A_is.sqrt[[i]] %*% R_i_inv %*% A_is.sqrt[[i]] %*% (X_is[[i]])
      }
      DWD <- Reduce("+", DWD_is)
      DWD_inv <- solve(DWD)
      for (i in 1:I) {
        S_is[[i]]  <-  solve(DWD) # Naive version
        resids_i <- residuals(unw.mod, "pearson")
        B_is[[i]]  <-  t(X_is[[i]])%*%A_is.sqrt[[i]]%*%R_i_inv%*%((resids_i[rdf$id==i])%*%t(resids_i[rdf$id==i]))%*%R_i_inv%*%A_is.sqrt[[i]]%*%X_is[[i]]
        SBS_is[[i]] <- S_is[[i]] %*% B_is[[i]] %*% S_is[[i]]
      }
      SL <- Reduce("+", SBS_is)
      SL[1,1]
    }

    if (corstructure == "equicorr") params_here_naive = optim(c(0), find_params_naive, method="Brent", lower=-0.24, upper=0.99)$par
    if (corstructure == "ar1") params_here_naive = optim(c(0), find_params_naive, method="Brent", lower=-0.99, upper=0.99)$par

    if (corstructure == "equicorr") cor.fixed = matrix(params_here_naive,n,n) + diag(1-params_here_naive,n)
    if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=params_here_naive,lag.max=n-1))
    zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)
    sl.naive.mod <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
    betas_sl_naive = append(betas_sl_naive, summary(sl.naive.mod)$coefficients[1,1] )

  }
  #
  betas_unw_ALL[[arc]] = betas_unw
  betas_gee_ALL[[arc]] = betas_gee
  betas_mle_ALL[[arc]] = betas_mle
  betas_sl_ALL[[arc]] = betas_sl
  betas_sl_naive_ALL[[arc]] = betas_sl_naive
  #
  MSE_unw[arc] = mean((betas_unw-1)^2)
  MSE_gee[arc] = mean((betas_gee-1)^2)
  MSE_mle[arc] = mean((betas_mle-1)^2)
  MSE_sl[arc] = mean((betas_sl-1)^2)
  MSE_sl_naive[arc] = mean((betas_sl_naive-1)^2)
  #
}


# Binomial simulation (misspecified covariance)
corstructure="ar1"
betas_unw_ALL = betas_gee_ALL = betas_mle_ALL = betas_sl_ALL = betas_sl_naive_ALL = list()
Is = c(5,10,seq(20,200,by=20));
MSE_unw = MSE_gee = MSE_mle = MSE_sl = MSE_sl_naive = rep(0,length(Is))
n=20
arc = 0
MMM = 400
sqrtCov1 = expm::sqrtm(matrix(0.5,n,n)+diag(0.5,n))
sqrtZcov = expm::sqrtm( toeplitz(ARMAacf(ar=c(0.4,0.5), ma=c(-0.9,0.4), lag.max=n-1)) )
for (I in Is) {
  arc=arc+1
  betas_unw = betas_gee = betas_mle = betas_sl = betas_sl_naive = numeric()
  SEEDj = 0
  while (length(betas_gee)<MMM) {
    SEEDj = SEEDj + 1
    {
      x = rep(0,n*I)
      for (i in 1:I) {
        x[(1+(i-1)*n):(i*n)] = sqrtCov1 %*% as.vector(rnorm(n))
      }
      y = rep(0,n*I)
      for (i in 1:I) {
        x_i = x[(1+(i-1)*n):(i*n)]
        z_i = sqrtZcov %*% as.vector(rnorm(n))
        U_i = pnorm(z_i)
        Y_i = as.numeric( U_i >= 1-exp(x_i)/(1+exp(x_i)) )
        y[(1+(i-1)*n):(i*n)] = Y_i
      }
      rdf = data.frame(x=x,y=y,id=rep(1:I,each=n))
    }

    betas_gee=append(betas_gee,   summary(geepack::geeglm(y~x-1, data=rdf, family=binomial("logit"), id=id, corstr="ar1"))$coefficients[1,1]   )

    {
      cor.fixed = matrix(0,n,n) + diag(1-0,n)
      zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)
      unw.mod <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
      betas_unw = append(betas_unw, summary(unw.mod)$coefficients[1,1] )
    }

    {
      gam.mle.finder <- function(gam) {

        if (corstructure == "equicorr") cor.fixed = matrix(gam,n,n) + diag(1-gam,n)
        if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
        zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)

        mod0 <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
        resids <- residuals(mod0, "pearson")
        MLE.obj <- 0
        phi <- 0
        if (corstructure == "equicorr") R_i <- matrix(gam,n,n) + diag(1-gam,n)
        if (corstructure == "ar1") R_i <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
        R_i_inv <- solve(R_i)
        for (i in 1:I) {
          r_i <- resids[rdf$id==i]
          phi <- phi + t(r_i)%*%R_i_inv%*%r_i
        }
        phi <- phi/(dim(rdf)[1])
        for (i in 1:I) {
          r_i <- resids[rdf$id==i]
          MLE.obj <- MLE.obj + (1/phi) * t(r_i)%*%R_i_inv%*%t(t(r_i)) + log(det(R_i))
        }
        MLE.obj.solve.for.zero = MLE.obj
        MLE.obj.solve.for.zero^2
      }
      if (corstructure == "equicorr")  gam.mle = optim(0, gam.mle.finder, method="Brent", lower=-0.001, upper=0.99)$par
      if (corstructure == "ar1")  gam.mle = optim(0, gam.mle.finder, method="Brent", lower=-0.99, upper=0.99)$par
      if (corstructure == "equicorr") cor.fixed = matrix(gam.mle,n,n) + diag(1-gam.mle,n)
      if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=gam.mle,lag.max=n-1))
      zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)
      mle.mod <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
      betas_mle = append(betas_mle, summary(mle.mod)$coefficients[1,1] )
    }

    gam.sl.finder <- function(gam) {
      X = model.matrix( ~ x-1, data = rdf)
      if (corstructure == "equicorr") cor.fixed <- matrix(gam,n,n) + diag(1-gam,n)
      if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=c(gam),lag.max=n-1))
      zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)

      mod0 <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
      beta_hat_gam <- summary(mod0)$coefficients[1,1]
      resids <- residuals(mod0, "pearson")
      SL.obj <- 0
      g_inv = function(x) log(x/(1-x))
      g_dif = function(mu) 1/(mu*(1-mu))
      v = function(mu) mu*(1-mu)

      X_is = A_is = A_is.sqrt = DWD_is = M_is = S_is = B_is = SBS_is = list()
      if (corstructure == "equicorr") {
        R_i <- matrix(gam,n,n) + diag(1-gam,n)
        R_i_inv <- solve(R_i)
        R_i_deriv <- matrix(1,n,n) - diag(1,n)
      }
      if (corstructure == "ar1") {
        R_i <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
        R_i_inv <- solve(R_i)
        R_i_deriv <- matrix(0,n,n)
        for (j in 1:n) {
          for (k in 1:n) {
            if (j!=k) R_i_deriv[j,k] = abs(j-k)*(gam)^(abs(j-k)-1)
          }
        }
      }

      for (i in 1:I) {
        X_i = X[ rdf$id==i , ]
        X_is[[i]] = X_i
        A_is[[i]] = diag(v(mod0$fitted.values[rdf$id==i]))
        A_is.sqrt[[i]] = diag(sqrt(v(mod0$fitted.values[rdf$id==i])))
        DWD_is[[i]] <- t(X_is[[i]]) %*% A_is.sqrt[[i]] %*% R_i_inv %*% A_is.sqrt[[i]] %*% (X_is[[i]])
      }
      DWD <- Reduce("+", DWD_is)
      DWD_inv <- solve(DWD)
      for (i in 1:I) {
        S_is[[i]]  <-  solve(DWD - DWD_is[[i]])
        resids_i <- residuals(mod0, "pearson")
        B_is[[i]]  <-  t(X_is[[i]])%*%A_is.sqrt[[i]]%*%R_i_inv%*%((resids_i[rdf$id==i])%*%t(resids_i[rdf$id==i]))%*%R_i_inv%*%A_is.sqrt[[i]]%*%X_is[[i]]
        SBS_is[[i]] <- S_is[[i]] %*% B_is[[i]] %*% S_is[[i]]
      }
      SL <- Reduce("+", SBS_is)
      SL[1,1]
    }
    if (corstructure == "equicorr")  gam.sl = optim(0, gam.sl.finder, method="Brent", lower=-0.01, upper=0.99)$par
    if (corstructure == "ar1")  gam.sl = optim(0, gam.sl.finder, method="Brent", lower=-0.99, upper=0.99)$par
    if (corstructure == "equicorr") cor.fixed = matrix(gam.sl,n,n) + diag(1-gam.mle,n)
    if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=gam.sl,lag.max=n-1))
    zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)
    sl.mod <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
    betas_sl = append(betas_sl, summary(sl.mod)$coefficients[1,1] )


    find_params_naive <- function(gam) {
      X = model.matrix( ~ x-1, data = rdf)

      beta_unweighted <- summary(unw.mod)$coefficients[1,1]
      resids <- residuals(unw.mod, "pearson")

      v = function(mu) mu*(1-mu)

      X_is = A_is = A_is.sqrt = DWD_is = M_is = S_is = B_is = SBS_is = list()
      if (corstructure == "equicorr") {
        R_i <- matrix(gam,n,n) + diag(1-gam,n)
        R_i_inv <- solve(R_i)
        R_i_deriv <- matrix(1,n,n) - diag(1,n)
      }
      if (corstructure == "ar1") {
        R_i <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
        R_i_inv <- solve(R_i)
        R_i_deriv <- matrix(0,n,n)
        for (j in 1:n) {
          for (k in 1:n) {
            if (j!=k) R_i_deriv[j,k] = abs(j-k)*(gam)^(abs(j-k)-1)
          }
        }
      }

      for (i in 1:I) {
        X_i = X[ rdf$id==i , ]
        X_is[[i]] = X_i
        A_is[[i]] = diag(v(unw.mod$fitted.values[rdf$id==i]))
        A_is.sqrt[[i]] = diag(sqrt(v(unw.mod$fitted.values[rdf$id==i])))
        DWD_is[[i]] <- t(X_is[[i]]) %*% A_is.sqrt[[i]] %*% R_i_inv %*% A_is.sqrt[[i]] %*% (X_is[[i]])
      }
      DWD <- Reduce("+", DWD_is)
      DWD_inv <- solve(DWD)
      for (i in 1:I) {
        S_is[[i]]  <-  solve(DWD) # Naive  version
        resids_i <- residuals(unw.mod, "pearson")
        B_is[[i]]  <-  t(X_is[[i]])%*%A_is.sqrt[[i]]%*%R_i_inv%*%((resids_i[rdf$id==i])%*%t(resids_i[rdf$id==i]))%*%R_i_inv%*%A_is.sqrt[[i]]%*%X_is[[i]]
        SBS_is[[i]] <- S_is[[i]] %*% B_is[[i]] %*% S_is[[i]]
      }
      SL <- Reduce("+", SBS_is)
      SL[1,1]
    }

    if (corstructure == "equicorr") params_here_naive = optim(c(0), find_params_naive, method="Brent", lower=-0.24, upper=0.99)$par
    if (corstructure == "ar1") params_here_naive = optim(c(0), find_params_naive, method="Brent", lower=-0.99, upper=0.99)$par

    if (corstructure == "equicorr") cor.fixed = matrix(params_here_naive,n,n) + diag(1-params_here_naive,n)
    if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=params_here_naive,lag.max=n-1))
    zcor <- rep(cor.fixed[lower.tri(cor.fixed)], I)
    sl.naive.mod <- geepack::geeglm(y ~ x-1, family = binomial("logit"), id=id, data=rdf, corstr="fixed", zcor=zcor)
    betas_sl_naive = append(betas_sl_naive, summary(sl.naive.mod)$coefficients[1,1] )

  }
  #
  betas_unw_ALL[[arc]] = betas_unw
  betas_gee_ALL[[arc]] = betas_gee
  betas_mle_ALL[[arc]] = betas_mle
  betas_sl_ALL[[arc]] = betas_sl
  betas_sl_naive_ALL[[arc]] = betas_sl_naive
  #
  MSE_unw[arc] = mean((betas_unw-1)^2)
  MSE_gee[arc] = mean((betas_gee-1)^2)
  MSE_mle[arc] = mean((betas_mle-1)^2)
  MSE_sl[arc] = mean((betas_sl-1)^2)
  MSE_sl_naive[arc] = mean((betas_sl_naive-1)^2)
  #
}




# Generate plots
{
  # For ggplots
  library(RColorBrewer)
  coul <- brewer.pal(3, "Dark2")
  coul_rgb <- col2rgb(coul)
  library(ggplot2)
  library(gridExtra)
  library(ggpubr)
  library(grid)
  library(rlang)
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  cairo_pdf("mse_sims.pdf", 9, 9)

  TRIM=0
  # Linear mixed effects model well specified
  {
    BETAS=readRDS("~/Downloads/sandwichreg_simulations/Datasets/para_proj_linear_mem_wellspec.Rdata")
    Is = c(5,10,seq(20,200,by=20))
    MSE_unw = MSE_gee = MSE_mle = MSE_sl = MSE_sl_naive = numeric(12)
    for (pp in 1:12) {
      MSE_unw[pp] = mean((BETAS[[pp]]$unw-1)^2,trim=TRIM)
      MSE_gee[pp] = mean((BETAS[[pp]]$gee-1)^2,trim=TRIM)
      MSE_mle[pp] = mean((BETAS[[pp]]$mle-1)^2,trim=TRIM)
      MSE_sl[pp] = mean((BETAS[[pp]]$sl-1)^2,trim=TRIM)
      MSE_sl_naive[pp] = mean((BETAS[[pp]]$sl_naive-1)^2,trim=TRIM)
    }
    MSE_gee_vs_unw = MSE_gee/MSE_unw
    MSE_mle_vs_unw = MSE_mle/MSE_unw
    MSE_sl_vs_unw = MSE_sl/MSE_unw
    MSE_sl_naive_vs_unw = MSE_sl_naive/MSE_unw

    MSE.DF.on.unw <- data.frame(Estimators=c(rep("GEE",length(Is)),rep("QML",length(Is)),rep("SL (finite sample)",length(Is)),rep("SL (large sample)",length(Is))), NumberOfClusters=rep(Is,4), MSE=c(MSE_gee_vs_unw,MSE_mle_vs_unw,MSE_sl_vs_unw,MSE_sl_naive_vs_unw))

    JITTER = 0
    plot_MEM_wellspec <- ggplot(MSE.DF.on.unw, aes(x=NumberOfClusters, y=MSE, group=Estimators, color=Estimators, linetype=Estimators)) +
      geom_line(position=position_dodge(width = JITTER)) +
      geom_point(position=position_dodge(width = JITTER)) +
      theme_bw() +
      labs(x = "Number of Clusters", y = "Mean Squared Error") +
      scale_color_manual(values=c(coul[1], coul[3], coul[2], coul[2])) +
      scale_linetype_manual(values=c(1,1,1,2)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(expression("Linear (mixed effects) well specified model:" ~ lambda==0)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(colour = guide_legend(override.aes = list(shape = NA)))
  }
  # Linear mixed effects model misspecified
  {
      BETAS=readRDS("~/Downloads/sandwichreg_simulations/Datasets/para_proj_linear_mem_misspec.Rdata")
    Is = c(5,10,seq(20,200,by=20))
    MSE_unw = MSE_gee = MSE_mle = MSE_sl = MSE_sl_naive = numeric(12)
    for (pp in 1:12) {
      MSE_unw[pp] = mean((BETAS[[pp]]$unw-1)^2,trim=TRIM)
      MSE_gee[pp] = mean((BETAS[[pp]]$gee-1)^2,trim=TRIM)
      MSE_mle[pp] = mean((BETAS[[pp]]$mle-1)^2,trim=TRIM)
      MSE_sl[pp] = mean((BETAS[[pp]]$sl-1)^2,trim=TRIM)
      MSE_sl_naive[pp] = mean((BETAS[[pp]]$sl_naive-1)^2,trim=TRIM)
    }
    MSE_gee_vs_unw = MSE_gee/MSE_unw
    MSE_mle_vs_unw = MSE_mle/MSE_unw
    MSE_sl_vs_unw = MSE_sl/MSE_unw
    MSE_sl_naive_vs_unw = MSE_sl_naive/MSE_unw

    MSE.DF.on.unw <- data.frame(Estimators=c(rep("GEE",length(Is)),rep("QML",length(Is)),rep("SL (finite sample)",length(Is)),rep("SL (large sample)",length(Is))), NumberOfClusters=rep(Is,4), MSE=c(MSE_gee_vs_unw,MSE_mle_vs_unw,MSE_sl_vs_unw,MSE_sl_naive_vs_unw))

    JITTER = 0
    plot_MEM_misspec <- ggplot(MSE.DF.on.unw, aes(x=NumberOfClusters, y=MSE, group=Estimators, color=Estimators, linetype=Estimators)) +
      geom_line(position=position_dodge(width = JITTER)) +
      geom_point(position=position_dodge(width = JITTER)) +
      theme_bw() +
      labs(x = "Number of Clusters", y = "Mean Squared Error") +
      scale_color_manual(values=c(coul[1], coul[3], coul[2], coul[2])) +
      scale_linetype_manual(values=c(1,1,1,2)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(expression("Linear (mixed effects) misspecified model:"~lambda==3)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_hline(yintercept=1,linetype=3) +
      guides(colour = guide_legend(override.aes = list(shape = NA)))
  }
  TRIM=0.01
  # Binomial well specified covariance
  {
    BETAS=readRDS("~/Downloads/sandwichreg_simulations/Datasets/para_proj_binom_wellpec.Rdata")
    Is = c(5,10,seq(20,200,by=20))
    MSE_unw = MSE_gee = MSE_mle = MSE_sl = MSE_sl_naive = numeric(12)
    for (pp in 1:12) {
      MSE_unw[pp] = mean((BETAS[[pp]]$unw-1)^2,trim=TRIM)
      MSE_gee[pp] = mean((BETAS[[pp]]$gee-1)^2,trim=TRIM)
      MSE_mle[pp] = mean((BETAS[[pp]]$mle-1)^2,trim=TRIM)
      MSE_sl[pp] = mean((BETAS[[pp]]$sl-1)^2,trim=TRIM)
      MSE_sl_naive[pp] = mean((BETAS[[pp]]$sl_naive-1)^2,trim=TRIM)
    }
    MSE_gee_vs_unw = MSE_gee/MSE_unw
    MSE_mle_vs_unw = MSE_mle/MSE_unw
    MSE_sl_vs_unw = MSE_sl/MSE_unw
    MSE_sl_naive_vs_unw = MSE_sl_naive/MSE_unw

    MSE.DF.on.unw <- data.frame(Estimators=c(rep("GEE",length(Is)),rep("QML",length(Is)),rep("SL (finite sample)",length(Is)),rep("SL (large sample)",length(Is))), NumberOfClusters=rep(Is,4), MSE=c(MSE_gee_vs_unw,MSE_mle_vs_unw,MSE_sl_vs_unw,MSE_sl_naive_vs_unw))

    JITTER = 0
    plot_binom_wellspec <- ggplot(MSE.DF.on.unw, aes(x=NumberOfClusters, y=MSE, group=Estimators, color=Estimators, linetype=Estimators)) +
      geom_line(position=position_dodge(width = JITTER)) +
      geom_point(position=position_dodge(width = JITTER)) +
      theme_bw() +
      labs(x = "Number of Clusters", y = "Mean Squared Error") +
      scale_color_manual(values=c(coul[1], coul[3], coul[2], coul[2])) +
      scale_linetype_manual(values=c(1,1,1,2)) +
      #scale_x_continuous(breaks = Is) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(expression("Binomial well specified covariance model")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(colour = guide_legend(override.aes = list(shape = NA)))
  }
  # Binomial misspecified covariance
  {
    BETAS=readRDS("~/Downloads/sandwichreg_simulations/Datasets/para_proj_binom_misspec.Rdata")
    Is = c(5,10,seq(20,200,by=20))
    MSE_unw = MSE_gee = MSE_mle = MSE_sl = MSE_sl_naive = numeric(12)
    for (pp in 1:12) {
      MSE_unw[pp] = mean((BETAS[[pp]]$unw-1)^2,trim=TRIM)
      MSE_gee[pp] = mean((BETAS[[pp]]$gee-1)^2,trim=TRIM)
      MSE_mle[pp] = mean((BETAS[[pp]]$mle-1)^2,trim=TRIM)
      MSE_sl[pp] = mean((BETAS[[pp]]$sl-1)^2,trim=TRIM)
      MSE_sl_naive[pp] = mean((BETAS[[pp]]$sl_naive-1)^2,trim=TRIM)
    }
    MSE_gee_vs_unw = MSE_gee/MSE_unw
    MSE_mle_vs_unw = MSE_mle/MSE_unw
    MSE_sl_vs_unw = MSE_sl/MSE_unw
    MSE_sl_naive_vs_unw = MSE_sl_naive/MSE_unw

    MSE.DF.on.unw <- data.frame(Estimators=c(rep("GEE",length(Is)),rep("QML",length(Is)),rep("SL (finite sample)",length(Is)),rep("SL (large sample)",length(Is))), NumberOfClusters=rep(Is,4), MSE=c(MSE_gee_vs_unw,MSE_mle_vs_unw,MSE_sl_vs_unw,MSE_sl_naive_vs_unw))

    JITTER = 0
    plot_binom_misspec <- ggplot(MSE.DF.on.unw, aes(x=NumberOfClusters, y=MSE, group=Estimators, color=Estimators, linetype=Estimators)) +
      geom_line(position=position_dodge(width = JITTER)) +
      geom_point(position=position_dodge(width = JITTER)) +
      theme_bw() +
      labs(x = "Number of Clusters", y = "Mean Squared Error") +
      scale_color_manual(values=c(coul[1], coul[3], coul[2], coul[2])) +
      scale_linetype_manual(values=c(1,1,1,2)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(expression("Binomial misspecified  covariance model")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_hline(yintercept=1,linetype=3) +
      guides(colour = guide_legend(override.aes = list(shape = NA)))
  }
  # Plot to generate legend
  {
    legend_netadata = MSE.DF.on.unw
    sub.df = MSE.DF.on.unw[MSE.DF.on.unw$Estimators=="GEE",]
    sub.df[,1] = rep("Unweighted",length(Is))
    legend_netadata = rbind(legend_netadata, sub.df)
    legend_plot = ggplot(legend_netadata, aes(x=NumberOfClusters, y=MSE, group=Estimators, color=Estimators, linetype=Estimators)) +
      geom_line(position=position_dodge(width = 0.001)) +
      geom_point(position=position_dodge(width = 0.001)) +
      theme_bw() +
      labs(x = expression(lambda), y = "Ratio of MSE over Oracle") +
      scale_color_manual(values=c(coul[1],coul[3],coul[2],coul[2],"black"), labels=c("GEE","QML","SL (finite sample)","SL (large sample)","Unweighted"), name="Estimators") +
      scale_linetype_manual(values=c(1,1,1,2,3)) +
      scale_y_continuous(trans='log2') +
      theme(legend.text=element_text(size=10)) +
      theme(legend.position = "right",legend.box = "horizontal") +
      guides(colour = guide_legend(override.aes = list(shape = NA), nrow = 2))
    legend_plot
  }

  plot_MEM_wellspec
  plot_MEM_misspec
  legend = get_legend(legend_plot)
  ggarrange(plot_MEM_wellspec, plot_MEM_misspec, plot_binom_wellspec, plot_binom_misspec, ncol=2, nrow=2, common.legend = TRUE, legend="bottom", legend.grob = legend)

  dev.off()
}
