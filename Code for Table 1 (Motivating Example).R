# Reproducibility code for the results of Table 1 in "Sandwich Regression: Accurate and robust inference in generalised linear multilevel and longitudinal models"

library(nlme)
MSE_beta = AIC = BIC = matrix(NA,3,3)
I = 10^6/50
{
  set.seed(1)
  n=50
  Sigma0 = matrix(1,n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i!=j) Sigma0[i,j] = exp(-(abs(i-j))^0.25/2) + 0.25*cos(abs(i-j))*exp(-abs(i-j)/20)
    }
  }
  sqrtSigma0 = expm::sqrtm(Sigma0)
  sqrtCovX = expm::sqrtm(matrix(0.9,n,n) + diag(0.1,n,n))
  beta = as.vector(c(1))
  Xs = matrix(data=NA, nrow=n*I, ncol=1)
  Ys = eps = numeric(0)
  for (i in seq_len(I)) {
    X = matrix(NA, nrow=n, ncol=1)
    X[,1] = sqrtCovX %*% as.vector(rnorm(n))
    epsilon = sqrtSigma0 %*% as.vector(rnorm(n))
    Y = X %*% beta + epsilon
    Xs[(n*(i-1)+1):(n*i), ] = X
    Ys = append(Ys, Y)
    eps = append(eps, epsilon)
  }
  X = Xs
  Y = Ys
  ep = eps
  rdf = data.frame(X=X,Y=Y,id=rep(1:I,each=n))
}
for (P in 0:2) {
  for (Q in 0:2) {
    if (P+Q>0) {
      print(paste0("(P,Q) = (",P,",",Q,")"))
      find_params = function(params) {
        grouping_vec = rep(0, length(unique(rdf$id)))
        XWX_all = XWY_all = list()
        q_pos = 0
        for (q in unique(rdf$id)){
          q_pos = q_pos + 1
          grouping_vec[q_pos] = dim(rdf[rdf$id==q,])[1]
          X_q = X[rdf$id==q,]
          Y_q = Y[rdf$id==q]
          if (P==0) W = solve(  toeplitz(ARMAacf(ar=c(),ma=c(params),lag.max=n-1))  )
          if (Q==0) W = solve(  toeplitz(ARMAacf(ar=c(params),ma=c(),lag.max=n-1))  )
          if (Q>0 & P>0) W = solve(  toeplitz(ARMAacf(ar=params[1:P],ma=c(params[(P+1):(P+Q)]),lag.max=n-1))  )
          XWX_all[[q_pos]] = t(X_q) %*% W %*% t(t(X_q))
          XWY_all[[q_pos]] = t(X_q) %*% W %*% (Y_q)
        }
        XWX = Reduce("+",XWX_all)
        XWY = Reduce("+",XWY_all)
        beta_hat = solve(XWX) %*% (XWY)
        XWX_inv = solve(XWX)

        B_all=0
        q_pos = 0
        for (q in unique(rdf$id)){
          q_pos = q_pos + 1
          X_q = X[rdf$id==q,]
          Y_q = Y[rdf$id==q]
          if (P==0) W = solve(  toeplitz(ARMAacf(ar=c(),ma=c(params),lag.max=n-1))  )
          if (Q==0) W = solve(  toeplitz(ARMAacf(ar=c(params),ma=c(),lag.max=n-1))  )
          if (Q>0 & P>0) W = solve(  toeplitz(ARMAacf(ar=params[1:P],ma=c(params[(P+1):(P+Q)]),lag.max=n-1))  )
          B = t(X_q) %*% W %*% (Y_q-t(t(X_q))%*%beta_hat) %*% t((Y_q-t(t(X_q))%*%beta_hat)) %*% W %*% X_q
          B_all = B_all + B
        }
        SL = XWX_inv %*% B_all %*% XWX_inv
        return(SL)
      }

      AA = nlme::gls(Y~X-1, correlation=corARMA(value=rep(0,P+Q),p=P,q=Q,form=~1|id), data = rdf )
      mle_pars_vec = coef(AA$modelStruct$corStruct, unconstrained = FALSE)

      MSE_beta[(P+1),(Q+1)] = find_params(mle_pars_vec)
      AIC[(P+1),(Q+1)] = summary(AA)$AIC
      BIC[(P+1),(Q+1)] = summary(AA)$BIC
    }
  }
}

options(digits = 14)
print(MSE_beta*(n*I))
print(AIC/(n*I))
print(BIC/(n*I))

