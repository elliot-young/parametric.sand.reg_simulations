# Reproducibility code for analysis of 1988 Bangladesh Fertility Survey (Section 4.2.2) in "Sandwich Regression: Accurate and robust inference in generalised linear multilevel and longitudinal models". 

library(epiDisplay)
data(Bang)
Bang_mine = Bang
Bang_mine$living.children = as.factor(Bang_mine$living.children)
levels(Bang_mine$living.children)=c("none","1 child", "2 children", ">2 children")
Bang_mine$district[Bang_mine$district>=55] = Bang_mine$district[Bang_mine$district>=55]-1
corstructure = "equicorr"

age_vals = seq(-12,12,by=0.1)
y.sl.vs.unw = y.gee.vs.unw = y.mle.vs.unw = numeric()
for (x in age_vals) {
    #print(paste0("Progress: Simulation for age = ", x))
    pred_vec = c(1,1,0,0,x,x^2,1)
    gam.sl.finder <- function(gam) {
      X = model.matrix( ~ living.children + age_mean + I(age_mean^2) + urban, data = Bang_mine)

      I = length(unique(Bang_mine$district))
      n_is = rep(0,I)
      for (i in 1:I) n_is[i] = dim(Bang_mine[Bang_mine$district==i,])[1]

      zcor <- numeric()
      for (i in 1:I) {
        n_i = n_is[i]
        if (corstructure == "equicorr") cor.fixed <- matrix(gam,n_i,n_i) + diag(1-gam,n_i)
        if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=c(gam),lag.max=n_i-1))
        zcor <- append( zcor, cor.fixed[lower.tri(cor.fixed)] )
      }
      mod0 <- geepack::geeglm(user ~ living.children + age_mean + I(age_mean^2) + urban, family = binomial("logit"), id=district, data=Bang_mine, corstr="fixed", zcor=zcor)
      resids <- residuals(mod0, "pearson")
      SL.obj <- 0
      g_inv = function(x) log(x/(1-x))
      g_dif = function(mu) 1/(mu*(1-mu))
      v = function(mu) mu*(1-mu)

      X_is = A_is = A_is.sqrt = DWD_is = M_is = S_is = B_is = SBS_is = list()
      for (i in 1:I) {
        n = n_is[i]
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
        X_i = X[ Bang_mine$district==i , ]
        X_is[[i]] = X_i
        A_is[[i]] = diag(v(mod0$fitted.values[Bang_mine$district==i]))
        A_is.sqrt[[i]] = diag(sqrt(v(mod0$fitted.values[Bang_mine$district==i])))
        DWD_is[[i]] <- t(X_is[[i]]) %*% A_is.sqrt[[i]] %*% R_i_inv %*% A_is.sqrt[[i]] %*% (X_is[[i]])
      }
      DWD <- Reduce("+", DWD_is)
      DWD_inv <- solve(DWD)
      for (i in 1:I) {
        n = n_is[i]
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
        S_is[[i]]  <-  solve(DWD - DWD_is[[i]])
        resids_i <- residuals(mod0, "pearson")
        B_is[[i]]  <-  t(X_is[[i]])%*%A_is.sqrt[[i]]%*%R_i_inv%*%((resids_i[Bang_mine$district==i])%*%t(resids_i[Bang_mine$district==i]))%*%R_i_inv%*%A_is.sqrt[[i]]%*%X_is[[i]]
        SBS_is[[i]] <- S_is[[i]] %*% B_is[[i]] %*% S_is[[i]]
      }
      SL <- Reduce("+", SBS_is)
      t(pred_vec) %*% SL %*% pred_vec
    }
    gam.sl = optim(0, gam.sl.finder, method="Brent", lower=-0.95, upper=0.95)$par
    sl.eval = gam.sl.finder(gam.sl)
    gee.eval = gam.sl.finder(as.numeric(summary( geepack::geeglm(user ~ living.children + age_mean + I(age_mean^2) + urban, family = binomial("logit"), id = district, data = Bang_mine, corstr = "exchangeable") )$corr[1]))
    unw.eval = gam.sl.finder(0)

    # mle method
    {
      gam.mle.finder <- function(gam) {
        I = length(unique(Bang_mine$district))
        n_is = rep(0,I)
        for (i in 1:I) n_is[i] = dim(Bang_mine[Bang_mine$district==i,])[1]

        zcor <- numeric()
        for (i in 1:I) {
          n_i = n_is[i]
          if (corstructure == "equicorr") cor.fixed <- matrix(0,n_i,n_i) + diag(1-0,n_i)
          if (corstructure == "ar1") cor.fixed <- toeplitz(ARMAacf(ar=c(0),lag.max=n_i-1))
          zcor <- append( zcor, cor.fixed[lower.tri(cor.fixed)] )
        }
        mod0 <- geepack::geeglm(user ~ living.children + age_mean + I(age_mean^2) + urban, family = binomial("logit"), id=district, data=Bang_mine, corstr="fixed", zcor=zcor)
        resids <- residuals(mod0, "pearson")
        MLE.obj <- 0
        phi <- 0
        for (i in 1:I) {
          n = n_is[i]
          if (corstructure == "equicorr") R_i <- matrix(gam,n,n) + diag(1-gam,n)
          if (corstructure == "ar1") R_i <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
          R_i_inv <- solve(R_i)
          r_i <- resids[Bang_mine$district==i]
          phi <- phi + t(r_i)%*%r_i
        }
        phi <- phi/(dim(Bang_mine)[1])
        for (i in 1:I) {
          n = n_is[i]
          if (corstructure == "equicorr") R_i <- matrix(gam,n,n) + diag(1-gam,n)
          if (corstructure == "ar1") R_i <- toeplitz(ARMAacf(ar=gam,lag.max=n-1))
          R_i_inv <- solve(R_i)
          r_i <- resids[Bang_mine$district==i]
          MLE.obj <- MLE.obj + (1/phi) * t(r_i)%*%R_i_inv%*%t(t(r_i)) + log(det(R_i)) #<- FIXED (not phi is independent of gam)
        }
        MLE.obj
      }
      gam.mle = optim(0, gam.mle.finder, method="Brent", lower=-0.001, upper=0.99)$par
      mle.eval = gam.sl.finder(gam.mle)
    }

    y.sl.vs.unw = append(y.sl.vs.unw, 100-100*sl.eval/unw.eval)
    y.gee.vs.unw = append(y.gee.vs.unw, 100-100*gee.eval/unw.eval)
    y.mle.vs.unw = append(y.mle.vs.unw, 100-100*mle.eval/unw.eval)
}


library(ggplot2)
coul <- RColorBrewer::brewer.pal(3, "Dark2")
data = data.frame(
  x = rep(age_vals, 4),
  y = c(rep(1,length(age_vals)), 1-y.sl.vs.unw/100, 1-y.gee.vs.unw/100, 1-y.mle.vs.unw/100),
  estimator = factor(rep(c("Unweighted", "Sandwich Regression", "GEE", "QML"), each = length(age_vals)))
)
cairo_pdf("bang.pdf", 6, 4)
ggplot(data, aes(x=x, y=y, color=estimator, linetype=estimator)) +
  geom_line(aes(size = estimator)) +
  scale_size_manual(values = c("Unweighted"=0.4, "Sandwich Regression"=0.8, "GEE"=0.8, "QML"=0.8)) +
  scale_color_manual(values = c("Unweighted"="black", "QML"=coul[3], "GEE"=coul[1], "Sandwich Regression"=coul[2])) +
  scale_linetype_manual(values = c("Unweighted"="dashed", "QML"="solid", "GEE"="solid", "Sandwich Regression"="solid")) +
  labs(x = "Age", y = "Variance (relative to unweighted estimator)") +
  coord_cartesian(xlim = c(-11,11)) +
  theme_minimal() +
  theme(legend.title = element_blank())
dev.off()
