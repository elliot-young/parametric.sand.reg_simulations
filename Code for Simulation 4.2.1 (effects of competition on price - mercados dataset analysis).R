# Reproducibility code for analysis of effect of competition on market price (Section 4.2.1) in "Sandwich Regression: Accurate and robust inference in generalised linear multilevel and longitudinal models".

{
  # Download dataset from https://www.openicpsr.org/openicpsr/project/113651/version/V1/view
  DAT1=as.data.frame(haven::read_dta("~/Downloads/113651-V1/data/1.-Data/Consumers_panel.dta"))
  DAT1=DAT1[,c("round","log_P_s","Q_index","id","RAS_time","mercado","hogar_id","Hgt0","Hfgt0","M_col_pre_treatment","province_1","province_2","province_3","province_4","province_5","province_6","province_7","province_8","province_9","M_total_benef_pret","M_percent_ben","HH_income_dis","DD_educ1","DD_educ2","DD_educ3","urbano")]
  for (id in unique(DAT1$hogar_id)) {
    if( length(DAT1$hogar_id[DAT1$hogar_id==id]) != 2 )   DAT1 = DAT1[DAT1$hogar_id!=id,]
  }
  DAT1 = cbind(L.log_P_s=rep(0,length(DAT1$round)) , L.Q_index=rep(0,length(DAT1$round)),DAT1)
  for (id in unique(DAT1$hogar_id)) {
    DAT1_id = DAT1[DAT1$hogar_id==id,]
    DAT1_id[2,"L.log_P_s"] = DAT1_id[1,"log_P_s"]
    DAT1_id[2,"L.Q_index"] = DAT1_id[1,"Q_index"]
    DAT1[DAT1$hogar_id==id,] = DAT1_id
  }
  DAT1=DAT1[DAT1$round==2,]
  DAT1=na.omit(DAT1)
  DAT1=subset(DAT1,select=-c(round))
  TARGET = 2
  X = cbind(rep(1,length(DAT1$mercado)), subset(DAT1,select=-c(Hfgt0,log_P_s,Q_index,id,mercado,hogar_id,province_9,DD_educ3,RAS_time))  )
  X = X[,c(1,4,2,3,(5:dim(X)[2]))]
  Y = as.vector(DAT1$log_P_s)
}

library(nlme)
vecc=rep(0,19); vecc[2]=1
#############################
{Q=3
print(Q)
find_params = function(params) {
  L = matrix(0,Q,Q)
  L[lower.tri(L,diag=TRUE)] = params
  LLt = L%*%t(L)

  grouping_vec = rep(0, length(unique(DAT1$mercado)))
  XWX_all = XWY_all = list()
  q_pos = 0
  for (q in unique(DAT1$mercado)){
    q_pos = q_pos + 1
    grouping_vec[q_pos] = dim(DAT1[DAT1$mercado==q,])[1]
    X_q = X[DAT1$mercado==q,]
    Y_q = Y[DAT1$mercado==q]
    Z_q = as.matrix(cbind(X_q[,c("ones")]))
    Qq = 0
    while (dim(Z_q)[2] < Q) {
      Qq=Qq+1
      Z_q = cbind(Z_q, (DAT1$L.log_P_s[DAT1$mercado==q])^Qq)
    }

    W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
    XWX_all[[q_pos]] = t(X_q) %*% W %*% t(t(X_q))
    XWY_all[[q_pos]] = t(X_q) %*% W %*% (Y_q)
  }
  XWX = Reduce("+",XWX_all)
  XWY = Reduce("+",XWY_all)
  beta_hat = solve(XWX) %*% (XWY)

  V = matrix(0,dim(X)[2],dim(X)[2])
  q_pos = 0
  for (q in unique(DAT1$mercado)){
    q_pos = q_pos + 1
    X_q = X[DAT1$mercado==q,]
    Y_q = Y[DAT1$mercado==q]
    Z_q = as.matrix(cbind(X_q[,c("ones")]))
    Qq = 0
    while (dim(Z_q)[2] < Q) {
      Qq=Qq+1
      Z_q = cbind(Z_q, (DAT1$L.log_P_s[DAT1$mercado==q])^Qq)
    }
    
    W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
    A = XWX - t(X_q) %*% W %*% t(t(X_q))
    B = t(X_q) %*% W %*% (Y_q-t(t(X_q))%*%beta_hat) %*% t((Y_q-t(t(X_q))%*%beta_hat)) %*% W %*% t(t(X_q))
    V = V + solve(A) %*% (B) %*% solve(A)
  }

  t(vecc)%*%V%*%t(t(vecc))

}

AA = lme4::lmer( log_P_s ~ Hgt0 + L.log_P_s + L.Q_index + M_col_pre_treatment + province_1 + province_2 + province_3 + province_4 + province_5 + province_6 + province_7 + province_8 + M_total_benef_pret + M_percent_ben + HH_income_dis + DD_educ1 + DD_educ2 + urbano + (1+L.log_P_s+I(L.log_P_s^2)|mercado), data = DAT1 )
summary(AA)
LLtt = as.matrix(VarCorr(AA)$mercado[1:Q,1:Q])
LLtt = LLtt/(sigma(AA)^2)
if (min(eigen(LLtt)$values)<=0) LLtt = LLtt+diag(1e-5,Q,Q)
mle_pars = t(chol(LLtt))
mle_pars_vec = as.vector(mle_pars)[as.vector(mle_pars)!=0]

parameters_here = optim(mle_pars_vec, find_params, method = "Nelder-Mead")

parameters_here_all = list()
i=0
for (qq in unique(DAT1$mercado)) {
  i=i+1
  X_mini = X[DAT1$mercado!=qq,]
  Y_mini = Y[DAT1$mercado!=qq]
  DAT1_mini = DAT1[DAT1$mercado!=qq,]
  find_params = function(params) {
    L = matrix(0,Q,Q)
    L[lower.tri(L,diag=TRUE)] = params
    LLt = L%*%t(L)

    grouping_vec = rep(0, length(unique(DAT1_mini$mercado)))
    XWX_all = XWY_all = list()
    q_pos = 0
    for (q in unique(DAT1_mini$mercado)){
      q_pos = q_pos + 1
      grouping_vec[q_pos] = dim(DAT1_mini[DAT1_mini$mercado==q,])[1]
      X_q = X_mini[DAT1_mini$mercado==q,]
      Y_q = Y_mini[DAT1_mini$mercado==q]
      Z_q = as.matrix(cbind(X_q[,c("ones")]))
      Qq = 0
      while (dim(Z_q)[2] < Q) {
        Qq=Qq+1
        Z_q = cbind(Z_q, (DAT1_mini$L.log_P_s[DAT1_mini$mercado==q])^Qq)
      }

      W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
      XWX_all[[q_pos]] = t(X_q) %*% W %*% t(t(X_q))
      XWY_all[[q_pos]] = t(X_q) %*% W %*% (Y_q)
    }
    XWX = Reduce("+",XWX_all)
    XWY = Reduce("+",XWY_all)
    beta_hat = solve(XWX) %*% (XWY)

    V = matrix(0,dim(X_mini)[2],dim(X_mini)[2])
    q_pos = 0
    for (q in unique(DAT1_mini$mercado)){
      q_pos = q_pos + 1
      X_q = X_mini[DAT1_mini$mercado==q,]
      Y_q = Y_mini[DAT1_mini$mercado==q]
      Z_q = as.matrix(cbind(X_q[,c("ones")]))
      Qq = 0
      while (dim(Z_q)[2] < Q) {
        Qq=Qq+1
        Z_q = cbind(Z_q, (DAT1_mini$L.log_P_s[DAT1_mini$mercado==q])^Qq)
      }

      W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
      A = XWX - t(X_q) %*% W %*% t(t(X_q))
      B = t(X_q) %*% W %*% (Y_q-t(t(X_q))%*%beta_hat) %*% t((Y_q-t(t(X_q))%*%beta_hat)) %*% W %*% t(t(X_q))
      V = V + solve(A) %*% (B) %*% solve(A)
    }

    t(vecc)%*%V%*%t(t(vecc))

  }

  parameters_here_all[[i]] = nlm(find_params, parameters_here$par, iterlim=ITERLIM)$estimate
  print(i)
}
beta_hat_all = list()
i=0
for (qq in unique(DAT1$mercado)) {
  i=i+1
  X_mini = X[DAT1$mercado!=qq,]
  Y_mini = Y[DAT1$mercado!=qq]
  DAT1_mini = DAT1[DAT1$mercado!=qq,]
  {
    L = matrix(0,Q,Q)
    L[lower.tri(L,diag=TRUE)] = parameters_here_all[[i]]
    LLt = L%*%t(L)

    grouping_vec = rep(0, length(unique(DAT1_mini$mercado)))
    XWX_all = XWY_all = list()
    q_pos = 0
    for (q in unique(DAT1_mini$mercado)){
      q_pos = q_pos + 1
      grouping_vec[q_pos] = dim(DAT1_mini[DAT1_mini$mercado==q,])[1]
      X_q = X_mini[DAT1_mini$mercado==q,]
      Y_q = Y_mini[DAT1_mini$mercado==q]
      Z_q = as.matrix(cbind(X_q[,c("ones")]))
      Qq = 0
      while (dim(Z_q)[2] < Q) {
        Qq=Qq+1
        Z_q = cbind(Z_q, (DAT1_mini$L.log_P_s[DAT1_mini$mercado==q])^Qq)
      }
      
      W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
      XWX_all[[q_pos]] = t(X_q) %*% W %*% t(t(X_q))
      XWY_all[[q_pos]] = t(X_q) %*% W %*% (Y_q)
    }
    XWX = Reduce("+",XWX_all)
    XWY = Reduce("+",XWY_all)
    beta_hat_all[[i]] = solve(XWX) %*% (XWY)
  }
}
beta_5_all=rep(0,length(unique(DAT1$mercado)))
i=0
for (qq in unique(DAT1$mercado)) {
  i=i+1
  beta_5_all[i] = t(vecc)%*%beta_hat_all[[i]]
}

V_all_AV = (length(unique(DAT1$mercado))-1)/(length(unique(DAT1$mercado)))*sum((beta_5_all-mean(beta_5_all))^2)
print(paste0("AV var: ", V_all_AV))
print(paste0("AV est: ", mean(beta_5_all)))
}

{
  beta_hat_all = list()
  i=0
  for (qq in unique(DAT1$mercado)) {
    i=i+1
    X_mini = X[DAT1$mercado!=qq,]
    Y_mini = Y[DAT1$mercado!=qq]
    DAT1_mini = DAT1[DAT1$mercado!=qq,]
    beta_hat_all[[i]] = coefficients(nlme::gls( log_P_s ~ Hgt0 + L.log_P_s + L.Q_index + M_col_pre_treatment + province_1 + province_2 + province_3 + province_4 + province_5 + province_6 + province_7 + province_8 + M_total_benef_pret + M_percent_ben + HH_income_dis + DD_educ1 + DD_educ2 + urbano, correlation=corAR1(value=0,form=~1|mercado,fixed=TRUE), data = DAT1_mini )
    )
  }
  beta_5_all=rep(0,length(unique(DAT1$mercado)))
  i=0
  for (qq in unique(DAT1$mercado)) {
    i=i+1
    beta_5_all[i] = t(vecc)%*%beta_hat_all[[i]]
  }
  V_all_AV = (length(unique(DAT1$mercado))-1)/(length(unique(DAT1$mercado)))*sum((beta_5_all-mean(beta_5_all))^2)
  print(paste0("Unweighted var: ",V_all_AV))
  print(paste0("Unweighted est: ",mean(beta_5_all)))
}
{
  beta_hat_all = list()
  i=0
  for (qq in unique(DAT1$mercado)) {
    i=i+1
    X_mini = X[DAT1$mercado!=qq,]
    Y_mini = Y[DAT1$mercado!=qq]
    DAT1_mini = DAT1[DAT1$mercado!=qq,]
    beta_hat_all[[i]] = coefficients(nlme::gls( log_P_s ~ Hgt0 + L.log_P_s + L.Q_index + M_col_pre_treatment + province_1 + province_2 + province_3 + province_4 + province_5 + province_6 + province_7 + province_8 + M_total_benef_pret + M_percent_ben + HH_income_dis + DD_educ1 + DD_educ2 + urbano, correlation=corAR1(value=0,form=~1|mercado,fixed=TRUE), data = DAT1_mini )
    )
  }
  beta_5_all=rep(0,length(unique(DAT1$mercado)))
  i=0
  for (qq in unique(DAT1$mercado)) {
    i=i+1
    beta_5_all[i] = t(vecc)%*%beta_hat_all[[i]]
  }
  V_all_AV = (length(unique(DAT1$mercado))-1)/(length(unique(DAT1$mercado)))*sum((beta_5_all-mean(beta_5_all))^2)
  print(paste0("Unweighted var: ",V_all_AV))
  print(paste0("Unweighted est: ",mean(beta_5_all)))
}
{Q=3
  print(Q)
  find_params_gee = function(params) {
    sigg = params[(1+Q*(Q+1)/2)]

    L = matrix(0,Q,Q)
    L[lower.tri(L,diag=TRUE)] = params[1:(Q*(Q+1)/2)]
    LLt = L%*%t(L)

    grouping_vec = rep(0, length(unique(DAT1$mercado)))
    XWX_all = XWY_all = list()
    q_pos = 0
    for (q in unique(DAT1$mercado)){
      q_pos = q_pos + 1
      grouping_vec[q_pos] = dim(DAT1[DAT1$mercado==q,])[1]
      X_q = X[DAT1$mercado==q,]
      Y_q = Y[DAT1$mercado==q]
      Z_q = as.matrix(cbind(X_q[,c("ones")]))
      Qq = 0
      while (dim(Z_q)[2] < Q) {
        Qq=Qq+1
        Z_q = cbind(Z_q, (DAT1$L.log_P_s[DAT1$mercado==q])^Qq)
      }

      W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
      XWX_all[[q_pos]] = t(X_q) %*% W %*% t(t(X_q))
      XWY_all[[q_pos]] = t(X_q) %*% W %*% (Y_q)
    }
    XWX = Reduce("+",XWX_all)
    XWY = Reduce("+",XWY_all)
    beta_hat = solve(XWX) %*% (XWY)

    Rsq=0
    q_pos = 0
    for (q in unique(DAT1$mercado)){
      q_pos = q_pos + 1
      X_q = X[DAT1$mercado==q,]
      Y_q = Y[DAT1$mercado==q]
      Z_q = as.matrix(cbind(X_q[,c("ones")]))
      Qq = 0
      while (dim(Z_q)[2] < Q) {
        Qq=Qq+1
        Z_q = cbind(Z_q, (DAT1$L.log_P_s[DAT1$mercado==q])^Qq)
      }
      
      W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
      mat_vals =  (Y_q-t(t(X_q))%*%beta_hat) %*% t((Y_q-t(t(X_q))%*%beta_hat)) - sigg*solve(W)
      b = sum(mat_vals^2)
      Rsq = Rsq + b
    }

    Rsq
  }

  AA = lme4::lmer( log_P_s ~ Hgt0 + L.log_P_s + L.Q_index + M_col_pre_treatment + province_1 + province_2 + province_3 + province_4 + province_5 + province_6 + province_7 + province_8 + M_total_benef_pret + M_percent_ben + HH_income_dis + DD_educ1 + DD_educ2 + urbano + (1+L.log_P_s+I(L.log_P_s^2)|mercado), data = DAT1 )
  summary(AA)
  LLtt = as.matrix(VarCorr(AA)$mercado[1:Q,1:Q])
  LLtt = LLtt/(sigma(AA)^2)
  if (min(eigen(LLtt)$values)<0) LLtt = LLtt+diag(1e-5,Q,Q)
  mle_pars = t(chol(LLtt))
  mle_pars_vec = as.vector(mle_pars)[as.vector(mle_pars)!=0]

  parameters_here = optim(c(rep(0,length(mle_pars_vec)),1), find_params_gee, method = "Nelder-Mead")

  parameters_here_all = list()
  i=0
  for (qq in unique(DAT1$mercado)) {
    i=i+1
    X_mini = X[DAT1$mercado!=qq,]
    Y_mini = Y[DAT1$mercado!=qq]
    DAT1_mini = DAT1[DAT1$mercado!=qq,]
    find_params = function(params) {
      sigg = params[(1+Q*(Q+1)/2)]

      L = matrix(0,Q,Q)
      L[lower.tri(L,diag=TRUE)] = params[1:(Q*(Q+1)/2)]
      LLt = L%*%t(L)

      grouping_vec = rep(0, length(unique(DAT1_mini$mercado)))
      XWX_all = XWY_all = list()
      q_pos = 0
      for (q in unique(DAT1_mini$mercado)){
        q_pos = q_pos + 1
        grouping_vec[q_pos] = dim(DAT1_mini[DAT1_mini$mercado==q,])[1]
        X_q = X_mini[DAT1_mini$mercado==q,]
        Y_q = Y_mini[DAT1_mini$mercado==q]
        Z_q = as.matrix(cbind(X_q[,c("ones")]))
        Qq = 0
        while (dim(Z_q)[2] < Q) {
          Qq=Qq+1
          Z_q = cbind(Z_q, (DAT1_mini$L.log_P_s[DAT1_mini$mercado==q])^Qq)
        }

        W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
        XWX_all[[q_pos]] = t(X_q) %*% W %*% t(t(X_q))
        XWY_all[[q_pos]] = t(X_q) %*% W %*% (Y_q)
      }
      XWX = Reduce("+",XWX_all)
      XWY = Reduce("+",XWY_all)
      beta_hat = solve(XWX) %*% (XWY)

      Rsq = 0
      q_pos = 0
      for (q in unique(DAT1_mini$mercado)){
        q_pos = q_pos + 1
        X_q = X_mini[DAT1_mini$mercado==q,]
        Y_q = Y_mini[DAT1_mini$mercado==q]
        Z_q = as.matrix(cbind(X_q[,c("ones")]))
        Qq = 0
        while (dim(Z_q)[2] < Q) {
          Qq=Qq+1
          Z_q = cbind(Z_q, (DAT1_mini$L.log_P_s[DAT1_mini$mercado==q])^Qq)
        }
        
        W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
        mat_vals =  (Y_q-t(t(X_q))%*%beta_hat) %*% t((Y_q-t(t(X_q))%*%beta_hat)) - sigg*solve(W)
        b = sum(mat_vals^2)
        Rsq = Rsq + b
      }

      Rsq
    }

    parameters_here_all[[i]] = nlm(find_params, parameters_here$par, iterlim=ITERLIM)$estimate
    print(i)
  }
  beta_hat_all = list()
  i=0
  for (qq in unique(DAT1$mercado)) {
    i=i+1
    X_mini = X[DAT1$mercado!=qq,]
    Y_mini = Y[DAT1$mercado!=qq]
    DAT1_mini = DAT1[DAT1$mercado!=qq,]
    {
      L = matrix(0,Q,Q)
      L[lower.tri(L,diag=TRUE)] = parameters_here_all[[i]]
      LLt = L%*%t(L)

      grouping_vec = rep(0, length(unique(DAT1_mini$mercado)))
      XWX_all = XWY_all = list()
      q_pos = 0
      for (q in unique(DAT1_mini$mercado)){
        q_pos = q_pos + 1
        grouping_vec[q_pos] = dim(DAT1_mini[DAT1_mini$mercado==q,])[1]
        X_q = X_mini[DAT1_mini$mercado==q,]
        Y_q = Y_mini[DAT1_mini$mercado==q]
        Z_q = as.matrix(cbind(X_q[,c("ones")]))
        Qq = 0
        while (dim(Z_q)[2] < Q) {
          Qq=Qq+1
          Z_q = cbind(Z_q, (DAT1_mini$L.log_P_s[DAT1_mini$mercado==q])^Qq)
        }

        W = solve(  diag(1,length(Y_q)) + Z_q%*%LLt%*%t(Z_q)  )
        XWX_all[[q_pos]] = t(X_q) %*% W %*% t(t(X_q))
        XWY_all[[q_pos]] = t(X_q) %*% W %*% (Y_q)
      }
      XWX = Reduce("+",XWX_all)
      XWY = Reduce("+",XWY_all)
      beta_hat_all[[i]] = solve(XWX) %*% (XWY)
    }
  }
  beta_5_all=rep(0,length(unique(DAT1$mercado)))
  i=0
  for (qq in unique(DAT1$mercado)) {
    i=i+1
    beta_5_all[i] = t(vecc)%*%beta_hat_all[[i]]
  }

  V_all_AV = (length(unique(DAT1$mercado))-1)/(length(unique(DAT1$mercado)))*sum((beta_5_all-mean(beta_5_all))^2)
  print(paste0("GEE: var",V_all_AV))
  print(paste0("GEE: est",mean(beta_5_all)))
}
