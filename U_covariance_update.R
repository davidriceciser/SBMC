




U_update<-function(U, solveV,solveU,data,psi_mat_U, rho_mat_U,phi_mat_U,W_mat_U){

  p<-dim(data)[1]
  q<-dim(data)[2]
  n<-dim(data)[3]

  #Inititializing values
  a<-b<-1/2

  C1<-10^4
  tausqU<-C1/(n*p^4)
  tausqV<-C1/(n*q^4)
  lambda = 1

  #Intitalizing return U
  ret_U<-U

  SV<-matrix(0,nrow = p,ncol = p)

  #Complexity = n*q^3
  for(i in 1:n){
    SV<-SV+(data[,,i])%*%solveV%*%t(data[,,i])
  }


  #Begin Column updates
  #Complexity = p*p^3 + p*p*[complexity of invSig_U11u1<-(SolveU11)%*%(u1)
  #p^2? + np^3


  #solveU[off_diag_U,off_diag_U]<-SolveU11+invSig_U11u1%*%t(invSig_U11u1)/v1]]
  for(j in 1:p){



    SolveU11<-solveU[-j,-j]-solveU[j,-j] %*% t(solveU[j,-j])/solveU[j,j]

    #Partition of SV
    SV11<-SV[-j,-j]
    SV12<-SV[j,-j]


    ##########################
    #Updating gamma (v)
    ##########################

    chiU1<-t(U[j,-j])%*%SolveU11%*%SV11%*%SolveU11%*%U[j,-j]
    chiU2<- (-2)*SV12%*%SolveU11%*%U[j,-j] + SV[j,j]
    chiU<-chiU1+chiU2

    #Full conditional go gamma given the others
    #set.seed(3)
    v1<-GIGrvg::rgig(n=1, lambda = (1-n*q/2),chi = chiU,psi = lambda )




    ########################
    #updating nu (u)
    ########################
    w1_u<-SolveU11%*%SV11%*%SolveU11
    w_u<-w1_u/v1 + diag(1/W_mat_U[j,-j])+lambda * SolveU11


    w_u_chol<-chol(w_u)
    w_u_chol_i<-solve(w_u_chol)
    mu_u_i<-(w_u_chol_i)%*%(t(w_u_chol_i)%*% SolveU11%*%SV12)/v1
    #set.seed(3)
    u1<-mu_u_i + w_u_chol_i%*%matrix(stats::rnorm(p-1))

    ret_U[-j,j]<-u1
    ret_U[j,-j]<-u1
    ret_U[j,j]<-v1 + t(u1)%*%SolveU11%*%(u1)



  off_diag_U<-(1:p)[-j]

  for(k in off_diag_U){

    #Updating phi
    chi_v<-(U[j,k])^2/tausqU
    if(chi_v < 10^-6){
      chi_v<- 10^-6
    }
    #set.seed(3)
    phi_mat_U[k,j]<-rgig(1,lambda = a-1/2,psi = 2*psi_mat_U[k,j],chi = chi_v )
    phi_mat_U[j,k]<-phi_mat_U[k,j]

    #updating psi based on new phi
    #set.seed(3)
    psi_mat_U[k,j]<-rgamma(n = 1,a+b,phi_mat_U[k,j]+1)
    psi_mat_U[j,k]<-psi_mat_U[k,j]

    #updating rho based on new phi
    rho_mat_U[k,j]<-phi_mat_U[k,j]/(1+phi_mat_U[k,j])
    rho_mat_U[j,k]<-rho_mat_U[k,j]

    #updating W_mat_U based on new rho

    W_mat_U[k,j]<-((rho_mat_U[k,j]*tausqU)/(1-rho_mat_U[k,j])) ###srqt removed
    W_mat_U[j,k]<-W_mat_U[k,j]


    #Updating sigma inverse
    invSig_U11u1<-(SolveU11)%*%(u1)
    solveU[off_diag_U,off_diag_U]<-SolveU11+invSig_U11u1%*%t(invSig_U11u1)/v1
    solveU[off_diag_U,j]<- -invSig_U11u1/v1
    solveU[j,off_diag_U]<- -invSig_U11u1/v1
    solveU[j,j]<-1/v1

  }
  }

  retlist<-list(ret_U,solveU,phi_mat_U,rho_mat_U,psi_mat_U,W_mat_U)
  names(retlist)<-c("U","solveU", "phi_mat_U","rho_mat_U","psi_mat_U","W_mat_U")
  return(retlist)

}




V_update<-function(V, solveU,solveV,data,psi_mat_V, rho_mat_V,phi_mat_V,W_mat_V){

  p<-dim(data)[1]
  q<-dim(data)[2]
  n<-dim(data)[3]
  #Inititializing values
  a<-b<-1/2

  C1<-10^4
  tausqU<-C1/(n*p^4)
  tausqV<-C1/(n*q^4)
  lambda = 1

  #Intitalizing return V
  ret_V<-V

  SU<-matrix(0,nrow = q,ncol = q)


  for(i in 1:n){

    SU<- SU + t(data[,,i])%*%solveU%*%(data[,,i])
  }





  for(j in 1:q){
    #solveV<-solve(V)
    #solveU<-solve(U)
    #Partition of SV
    SU11<-SU[-j,-j]
    SU12<-SU[j,-j]


    #SolveV11<-solveV[-j,-j]

    #SolveV12<-solveV[j,-j]
    #try this
    #SolveV11<-solve(V[-j,-j])
    SolveV11<-solveV[-j,-j]-solveV[j,-j] %*% t(solveV[j,-j])/solveV[j,j]





    ##########################
    #Updating gamma (v)
    #Note: We are updating V, must use the SU11, SU12
    ##########################


    chiV1<-(V[j,-j])%*%SolveV11%*%SU11%*%SolveV11%*%V[j,-j]
    chiV2<- (-2)*SU12%*%SolveV11%*%V[j,-j] + SU[j,j]
    chiV<-chiV1+chiV2

    #Full conditional go gamma given the others
    v2<-rgig(n=1, lambda = (1-n*p/2),chi = chiV,psi = lambda )
    #This works

   ###########
    #update u
    #######

    #trying again
    w1_v<-SolveV11%*%SU11%*%SolveV11
    w_v<-w1_v/v2 + diag(1/W_mat_V[j,-j])+lambda * SolveV11

    w_v_chol<-chol(w_v)
    w_v_chol_i<-solve(w_v_chol)
    mu_v_i<-(w_v_chol_i)%*%(t(w_v_chol_i)%*% SolveV11%*%SU12)/v2
    #set.seed(3)
    u2<-mu_v_i + w_v_chol_i%*%matrix(stats::rnorm(q-1))






    #NOW WE UPDATE COVARIANCE MATRIX
    ret_V[-j,j]<-u2
    ret_V[j,-j]<-u2
    ret_V[j,j]<-v2 + t(u2)%*%SolveV11%*%(u2)

    #NOW WE UPDATE PHI, RHO WITHIN THIS LOOP
    #we need to update every phi except for the jth one,
    #in our loop

    #off diagonals of V

    #Updating the phi vector
    #####################################
    # m <- p*q
    # my_n <- 1
    # my_rho = 0.5
    # my_sigma = ARgenerate(n = m,rho = my_rho)
    # my_z <- mvrnorm(n,mu=rep(0, m),Sigma=sigma,empirical=T)

    # #Get cumulative distributions
    # my_u <- pnorm(my_z)
    # my_formatted_u = matrix(my_u, dim = (length(off_diag_V),q))



    #######################################


    off_diag_V<-(1:q)[-j]

    for(k in off_diag_V){

      #Updating phi
      chi_v<-(V[j,k])^2/tausqV
      if(chi_v < 10^-6){
        chi_v<- 10^-6
      }
      #set.seed(3)
      #Instead of sampling from Rgig, just do a marginal sample using these variables.
      phi_mat_V[k,j]<-rgig(1,lambda = a-1/2,psi = 2*psi_mat_V[k,j],chi = chi_v )
      phi_mat_V[j,k]<-phi_mat_V[k,j]

      #phi_mat_V[k,j] = HyperbolicDist::qgig(u[,k*p + j],Theta = c(a-1/2,2*psi_mat_V[k,j],chi_v))

      #How do we update this with everything at once?




      #updating psi based on new phi
      #set.seed(3)
      psi_mat_V[k,j]<-rgamma(n = 1,a+b,phi_mat_V[k,j]+1)
      psi_mat_V[j,k]<-psi_mat_V[k,j]

      #updating rho based on new phi
      rho_mat_V[k,j]<-phi_mat_V[k,j]/(1+phi_mat_V[k,j])
      rho_mat_V[j,k]<-rho_mat_V[k,j]

      #updating W_mat_V based on new rho

      W_mat_V[k,j]<-((rho_mat_V[k,j]*tausqV)/(1-rho_mat_V[k,j])) ##sqrt removed
      W_mat_V[j,k]<-W_mat_V[k,j]
    }



    invSig_V11u2<-(SolveV11)%*%(u2)
    solveV[off_diag_V,off_diag_V]<-SolveV11+invSig_V11u2%*%t(invSig_V11u2)/v2
    solveV[off_diag_V,j]<- -invSig_V11u2/v2
    solveV[j,off_diag_V]<- -invSig_V11u2/v2
    solveV[j,j]<-1/v2



  }
  retlist<-list(ret_V,solveV,phi_mat_V,rho_mat_V,psi_mat_V,W_mat_V)
  names(retlist)<-c("V","solveV", "phi_mat_V","rho_mat_V","psi_mat_V","W_mat_V")
  return(retlist)

}
















