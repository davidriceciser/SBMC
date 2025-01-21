#Creating 200 500 by 500 matricies





dauls_method<-function(myMCMC, niter, burn_in){

  n<-length(myMCMC[[1]]$zout)
  daul_array<-array(data = NA, dim = c(n,n,niter - burn_in))


  #Calculating array
  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:(niter - burn_in )){
        daul_array[i,j,k]<-myMCMC[[k + burn_in]]$zout[i]==myMCMC[[k + burn_in]]$zout[j]
      }
    }
  }

  #Calculating Abar

  mean_matrix<-matrix(data=  0, nrow=n  , ncol =n)
  for(i in 1:(niter - burn_in)){
    mean_matrix<-mean_matrix + daul_array[,,i]
  }
  mean_matrix<-mean_matrix/(niter-burn_in)
  #Finding the one that is closest (just want index)
  ####################################################


  #Identifying closest one

  dist_vec<-c()


  for(i in 1:(niter-burn_in)){
    dist_vec<-c(dist_vec,sum((mean_matrix - daul_array[,,i])^2))

  }
  close_index<-which(dist_vec == min(dist_vec))[1]

  return(close_index + burn_in)

}


