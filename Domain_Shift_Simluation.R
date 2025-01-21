#install.packages("plot.matrix")
library(plot.matrix)
#install.packages("abind")
library(abind)
source("MFM_covariance_edited.R")
source("U_covariance_update.R")
source("Dauls_method.R")

mean1 = matrix(data=  0, nrow = 16, ncol = 16)
mean2 = matrix(data=  0, nrow = 16, ncol = 16)
mean3 = matrix(data=  0, nrow = 16, ncol = 16)

mean1[1:4,1:4] = 1
plot(mean1)
mean2[8:14,8:14] = 1
plot(mean2)
mean3[6:16,6:16] = 1
plot(mean3)


#Creating matrix observations function
create_matrix_observations<-function(signal_strength, U_true, V_true, observations_each,M_array){
  rows<-dim(U_true)[1]
  cols<-dim(V_true)[1]


  Cluster_start<-array(NA, dim = c(rows,cols,1))
  #Clusters<-array(data = NA, dim = c(10,6,sum(observations_each)))
  for(i in 1:length(observations_each)){
    Cluster_new<-mniw::rMNorm(observations_each[i],
                              Lambda = M_array[,,i]*signal_strength,

                              SigmaR = U_true, SigmaC = V_true)
    #plot(M[,,i])
    Cluster_start<-abind(Cluster_start,Cluster_new,along = 3)
  }
  Clusters<-Cluster_start[,,-1]
  return(Clusters)

}
#install.packages("MixMatrix")
library(MixMatrix)
U_true = ARgenerate(16,rho = 0.1)
V_true=  ARgenerate(16,rho = 0.1)


M_array = array(data = NA, dim = c(16,16,3))
M_array[1:16,1:16,1] = mean1
M_array[1:16,1:16,2] = mean2
M_array[1:16,1:16,3] = mean3


obs = create_matrix_observations(signal_strength = 1, U_true = U_true, V_true = V_true, observations_each = c(100, 50, 10), M_array = M_array)


gamma = 1
cluster_num = 20
GAMMA = 1
MCMC.total = 1000


domain_shift = c()
for(i in 1:10){
  #Independent prior
  obs = create_matrix_observations(signal_strength = 1.5, U_true = U_true, V_true = V_true, observations_each = c(100, 50, 10), M_array = M_array)
  domain_shift[i]<-MFM_MxN_equal_cov_sparse(data=obs,
                                                      niterations=MCMC.total,
                                                      alpha=(1+dim(obs)[1])/2, beta=diag(rep(0.5,dim(obs)[1])),
                                                      psi=(1+dim(obs)[2])/2, rho=diag(rep(0.5,dim(obs)[2])),
                                                      GAMMA=GAMMA,
                                                      M0=(apply(obs, c(1,2), max) + apply(obs, c(1,2), min))/2,
                                                      Sigma0=diag( ((apply(obs, c(1), max)-apply(obs, c(1), min))/4)^2 ),
                                                      Omega0=diag( ((apply(obs, c(2), max)-apply(obs, c(2), min))/4)^2 ),
                                                      initNClusters=cluster_num,
                                                      MLE.initial=T)
}







#############################################################

generate_statistics<-function(clustering_object,niter,nreps,theoretical_groups,daul_method = F){

  # clustering_object = SBMC_lowdim_lowsig
  # niter = 200
  # nreps = 50
  # theoretical_groups =theoretical_groups
  # daul_method = T
  #
  #

  #Generate randint values
  if(daul_method == T){
    randint_vec<-c()

    for(i in 1:nreps){
      #clustering_object = SBMC_lowdim_lowsig
      #niter = 200
      #burn_in = floor(niter/2)
      daul_num<-dauls_method(clustering_object[[i]],niter = niter, burn_in = floor(niter/2))
      randint_vec[i]<-fossil::adj.rand.index(clustering_object[[i]][[daul_num]]$zout,theoretical_groups )
      print(daul_num)
    }
  }
  ret_vals<-data.frame("Mean adj randint" = mean(randint_vec), "SD" =sd(randint_vec))
  return(ret_vals)
}


theoretical_groups = rep(c(1,2,3), times = c(100, 50, 10))
generate_statistics(clustering_object = domain_shift, niter = 1000,nreps = 10,theoretical_groups = theoretical_groups,daul_method = T)






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



SBMC_highdim_lowsig_500<-generate_statistics(clustering_object = data,niter = 1000,
                                             nreps = 10, theoretical_groups =c(100, 50, 10),daul_method = T)

#############################################################
#saveRDS(domain_shift,"domain_shift.rds")
domain_shift = readRDS("domain_shift.rds")

par(mfrow = c(1,3))
plot(mean1,main = "Cluster 1",key = NULL)
plot(mean2, main = "Cluster 2", key = NULL)
plot(mean3, main = "Cluster 3", key = NULL)


par(mfrow = c(1,1))
plot(domain_shift[[5]][501][[1]]$Mout[,,3])



par(mfrow = c(1,3))
plot(domain_shift[[5]][501][[1]]$Mout[,,3],main = "Cluster 1",key = NULL)
plot(domain_shift[[5]][501][[1]]$Mout[,,2],main = "Cluster 2", key = NULL)
plot(domain_shift[[5]][501][[1]]$Mout[,,1],main = "Cluster 3", key = NULL)



