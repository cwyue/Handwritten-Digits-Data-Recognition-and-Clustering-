##Read in data:
library(mvtnorm)
myData=read.csv("semeion.csv",header=FALSE)
myX=data.matrix(myData[,1:256])
# "myLabel" records the true value of each data point.
myLabel=data.matrix(apply(myData[,257:266],1,function(xx){
  return(which(xx=="1")-1)
}))
d=256
data_number=1593
Iteration=10
likelihood=matrix(0,Iteration,4)
AICs=array(0,dim=4)
final_gamma=array(0,dim=c(data_number,k_num,4))
final_mean=array(0,dim=c(k_num,d,4))
final_sigma=array(0,dim=c(d,d,k_num,4))
final_cluster=array(0,dim=data_number)

##1.Initialize 10 clusters
k_num=10
pre_cluster=kmeans(myX, k_num, nstart = 10, iter.max =10)
ini_label=pre_cluster$cluster  #initial assignment of data points

#2.Convergence
for(q in c(0,2,4,6)){
  #Initialize gamma, calculate the value of mean, pi,N
  gamma=matrix(0,data_number,10)
  for (i in 1:data_number){
    gamma[i,ini_label[i]]=1
  }
  N=apply(gamma,2,sum)
  mean=t(gamma)%*%myX
  for(i in 1:k_num){
    mean[i,]=mean[i,]/N[i]
  }
  Pi=N/data_number
  #calculate sigma
  sigma=array(0,dim=c(d,d,k_num)) 
  for(k in 1:k_num){
    Dec=matrix(0,d,d)
    for (i in 1:data_number){
      Dec=Dec+gamma[i,k]*(myX[i,]-mean[k,])%*%(t(myX[i,]-mean[k,]))
    }
    Dec=Dec/N[k]
    myEig=eigen(Dec,symmetric = TRUE)
    Eig_val=myEig$value
    var=sum(Eig_val[(q+1):d])/(d-q)
    if(q==0){
      sigma[,,k]=var*diag(d)
    }else{
      Vq=myEig$vectors[,1:q]
      diag_matrix=diag(q)
      for (i in 1:q){
        diag_matrix[i,i]=sqrt(Eig_val[i]-var)
      }
      Wq=Vq%*%diag_matrix
      sigma[,,k]=Wq%*%t(Wq)+var*diag(d)
    }
    
  }
  
  #iterations start
  for (iter_count in 1:Iteration){
    #E step: update gamma matrix
    prob=matrix(0,data_number,k_num)
    for (k in 1:k_num){
      prob[,k]=dmvnorm(myX,mean[k,],sigma[,,k],log = FALSE)*Pi[k]
    }
    
    for(i in 1:data_number){
      for(k in 1:k_num){
        gamma[i,k]=prob[i,k]/sum(prob[i,])
      }
      #Calculate log-likelihood
      likelihood[iter_count,q/2+1]=likelihood[iter_count,q/2+1]+log(sum(prob[i,]))
    }
    
    #M step
    N=apply(gamma,2,sum)
    mean=t(gamma)%*%myX
    for(i in 1:k_num){
      mean[i,]=mean[i,]/N[i]
    }
    Pi=N/data_number
    #calculate sigma
    sigma=array(0,dim=c(d,d,k_num)) 
    for(k in 1:k_num){
      Dec=matrix(0,d,d)
      for (i in 1:data_number){
        Dec=Dec+gamma[i,k]*(myX[i,]-mean[k,])%*%(t(myX[i,]-mean[k,]))
      }
      Dec=Dec/N[k]
      myEig=eigen(Dec,symmetric = TRUE)
      Eig_val=myEig$value
      var=sum(Eig_val[(q+1):d])/(d-q)
      if(q==0){
        sigma[,,k]=var*diag(d)
      }else{
        Vq=myEig$vectors[,1:q]
        diag_matrix=diag(q)
        for (i in 1:q){
          diag_matrix[i,i]=sqrt(Eig_val[i]-var)
        }
        Wq=Vq%*%diag_matrix
        sigma[,,k]=Wq%*%t(Wq)+var*diag(d)
      }
    }
  }
  
  #calculate AIC
  AICs[q/2+1]=-2*likelihood[Iteration,q/2+1]+2*(d*q+1-q*(q-1)/2)
  #record the final value of mean and sigma matrix
  final_mean[,,q/2+1]=mean
  final_sigma[,,,q/2+1]=sigma
  final_gamma[,,q/2+1]=gamma
}


# 3.choose optimial number of principle components q:
q_op=2*(which.min(AICs)-1)

#4.Visualization of Clusters
#Visualizaiton of Clusters
dev.new(width=7,height=3.5)
par(mai=c(0.05,0.05,0.05,0.05),mfrow=c(10,6))
for (k in 1:k_num){
  image(t(matrix(final_mean[k,,(q_op/2+1)],byrow=TRUE,16,16)[16:1,]),col=gray(seq(0,1,length=d)),axes=FALSE)
  for(j in 1:5){
    image(t(matrix(rmvnorm(1,mean=final_mean[k,,(q_op/2+1)],sigma=final_sigma[,,k,(q_op/2+1)]),byrow=TRUE,16,16)[16:1,]),col=gray(seq(0,1,length=d)),axes=FALSE)
  }
}

#Accuracy
accuracy=0
for (i in 1:data_number){
  final_cluster[i]=which.max(final_gamma[i,,q_op/2+1])
}
#find the true diginal number of each cluster
true_k=array(-1,dim=k_num)
#myLabel[which(final_cluster==1)] 
for(k in 1:k_num){
  true_k[k]=strtoi(names(sort(table(myLabel[which(final_cluster==k)]),decreasing = TRUE))[1])
}
mis_count=array(0,dim=k_num)
mis_rate=array(0,dim=k_num)
for(k in 1:k_num){
  mis_count[k]=length(which(myLabel[which(final_cluster==k)]!=true_k[k]))
  mis_rate[k]=mis_count[k]/length(myLabel[which(final_cluster==k)])
}
total_mis=sum(mis_count)/data_number


PlotPic=matrix(1:Iteration,Iteration,1)
dev.new(width=8,height=8)
plot(PlotPic,likelihood[,1],main="q=0",xlab="iteration")

PlotPic=matrix(1:Iteration,Iteration,1)
dev.new(width=8,height=8)
plot(PlotPic,likelihood[,2])
