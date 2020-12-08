###############################################################################
# The uncommented code can be found at
# http://www2.stat.duke.edu/~pdh10/FCBS/Replication/chapter7.R
###############################################################################
# our prior
rmvnorm = function(n,mu,Sigma) {
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 ) {
  E<-matrix(rnorm(n*p),n,p)
  res<-t(  t(E%*%chol(Sigma)) +c(mu))
                   }
  res
}

## prior parameters
#rows
n = dim(Y)[1]
#columns
p = dim(Y)[2]

# gets the mean of every column
mu0 = apply(Y, 2, mean, na.rm = TRUE)

#takes mu0 and divides it by 2
sd0 = (mu0/2)

#creates a matrix placeholder where all values are 0.1
L0 = matrix(.1,p,p)

#sets diagonals to 1
diag(L0) = 1 

#turn the product of mu0/2 into a square matrix that we can multiply by L0
L0 = L0*outer(sd0,sd0)

#set some more parameters and init our sigma matrix S0 for our approximation using invWish 
nu0 = p+2
S0 = Sigma = L0

#reset Y.full
Y.full = Y

# a matrix where 1 includes values and 0 are missing values
O = 1 * (!is.na(Y))

#initiate our values to be imputed by setting them to the mean
for(j in 1:p)
{
  Y.full[is.na(Y.full[,j]),j] = mean(Y.full[,j],na.rm = TRUE)
}
###


### Gibbs sampler
# initiate our parameter space and the matrix Y.MISS where we will store our imputations 
THETA = SIGMA = Y.MISS = NULL
#set.seed(1)


#start the loop
for(s in 1:iterations)
{
  
  ###update theta
  # we set the mean again - this is updated each iteration so it is not mu0
  ybar = apply(Y.full,2,mean)
  
  # this is our initial precision matrix plus rows times the inverse of our parameter Sigma
  # we use it as the variance in our theta draw from the multivariate normal
  # note that Ln is a pxp matrix
  Ln = solve( solve(L0) + n*solve(Sigma) )
  
  #this updates the means - Ln is a pxp matrix and the remainder is a px1 matrix
  # note we are multiplying our pxp matrices by mu0 - the initial means, and adding n*inv(Sigma)*ybar
  # where ybar is going to be our dynamically updated new values
  # to me this looks like we are adding fractoinal values to the mean at each iteration?
  mun = Ln%*%( solve(L0)%*%mu0 + n*solve(Sigma)%*%ybar )
  
  #get a single draw for theta where the matrix for theta is size 1xp
  theta = rmvnorm(1,mun,Ln)
  ###
  
  ###update Sigma
  # note that S0 is the initialized covariance matrix - it is not updated
  t(Y.full[1:8,])-c(theta)
  c(theta)
  t(Y.full[1:8,])
  #then t(Y.full)-c(theta) subtracts each theta value from the value at Y.full. You multiply its transpose. This gives us a pxp matrix
  #this generates our Sigma for our inverse wishart
  Sn = S0 + ( t(Y.full)-c(theta) )%*%t( t(Y.full)-c(theta) )
  
  # nu0+n does not change between runs 
  #Sn is dependent on updating from the subtraction of values from their approximated values with theta
  Sigma = solve( rwish(solve(Sn), nu0+n) )
  ###
  ###update missing data doing this once for each row (1:n times)
  for(i in 1:n)
  { 
    i=3
    # creates FALSE if value is not missing
    # creates TRUE if value is missing
    b = ( O[i,] == 0 )
    
    #creates TRUE if value is not missing
    #creates FALSE if value is missing
    a = ( O[i,] == 1 )
    Sigma
    Sigma[a,a]
    solve(Sigma[a,a])
    #this generates a covariance matrix based only on values that exist in the dataset
    iSa = solve(Sigma[a,a])
    
    
    # beta.j is a (#of missing values)x(p-#of missing values) matrix
    # Sigma[b,a] is the rows of the covariance matrix for values that are missing but does not include those values themselves
    # eg the table for the missing values of SkinThickness and insulin for the 3rd row in the kaggle diabetes csv look like:
    #                         Glucose    BloodPressure       BMI              DiabetesPedigreeFunction       Age
    #        SkinThickness   85.33774      37.19539        54.94596                0.2156516             20.20017
    #        Insulin       2285.08519     218.14113       224.21611                4.8971403             312.56007      
    beta.j = Sigma[b,a]%*%iSa
    
    # s2.j is a (# of missing values)x(# of missing values) size matrix
    # where sigma[b,b] is the covariance matrix of missing values only
    # and the second term is essentially Sigma[b,a]*isa*Sigma[b,a]^T where the last term is a transpose...
    s2.j  = Sigma[b,b] - Sigma[b,a]%*%iSa%*%Sigma[a,b]
    
    # theta.j is a 1x(# of missing values) that contains our estimates
    # theta[b] returns the missing values in this row for thetas - this is a (#of missing values)x1 matrix
    # add that to our beta.j multiplied by 
    # the transpose(data row without missing values) minus the new theta values of values that are not missing 
    theta.j =  theta[b] + beta.j%*%(t(Y.full[i,a])-theta[a])
    
    # this updates our row with values that are missing using theta.j and s2.j as the parmaeters to update the prediction
    Y.full[i,b]  =  rmvnorm(1,theta.j,s2.j )
  }
  
  ### save results
  THETA = rbind(THETA,theta) ; SIGMA = rbind(SIGMA,c(Sigma))
  
  #add results from each iteration for each missing value
  Y.MISS = rbind(Y.MISS, Y.full[O == 0] )

}

#this gives us mean estimates for our values, but do we really care about this?
results = rbind(quantile(THETA[,1], c(0.025, 0.5, 0.975)),
                quantile(THETA[,2], c(0.025, 0.5, 0.975)),
                quantile(THETA[,3], c(0.025, 0.5, 0.975)),
                quantile(THETA[,4], c(0.025, 0.5, 0.975)),
                quantile(THETA[,5], c(0.025, 0.5, 0.975)),
                quantile(THETA[,6], c(0.025, 0.5, 0.975)),
                quantile(THETA[,7], c(0.025, 0.5, 0.975)))

results = as.data.frame(results, row.names = c("Glucose", "Blood Pressure", 
                                               "Skin Thickness", "Insulin",
                                               "BMI", "Pedigree Function",
                                               "Age"))

knitr::kable(results)
