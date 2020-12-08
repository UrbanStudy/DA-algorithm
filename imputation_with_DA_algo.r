#######################################################################################
# This code is designed to compare the Data Augmentation algorithm to the MICE 
# algorithm for data imputation. How it works:
# 1. Import a dataset with no missing values and subset to only include continuous values
#    as it will likely do a poor job (or fail if nonnumeric) with categorical values
# 2. Either use bayes_v_mice() or bayes_v_mice_JEFFREYS() to compare this with MICE.
#    Both perform best with normally distributed data, which is largely what we've
#    used in our experiment below.
# 3. For either function specify the dataset, the count of NAs you want to make in
#    the dataset, the number of iterations you want the Bayes code to run
#    and if you want to just do one overall removal of values or if you want to 
#    perform a single removal of several values.
#    An iterative removal of values will result in a successive (i.e. 1, 2, 3, ...)
#    number of imputations and simulations... this can take a very long time
#    An example of the iterative removal that would do 30 successive imputations:
#    bayes_v_mice(
#        Y.fish,
#        CountOfNewNA = 30,
#        DAIterations = 1000,
#        JustDoOne = FALSE
#    )
#   and this would be an example of just removing the values and doing one imputation
#   of 400 values
#    bayes_v_mice(
#        Y.fish,
#        CountOfNewNA = 400,
#        DAIterations = 1000,
#        JustDoOne = TRUE
#    )
# 4. Visualizing the results can be done with plot_my_bayes() - the JustDoOne flag
#    must be identical to what was used in your experiment above. See the section
#    on plotting below to find out how to best use these.
# 
# You can run this full code and it will take several hours or you can use the functions
# within to do your own simulations. If you want to use the functions below just delete
# everything below line 650. Note that we made some tables at the end you might also
# find useful at line 862.
# YOUR DATASET MUST HAVE NO MISSING VALUES BEFORE USING THESE FUNCTIONS, OR ELSE
# MODIFICATIONS WILL HAVE TO BE MADE TO THE FUNCTIONS BELOW!
#######################################################################################
library(mice) #for imputation comparisons
library(MASS) #for multivariate normal draws
library(covreg) #for Wishart draws
library(ggplot2) #for plotting
library(GGally) #for pairs plots
library(knitr) #for nice tables
library(mice) #its so nice we loaded it twice
library(reshape2) #for melt

#######################################################################################
# This is our Bayes function we will use for imputation - a slightly modified version 
# of Hoff's code located at http://www2.stat.duke.edu/~pdh10/FCBS/Replication/chapter7.R
# There were some slight changes in Hoff's code to ensure it ran with all datasets.
# this function "double dips" into the data using means, and it has been suggested
# this is bad practice. Below we will use a version using Jeffreys prior 
# to copmare an "uninformative" prior vs this "mean" prior that 
# is significantly more dependent on our sample
#######################################################################################
daa_for_bay_ay_ayes = function(MissingValuesMatrix, iterations) {
    set.seed(1)
    
    # set up the multivariate normal function
    rmvnorm = function(n, mu, Sigma) {
        p <- length(mu)
        res <- matrix(0, nrow = n, ncol = p)
        if (n > 0 & p > 0) {
            E <- matrix(rnorm(n * p), n, p)
            res <- t(t(E %*% chol(Sigma)) + c(mu))
        }
        res
    }
    
    ## prior parameters
    n = dim(MissingValuesMatrix)[1]
    p = dim(MissingValuesMatrix)[2]
    
    mu0 = apply(MissingValuesMatrix, 2, mean, na.rm = TRUE)
    sd0 = (mu0 / 2)
    L0 = matrix(.1, p, p)
    diag(L0) = 1
    L0 = L0 * outer(sd0, sd0)
    nu0 = p + 2
    S0 = L0
    ###
    
    ### starting values
    Sigma = S0
    MissingValuesMatrix.full = MissingValuesMatrix
    

    ### this O notation is super painful and confusing, but it is Hoff's so we will go with it.
    O = 1 * (!is.na(MissingValuesMatrix))
    for (j in 1:p)
    {
        MissingValuesMatrix.full[is.na(MissingValuesMatrix.full[, j]), j] = mean(MissingValuesMatrix.full[, j], na.rm = TRUE)
    }
    
    ### Gibbs sampler
    THETA = SIGMA = MissingValuesMatrix.MISS = NULL
    #set.seed(1)
    
    for (s in 1:iterations)
    {
        ###update theta
        ybar = apply(MissingValuesMatrix.full, 2, mean)
        Ln = solve(solve(L0) + n * solve(Sigma))
        mun = Ln %*% (solve(L0) %*% mu0 + n * solve(Sigma) %*% ybar)
        theta = rmvnorm(1, mun, Ln)
        ###
        
        ###update Sigma
        Sn = S0 + (t(MissingValuesMatrix.full) - c(theta)) %*% t(t(MissingValuesMatrix.full) -
                                                                     c(theta))
        Sigma = solve(rwish(solve(Sn), nu0 + n))
        ###
        
        ###update missing data
        for (i in 1:n)
        {
            b = (O[i,] == 0)
            a = (O[i,] == 1)
            iSa = solve(Sigma[a, a])
            beta.j = Sigma[b, a] %*% iSa
            s2.j  = Sigma[b, b] - Sigma[b, a] %*% iSa %*% Sigma[a, b]
            theta.j =  theta[b] + beta.j %*% (t(MissingValuesMatrix.full[i, a]) -
                                                  theta[a])
            
            ### ALLSO IS THIS SUPPOSED TO BE COMMENTED OUT?!?!
            MissingValuesMatrix.full[i, b]  =  rmvnorm(1, theta.j, s2.j)
            
            imputed_value = matrix(-1, nrow = 1, ncol = nrow(theta.j))
            ####### HERE IS THE IF STATEMENT THAT I ADDED IN
            if (!is.na(imputed_value[1])) {
                while (all(imputed_value < 0)) {
                    imputed_value = rmvnorm(1, theta.j, s2.j)
                }
                MissingValuesMatrix.full[i, b]  =  imputed_value
            }
        }
        
        ### save results
        THETA = rbind(THETA, theta)
        SIGMA = rbind(SIGMA, c(Sigma))
        
        MissingValuesMatrix.MISS = rbind(MissingValuesMatrix.MISS, MissingValuesMatrix.full[O == 0])
        ###
        
        #cat(s,theta,"\n")
    }
    return(list(
        MissingValuesMatrix.MISS,
        MissingValuesMatrix.full,
        THETA,
        SIGMA
    ))
}


#######################################################################################
# This simply runs things through MICE using its most basic settings.
# Until the last decade or so it was widely accepted that 5 iterations were enough to solve this problem
# So this is the default of MICE
# Tuning MICE is beyond the scope of our project
# to call this object you should set a variable equal to miceRun(your_matrix)
# EXAMPLE USE:
# my_mice_run = mice_run(A.matrix.with.missing.values)
# completedData <- complete(my_mice_run,1)
#######################################################################################
mice_run = function(MissingValuesMatrix) {
    # use all mice defaults
    # the method mice typically uses for continuous variables is predictive mean matching or pmm
    # it is a popular algorithm
    # It calculates the predicted value using a regression model and picks the 5 closest elements to the predicted value (by Euclidean distance).
    # These chosen elements are called the donor pool and the final value is chosen at random from this donor pool.
    # Note that the row is the target variable and the column the predictors
    imputed = mice(MissingValuesMatrix, printFlag = FALSE, silent=TRUE, ridge=0.0001)
    completed_data <- complete(imputed,1)
    return(completed_data)
}


#######################################################################################
# This plots our data across a 45 degree line where the x-axis will be the values
# removed using our random removal process and the y-axis will be the imputed values
# from Bayes and MICE. 
# The first function get_plottable_data adds all values into a single dataframe 
# to make plotting easier.
#######################################################################################

#plotting
get_plottable_data = function(AnalyzeDF,  CountOfNA)
{
    df = as.data.frame(AnalyzeDF[[CountOfNA]][1], col.names = "Winner")
    df$Actual = AnalyzeDF[[CountOfNA]][3][[1]]
    df$Bayes = AnalyzeDF[[CountOfNA]][4][[1]]
    df$Mice = AnalyzeDF[[CountOfNA]][5][[1]]
    return(df)
}


#MUST SET JustDoOne=TRUE IF YOU DID IT PRIOR
#LEAVE CountOfNA=1 IF JustDoOne
plot_my_bayes = function(df_plot,
                         CountOfNA  = 5,
                         JustDoOne = FALSE,
                         SubsetByVals = FALSE,
                         subset_val = 100) {
    if (JustDoOne == TRUE) {
        plot_data = get_plottable_data(list(df_plot), 1)
    } else {
        plot_data = get_plottable_data(df_plot, CountOfNA)
    }
    if (SubsetByVals == TRUE) {
        
        plot_data_small = subset(plot_data, Actual < subset_val)
       
        plot_data_large = subset(plot_data, Actual > subset_val)
        
        dataM_small <- melt(plot_data_small[2:4], id.vars = "Actual")
        dataM_large <- melt(plot_data_large[2:4], id.vars = "Actual")
        print(
            ggplot(dataM_small[order(dataM_small$variable, decreasing = T), ], aes(
                x = Actual,
                y = value,
                colour = variable
            )) +
                geom_point(size = 3, alpha = 0.5) + geom_abline(h = 45)
        )
        print(
            ggplot(dataM_large[order(dataM_large$variable, decreasing = T), ], aes(
                x = Actual,
                y = value,
                colour = variable
            )) +
                geom_point(size = 2, alpha = 0.5) + geom_abline(h = 45)
        )
    } else {
        dataM <- melt(plot_data[2:4], id.vars = "Actual")
        ggplot(dataM, aes(
            x = Actual,
            y = value,
            colour = variable
        )) +
            geom_point(size = 3, alpha = 0.5) + geom_abline(h = 45)
    }
    
}

#######################################################################################
# Plot a barchart that shows the percentage of Bayes values that are closer to the
# deleted values when compared to MICE
#######################################################################################

generate_winning_plots = function(df) {
    list_to_fill = list()
    percent_of_bayes = list()
    for (jj in 1:length(df)) {
        my_results <- df[[jj]][[1]]
        my_results[is.na(my_results)] = 1
        list_to_fill = c(df, sum(my_results))
        percent_of_bayes = c(percent_of_bayes, sum(my_results) / jj)
    }
    unlist(percent_of_bayes)
    bayes_percent = data.frame("Percent" = unlist(percent_of_bayes),
                               "NAs" = seq(1:length(df)))
    ggplot(data = bayes_percent, aes(x = NAs, y = Percent)) + geom_bar(stat =
                                                                           "identity", fill = "#F8766D")
}

#######################################################################################
# This function chooses a winner (loser) and assigns them based on the value that is
# the smallest absolute distance away from the actual removed value. 
# Returns a vector of 0s and 1s - if MICE is closer it is 0 and if Bayes is closer
# it returns a 1. 
#######################################################################################
soy_un_perderdor = function(WhyDidIMakeThisListSoLong, DAIterations, N) {
    #get the missing values
    removed_values = WhyDidIMakeThisListSoLong[[3]]
    
    #get the indices to examine because we need to do the thing
    indices_to_examine = WhyDidIMakeThisListSoLong[[4]]
    
    #get the mice imputed values
    mice_return = WhyDidIMakeThisListSoLong[[1]]
    
    
    ########
    # do some bayes processing
    ########
    
    #get the bayes list that we had returned
    bayes_return = WhyDidIMakeThisListSoLong[[2]]
    
    #pull out our bayes variables
    #a column of all imputations
    bayes_return.MISS = bayes_return[[1]]
    #a dataframe with the last updated value
    bayes_return.FULL = bayes_return[[2]]
    
    #get the final values that are correct
    if(N==1){
        bayes_in_order = bayes_return.FULL[indices_to_examine[[1]],indices_to_examine[[2]]]
    } else{
        bayes_in_order = bayes_return.FULL[indices_to_examine]
    }
    
    ###
    
    
    #get each run for every NA value to average it - as.numeric returns it as a vector
    bayes_out_of_order = as.numeric(tail(bayes_return.MISS, 1))

    #order this... painfully
    vector=c(0)
    for (jj in 1:length(bayes_out_of_order)) {
        vector[[jj]] = match(bayes_out_of_order[[jj]], bayes_in_order)
    }
    # get the column means
    ninety_percent_of_DAIterations = DAIterations*.9
    burn_in_df = tail(bayes_return.MISS, ninety_percent_of_DAIterations)
    unordered_column_means = colMeans(burn_in_df)
    #order them by value
    ordered_mean_values = unordered_column_means[order(vector)]
    #########
    # now we can get the differences
    ########
    mice_diff = mice_return - removed_values
    bayes_diff =  ordered_mean_values - removed_values
    #this checks to see which was furthest away...
    absolute_diff = abs(mice_diff) - abs(bayes_diff)
    
    #this creates a vector of 0s and 1s, where 1s are where the bayes outperformed mice and vice-versa
    if_bayes_is_closer_then_ONES = ifelse(absolute_diff > 0, 1, 0)
    
    return(
        list(
            if_bayes_is_closer_then_ONES,
            absolute_diff,
            removed_values,
            ordered_mean_values,
            mice_return
        #    mice_diff,
        #    bayes_diff,
        )
    )
}


#######################################################################################
# This function randomly removes N items from a matrix and replaces them with NAs
# shoutout to this edited thing on stackexchange I had to dig deep for
# https://stackoverflow.com/posts/21372274/revisions
# I was having difficulty with it working at N=1, so I updated the function accordingly
#######################################################################################


replace_with_NAs = function(NoMissingValuesMatrix, N) {
    set.seed(1)
    # Generate a list of i, j coordinates
    indices <-
        as.matrix(expand.grid(
            1:nrow(NoMissingValuesMatrix),
            1:ncol(NoMissingValuesMatrix)
        ))
    
    # Get the random ij_th value in the matrix
    sampled_indices <- indices[sample(nrow(indices), N), ]
    #subset the matrix with only those values
    #for some reason have to handle N=1 differently
    if (N == 1) {
        little_x = sampled_indices[[1]]
        little_y = sampled_indices[[2]]
        littles = c(little_x, little_y)
        values = NoMissingValuesMatrix[little_x, little_y]
        NoMissingValuesMatrix[little_x, little_y] = NA
    } else {
        values = NoMissingValuesMatrix[sampled_indices]
        NoMissingValuesMatrix[sampled_indices] <- NA
    }
    
    # This is accessed as TABLE_NAME[[1]] for results and TABLE_NAME[[2]] for values
    # and TABLE_NAME[[3]] for the indices that can plugged in to get the original vs imputed values!
    return(list(NoMissingValuesMatrix, values, sampled_indices, littles))
}




#######################################################################################
# Actual Function to call that simulates Bayes V MICE
#######################################################################################

bayes_v_mice = function(SomeMatrix,
                        CountOfNewNA = 10,
                        DAIterations = 100,
                        JustDoOne = FALSE) {
    set.seed(1)
    # just do 1
    # handle the addition of variables - do we just want one run?
    our_returnable_list = list()
    if (sum(is.na(SomeMatrix)) == 0 && JustDoOne == TRUE) {
        # call the function to update the matrix
        replaced_matrix = replace_with_NAs(SomeMatrix, CountOfNewNA)
        # set each variable - the first being the new matrix with NAs
        MissingValuesMatrix = replaced_matrix[[1]]
        # actual values that were removed from the matrix
        removed_values = our_returnable_list[[3]] = replaced_matrix[[2]]
        # the indices taken from the amtrix
        sampled_indices = our_returnable_list[[4]] = replaced_matrix[[3]]
        #add these things to our lists
        #first mice
        mice_df = mice_run(MissingValuesMatrix)
        #get the mice_df missing values
        mice_missing_vals = our_returnable_list[[1]] = mice_df[sampled_indices]
        #second bayes
        all_the_bayes = our_returnable_list[[2]] = daa_for_bay_ay_ayes(MissingValuesMatrix, DAIterations)
        #get the bayes df
        bayes_df = all_the_bayes[[2]]
        
        #check if bayes won
        bayes_v_mice_results = soy_un_perderdor(our_returnable_list, DAIterations, CountOfNewNA)
    }
    ## do this thing when you want to see several values of iterations and how it performs over time
    else if (sum(is.na(SomeMatrix)) == 0 && JustDoOne == FALSE) {
        bayes_v_mice_results = list()
        for (ii in 1:CountOfNewNA) {
            print(ii)
            our_returnable_list = list()
            replaced_matrix = replace_with_NAs(SomeMatrix, ii)
            # set each variable - the first being the new matrix with NAs
            MissingValuesMatrix = replaced_matrix[[1]]
            # actual values that were removed from the matrix
            removed_values = our_returnable_list[[3]] = replaced_matrix[[2]]
            
            #first mice
            mice_df = mice_run(MissingValuesMatrix)
            
            # the indices taken from the matrix
            if (ii == 1) {
                sampled_indices = our_returnable_list[[4]] = replaced_matrix[[4]]
                mice_missing_vals = our_returnable_list[[1]] = mice_df[sampled_indices[[1]], sampled_indices[[2]]]
                
            } else{
                sampled_indices = our_returnable_list[[4]] = replaced_matrix[[3]]
                mice_missing_vals = our_returnable_list[[1]] = mice_df[sampled_indices]
            }
            
            #second bayes
            all_the_bayes = our_returnable_list[[2]] = daa_for_bay_ay_ayes(MissingValuesMatrix, DAIterations)
            #get the bayes df

            #check if bayes won
            bayes_v_mice_results[[ii]] = soy_un_perderdor(our_returnable_list, DAIterations, ii)
        }
    }
    if (sum(is.na(SomeMatrix)) > 0) {
        print('that matrix already works, fool')
    }
    return(bayes_v_mice_results)
    
}


#######################################################################################
# Impute values using Jeffreys Prior
#######################################################################################

daa_for_bay_ay_ayes_JEFFREY = function(MissingValuesMatrix, iterations) {
    set.seed(1)
   
    # set up the multivariate normal function
    rmvnorm = function(n, mu, Sigma) {
        p <- length(mu)
        res <- matrix(0, nrow = n, ncol = p)
        if (n > 0 & p > 0) {
            E <- matrix(rnorm(n * p), n, p)
            res <- t(t(E %*% chol(Sigma)) + c(mu))
        }
        res
    }
   
    ## prior parameters
    n = dim(MissingValuesMatrix)[1]
    p = dim(MissingValuesMatrix)[2]
   
    #mu0 = apply(MissingValuesMatrix, 2, mean, na.rm = TRUE)
    # sd0 = (mu0 / 2)
    # L0 = matrix(.1, p, p)
    # diag(L0) = 1
    # L0 = L0 * outer(sd0, sd0)
    # nu0 = p + 2
    # S0 = L0
   
    ###
   
    ### starting values
    #Sigma = S0
    MissingValuesMatrix.full = MissingValuesMatrix
   
    O = 1 * (!is.na(MissingValuesMatrix))
    for (j in 1:p)
    {
        MissingValuesMatrix.full[is.na(MissingValuesMatrix.full[, j]), j] = mean(MissingValuesMatrix.full[, j], na.rm = TRUE)
    }
   
    #get starter Gibbs Sigma value
    Sigma = cov(MissingValuesMatrix.full)
    ###
   
   
    ### Gibbs sampler
    THETA = SIGMA = MissingValuesMatrix.MISS = NULL
    #set.seed(1)
   
    for (s in 1:iterations)
    {
        ###update theta
        ybar = apply(MissingValuesMatrix.full, 2, mean)
        # Ln = solve(solve(L0) + n * solve(Sigma))
        # mun = Ln %*% (solve(L0) %*% mu0 + n * solve(Sigma) %*% ybar)
        # theta = rmvnorm(1, mun, Ln)
        theta = rmvnorm(1, ybar, Sigma*(1/n))
        ###
       
        ###update Sn (the sample covariance matrix)
        Sn = cov(MissingValuesMatrix.full)
       
        ###update Sigma
        # Sn = S0 + (t(MissingValuesMatrix.full) - c(theta)) %*% t(t(MissingValuesMatrix.full) -
        #                                                              c(theta))
        # Sigma = solve(rwish(solve(Sn), nu0 + n))
        Sigma = solve(rwish(solve(n*Sn), n))
        ###
       
        ###update missing data
        for (i in 1:n)
        {
            b = (O[i,] == 0)
            a = (O[i,] == 1)
            iSa = solve(Sigma[a, a])
            beta.j = Sigma[b, a] %*% iSa
            s2.j  = Sigma[b, b] - Sigma[b, a] %*% iSa %*% Sigma[a, b]
            theta.j =  theta[b] + beta.j %*% (t(MissingValuesMatrix.full[i, a]) -
                                                  theta[a])
           
            ### ALLSO IS THIS SUPPOSED TO BE COMMENTED OUT?!?!
            MissingValuesMatrix.full[i, b]  =  rmvnorm(1, theta.j, s2.j)
           
            imputed_value = matrix(-1, nrow = 1, ncol = nrow(theta.j))
            ####### HERE IS THE IF STATEMENT THAT I ADDED IN
            if (!is.na(imputed_value[1])) {
                while (all(imputed_value < 0)) {
                    imputed_value = rmvnorm(1, theta.j, s2.j)
                }
                MissingValuesMatrix.full[i, b]  =  imputed_value
            }
        }
       
        ### save results
        THETA = rbind(THETA, theta)
        SIGMA = rbind(SIGMA, c(Sigma))
       
        MissingValuesMatrix.MISS = rbind(MissingValuesMatrix.MISS, MissingValuesMatrix.full[O == 0])
        ###
       
        #cat(s,theta,"\n")
    }
    return(list(
        MissingValuesMatrix.MISS,
        MissingValuesMatrix.full,
        THETA,
        SIGMA
    ))
}


#######################################################################################
# use Jeffreys prior to simulate rather than MVN using means from dataset
#######################################################################################

bayes_v_mice_JEFFREYS = function(SomeMatrix,
                        CountOfNewNA = 10,
                        DAIterations = 100,
                        JustDoOne = FALSE) {
    set.seed(1)
    # just do 1
    # handle the addition of variables - do we just want one run?
    our_returnable_list = list()
    if (sum(is.na(SomeMatrix)) == 0 && JustDoOne == TRUE) {
        # call the function to update the matrix
        replaced_matrix = replace_with_NAs(SomeMatrix, CountOfNewNA)
        # set each variable - the first being the new matrix with NAs
        MissingValuesMatrix = replaced_matrix[[1]]
        # actual values that were removed from the matrix
        removed_values = our_returnable_list[[3]] = replaced_matrix[[2]]
        # the indices taken from the amtrix
        sampled_indices = our_returnable_list[[4]] = replaced_matrix[[3]]
        #add these things to our lists
        #first mice
        mice_df = mice_run(MissingValuesMatrix)
        #get the mice_df missing values
        mice_missing_vals = our_returnable_list[[1]] = mice_df[sampled_indices]
        #second bayes
        all_the_bayes = our_returnable_list[[2]] = daa_for_bay_ay_ayes(MissingValuesMatrix, DAIterations)
        #get the bayes df
        bayes_df = all_the_bayes[[2]]
        
        #check if bayes won
        bayes_v_mice_results = soy_un_perderdor(our_returnable_list, DAIterations, CountOfNewNA)
    }
    ## do this thing when you want to see several values of iterations and how it performs over time
    else if (sum(is.na(SomeMatrix)) == 0 && JustDoOne == FALSE) {
        bayes_v_mice_results = list()
        for (ii in 1:CountOfNewNA) {
            print(ii)
            our_returnable_list = list()
            replaced_matrix = replace_with_NAs(SomeMatrix, ii)
            # set each variable - the first being the new matrix with NAs
            MissingValuesMatrix = replaced_matrix[[1]]
            # actual values that were removed from the matrix
            removed_values = our_returnable_list[[3]] = replaced_matrix[[2]]
            
            #first mice
            mice_df = mice_run(MissingValuesMatrix)
            
            # the indices taken from the matrix
            if (ii == 1) {
                sampled_indices = our_returnable_list[[4]] = replaced_matrix[[4]]
                mice_missing_vals = our_returnable_list[[1]] = mice_df[sampled_indices[[1]], sampled_indices[[2]]]
                
            } else{
                sampled_indices = our_returnable_list[[4]] = replaced_matrix[[3]]
                mice_missing_vals = our_returnable_list[[1]] = mice_df[sampled_indices]
            }
            
            #second bayes
            all_the_bayes = our_returnable_list[[2]] = daa_for_bay_ay_ayes_JEFFREY(MissingValuesMatrix, DAIterations)
            #get the bayes df

            #check if bayes won
            bayes_v_mice_results[[ii]] = soy_un_perderdor(our_returnable_list, DAIterations, ii)
        }
    }
    if (sum(is.na(SomeMatrix)) > 0) {
        print('that matrix already works, fool')
    }
    return(bayes_v_mice_results)
    
}



#######################################################################################
# LOAD DATA
# to repeat our experiments we used the below datasets to impute values
#######################################################################################

# https://www.kaggle.com/harshit831/fish-weight
file.fish = "Z:\\Classes\\Stat 572 - Bayesian Statistics\\Project\\fish.csv"
data.fish = read.csv(file.fish)
Y.fish = data.fish[,2:7]

# https://www.kaggle.com/yersever/500-person-gender-height-weight-bodymassindex?select=500_Person_Gender_Height_Weight_Index.csv
file.hw = "Z:\\Classes\\Stat 572 - Bayesian Statistics\\Project\\height_weight_index.csv"
data.hw = read.csv(file.hw)
Y.hw = data.hw[2:3]

#### do just one big removal because it will take forever
# https://www.kaggle.com/piyushgoyal443/red-wine-dataset?select=wineQualityReds.csv
file.wine = "Z:\\Classes\\Stat 572 - Bayesian Statistics\\Project\\wine.csv"
data.wine = read.csv(file.wine)
Y.wine = data.wine[,-1]

#### do just one big removal because it will take forever
# https://www.kaggle.com/justinas/nba-players-data?select=all_seasons.csv
file.nba = "Z:\\Classes\\Stat 572 - Bayesian Statistics\\Project\\nba.csv"
data.nba = read.csv(file.nba)
Y.nba = data.nba[,-which(sapply(data.nba, class) == "factor")]
Y.nba = Y.nba[complete.cases(Y.nba), ]
summary(Y.nba)

#### do just one big removal because it will take forever
# https://www.kaggle.com/jinxbe/wnba-player-stats-2017?select=WNBA+Stats.csv
file.wnba = "Z:\\Classes\\Stat 572 - Bayesian Statistics\\Project\\wnba.csv"
data.wnba = read.csv(file.wnba)
Y.wnba = data.wnba[,-which(sapply(data.wnba, class) == "factor")]
Y.wnba = Y.wnba[complete.cases(Y.wnba), ]



#######################################################################################
# Run our experiments using the "mean" prior
#######################################################################################

fish_experiment = bayes_v_mice(
    Y.fish,
    CountOfNewNA = 30,
    DAIterations = 1000,
    JustDoOne = FALSE
)

hw_experiment = bayes_v_mice(
    Y.hw,
    CountOfNewNA = 15,
    DAIterations = 1000,
    JustDoOne = FALSE
)

wine_experiment_100 = bayes_v_mice(
    Y.wine,
    CountOfNewNA = 100,
    DAIterations = 1000,
    JustDoOne = TRUE
)
wine_experiment_200 = bayes_v_mice(
    Y.wine,
    CountOfNewNA = 200,
    DAIterations = 1000,
    JustDoOne = TRUE
)
wine_experiment_300 = bayes_v_mice(
    Y.wine,
    CountOfNewNA = 300,
    DAIterations = 1000,
    JustDoOne = TRUE
)
wine_experiment_400 = bayes_v_mice(
    Y.wine,
    CountOfNewNA = 400,
    DAIterations = 1000,
    JustDoOne = TRUE
)
wine_experiment_500 = bayes_v_mice(
    Y.wine,
    CountOfNewNA = 500,
    DAIterations = 1000,
    JustDoOne = TRUE
)

nba_experiment = bayes_v_mice(
    Y.nba,
    CountOfNewNA = 30,
    DAIterations = 1000,
    JustDoOne = FALSE
)

wnba_experiment = bayes_v_mice(
    Y.wnba,
    CountOfNewNA = 400,
    DAIterations = 1000,
    JustDoOne = TRUE
)



#######################################################################################
# Plot the experiments using the "mean" prior
# REMEMBER TO SET JustDoOne=TRUE if you have JustDoOne=TRUE for *_experiment
# I was able to generate a gif using imgflip and the plots for the wine_experiments
# and then we compared them to those using Jeffreys prior
#######################################################################################

plot_my_bayes(fish_experiment, 30, SubsetByVals = TRUE)
plot_my_bayes(hw_experiment,15)
plot_my_bayes(wine_experiment, JustDoOne=TRUE)
plot_my_bayes(nba_experiment,30)
plot_my_bayes(nba_experiment,5)
plot_my_bayes(wnba_experiment,JustDoOne=TRUE)
plot_my_bayes(wine_experiment_100, JustDoOne=TRUE)
plot_my_bayes(wine_experiment_200, JustDoOne=TRUE)
plot_my_bayes(wine_experiment_300, JustDoOne=TRUE)
plot_my_bayes(wine_experiment_400, JustDoOne=TRUE)
plot_my_bayes(wine_experiment_500, JustDoOne=TRUE)


#######################################################################################
# Run Jeffreys Prior experiments
#######################################################################################
fish_experiment_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.fish,
    CountOfNewNA = 30,
    DAIterations = 1000,
    JustDoOne = FALSE
)

hw_experiment_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.hw,
    CountOfNewNA = 15,
    DAIterations = 1000,
    JustDoOne = FALSE
)

wine_experiment_100_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.wine,
    CountOfNewNA = 100,
    DAIterations = 1000,
    JustDoOne = TRUE
)
wine_experiment_200_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.wine,
    CountOfNewNA = 200,
    DAIterations = 1000,
    JustDoOne = TRUE
)
wine_experiment_300_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.wine,
    CountOfNewNA = 300,
    DAIterations = 1000,
    JustDoOne = TRUE
)
wine_experiment_400_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.wine,
    CountOfNewNA = 400,
    DAIterations = 1000,
    JustDoOne = TRUE
)
wine_experiment_500_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.wine,
    CountOfNewNA = 500,
    DAIterations = 1000,
    JustDoOne = TRUE
)

nba_experiment_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.nba,
    CountOfNewNA = 30,
    DAIterations = 1000,
    JustDoOne = FALSE
)

wnba_experiment_JEFFREYS = bayes_v_mice_JEFFREYS(
    Y.wnba,
    CountOfNewNA = 400,
    DAIterations = 1000,
    JustDoOne = TRUE
)
plot_my_bayes(wnba_experiment_JEFFREYS,JustDoOne=TRUE)
sum(wnba_experiment_JEFFREYS[[1]])


#######################################################################################
# Plot the experiments using Jeffreys prior
# REMEMBER TO SET JustDoOne=TRUE if you have JustDoOne=TRUE for *_experiment
#######################################################################################


plot_my_bayes(fish_experiment_JEFFREYS, 30, SubsetByVals = FALSE)
plot_my_bayes(hw_experiment_JEFFREYS,15)
plot_my_bayes(wine_experiment_JEFFREYS, JustDoOne=TRUE)
plot_my_bayes(nba_experiment_JEFFREYS,30)
plot_my_bayes(wnba_experiment,JustDoOne=TRUE)
plot_my_bayes(wine_experiment_100_JEFFREYS, JustDoOne=TRUE)
plot_my_bayes(wine_experiment_200, JustDoOne=TRUE)
plot_my_bayes(wine_experiment_300, JustDoOne=TRUE)
plot_my_bayes(wine_experiment_400, JustDoOne=TRUE)
plot_my_bayes(wine_experiment_500, JustDoOne=TRUE)





#######################################################################################
# Create a list for the percentage of Bayes winners when not using the iterative
# version - aka JustDoOne = FALSE - but still imputing across multiple sets of NAs.
# The wine_list and wine_JEFFREYS hold lists of percentages used below. 
# Note the hardcoded values for wine_JEFFREYS is because we ran the data on multiple
# machines because it takes some time to run the simulations.
# Because the results are very similar for our imputations we created a table to 
# explicitly call out the values that were close to one another.
# For example, the NA=400 values are .61 for the "mean" prior & .63 for Jeffreys prior,
# so we can show both plots and clarify the differences if they are visually
# dificult to separate.
#######################################################################################
wine_list = c(sum(wine_experiment_100[[1]])/length(wine_experiment_100[[1]]), sum(wine_experiment_200[[1]])/length(wine_experiment_200[[1]]), sum(wine_experiment_300[[1]])/length(wine_experiment_300[[1]]), sum(wine_experiment_400[[1]])/length(wine_experiment_400[[1]]), sum(wine_experiment_500[[1]])/length(wine_experiment_500[[1]]))
bayes_percent = data.frame("Percent"=unlist(wine_list), "NAs"=c("100","200","300","400","500"))
ggplot(data=bayes_percent, aes(x=NAs,y=Percent)) +geom_bar(stat="identity", fill="#F8766D") 


wine_JEFFREYS = c(sum(wine_experiment_100_JEFFREYS[[1]])/length(wine_experiment_100_JEFFREYS[[1]]), 0.61,0.59,.63,.626)
bayes_percent_JEFFREYS = data.frame("Percent"=unlist(wine_JEFFREYS), "NAs"=c("100","200","300","400","500"))
ggplot(data=bayes_percent_JEFFREYS, aes(x=NAs,y=Percent)) +geom_bar(stat="identity", fill="#F8766D") 


results = rbind(wine_list,wine_JEFFREYS)
rownames(results) <- c('MVN with Mean Prior', 'Jeffreys Prior')
t(results)
colnames(results) <- c("100 Missing", "200 Missing", "300 Missing", "400 Missing", "500 Missing")
knitr::kable(t(results), digits=2)

