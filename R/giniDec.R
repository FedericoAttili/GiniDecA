#prova

setwd(getwd())

library(Rcpp)
#Rcpp.package.skeleton("GiniDecA")
#sourceCpp('https://github.com/FedericoAttili/GiniDecA/src/ginicpp.cpp') # Replace with the actual path to your C++ file


########################
#useful functions

mapping=function(groups){#mapping of factor labels for conversion to numeric

  # Convert groups label to factor
  group_factor <- factor(groups)

  # Convert factor to numeric
  groups <- as.numeric(group_factor)

  # Convert numeric back to original categories
  key <- data.frame(label=levels(group_factor)[unique(groups)],id_group=as.numeric(unique(groups))-1)
  mapping <- setNames(key$label, key$id_group)
}


# Function to calculate the greatest common divisor (GCD) of two numbers
gcd <- function(a, b) {
  if(b == 0) {
    return(abs(a))
  } else {
    return(gcd(b, a %% b))
  }
}

# Function to calculate the least common multiple (LCM) of two numbers
lcm_pair <- function(a, b) {
  return(abs(a * b) / gcd(a, b))
}

# Function to calculate the LCM of a vector of numbers
lcm <- function(x) {
  # Check if the vector has at least one element
  if(length(x) < 1) {
    stop("Input vector must have at least one element.")
  }

  # Initialize LCM with the first element
  result_lcm <- x[1]

  # Iteratively calculate LCM for the rest of the elements
  for(i in x[-1]) {
    result_lcm <- lcm_pair(result_lcm, i)
  }

  return(result_lcm)
}

# Example usage
numbers <- c(3, 7, 8)
lcm(numbers) # Should return 24, which is the LCM of 4, 6, and 8


group_size_recalib <- function(x, z) {
  # The function takes an income vector x and a grouping vector z, then recalibrates
  # the data such that each group has the same number of observations based on the
  # least common multiple of the original group sizes. It scales the data using
  # inverse proportions to maintain the overall distribution while ensuring
  # uniform group sizes. The recalibration could be used for certain types of
  # analysis where uniform group sizes are required for comparison or aggregation
  # purposes.

  # Check if x and z are of the same length
  if (length(x) != length(z)) {
    stop("Lengths of 'x' and 'z' must be the same")
  }

  # Calculate the number of unique groups in z
  K <- length(table(z))

  # Initialize vector to store the size of each group
  n_k <- numeric(K)

  # Populate n_k with the size of each group
  for (i in 1:K) {
    n_k[i] <- table(z)[[i]]
  }

  # Calculate the least common multiple (LCM) of group sizes
  n <- lcm(n_k)

  # Calculate the inverse proportion for each group size relative to LCM
  p_k_inv <- n / n_k

  # Combine income data with group labels and sort, by group and income
  x_old <- cbind(x, z)
  x_old <- x_old[order(x_old[, 2], x_old[, 1]), ]

  # Initialize new matrix for recalibrated data
  x_new <- matrix(nrow = n * K, ncol = 3)

  # Recalibrate data for each group
  for (i in 1:K) {
    start_index <- sum(n_k[1:i]) - n_k[i] + 1
    end_index <- sum(n_k[1:i])
    group_data <- x_old[start_index:end_index, 1]

    # Replicate data to fill new group size and assign group label and inverse proportion
    x_new[(n * (i - 1) + 1):(n * i), 1] <- rep(group_data, length.out = n)
    x_new[(n * (i - 1) + 1):(n * i), 2] <- rep(i, n)
    x_new[(n * (i - 1) + 1):(n * i), 3] <- rep(p_k_inv[i], n)
  }

  # Sort the recalibrated data by group and income
  x_new <- x_new[order(x_new[, 2], x_new[, 1]), ]

  # Return the recalibrated data
  return(x_new)
}

library(Hmisc)
quant <- function(gruppi, indice_gruppo, m, group_sizes, type=7, weights=F) {
  quantili <- seq(0, 1, 1 / (m - 1))  # Sequence of quantiles to calculate
  valori_quantili <- c()  # Initialize vector to store quantile values
  g <- c()  # Initialize vector to store group identifiers

  if (!is.numeric(weights[1])) {
    # If weights are not provided
    for (i in group_sizes) {
      # Calculate quantiles for each group and append to valori_quantili
      valori_quantili <- append(valori_quantili, as.vector(quantile(gruppi[indice_gruppo == i], quantili, type = type)))
      # Append group identifier repeated m times
      g <- append(g, rep.int(i, m))
    }
  } else {
    # If weights are provided
    for (i in group_sizes) {
      # Calculate weighted quantiles for each group and append to valori_quantili
      valori_quantili <- append(valori_quantili, as.vector(wtd.quantile(gruppi[indice_gruppo == i], weights[indice_gruppo == i], quantili)))
      # Append group identifier repeated m times
      g <- append(g, rep.int(i, m))
    }
  }

  return(list(gruppi = valori_quantili, indice_gruppo = g))
}


########################

#CORE FUNCTION

library(dineq)

giniDec=function(x,z,w=F,m=NULL,n_equalizer='quant',type=7,contrib=F){
  #x is the vector with the input variable in the population
  #z is the vector of group belonging; it defines the partition in groups
  #w represents the vector of individual weights
  # m tra 0 e 1 indica il livello del quantile, or if >1 directly the value to employ,  if NULL the value suggested in Attili (2020) is employed
  #type indica la def di quantile utilizzata
  if(!is.numeric(w[1])){data=data.frame(x,z,1)}else{data=data.frame(x,z,w)}
  colnames(data)=c('x','groups','w')
  data_store=data
  #force the three columns to be numeric,
  data$x=as.numeric(data$x)
  data$w=as.numeric(data$w)
  data$groups=as.numeric(as.factor(data$groups))
  Gini=gini.wtd(data$x,data$w)

  data=data[order(data$groups),]
  x=data$x
  w=data$w
  z=data$groups


  #recalibration
  #nb if n_equalizer=="ricalibrazione", data should not contain weights

  n_k=as.vector(table(z))
  K=length(n_k)
  kapp=unique(z)

  N_orig=sum(w)
  mu=sum(x*w)/N_orig
  w_k=as.vector(tapply(w, z, sum))
  if(is.null(m)){
    mm=0
    for (i in 1:K){mm=mm+n_k[i]*(w_k[i]/sum(w_k))}
    q=seq(0,1,0.005)
    qq=quantile(n_k,q)
    dif=abs(qq-mm)
    m=q[which.min(dif)]
  }

  if(n_equalizer=='quant'){
    if(m<=1){m=as.integer(quantile(n_k,m))}
    quantili=quant(x, z, m, kapp,type, w)
    z=quantili$indice_gruppo
    x=quantili$gruppi
    numerosita=c(matrix(rep(w_k,m),m,K,byrow = T))
    p_k_inv=m/numerosita
    x=cbind(x,z,p_k_inv)
    rm(z,p_k_inv)

  }

  if(n_equalizer=="ricalibrazione") {
    x=ricalibrazione_gruppi(x,z)
    m=mcm(n_k)    #restituisce x=cbind(x,z,p_k_inv) con numerosit? apparate e ordinato prima su z e poi su x
  }

  x[,3]=1/x[,3]

  data=as.matrix(x)

  output_giniDecomposition <- giniDecomposition(data,contrib=contrib)#contrib=F by default
  #this function is written in c++. It takes data and return
  #the value of the two components (g_w and g_b). If contrib=T,
  #also two matrices are provided with the contribution of each
  #pair of units to within and between component (eq. ...)
  #organise results
  G_dot=output_giniDecomposition$g_w+output_giniDecomposition$g_b
  #evaluate the components
  G_w=output_giniDecomposition$g_w*Gini/G_dot
  G_b=output_giniDecomposition$g_b*Gini/G_dot
  #evaluate the component shares
  within_share=G_w/Gini
  between_share=G_b/Gini

  if(contrib==T){
    map=mapping(data_store$groups)
    #evaluate (absolute) contribution of groups to within inequality
    group_w_contrib=(tapply(output_giniDecomposition$w_contributions$value,output_giniDecomposition$w_contributions$k,sum)+tapply(output_giniDecomposition$w_contributions$value,output_giniDecomposition$w_contributions$h,sum))/2*Gini/G_dot
    names(group_w_contrib) <- map[names(group_w_contrib)]


    #evaluate (absolute) contribution of groups to between inequality
    group_b_contrib=(tapply(output_giniDecomposition$b_contributions$value,output_giniDecomposition$b_contributions$k,sum)+tapply(output_giniDecomposition$b_contributions$value,output_giniDecomposition$b_contributions$h,sum))/2*Gini/G_dot
    # Replace the names in group_b_contrib with corresponding labels
    names(group_b_contrib) <- map[names(group_b_contrib)]


    #evaluate (absolute) pairwise contribution of groups to between inequality
    group_pairwise_contr=tapply(output_giniDecomposition$b_contributions$value,list(output_giniDecomposition$b_contributions$k,output_giniDecomposition$b_contributions$h),sum)*Gini/G_dot
    group_pairwise_contr=group_pairwise_contr[upper.tri(group_pairwise_contr)]+t(group_pairwise_contr)[upper.tri(t(group_pairwise_contr))]
    #set names
    names(group_pairwise_contr)=apply(combn(names(group_b_contrib),2),2,paste,collapse='-')


    #ADD
    #evaluate (absolute) contribution of ranks to within and between

  }

  results <- list(
    Gini = Gini,  # Assuming Gini is calculated somewhere in your function
    G_w = G_w,
    G_b = G_b,
    within_share = within_share,
    between_share = between_share
  )

  if (contrib) {
    # Additional results when contrib is TRUE
    results$group_w_contrib = group_w_contrib
    results$group_b_contrib = group_b_contrib
    results$group_pairwise_contr = group_pairwise_contr
    # ... [any other additional results you want to include]
  }

  return(results)
}
