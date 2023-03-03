# Adaptation of `fa.parallel` to allow {Turbofuns} correlations
# Drastically speeds up processing
# Much cleaner code
fa.parallel.turbo <- function(
    x, fm = "minres", fa = "both",
    n.iter = 20, SMC = FALSE, quant = 0.95,
    cor = "cor", use = "pairwise"
)
{
  
  # Require {Turboruns} package
  require(Turbofuns)
  
  # Ensure data is matrix
  data <- as.matrix(x)
  
  # Obtain number of variables
  nvariables <- ncol(data)
  
  # Obtain number of observations
  nsub <- nrow(data)
  
  # Obtain correlations
  rx <- switch(
    cor,
    "tet" = Turbofuns::PolychoricRM(data, IMissing = 1, NCore = 1)$correlation,
    "poly" = Turbofuns::PolychoricRM(data, IMissing = 1, NCore = 1)$correlation,
    "cor" = cor(data, use = use)
  )
  
  # Obtain eigenvalues
  valuesx  <- eigen(rx)$values #these are the PC values
  
  # Check for SMC
  if(isTRUE(SMC)){
    diag(rx) <- smc(rx)
    fa.valuesx <- eigen(rx)$values
  }else{
    fa.valuesx  <- suppressWarnings(
      fa(rx, nfactors = 1, rotate = "none", fm = fm, warnings = FALSE)$values
    )
  }  # t
  
  # Set up temporary list
  temp <- list(
    samp = vector("list", n.iter),
    samp.fa = vector("list", n.iter),
    sim = vector("list", n.iter),
    sim.fa = vector("list", n.iter)
  )
  
  # Perform resampling
  templist <- lapply(1:n.iter, function(XX){
    
    # Check for bad data
    bad_data <- TRUE
    
    # Get good data
    while(bad_data){
      
      # Generate data
      sampledata <- matrix(
        apply(data, 2, function(y){
          sample(y, nsub, replace = TRUE)
        }), ncol = nvariables
      )
      
      # Obtain correlations
      C <- switch(
        cor,
        "tet" = Turbofuns::PolychoricRM(sampledata, IMissing = 1, NCore = 1)$correlation,
        "poly" = Turbofuns::PolychoricRM(sampledata, IMissing = 1, NCore = 1)$correlation,
        "cor" = cor(sampledata, use = use)
      )
  
      # Check for NA correlations
      bad_data <- any(is.na(C))
  
    }
    
    # Try resampling until we get a correlation matrix that works                    
    values.samp <- eigen(C)$values
    
    # Check for SMC
    if(isTRUE(SMC)){
      sampler <- C 
      diag(sampler) <- smc(sampler)
      samp.fa <- eigen(sampler)$values
    }else{
      samp.fa <- suppressWarnings(fa(C, fm = fm, nfactors = 1, SMC = FALSE, warnings = FALSE)$values)
    }
    
    # Set up return list
    temp <- list(
      samp = values.samp,
      samp.fa = samp.fa
    )
  
    # Return list
    return(temp)
    
  })
  
  # Convert to matrix
  values <- t(matrix(unlist(templist), ncol = n.iter))
  
  # Obtain means
  values.sim.mean <- colMeans(values, na.rm = TRUE)

  # Obtain quantiles
  values.ci = apply(values, 2 ,function(x){quantile(x, quant)}) #always apply quant
  
  # Obtain PCA
  sim.pcr <- values.sim.mean[1:nvariables] 
  sim.pcr.ci <- values.ci[1:nvariables]
  pc.test <- which(!(valuesx > sim.pcr.ci))[1]-1
  
  # Obtain FA
  sim.far <- values.sim.mean[(nvariables+1):(2*nvariables)]
  sim.far.ci <- values.ci[(nvariables+1):(2*nvariables)]
  fa.test <- which(!(fa.valuesx > sim.far.ci))[1]-1
  
  # Set up results
  results <- list(
    nfact = fa.test,
    ncomp = pc.test,
    values = values,
    values_mean = values.sim.mean,
    values_quant = values.ci,
    pc.values = valuesx,
    fa.values = fa.valuesx
  )
  
  # Return results
  return(results)
  
}
