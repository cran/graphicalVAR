# # Augmented graphical VAR:
# 
# 
# AGV <-
#   function(
#     data, # A n by p data frame containing repeated measures
#     nLambda = 50, # Either single value or vector of two corresponding to c(kappa, beta)
#     verbose = TRUE,
#     gamma = 0.5,
#     scale = TRUE,
#     time, # Vector with response times
#     day, # Vector with days
#     include = c("lags","change","interactions"), # Add day!
#     lags = 1,
#     lambda_beta,
#     lambda_kappa, maxit.in = 100, maxit.out = 100,
#     penalize.diagonal
#   ){
#     
#     # Check input:
#     if (is.data.frame(data)){
#       data <- as.matrix(data)
#     }
#     
#     stopifnot(is.matrix(data))
#     
#     
#     Nvar <- ncol(data)
#     Ntime <- nrow(data)
#     
#     # Center data:
#  
#     if (is.matrix(data)){
#       data <- as.data.frame(data)
#     }
#     if (is.null(names(data))){
#       names(data) <- paste0("V",seq_len(ncol(data)))
#     }
#     Scale <- function(x,center=TRUE,scale=TRUE){
#       if (center){
#         y <- x - mean(x,na.rm=TRUE)
#       } else y <- x
#       
#       if (scale){
#         if(sd(x,na.rm=TRUE)==0){
#           return(0) 
#         }
#         y <- y / sd(x,na.rm=TRUE)
#       }
#       return(y)
#     }
#     
#     
#     data <- data %>% mutate_each(funs(Scale(.,TRUE,scale)))
#     varNames <- names(data)
#     augData <- data
#     names(augData) <- paste0("Y_@_",names(data))
#     
#     # Changescores:
#     changescore <- function(x)c(NA,diff(x))
#     changeData <- data %>% mutate_each(funs(changescore)) %>% mutate_each(funs(Scale(.,TRUE,scale)))
#       
#     # Lag 1 correlations:
#     if (!"lags" %in% include){
#       stop("Lag effects must be included")
#     }
#     for (l in lags){
#       lag <- do.call(rbind,c(rep(list(NA),l),list(data[1:(nrow(data)-l),])))
#       names(lag) <- paste0("L",l,"_@_",names(data))
#       augData <- cbind(augData,lag)
# 
#       # Changescore:
#       if ("change" %in% include){
#         lc <- l
#         change <-  do.call(rbind,c(rep(list(NA),lc),list(changeData[1:(nrow(changeData)-lc),])))
#         names(change) <- paste0("C",l,"_@_",names(data))
#         augData <- cbind(augData,change)
#       }
#     }
# 
#     # Remove all rows with missing data:
#     augDataNoMissing <- augData[rowSums(is.na(augData))==0,]
#     whichY <- grepl("^Y_@_",names(augDataNoMissing))
#     X <- as.matrix(augDataNoMissing[!whichY])
#     Y <- as.matrix(augDataNoMissing[whichY])
#  
#     # Generate lambdas (from SparseTSCGM package):
#     if (missing(lambda_beta) | missing(lambda_kappa)){
#       lams <- SparseTSCGM_lambdas(X, Y, nLambda)
#       if (missing(lambda_beta)){
#         lambda_beta <- lams$lambda_beta
#       }
#       if (missing(lambda_kappa)){
#         lambda_kappa <- lams$lambda_kappa
#       }
#     }
#     
#     Nlambda_beta <- length(lambda_beta)
#     Nlambda_kappa <- length(lambda_kappa)
#     
#     
#     # Expand lambda grid:
#     lambdas <- expand.grid(kappa = lambda_kappa, beta = lambda_beta)
#     Estimates <- vector("list", nrow(lambdas))
#     
#     ### Algorithm 2 of Rothmana, Levinaa & Ji Zhua
#     if (verbose){
#       pb <- txtProgressBar(0, nrow(lambdas), style = 3) 
#     }
#     for (i in seq_len(nrow(lambdas))){
#       Estimates[[i]] <- Rothmana(X, Y, lambdas$beta[i],lambdas$kappa[i], gamma=gamma,maxit.in=maxit.in, maxit.out = maxit.out, penalize.diagonal=penalize.diagonal)
#       if (verbose){
#         setTxtProgressBar(pb, i)
#       } 
#     }
#     if (verbose){
#       close(pb)
#     }
#     
#     #   
#     #   logandbic <- LogLik_and_BIC(data_l, data_c, Estimates)
#     #   lambdas$bic <- logandbic$BIC
#     #   lambdas$loglik <- logandbic$logLik
#     lambdas$ebic <- sapply(Estimates,'[[','EBIC')
#     # Which minimal BIC:
#     min <- which.min(lambdas$ebic)
#     Results <- Estimates[[min]]
#     
#     # Standardize matrices (Wild et al. 2010)
#     # partial contemporaneous correlation (PCC) 
#     # Results$PCC <- computePCC(Results$kappa)
#     # Results$PDC <- computePDC(Results$beta, Results$kappa)  
#     
#     Results$path <- lambdas
#     Results$labels <- colnames(data)
#     
#     # colnames(Results$beta) <- rownames(Results$beta) <- colnames(Results$kappa) <- rownames(Results$kappa) <-
#       # colnames(Results$PCC) <- rownames(Results$PCC) <- colnames(Results$PDC) <- rownames(Results$PDC) <-
#       # Results$labels
#     Results$gamma <- gamma
#     
#     class(Results) <- "AGV"
#     
#     return(Results)
#   }
