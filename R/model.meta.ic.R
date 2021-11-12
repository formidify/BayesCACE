#' This function generates the model code for meta-analysis 
#' when the dataset has incomplete compliance information for all studies, 
#' as described in the paper Section 2.2.2, the Bayesian hierarchical model.
#' @title Bayesian hierarchical model code for CACE meta-analysis with complete compliance data
#' @param prior.type type of prior information. Default to "default".
#' @param random.effects a list of logical values indicating whether random effects are included in the model.
#' The list should contain the assignment for these parameters only: \code{delta.n} (\eqn{\delta_{in}}), 
#' \code{delta.a} (\eqn{\delta_{ia}}), \code{delta.u} (\eqn{\delta_{iu}}), \code{delta.v} (\eqn{\delta_{iv}}), 
#' \code{delta.s} (\eqn{\delta_{is}}), \code{delta.b} (\eqn{\delta_{ib}}), \code{cor}. The list should be in the
#' form of \code{list(delta.a = FALSE, cor = FALSE, ...)}. By default, this
#' is an empty list, and all parameters are default to \code{TRUE}. Parameters that are not listed in the list
#' are assumed to be \code{TRUE}. Note that \eqn{\rho} (\code{cor}) can only be included when both \eqn{\delta_{in}} 
#' (\code{delta.n}) and \eqn{\delta_{ia}} (\code{delta.a}) are set to \code{TRUE}. Otherwise, a warning 
#' occurs and the model continues running by forcing \code{delta.n = TRUE} and \code{delta.a = TRUE}. 
#' @return It returns a model string
#' @export
#' @examples
#' # use default settings
#' model.string <- model.meta.ic()
model.meta.ic <- function(prior.type="default", random.effects = list()){

    delta.n <- delta.a <- delta.u <- delta.v <- delta.s <- delta.b <- cor <- TRUE
    if ("delta.n" %in% names(random.effects)) {delta.n <- random.effects[['delta.n']]}
    if ("delta.a" %in% names(random.effects)) {delta.a <- random.effects[['delta.a']]}
    if ("delta.u" %in% names(random.effects)) {delta.u <- random.effects[['delta.u']]}
    if ("delta.v" %in% names(random.effects)) {delta.v <- random.effects[['delta.v']]}
    if ("delta.s" %in% names(random.effects)) {delta.s <- random.effects[['delta.s']]}
    if ("delta.b" %in% names(random.effects)) {delta.b <- random.effects[['delta.b']]}
    if ("cor" %in% names(random.effects)) {delta.n <- random.effects[['cor']]}
    
    if ((!(delta.n & delta.a)) & cor){
      warning("'cor' can be assigned as TRUE only if both delta.n and delta.a are TRUE.\n
              the model is continued by forcing delta.n=TRUE and delta.a=TRUE")
      delta.n <- TRUE
      delta.a <- TRUE
    }
    
    Ind <- rep(1, 7)
    if (!delta.n) Ind[1] <- 0
    if (!delta.a) Ind[2] <- 0
    if (!delta.u) Ind[3] <- 0
    if (!delta.v) Ind[4] <- 0
    if (!delta.s) Ind[5] <- 0
    if (!delta.b) Ind[6] <- 0
    if (!cor) Ind[7] <- 0

    string1 <-
    "model{
      for (i in 1:n1) {

      R0_4[i, ] ~ dmulti(prob0[i, ], N0_4[i])
      prob0[i, 1] <- pi.n[i]*(1-s1[i]) + pi.c[i]*(1-v1[i])
      prob0[i, 2] <- pi.n[i]*s1[i] + pi.c[i]*v1[i]
      prob0[i, 3] <- pi.a[i]*(1-b1[i])
      prob0[i, 4] <- pi.a[i]*b1[i]
      
      R1_4[i, ] ~ dmulti(prob1[i, ], N1_4[i])
      prob1[i, 1] <- pi.n[i]*(1-s1[i])
      prob1[i, 2] <- pi.n[i]*s1[i]
      prob1[i, 3] <- pi.c[i]*(1-u1[i])+pi.a[i]*(1-b1[i])
      prob1[i, 4] <- pi.c[i]*u1[i]+pi.a[i]*b1[i]
    } 
      
      for (i in (n1+1):(n1+n2)) {
      
      R0_4[i, ] ~ dmulti(prob0[i, ], N0_4[i])
      prob0[i, 1] <- pi.n[i]*(1-s1[i]) + pi.c[i]*(1-v1[i])
      prob0[i, 2] <- pi.n[i]*s1[i] + pi.c[i]*v1[i]
      prob0[i, 3] <- pi.a[i]*(1-b1[i])
      prob0[i, 4] <- pi.a[i]*b1[i]
      
      R1_2[i, 2] ~ dbin(p1[i], N1_2[i])
      p1[i] <- pi.n[i]*s1[i] + pi.c[i]*u1[i] + pi.a[i]*b1[i]
      } 
      
      for (i in (n1+n2+1):(n1+n2+n3)) {
      
      R0_2[i, 2] ~ dbin(p0[i], N0_2[i])
      p0[i] <-  pi.n[i]*s1[i] + pi.c[i]*v1[i] + pi.a[i]*b1[i]
      
      R1_4[i, ] ~ dmulti(prob1[i, ], N1_4[i])
      prob1[i, 1] <- pi.n[i]*(1-s1[i])
      prob1[i, 2] <- pi.n[i]*s1[i]
      prob1[i, 3] <- pi.c[i]*(1-u1[i])+pi.a[i]*(1-b1[i])
      prob1[i, 4] <- pi.c[i]*u1[i]+pi.a[i]*b1[i]
      } 
      
      for (i in (n1+n2+n3+1):(n1+n2+n3+n4)) {
      
      R0_2[i, 2] ~ dbin(p0[i], N0_2[i])
      p0[i] <-  pi.n[i]*s1[i] + pi.c[i]*v1[i] + pi.a[i]*b1[i]
      
      R1_2[i, 2] ~ dbin(p1[i], N1_2[i])
      p1[i] <- pi.n[i]*s1[i] + pi.c[i]*u1[i] + pi.a[i]*b1[i]
      } 
      
      for (i in 1:(n1+n2+n3+n4)) {

      n[i] <- alpha.n + Ind[1]*delta.n[i]
      a[i] <- alpha.a + Ind[2]*delta.a[i]
      pi.n[i] <- exp(n[i])/(1+exp(n[i])+exp(a[i]))
      pi.a[i] <- exp(a[i])/(1+exp(n[i])+exp(a[i]))
      pi.c[i] <- 1-pi.a[i]-pi.n[i]
      probit(u1[i]) <- alpha.u + Ind[3]*delta.u[i]
      probit(v1[i]) <- alpha.v + Ind[4]*delta.v[i]
      logit(s1[i]) <- alpha.s + Ind[5]*delta.s[i]
      logit(b1[i]) <- alpha.b + Ind[6]*delta.b[i]

      cacei[i] <- u1[i]-v1[i]
    "
    string3 <- 
    "pin <- exp(alpha.n)/(1+exp(alpha.n)+exp(alpha.a))
    pia <- exp(alpha.a)/(1+exp(alpha.n)+exp(alpha.a))
    pic <- 1-pia-pin

    # priors
    alpha.n ~  dnorm(0, 0.16)
    alpha.a ~ dnorm(0, 0.16)
    alpha.s ~  dnorm(0, 0.25)
    alpha.b ~  dnorm(0, 0.25)
    alpha.u ~  dnorm(0, 0.25)
    alpha.v ~  dnorm(0, 0.25)
    "
    if(!is.element(prior.type,c("default"))){
      stop("specified prior type is wrong.")
    }

    if(Ind[7]==1){
      string2_0 <- 
    "delta.n[i] <- delta.rho[1, i]
    delta.a[i] <- delta.rho[2, i]
    delta.rho[1:2, i] ~ dmnorm(c(0, 0), Omega.rho)"
      string4_0 <- 
    "II[1,1] <- 1
    II[2,2] <- 1
    II[1,2] <- 0
    II[2,1] <- 0
    Omega.rho ~  dwish (II[,], 3)
    Sigma.rho <- inverse(Omega.rho)
    sigma.n <- Sigma.rho[1, 1]
    sigma.a <- Sigma.rho[2, 2]
    rho <- Sigma.rho[1, 2]"
    string2_1 <- string2_2 <- ""
    string4_1 <- string4_2 <- ""
    }
      
    else if (Ind[7]==0){
      string2_0 <- string4_0 <- ""
      
      if(Ind[1]==1){
        string2_1 <- "delta.n[i] ~ dnorm(0, tau.n)"
        string4_1 <- 
      "tau.n ~ dgamma(2, 2)
    sigma.n <- 1/sqrt(tau.n)"
      }
      else if (Ind[1]==0){
        string2_1 <- ""
        string4_1 <- ""
      }
      
      if(Ind[2]==1){
        string2_2 <- "delta.a[i] ~ dnorm(0, tau.a)"
        string4_2 <- 
      "tau.a ~ dgamma(2, 2)
    sigma.a <- 1/sqrt(tau.a)"
      }
      else if (Ind[2]==0){
        string2_2 <- ""
        string4_2 <- ""
      }
    }  
      
    if(Ind[3]==1){
      string2_3 <- "delta.u[i] ~ dnorm(0, tau.u)"
      string4_3 <- 
    "u1out <- phi(alpha.u/sqrt(1+sigma.u^2))
    tau.u ~ dgamma(2, 2)
    sigma.u <- 1/sqrt(tau.u)"
    }
    else if (Ind[3]==0){
      string2_3 <- ""
      string4_3 <- "u1out <- phi(alpha.u)"
    }
      
    if(Ind[4]==1){
      string2_4 <- "delta.v[i] ~ dnorm(0, tau.v)"
      string4_4 <- 
    "v1out <- phi(alpha.v/sqrt(1+sigma.v^2))
    CACE <- u1out-v1out
    tau.v ~ dgamma(2, 2)
    sigma.v <- 1/sqrt(tau.v)"
    }
    else if (Ind[4]==0){
      string2_4 <- ""
      string4_4 <- 
    "v1out <- phi(alpha.v)
    CACE <- u1out-v1out"
    }
      
    if(Ind[5]==1){
      string2_5 <- "delta.s[i] ~ dnorm(0, tau.s)"
      string4_5 <- 
    "s1out <- ilogit(alpha.s/sqrt(1 + (16^2*3/(15^2*pi^2))*sigma.s^2))
    tau.s ~ dgamma(2, 2)
    sigma.s <- 1/sqrt(tau.s)"
    }
    else if (Ind[5]==0){
      string2_5 <- ""
      string4_5 <- "s1out <- ilogit(alpha.s)"
    }  
      
    if(Ind[6]==1){
      string2_6 <- "delta.b[i] ~ dnorm(0, tau.b)
    }"
      string4_6 <- 
    "b1out <- ilogit(alpha.b/sqrt(1 + (16^2*3/(15^2*pi^2))*sigma.b^2))
    tau.b ~ dgamma(2, 2)
    sigma.b <- 1/sqrt(tau.b)
    }"
    }
    else if (Ind[6]==0){
      string2_6 <- "
    }"
      string4_6 <- "b1out <- ilogit(alpha.b)
    }"
    }   

      
    modelstring <- paste(string1, string2_0, string2_1, string2_2, string2_3, string2_4, string2_5, string2_6,
                         string3, string4_0, string4_1, string4_2, string4_3, string4_4, string4_5, string4_6, sep="\n")

    return(modelstring)
}
