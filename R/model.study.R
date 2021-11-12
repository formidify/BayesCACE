#' This function generates the model code for a single study using the 
#' likelihood and model specified in the paper Section 2.1, or a two-step 
#' approach for meta-analysis with complete compliance information as 
#' described in the paper Section 2.2 "The two-step approach".
#' @title Model code of CACE analysis for a single study, or a two-step approach for meta-analysis 
#' with complete complice information
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
#' model.string <- model.study()
model.study <- function(prior.type="default", random.effects = list()){

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

    if(prior.type == "default"){ 
    string1 <- "model{
      prob[1] <- (pi.n*(1-s1) + pi.c*(1-v1))
      prob[2] <- (pi.n*s1 + pi.c*v1)
      prob[3] <- (pi.a*(1-b1))
      prob[4] <- (pi.a*b1)
      prob[5] <- (pi.n*(1-s1))
      prob[6] <- (pi.n*s1)
      prob[7] <- (pi.c*(1-u1)+pi.a*(1-b1))
      prob[8] <- (pi.c*u1+pi.a*b1)
      
      R[1:4] ~ dmulti(prob[1:4], N0)
      R[5:8] ~ dmulti(prob[5:8], N1)
      
      pi.n <- exp(n)/(1+exp(n)+exp(a))
      pi.a <- exp(a)/(1+exp(n)+exp(a))
      pi.c <- 1-pi.a-pi.n
      probit(u1) <- alpha.u
      probit(v1) <- alpha.v
      logit(s1) <- alpha.s
      logit(b1) <- alpha.b
      
      CACE <- u1-v1
    "
    string2 <-
      "# priors
      n ~ dnorm(0, 0.16)
      a ~ dnorm(0, 0.16)
      alpha.s ~ dnorm(0, 0.25)
      alpha.b ~ dnorm(0, 0.25)
      alpha.u ~ dnorm(0, 0.25)
      alpha.v ~ dnorm(0, 0.25)
      }
      "
    modelstring <- paste(string1, string2, sep="\n")
    }

      
    else if (prior.type == "custom"){
      string2 <- prior.study(prior.type)
      modelstring <- paste(string1, string2, sep="\n")
    }
      
    if(!is.element(prior.type,c("default", "custom"))){
      stop("specified prior type should be either 'default' or 'custom'.")
    }
      
    return(modelstring)
}
