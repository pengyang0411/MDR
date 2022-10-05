#' @importFrom stats dbinom pbinom
#'
##----------------------##
## Define Joint func    ##
##----------------------##

## Complementary of the Joint probability of low tox and high eff
## given the true tox and eff
Joint_LH <- function(mT, mE, n1,
                     phi_T, phi_E, p){

  return((1 - pbinom(mT, n1, phi_T) * (1 - pbinom(mE - 1, n1, phi_E)))^p)

}

## Joint probability of low tox and high eff given the true tox and eff
Joint_LH_power <- function(mT, mE, n1,
                           phi_T, phi_E, p){
  return((pbinom(mT, n1, phi_T) * (1 - pbinom(mE - 1, n1, phi_E)))^p)
}


##----------------------##
## Local Type I error   ##
##----------------------##

## Local type I error function give s and k
Prob_reject_H0sk <- function(mT, mE, n1, J,
                             s, k,
                             phi_TL, phi_TH,
                             phi_EL, phi_EH){

  res_LL <- Joint_LH(mT, mE, n1, phi_TL, phi_EL, s)

  res_HL <- Joint_LH(mT, mE, n1, phi_TH, phi_EL, k - s)

  res_HH <- Joint_LH(mT, mE, n1, phi_TH, phi_EH, J - k)

  return(1 - res_LL*res_HL*res_HH)

}


##----------------------##
## Global Type I error  ##
##----------------------##

## Global type I error function iterating over s and k
TypeI_Global <- function(mT, mE, n1, J,
                         phi_TL, phi_TH,
                         phi_EL, phi_EH){

  res <- c()

  for(s in 0:J){

    for(k in s:J){

      res <- rbind(res, Prob_reject_H0sk(mT, mE, n1, J,
                                         s, k,
                                         phi_TL, phi_TH,
                                         phi_EL, phi_EH))

    }
  }

  return(apply(res, 2, max))

}

##----------------------##
## Power function       ##
##----------------------##

## Define the general power function under LFC condition
power_LFC <- function(j, mT, mE, n1, delta, J = 3,
                      phi_TL, phi_TH,
                      phi_EL, phi_EH){


  res <- 0

  ## Under LFC condition where OBD is not the last dose
  if(j != J){

    for(y_J in 0:n1){

      ## First j - 1 doses are with low efficacy
      res <- res + pbinom(mE - 1, n1, phi_EL)^(j-1) *

      ## Dose j with low toxicity and high efficacy and meet delta
        Joint_LH_power(mT, max(mE, floor(y_J - delta*n1)), n1, phi_TL, phi_EH, 1) *

      ## Last J - j doses are with high toxicity
        (1 - pbinom(mT, n1, phi_TH))^(J - j) *

        dbinom(y_J, n1, phi_EH)

    }

  }else{

    ## Under LFC condition where OBD is not the last dose
    ## First J - 1 doses are with low efficacy
    res <- pbinom(mE - 1, n1, phi_EL)^(j-1) *
    ## Dose J with low toxicity and high efficacy
      Joint_LH_power(mT, mE, n1, phi_TL, phi_EH, 1)

  }

  return(res)

}

##-------------------------------##
##       Main function           ##
##-------------------------------##

#' @title Find the minimized sample size and decision rules
#'
#' @description This function is designed to estimate the minimized sample size
#' and decision boundaries of toxicity and efficacy endpoints for trial
#' enrollment and decision-making.
#'
#' @param J Number of doses considered in the trial
#' @param phi_TL True tolerable toxicity rate for each dose
#' @param phi_TH True insuffurable toxicity rate for each dose
#' @param phi_EL True insignificant efficacy rate for each dose
#' @param phi_EH True desirable efficacy rate for each dose
#' @param alpha Targeted type I error rate
#' @param beta Targeted general power
#' @param delta Targeted plateau criteria
#' @param Nmin Minimal number of patients considered
#' @param Nmax Maximum number of patients considered
#' @param if.verbose Display extended information
#'
#' @return A list contains the estimated design parameters.
#' \itemize{
#'     \item best_N - Estimated minimal sample size.
#'     \item best_mT - Estimated decision rule for toxicity endpoints (counts).
#'     \item best_mTr - Estimated decision rule for toxicity endpoints (rate).
#'     \item best_mE - Estimated decision rule for efficacy endpoints (counts).
#'     \item best_mEr - Estimated decision rule for efficacy endpoints (rate).
#'     \item best_Ind - Index of the estimated design parameter in Metadata.
#'     \item meta.data - Metadata for the estimation.
#' }
#'
#' @export
#'
#' @examples
#' res <- FindSampleSize(J = 3, delta = 0.2,
#'                       phi_TL = 0.2, phi_TH = 0.4,
#'                       phi_EL = 0.2, phi_EH = 0.4,
#'                       alpha = 0.1, beta = 0.8,
#'                       Nmin = 10, Nmax = 100, if.verbose = FALSE)
#'
FindSampleSize <- function(J,                ## Number of dose in design
                           phi_TL,           ## Tolerable Toxicity
                           phi_TH,           ## Insuffurable Toxicity
                           phi_EL,           ## Non-significant Efficacy
                           phi_EH,           ## Desirable Efficacy
                           alpha = 0.1,      ## Size for global test
                           beta = 0.8,       ## Targeted power for stage I
                           delta  = 0.2,     ## Plateau criteria
                           Nmin,             ## Min number of samples considered
                           Nmax,             ## Max number of samples considered
                           if.verbose = FALSE
){

  ##--------Check Input data-----------##

  if(J != round(J)){
    stop('Number of dose must be an integer!')
  }

  if(phi_TL >= phi_TH){
    stop('phi_TL must be smaller than phi_TH')
  }

  if(phi_EL >= phi_EH){
    stop('phi_EL must be smaller than phi_EH')
  }

  # Save meta data
  meta.data <- c()

  ##---Iteratively finding sample size---##

  best_N <- best_mT <- best_mE <- NULL

  ## Searching the sample size
  for(n in Nmin:Nmax){

    ## Searching sapce for fixed cutoffs
    m_grid <- 1:n

    ## Searching cutoff for Efficacy
    for(mE in m_grid){

      # ## Evaluating Global Type I error
      Alpha_Global <- TypeI_Global(mT = m_grid, mE, n1 = n, J,
                                   phi_TL, phi_TH,
                                   phi_EL, phi_EH) #/ J

      ## For the efficacy meets type I error
      Ind_eff <- which(Alpha_Global <= alpha)

      if(length(Ind_eff) == 0) next

      ## Evaluating Power under LFC condition
      ## Selecting the smallest power
      Power <- do.call(cbind,
                       lapply(X = 1:J, FUN = power_LFC,
                              mT = m_grid[Ind_eff], mE, n1 = n,
                              delta = delta, J = J,
                              phi_TL = phi_TL, phi_TH = phi_TH,
                              phi_EL = phi_EL, phi_EH = phi_EH))

      ## GLFC
      GLFC <- apply(Power, 1, FUN = function(x) min(which(x == min(x))))

      ## Obtian the minimized power acorss all LFC conditions
      Power <- apply(Power, 1, min)


      ## Tmp dataframe for meta.data
      n.tmp <- length(Ind_eff)
      tmp   <- data.frame(n = rep(n, n.tmp),
                          mT = m_grid[Ind_eff],
                          cutT = round(m_grid[Ind_eff]/n,3),
                          mE = rep(mE, n.tmp),
                          cutE = round(mE/n,3),
                          sizeG = round(Alpha_Global[Ind_eff],2),
                          Power = Power,
                          GLFC = GLFC)

      ## Print verbose
      if(if.verbose){

        print(tmp)

      }


      ## Save Meta data
      meta.data = rbind(meta.data, tmp)

      if(max(Power) >= beta){

        ## Find the best sample size!
        bestN <- n

        cat('Best sample size :', bestN, '\n')

        break

      }

    }

    if(max(Power) >= beta) break

  }

  ## Name meta.data

  colnames(meta.data) <- c('n', 'mT', 'mT-rate', 'mE', 'mE-rate',
                           'alpha', 'power', 'GLFC')
  #
  # meta.data <- as.data.frame(meta.data)

  ## Find the best power
  Index <- which(meta.data$power >= beta & meta.data$alpha <= alpha)

  if(length(Index) > 0){
    best_N   <- meta.data[Index[1], 1]
    best_mT  <- meta.data[Index[1], 2]
    best_mE  <- meta.data[Index[1], 3]
  }else{
    stop('No minimized sample size obtianed, try to enlarge the searching
         bounder Nmax')
  }

  return(list(best_N    = meta.data[Index[1], 1],
              best_mT   = meta.data[Index[1], 2],
              best_mTr  = meta.data[Index[1], 3],
              best_mE   = meta.data[Index[1], 4],
              best_mEr  = meta.data[Index[1], 5],
              best_Ind  = Index,
              meta.data = meta.data))

}

