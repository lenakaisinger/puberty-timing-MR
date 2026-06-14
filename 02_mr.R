# =============================================================================
# Title: Likely Causal Effects of Circulating Hormones on Puberty Timing:
#        a Mendelian Randomisation Study
# Description: Code to run Mendelian Randomisation. 
# Author: Lena Kaisinger (MR originally written by Felix R. Day)
# =============================================================================

# Dependencies
library(tidyverse)

# ---------------------------------------------------------------------------
# Load and prepare input data
# ---------------------------------------------------------------------------
# Input: harmonised data file produced by the data preparation script

input_data <- read_csv("/your/file/path/formatted_data.csv", col_names = TRUE) %>%
  transmute(
    SNP            = SNP,
    Beta.exposure  = beta.exposure,
    SE.exposure    = se.exposure,
    Beta.outcome   = beta.outcome,
    SE.outcome     = se.outcome,
    Beta.confounder = beta.confounder, # drop these if you have no covariates
    SE.confounder  = se.confounder     # drop these if you have no covariates
  )

# ---------------------------------------------------------------------------
# Analysis settings
# ---------------------------------------------------------------------------

# Optional analyses (set to 1 to enable, 0 to disable)
Radial <- 1                 # Radial MR
Multi <- 1                  # Multivariable MR
Steiger <- 1                # Steiger filtering

inputs <- c("exposure")     # exposure variable name
outputs <- c("outcome")     # outcome variable name
covars <- c("confounder")   # only use if Multi == 1
save <- "Results"           # results file prefix

# Seed for bootstrapped standard errors used in WM and PWM models
set.seed(30)

if (Steiger == 1) {
  input_data <- input_data %>%
    filter(steiger_dir != FALSE)
}

# ---------------------------------------------------------------------------
# MR Analysis pipeline (developed by Felix R. Day)
# ---------------------------------------------------------------------------

# Install packages if not already available
if(require("meta")){
  print("meta is loaded correctly")
} else {
  print("Trying to install meta")
  install.packages("meta")
  if(require(meta)){
    print("meta installed and loaded")
  } else {
    stop("Could not install meta")
  }
}

if(require("ggplot2")){
  print("ggplot2 is loaded correctly")
} else {
  print("Trying to install ggplot2")
  install.packages("ggplot2")
  if(require(ggplot2)){
    print("ggplot2 installed and loaded")
  } else {
    stop("Could not install ggplot2")
  }
}

if(require("RadialMR")){
  print("RadialMR is loaded correctly")
} else {
  print("Trying to install RadialMR")
  library(devtools)
  install_github("WSpiller/RadialMR")
  if(require(RadialMR)){
    print("RadialMR installed and loaded")
  } else {
    stop("Could not install RadialMR")
  }
}

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Function for weighted medians
weighted.median <- function(betaIV.in, weights.in) {
  betaIV.order = betaIV.in[order(betaIV.in)]
  weights.order = weights.in[order(betaIV.in)]
  weights.sum = cumsum(weights.order)-0.5*weights.order
  weights.sum = weights.sum/sum(weights.order)
  below = max(which(weights.sum<0.5))
  weighted.est = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
    (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
  return(weighted.est) }

# Bootstrapping the weighted median to get SEs
weighted.median.boot <- function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in) {
  med = NULL
  for(i in 1:1000){
    betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
    betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
    betaIV.boot = betaYG.boot/betaXG.boot
    med[i] = weighted.median(betaIV.boot, weights.in)
  }
  return(sd(med)) 
}

# ---------------------------------------------------------------------------
# Initialise results storage
# ---------------------------------------------------------------------------

n1 <- length(inputs)
n2 <- length(outputs)
n <- n1*n2
l <- vector("list", n)
k = 1

# ---------------------------------------------------------------------------
# Main analysis loop
# ---------------------------------------------------------------------------

for (j in inputs) {
  
  # Defining exposure variable names
  BetaIn <- eval(parse(text=paste("input_data$Beta.",j, sep="")))
  SEIn <- eval(parse(text=paste("input_data$SE.",j, sep="")))
  savefile <- paste(save, ".csv", sep="")
  
  for (i in outputs) {
    
    # Defining outcome variable names
    BetaOut <- eval(parse(text=paste("input_data$Beta.",i, sep="")))
    SEOut <- eval(parse(text=paste("input_data$SE.",i, sep="")))
    ColOut <- (parse(text=paste(j, ".", i, sep="")))
    l[k] <- paste(j, ".", i, sep="")
    k = k+1
    
    # Defining graph output options
    graph1 <- paste(j,"_", i, "_funnel.png", sep="")
    graph2 <- paste(j,"_", i, "_forest.png", sep="")
    graph3 <- paste(j,"_", i, "_dosage.png", sep="")
    graph_title <- paste("Effect of SNPs for", j, "on", i, sep=" ")
    x_title <- paste("Effect on", j, sep=" ")
    y_title <- paste("Effect on", i, sep=" ")
    
    # Getting number of SNPs
    n_SNPs <- length(BetaOut[!is.na(BetaOut)])
    
    # IVW analysis 
    betaIVW = summary(lm(BetaOut~BetaIn-1, weights=SEOut^-2))$coef[1,1]
    sebetaIVW = summary(lm(BetaOut~BetaIn-1, weights=SEOut^-2))$coef[1,2]/
      min(summary(lm(BetaOut~BetaIn-1, weights=SEOut^-2))$sigma, 1)
    
    # MR Egger analysis, including Egger intercept
    betaEGGER = summary(lm(BetaOut~BetaIn, weights=SEOut^-2))$coef[2,1]
    sebetaEGGER = summary(lm(BetaOut~BetaIn, weights=SEOut^-2))$coef[2,2]/
      min(summary(lm(BetaOut~BetaIn, weights=SEOut^-2))$sigma, 1)
    interEGGER = summary(lm(BetaOut~BetaIn, weights=SEOut^-2))$coef[1,1]
    seinterEGGER = summary(lm(BetaOut~BetaIn, weights=SEOut^-2))$coef[1,2]
    
    # Weighted median and penalised weighted median analysis
    for_median <- data.frame(BetaIn, SEIn, BetaOut, SEOut)
    for_median <- na.omit(for_median)
    for_median$betaIV = for_median$BetaOut/for_median$BetaIn 
    for_median$weights = (for_median$SEOut/for_median$BetaIn)^-2
    for_median$penalty = pchisq(for_median$weights*(for_median$betaIV-betaIVW)^2, df=1, lower.tail=FALSE)
    for_median$pen.weights = for_median$weights*pmin(1, for_median$penalty*20)
    betaWM = weighted.median(for_median$betaIV, for_median$weights)
    sebetaWM = weighted.median.boot(for_median$BetaIn, for_median$BetaOut,
                                    for_median$SEIn, for_median$SEOut, for_median$weights)
    betaPWM = weighted.median(for_median$betaIV, for_median$pen.weights)
    sebetaPWM = weighted.median.boot(for_median$BetaIn, for_median$BetaOut,
                                     for_median$SEIn, for_median$SEOut, for_median$pen.weights)
    
    # Heterogeneity statistics
    metagen(BetaOut/BetaIn, SEOut/BetaIn)
    CochQp = 1-pchisq(metagen(BetaOut/BetaIn, SEOut/BetaIn)$Q,
                      metagen(BetaOut/BetaIn, SEOut/BetaIn)$df.Q)
    Q <- (metagen(BetaOut/BetaIn, SEOut/BetaIn)$Q)
    Q.df <- (metagen(BetaOut/BetaIn, SEOut/BetaIn)$df.Q)
    Isq <- max(0, ((Q - Q.df)/Q)*100)
    
    # -------------------------------------------------------------------------
    # Radial MR (optional)
    # -------------------------------------------------------------------------
    
    if(Radial == 1) {
      
      for_radial <- data.frame(BetaIn, BetaOut, 
                               SEIn, SEOut, input_data$SNP)
      for_radial <- na.omit(for_radial)
      radial_format <- format_radial(BXG = for_radial$BetaIn, BYG = for_radial$BetaOut, 
                                     seBXG = for_radial$SEIn, seBYG = for_radial$SEOut, 
                                     RSID = for_radial$input_data.SNP)
      radial_res <- ivw_radial(radial_format)	
      radial_data <- merge(input_data, radial_res$data, by.x="SNP", by.y="SNP")
      
      # Save per-SNP radial output
      rad_out <- data.frame(
        SNP = radial_res$data$SNP, 
        Qj_Chi = radial_res$data$Qj_Chi, 
        Outliers = radial_res$data$Outliers)
      
      write.table(rad_out, file = paste0(save, ".", ColOut, ".csv"),
                  sep = ",", row.names = FALSE)
      
      # Run MR analysis without outliers
      radial_data <- radial_data[which(radial_data$Outliers == "Variant"), ]
      
      RBetaIn <- eval(parse(text=paste("radial_data$Beta.",j, sep="")))
      RSEIn <- eval(parse(text=paste("radial_data$SE.",j, sep="")))
      RBetaOut <- eval(parse(text=paste("radial_data$Beta.",i, sep="")))
      RSEOut <- eval(parse(text=paste("radial_data$SE.",i, sep="")))		
      
      Rn_SNPs <- length(RBetaIn)
      
      # IVW
      RbetaIVW = summary(lm(RBetaOut~RBetaIn-1, weights=RSEOut^-2))$coef[1,1]
      RsebetaIVW = summary(lm(RBetaOut~RBetaIn-1, weights=RSEOut^-2))$coef[1,2]/
        min(summary(lm(RBetaOut~RBetaIn-1, weights=RSEOut^-2))$sigma, 1)
      
      # MR Egger
      RbetaEGGER = summary(lm(RBetaOut~RBetaIn, weights=RSEOut^-2))$coef[2,1]
      RsebetaEGGER = summary(lm(RBetaOut~RBetaIn, weights=RSEOut^-2))$coef[2,2]/
        min(summary(lm(RBetaOut~RBetaIn, weights=RSEOut^-2))$sigma, 1)
      RinterEGGER = summary(lm(RBetaOut~RBetaIn, weights=RSEOut^-2))$coef[1,1]
      RseinterEGGER = summary(lm(RBetaOut~RBetaIn, weights=RSEOut^-2))$coef[1,2]
      
      # Weighted median and penalised weighted median analysis
      Rfor_median <- data.frame(RBetaIn, RSEIn, RBetaOut, RSEOut)
      Rfor_median <- na.omit(Rfor_median)
      RbetaIV = Rfor_median$RBetaOut/Rfor_median$RBetaIn	
      Rweights = (Rfor_median$RSEOut/Rfor_median$RBetaIn)^-2
      Rpenalty = pchisq(Rweights*(RbetaIV-RbetaIVW)^2, df=1, lower.tail=FALSE)
      Rpen.weights = Rweights*pmin(1, Rpenalty*20)
      RbetaWM = weighted.median(RbetaIV, Rweights)
      RsebetaWM = weighted.median.boot(Rfor_median$RBetaIn, Rfor_median$RBetaOut, 
                                       Rfor_median$RSEIn, Rfor_median$RSEOut, Rweights)
      RbetaPWM = weighted.median(RbetaIV, Rpen.weights)
      RsebetaPWM = weighted.median.boot(Rfor_median$RBetaIn, Rfor_median$RBetaOut, 
                                        Rfor_median$RSEIn, Rfor_median$RSEOut, Rpen.weights)
      
      # Heterogeneity statistics
      metagen(RBetaOut/RBetaIn, RSEOut/RBetaIn)
      RCochQp = 1-pchisq(metagen(RBetaOut/RBetaIn, RSEOut/RBetaIn)$Q,
                         metagen(RBetaOut/RBetaIn, RSEOut/RBetaIn)$df.Q)
      RQ <- (metagen(RBetaOut/RBetaIn, RSEOut/RBetaIn)$Q)
      RQ.df <- (metagen(RBetaOut/RBetaIn, RSEOut/RBetaIn)$df.Q)
      RIsq <- max(0, ((RQ - RQ.df)/RQ)*100)
      
    }
    
    # -------------------------------------------------------------------------
    # Multivariable MR (optional)
    # -------------------------------------------------------------------------
    
    if(Multi == 1){
      
      # Setting up the multivariate formula
      form <- ""
      f <- length(covars)
      for(b in 1:length(covars)){
        f[b] <- paste("+ input_data$Beta.", paste(covars[b]), sep="") 
        form <- paste(form, paste(f[b]))
      }
      formula_IVW = paste(paste("BetaOut ~ BetaIn"), form, paste("- 1"))
      formula_EGGER = paste(paste("BetaOut ~ BetaIn"), form)
      
      # IVW  
      MbetaIVW = summary(lm(paste(formula_IVW), weights=SEOut^-2), data=input_data)$coef[1,1]
      MsebetaIVW = summary(lm(paste(formula_IVW), weights=SEOut^-2), data=input_data)$coef[1,2]/
        min(summary(lm(formula_IVW, weights=SEOut^-2))$sigma, 1)
      
      Mn_SNPs <- length(lm(formula_IVW, weights=SEOut^-2, data=input_data)$residuals)	
      
      # MR Egger
      MbetaEGGER = summary(lm(paste(formula_EGGER), weights=SEOut^-2), data=input_data)$coef[2,1]
      MsebetaEGGER = summary(lm(paste(formula_EGGER), weights=SEOut^-2), data=input_data)$coef[2,2]/
        min(summary(lm(formula_EGGER, weights=SEOut^-2))$sigma, 1)
      MinterEGGER = summary(lm(formula_EGGER, weights=SEOut^-2), data=input_data)$coef[1,1]
      MseinterEGGER = summary(lm(paste(formula_EGGER), weights=SEOut^-2), data=input_data)$coef[1,2]
      
      if(Radial==1){
        
        rform <- ""
        rf <- length(covars)
        for(b in 1:length(covars)){
          rf[b] <- paste("+ radial_data$Beta.", paste(covars[b]), sep="") 
          rform <- paste(rform, paste(rf[b]))
        }
        Rformula_IVW = paste(paste("RBetaOut ~ RBetaIn"), rform, paste("- 1"))
        Rformula_EGGER = paste(paste("RBetaOut ~ RBetaIn"), rform)
        
        # IVW 
        RMbetaIVW = summary(lm(paste(Rformula_IVW), weights=RSEOut^-2), data=radial_data)$coef[1,1]
        RMsebetaIVW = summary(lm(paste(Rformula_IVW), weights=RSEOut^-2), data=radial_data)$coef[1,2]/
          min(summary(lm(Rformula_IVW, weights=RSEOut^-2))$sigma, 1)
        
        RMn_SNPs <- length(lm(Rformula_IVW, weights=RSEOut^-2, data=radial_data)$residuals)
        
        # MR Egger
        RMbetaEGGER = summary(lm(paste(Rformula_EGGER), weights=RSEOut^-2), data=radial_data)$coef[2,1]
        RMsebetaEGGER = summary(lm(paste(Rformula_EGGER), weights=RSEOut^-2), data=radial_data)$coef[2,2]/
          min(summary(lm(Rformula_EGGER, weights=RSEOut^-2))$sigma, 1)
        RMinterEGGER = summary(lm(Rformula_EGGER, weights=RSEOut^-2), data=radial_data)$coef[1,1]
        RMseinterEGGER = summary(lm(paste(Rformula_EGGER), weights=RSEOut^-2), data=radial_data)$coef[1,2]
      }	
    }			
    
    # -------------------------------------------------------------------------
    # Plots
    # -------------------------------------------------------------------------
    
    # Funnel plots
    png(file=graph1)
    plot(BetaOut/BetaIn, 1/(SEOut/BetaIn), 
         xlab=y_title, ylab="Instumental Variable Strength")
    for (b in 1:length(BetaIn)) {
      lines(c((BetaOut[b]-1.96*SEOut[b])/BetaIn[b], (BetaOut[b]+1.96*SEOut[b])/BetaIn[b]),
            c(BetaIn[b]/SEOut[b], BetaIn[b]/SEOut[b]))
    }
    abline(h=0, lwd=1); abline(v=0, lwd=1)
    abline(v=betaIVW, lwd=1, lty=2)
    dev.off()
    
    # Dosage plots
    dosage_data <- data.frame(BetaIn, SEIn, BetaOut, SEOut)
    dosage_data$inLCI <- dosage_data$BetaIn+(1.96*dosage_data$SEIn)
    dosage_data$inUCI <- dosage_data$BetaIn-(1.96*dosage_data$SEIn)
    dosage_data$outLCI <- BetaOut-(1.96*SEOut)
    dosage_data$outUCI <- BetaOut+(1.96*SEOut)
    
    cols <- c("IVW"="red", "EGGER"="blue", "WM"="yellow", "PWM"="orange")
    dose_graph <- ggplot(data = dosage_data, aes(x = BetaIn, y = BetaOut, show.legend=TRUE)) + 
      geom_point() + 
      theme_bw() +
      geom_errorbar(aes(ymin = outLCI, ymax = outUCI)) + 
      geom_errorbarh(aes(xmin = inLCI, xmax = inUCI)) +
      geom_abline(intercept = 0, slope = 0, size=1) +
      geom_abline(intercept = 0, slope = betaIVW, size=1, colour = "red", show.legende=FALSE) +
      geom_abline(intercept = interEGGER, slope = betaEGGER, size=1, colour = "blue", show.legend=FALSE) +
      geom_abline(intercept = 0, slope = betaWM, size=1, colour = "yellow", show.legend=FALSE) +
      geom_abline(intercept = 0, slope = betaPWM, size=1, colour = "orange", show.legend=FALSE) +
      geom_vline(xintercept = 0, size=1) +
      ggtitle(graph_title) + xlab(x_title) + ylab(y_title) +
      scale_colour_manual(name="MR estimates", values=cols)
    
    ggsave(filename=graph3, plot=dose_graph)
    rm(dosage_data)
    
    # -------------------------------------------------------------------------
    # P-values
    # -------------------------------------------------------------------------
    
    pIVW <- 2*pt(-abs(betaIVW/sebetaIVW), df=n_SNPs-1)
    pEGGER <- 2*pt(-abs(betaEGGER/sebetaEGGER), df=n_SNPs-1)
    pinterEGGER <- 2*pt(-abs(interEGGER/seinterEGGER), df=n_SNPs-1)
    pWM <- 2*pt(-abs(betaWM/sebetaWM), df=n_SNPs-1)
    pPWM <- 2*pt(-abs(betaPWM/sebetaPWM), df=n_SNPs-1)
    
    # -------------------------------------------------------------------------
    # Compile results
    # -------------------------------------------------------------------------
    
    base_row <- c(paste(ColOut), n_SNPs,
                betaIVW, sebetaIVW, pIVW, CochQp, Isq,
                betaEGGER, sebetaEGGER, pEGGER, interEGGER, seinterEGGER, pinterEGGER,
                betaWM, sebetaWM, pWM, 
                betaPWM, sebetaPWM, pPWM)
    
    outout <- base_row
    
    if(Radial == 1){
      RpIVW <- 2*pt(-abs(RbetaIVW/RsebetaIVW), df=Rn_SNPs-1)
      RpEGGER <- 2*pt(-abs(RbetaEGGER/RsebetaEGGER), df=Rn_SNPs-1)
      RpinterEGGER <- 2*pt(-abs(RinterEGGER/RseinterEGGER), df=Rn_SNPs-1)
      RpWM <- 2*pt(-abs(RbetaWM/RsebetaWM), df=Rn_SNPs-1)
      RpPWM <- 2*pt(-abs(RbetaPWM/RsebetaPWM), df=Rn_SNPs-1)
      
      radial_row <- c(paste("Radial filtered"), Rn_SNPs, RbetaIVW, RsebetaIVW, RpIVW, RCochQp, RIsq,
                        RbetaEGGER, RsebetaEGGER, RpEGGER, RinterEGGER, RseinterEGGER, RpinterEGGER,
                        RbetaWM, RsebetaWM, RpWM, RbetaPWM, RsebetaPWM, RpPWM)
      
      output <- rbind(base_row, radial_row)
    }	
    
    if(Multi == 1){
      MpIVW <- 2*pt(-abs(MbetaIVW/MsebetaIVW), df=Mn_SNPs-(1 + length(covars)))
      MpEGGER <- 2*pt(-abs(MbetaEGGER/MsebetaEGGER), df=Mn_SNPs-(1 + length(covars)))
      MpinterEGGER <- 2*pt(-abs(MinterEGGER/MseinterEGGER), df=Mn_SNPs-(1 + length(covars)))
      
      multi_row <- c(paste("Multivariate"), Mn_SNPs, MbetaIVW, MsebetaIVW, MpIVW, "NA", "NA",
                        MbetaEGGER, MsebetaEGGER, MpEGGER, MinterEGGER, MseinterEGGER, MpinterEGGER,
                        "NA", "NA", "NA", "NA", "NA")
      
      output <- rbind(base_row, multi_row)
      
      if(Radial==1){
        RMpIVW <- 2*pt(-abs(RMbetaIVW/RMsebetaIVW), df=RMn_SNPs-(1 + length(covars)))
        RMpEGGER <- 2*pt(-abs(RMbetaEGGER/RMsebetaEGGER), df=RMn_SNPs-(1 + length(covars)))
        RMpinterEGGER <- 2*pt(-abs(RMinterEGGER/RMseinterEGGER), df=RMn_SNPs-(1 + length(covars)))
        
        # Setting up vector of with-in loop results
        multi_radial_row <- c(paste("Multivariate and Radial filtered"), RMn_SNPs, RMbetaIVW, RMsebetaIVW, RMpIVW, "NA", "NA",
                          RMbetaEGGER, RMsebetaEGGER, RMpEGGER, RMinterEGGER, RMseinterEGGER, RMpinterEGGER,
                          "NA", "NA", "NA", "NA", "NA")
        
        output <- rbind(base_row, radial_row, multi_row, multi_radial_row)
      }				
    }	
    
    assign(paste(ColOut), output)
    
  }
}

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

# Formatting results table
results <- c("IV", "n_SNPs", "betaIVW", "sebetaIVW", "pIVW", "CochQp", "Isq", 
             "betaEGGER", "sebetaEGGER", "pEGGER", "interEGGER", "seinterEGGER", "pinterEGGER",
             "betaWM", "sebetaWM", "pWM", "betaPWM", "sebetaPWM", "pPWM")

for (p in l) { 
  q <- eval(parse(text=paste(p)))
  results <- rbind(results, q)
}

write.table(results, file=savefile, sep=",", row.names=FALSE, col.names=FALSE)

