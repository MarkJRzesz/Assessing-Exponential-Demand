# Libararies ----
library(beezdemand)
library(nls.multstart)
library(nlme)
library(emmeans)
library(psych)
library(weights)
library(VGAM)
library(tidyverse)
library(data.table)
library(AICcmodavg)

citation()
version

# Handy Functions ----
rSquarer <- function(beep, data){
  if(is.null(beep) || any(is.na(beep))) NA 
  else
    1 - sum(((predict(beep)) - na.omit(data$y))^2) / sum((mean(na.omit(data$y)) - na.omit(data$y))^2)
}

MAEmaker <- function(beep, data){
  if(is.null(beep) || any(is.na(beep))) NA 
  else
    mean(abs(predict(beep) - na.omit(data$y)))
}

AICcmaker <- function(beep, data){
  if(is.null(beep) || any(is.na(beep))) NA 
  else AICc(beep)
}


tickFunc <- function(from = .0002, to = .001, jumps = 6){
  tickHolder <- vector("numeric")
  for(a in 1:jumps){
    ticks <- seq(from*(10^a), to*(10^a), by = to*(10^a)/10)
    tickHolder<- append(tickHolder, ticks)
  }
  return(tickHolder)
}

Pmaxxer <- function(q0, alpha){
  1/(10^alpha*10^q0)
}

Omaxxer <- function(q0, alpha){
  Pmax <- Pmaxxer(q0, alpha)
  Pmax * (10^q0 * exp(-10^alpha * 10^q0 * Pmax))
}

tryCoef <- function(x){
    trying <- tryCatch(coef(x), error = function(e) c(q0= NA, alpha = NA))
    return(trying)
}

EXPDToSND <- function(k, base = 10){
    - 1/(log(1 - (1/(k*log(base)))))
}


# Consumption Expenditure K examples  ----
svg("expenditure.svg", 7.5, 7.5)
par(mfrow = c(2,2), pty = "s", mar = c(4.5, 4, 4, 2))
q0 <- 20
k <- 2
alpha <- .005
Pmax <- -VGAM::lambertW(-1/log(10^(k)))/(alpha*q0)
Pmax

curve(q0 * 10^(k*(exp(-alpha*q0*x)- 1)), ylim = c(0, 25), xlim = c(.001, 1000), col = "orange", 
      log = "x", main = expression(bold(atop("Consumption Curves,","Fixed "*Q[0]*" & "*alpha*", Different k"))), cex.main = 1.5,
      ylab = NA, xlab = NA, axes = FALSE, 
      from = .0001, to = 10000, lty = 1, lwd = 1.5)
title(ylab = "Units Consumed", xlab = "Cost per Unit", line = 2.5)
legend("topright", legend = c("k = 2", "k = 1.5", "k = 1.15", "SND"), ncol = 1, xpd = NA,
       bty = "n", lty = c(1, 5, 3, 5), col = c(rep("orange",3), "blue"),
       xjust = 0)

lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .8)
curve(q0 * 10^(1.15*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 3, lwd = 1.5)
curve(q0 * 10^(1.5*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 5, lwd = 1.5)
curve(q0 *(exp(-EXPDToSND(k)*alpha*q0*x)), col = "blue",add = TRUE, 
      from = .0001, to = 10000, lty = 5, lwd = 1.5)


curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), ylim = c(0, 50), xlim = c(.001, 1000), col = "orange", 
      log = "x", main = expression(bold(atop("Expenditure Curves,","Fixed "*Q[0]*" & "*alpha*", Different k"))), cex.main = 1.5,
      ylab = NA, xlab = NA, axes = FALSE, 
      from = .0001, to = 10000, lty = 1, lwd = 1.5)
title(ylab = "Consumption * Cost", xlab = "Cost per Unit", line = 2.5)
legend("topleft", legend = c("k = 2", "k = 1.5", "k = 1.15", "SND"), ncol = 1, xpd = NA,
       bty = "n", lty = c(1, 5, 3, 5), col = c(rep("orange",3), "blue"),
       xjust = 0)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .8)
curve(x * q0 * 10^(1.15*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 3, lwd = 1.5)
curve(x * q0 * 10^(1.5*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 5, lwd = 1.5)
curve(x * q0 *(exp(-EXPDToSND(k)*alpha*q0*x)), col = "blue",add = TRUE, 
      from = .0001, to = 10000, lty = 5, lwd = 1.5)

curve(x * q0 *(exp(-EXPDToSND(k)*alpha*q0*x)), ylim = c(0, 50), xlim = c(.001, 1000), col = "blue", 
      log = "x", main = expression(bold(atop("Expenditure Curves,","Fixed "*alpha*" & k, Different "*Q[0]))), cex.main = 1.5,
      ylab = NA, xlab = NA, axes = FALSE, 
      from = .0001, to = 10000, lty = 2, lwd = 1.5)
title(ylab = "Consumption * Cost", xlab = "Cost per Unit", line = 2.5)
legend("topleft", 
       legend = c(expression(Q[0]*" = 1"), 
                  expression(Q[0]*" = 20"), 
                  expression(Q[0]*" = 1000"), "HS", "SND"), ncol = 1, xpd = NA,
       bty = "n", lty = c(1, 2, 3, 1, 1), col = c(rep("black",3), "orange", "blue"),
       xjust = 0)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .8)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 2, lwd = 1.5)
q0 <- 1
curve(x * q0 *(exp(-EXPDToSND(k)*alpha*q0*x)), col = "blue",add = TRUE, 
      from = .0001, to = 10000, lty = 1, lwd = 1.5)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 1, lwd = 1.5)
q0 <- 1000
curve(x * q0 *(exp(-EXPDToSND(k)*alpha*q0*x)), col = "blue",add = TRUE, 
      from = .0001, to = 10000, lty = 3, lwd = 1.5)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 3, lwd = 1.5)

q0 <- 20
curve(x * q0 *(exp(-EXPDToSND(k)*alpha*q0*x)), ylim = c(0, 50), xlim = c(.001, 1000), col = "blue", 
      log = "x", main = expression(bold(atop("Expenditure Curves,","Fixed "*Q[0]*" & k, Different "*alpha))), cex.main = 1.5,
      ylab = NA, xlab = NA, axes = FALSE, 
      from = .0001, to = 10000, lty = 2, lwd = 1.5)
title(ylab = "Consumption * Cost", xlab = "Cost per Unit", line = 2.5)
legend("topleft", 
       legend = c(expression(alpha[HS]*" = 0.02"), 
                  expression(alpha[HS]*" = 0.005"), 
                  expression(alpha[HS]*" = 0.0025"), "HS", "SND"), ncol = 1, xpd = NA,
       bty = "n", lty = c(1, 2, 3, 1, 1), col = c(rep("black",3), "orange", "blue"),
       xjust = 0)
lines(x = exp(log(Pmax) +c(-log(100000), log(100000))), y = exp(log(q0 * 10^(k*(exp(-alpha*q0*Pmax)- 1)))+ c(log(100000), -log(100000))), lty = 2)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(2)
axis(1, at = 10^(-3:3), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000)), cex.axis = .8)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 2, lwd = 1.5)
alpha <- .0025
curve(x * q0 *(exp(-EXPDToSND(k)*alpha*q0*x)), col = "blue",add = TRUE, 
      from = .0001, to = 10000, lty = 1, lwd = 1.5)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 1, lwd = 1.5)
alpha <- .02
curve(x * q0 *(exp(-EXPDToSND(k)*alpha*q0*x)), col = "blue",add = TRUE, 
      from = .0001, to = 10000, lty = 3, lwd = 1.5)
curve(x * q0 * 10^(k*(exp(-alpha*q0*x)- 1)), col = "orange",add = TRUE, 
      from = .0001, to = 10000, lty = 3, lwd = 1.5)


1/(EXPDToSND(k)*.0025 * exp(1))
1/(EXPDToSND(k)*.005 * exp(1))
1/(EXPDToSND(k)*.02 * exp(1))

EXPDToSND(k)*.0025
EXPDToSND(k)*.005
EXPDToSND(k)*.02

dev.off()




# Ko et al. (2002) ----
ko2002Data <- read.csv("ko2002Organized.csv")
ko2002Data$combDrug <- paste(ko2002Data$drug, ko2002Data$dose)
ko2002Data$dose <- as.factor(ko2002Data$dose)
ko2002Data$aloneDose <- rep(c("high", "mid", "low"), each = 15)
ko2002Data$xy <- ko2002Data$x * ko2002Data$y
ko2002Data$k <- GetK(ko2002Data)

ko2002Exp <- function(ind){
  est<- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = data.frame(ind),
                      iter = c(10, 10),
                      start_lower = c(q0 = 1, alpha = -5),
                      start_upper = c(q0 = 3, alpha = 0))
  return(est)
}

ko2002FreeK <- function(ind){
  est<- nls_multstart(y ~ 10^(q0) * 10^(ks * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = data.frame(ind),
                      iter = c(7, 7, 7),
                      start_lower = c(q0 = 1, ks = 1, alpha = -5),
                      start_upper = c(q0 = 3, ks = 3, alpha = 0),
                      algorithm = "LM",
                      lower = c(-Inf, 0, -Inf),
                      upper = c(Inf, 10, Inf))
  return(est)
}

ko2002SND <- function(ind){
  est<- nls_multstart(y ~ 10^q0 * exp(-10^alpha *10^q0 * x), data = data.frame(ind),
                      iter = c(10, 10),
                      start_lower = c(q0 = 1, alpha = -5),
                      start_upper = c(q0 = 3, alpha = 0))
  return(est)
}

ko2002Fits <- ko2002Data |> 
  group_by(monkey, combDrug) |> 
  nest() |> 
  mutate(Exp = map(data, ko2002Exp), 
         FreeK = map(data, ko2002FreeK),
         SND = map(data, ko2002SND))

saveRDS(ko2002Fits, "ko2002fits.rds")
ko2002Fits <- readRDS("ko2002fits.rds")


ko2002Outputs <- ko2002Fits |> mutate(ExpCoef = map(Exp, coef),
                                      FreeKCoef = map(FreeK, coef),
                                      SNDCoef = map(SND, coef))

ko2002AlphaEXPD <- ko2002Outputs |> unnest(c(ExpCoef)) |> filter(ExpCoef < 0) |> pull(ExpCoef)
ko2002Q0EXPD <- ko2002Outputs |> unnest(c(ExpCoef)) |> filter(ExpCoef > 0) |> pull(ExpCoef)
ko2002AlphaFreeK <- ko2002Outputs |> unnest(c(FreeKCoef)) |> filter(FreeKCoef < 0) |> pull(FreeKCoef)
ko2002FreeValues <- ko2002Outputs |> unnest(c(FreeKCoef)) |> filter(FreeKCoef > 0) |> pull(FreeKCoef)
ko2002Q0FreeK <- ko2002FreeValues[seq(1,53, by = 2)]
ko2002KFreeK <- ko2002FreeValues[seq(2,54, by = 2)]
ko2002AlphaSND <- ko2002Outputs |> unnest(c(SNDCoef)) |> filter(SNDCoef < 0) |> pull(SNDCoef)
ko2002Q0SND <- ko2002Outputs |> unnest(c(SNDCoef)) |> filter(SNDCoef > 0) |> pull(SNDCoef)


# ko2002Pairs ----
svg("ko2002pairs.svg", 7.5, 7.5)
pairs.panels(data.frame(
  Q0EXPD = ko2002Q0EXPD,
  FreeQ0 = ko2002Q0FreeK,
  SNDQ0 = ko2002Q0SND,
  AlphaEXPD = ko2002AlphaEXPD,
  AlphaFree = ko2002AlphaFreeK,
  SNDAlpha = ko2002AlphaSND), 
  labels = c(expression(log[10](italic(Q)[0*EXPD])), 
             expression(log[10](italic(Q)[0*FreeK])), 
             expression(log[10](italic(Q)[0*SND])), 
             expression(log[10](alpha[EXPD])), 
             expression(log[10](alpha[FreeK])), 
             expression(log[10](alpha[SND]))),
  col = "black", hist.col = "lightblue1",
  ellipses = FALSE, lm = TRUE, main = "Ko et al. (2002) SND and EXPD Correlations", scale =FALSE)
dev.off()

# Ko 2002 GNLS----
ko2002SNDNLS <- nls_multstart(y ~ (10^q0) * exp(-(10^alpha) * (10^q0) * x),
                              data = ko2002Data, iter = c(10, 10),
                              start_lower = c(q0 = 1, alpha = -5),
                              start_upper = c(q0 = 3, alpha = 1))
ko2002SNDGNLS <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                      data = ko2002Data,
                      params = list(q0 ~ drug*aloneDose,  
                                    alpha ~ drug*aloneDose
                      ),
                      start = list(fixed = c(coef(ko2002SNDNLS)[1], rep(0, 8), coef(ko2002SNDNLS)[2], rep(0, 8))),
                      verbose = 2,
                      control = list(msMaxIter = 50000,
                                     niterEM = 5000,
                                     maxIter = 5000,
                                     pnlsTol = .01,
                                     tolerance = .04,
                                     apVar = T,
                                     minScale = .0000001,
                                     opt = "optim"), na.action = na.omit)

ko2002SNDIndGNLS <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                         data = ko2002Data,
                         params = list(q0 ~ combDrug,  
                                       alpha ~ combDrug
                         ),
                         start = list(fixed = c(coef(ko2002SNDNLS)[1], rep(0, 8), coef(ko2002SNDNLS)[2], rep(0, 8))),
                         verbose = 2,
                         control = list(msMaxIter = 50000,
                                        niterEM = 5000,
                                        maxIter = 5000,
                                        pnlsTol = .01,
                                        tolerance = .01,
                                        apVar = T,
                                        minScale = .0000001,
                                        opt = "optim"), na.action = na.omit)

anova(ko2002SNDGNLS, ko2002SNDIndGNLS)
anova(ko2002SNDGNLS)
coef(ko2002SNDGNLS)

emmeans(ko2002SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(ko2002SNDGNLS, param = "q0", specs =  pairwise~aloneDose|drug)

ko2002EXPDNLS <- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                               data = ko2002Data, iter = c(10, 10),
                               start_lower = c(q0 = 1, alpha = -5),
                               start_upper = c(q0 = 3, alpha = 1))
ko2002EXPDGNLS <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                       data = ko2002Data,
                       params = list(q0 ~ drug*aloneDose,  
                                     alpha ~ drug*aloneDose
                       ),
                       start = list(fixed = c(coef(ko2002EXPDNLS)[1], rep(0, 8), coef(ko2002EXPDNLS)[2], rep(0, 8))),
                       verbose = 2,
                       control = list(msMaxIter = 50000,
                                      niterEM = 5000,
                                      maxIter = 5000,
                                      pnlsTol = .01,
                                      tolerance = .01,
                                      apVar = T,
                                      minScale = .0000001,
                                      opt = "optim"), na.action = na.omit)

ko2002EXPDIndGNLS <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                          data = ko2002Data,
                          params = list(q0 ~ combDrug,  
                                        alpha ~ combDrug
                          ),
                          start = list(fixed = c(coef(ko2002EXPDNLS)[1], rep(0, 8), coef(ko2002EXPDNLS)[2], rep(0, 8))),
                          verbose = 2,
                          control = list(msMaxIter = 50000,
                                         niterEM = 5000,
                                         maxIter = 5000,
                                         pnlsTol = .01,
                                         tolerance = .01,
                                         apVar = T,
                                         minScale = .0000001,
                                         opt = "optim"), na.action = na.omit)

summary(ko2002EXPDGNLS)
anova(ko2002EXPDGNLS)
coef(ko2002EXPDGNLS)

emmeans(ko2002EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(ko2002EXPDGNLS, param = "q0", specs = pairwise~aloneDose|drug, data = ko2002Data)

anova(ko2002EXPDGNLS)
anova(ko2002SNDGNLS)
anova(ko2002EXPDGNLS, ko2002SNDGNLS)

1 - sum((predict(ko2002SNDGNLS) - na.omit(ko2002Data)$y)^2)/sum((mean(ko2002Data$y, na.rm = TRUE) - na.omit(ko2002Data)$y)^2)
1 - sum((predict(ko2002EXPDGNLS) - na.omit(ko2002Data)$y)^2)/sum((mean(ko2002Data$y, na.rm = TRUE) - na.omit(ko2002Data)$y)^2)

emmeans(ko2002EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(ko2002SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)


# Ko 2002 Figures ----
ko2002k <- GetK(ko2002Data)
ko2002xVals <- c(10, 32, 100, 320, 1000) 
svg("ko2002figures.svg", width = 9, height = 7.25)

layout(matrix(c(rep(1, 4), 2:13), byrow = TRUE, ncol = 4), heights = c(.06, rep(.3125, 3)))
par(mar = rep(0,4))
plot.new()
text(.5, .5, "Ko et al. (2002)", cex = 3)
legend(-.025, 1.1, c("Low Dose", "Medium Dose", "High Dose", "EXPD", "SND", NA), ncol = 2, xpd = NA,
      bty = "n", pch = c(21, 22, 23, NA, NA, NA), lty = c(NA, NA, NA, 1, 5, NA), col = c(1, 1, 1, "orange", "blue", NA), pt.cex = 1.5,
      xjust = 0)
par(mar = c(5, 4, 2, .5), pty = "m")


curve(10^(coef(ko2002EXPDIndGNLS)[1]) * 10^(ko2002k * 
                                              (exp(-(10^coef(ko2002EXPDIndGNLS)[10]) * 10^(coef(ko2002EXPDIndGNLS)[1]) * x) - 1)), log = "x",
      xlim = c(5, 1200), ylim = c(0,500), from = .01, to = 1200,
      col = "orange", ylab = "Mean Consumption", lwd = 1.5, main = "Alfentanil", cex.main = 2, cex.lab = 1.25, axes = FALSE)
curve(10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[2]) * 10^(ko2002k * 
                                                                           (exp(-(10^(coef(ko2002EXPDIndGNLS)[10] + coef(ko2002EXPDIndGNLS)[11])) * 10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[2]) 
                                                                                * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[3]) * 10^(ko2002k * 
                                                                           (exp(-(10^(coef(ko2002EXPDIndGNLS)[10] + coef(ko2002EXPDIndGNLS)[12])) * 10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[3]) 
                                                                                * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^coef(ko2002SNDIndGNLS)[1] * exp(-(10^coef(ko2002SNDIndGNLS)[10])*10^coef(ko2002SNDIndGNLS)[1]*x), 
      add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[2]) * 
        exp(-(10^(coef(ko2002SNDIndGNLS)[10] + coef(ko2002SNDIndGNLS)[11])*
                10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[2])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[3]) * 
        exp(-(10^(coef(ko2002SNDIndGNLS)[10] + coef(ko2002SNDIndGNLS)[12])*
                10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[3])*x)), add = TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
points(x = rep(ko2002xVals, 3), 
       y =aggregate(y ~ x + dose, subset(ko2002Data, drug == "alfentanil"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 5), cex = 1.5)

arrows(x0 = rep(ko2002xVals), x1 = rep(ko2002xVals), 
       y0= aggregate(y ~ x + dose, subset(ko2002Data, drug == "alfentanil"),mean, na.rm = TRUE)$y - 
         aggregate(y ~ x + dose, subset(ko2002Data, drug == "alfentanil"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(ko2002Data, drug == "alfentanil"),mean, na.rm = TRUE)$y + 
         aggregate(y ~ x + dose, subset(ko2002Data, drug == "alfentanil"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300, 400, 500))


for(a in 1:3){
  curve(10^(ko2002Q0EXPD)[a] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a] * 10^(ko2002Q0EXPD)[a] * x) - 1)), 
        xlim = c(5, 1200), ylim = c(0, 500), log = "x", from = .01, to = 1200, col = "orange", axes = FALSE, cex.main = 1.25,
        xlab = "FR", ylab = "Consumption", main = paste("Monkey", unique(ko2002Data$monkey)[a]), lwd = 1.25, cex.lab = 1.25)
  curve(10^(ko2002Q0EXPD)[a+3] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a+3] * 10^(ko2002Q0EXPD)[a+3] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^(ko2002Q0EXPD)[a+3*2] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a+3*2] * 10^(ko2002Q0EXPD)[a+3*2] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^ko2002Q0SND[a] * exp(- 10^ko2002AlphaSND[a] * 10^ko2002Q0SND[a] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(10^ko2002Q0SND[a+3] * exp(- 10^ko2002AlphaSND[a+3] * 10^ko2002Q0SND[a+3] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(10^ko2002Q0SND[a+3*2] * exp(- 10^ko2002AlphaSND[a+3*2] * 10^ko2002Q0SND[a+3*2] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  points(subset(ko2002Data, monkey == unique(ko2002Data$monkey)[a] & 
                  drug == unique(ko2002Data$drug)[1], select = c(x, y)), 
         pch =rep(c(23, 22, 21), each = 5)) 
  axis(1, at = tickFunc(), tcl = -.25, labels = NA)
  axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
  axis(2, at = c(0, 100, 200, 300, 400, 500))
}


curve(10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[7]) * 10^(ko2002k * 
         (exp(-(10^(coef(ko2002EXPDIndGNLS)[10] + coef(ko2002EXPDIndGNLS)[16])) * 10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[7]) 
              * x) - 1)), log = "x",
      xlim = c(5, 1200), ylim = c(0,500), from = .01, to = 1200,
      col = "orange", ylab = "Mean Consumption", xlab = "FR", lwd = 1.5, main = "Remifentanil", cex.main = 2, cex.lab = 1.25, axes = FALSE)
curve(10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[8]) * 10^(ko2002k * 
                                                                           (exp(-(10^(coef(ko2002EXPDIndGNLS)[10] + coef(ko2002EXPDIndGNLS)[17])) * 10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[8]) 
                                                                                * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[9]) * 10^(ko2002k * 
                                                                           (exp(-(10^(coef(ko2002EXPDIndGNLS)[10] + coef(ko2002EXPDIndGNLS)[18])) * 10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[9]) 
                                                                                * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[7]) * 
        exp(-(10^(coef(ko2002SNDIndGNLS)[10] + coef(ko2002SNDIndGNLS)[16])*
                10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[7])*x)), add  =TRUE,from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[8]) * 
        exp(-(10^(coef(ko2002SNDIndGNLS)[10] + coef(ko2002SNDIndGNLS)[17])*
                10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[8])*x)), add  =TRUE,from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[9]) * 
        exp(-(10^(coef(ko2002SNDIndGNLS)[10] + coef(ko2002SNDIndGNLS)[18])*
                10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[9])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
points(x = rep(ko2002xVals, 3), 
       y =aggregate(y ~ x + dose, subset(ko2002Data, drug == "remifentanil"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 6), cex = 1.5)
arrows(x0 = rep(ko2002xVals), x1 = rep(ko2002xVals), 
       y0= aggregate(y ~ x + dose, subset(ko2002Data, drug == "remifentanil"),mean, na.rm = TRUE)$y - 
         aggregate(y ~ x + dose, subset(ko2002Data, drug == "remifentanil"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(ko2002Data, drug == "remifentanil"),mean, na.rm = TRUE)$y + 
         aggregate(y ~ x + dose, subset(ko2002Data, drug == "remifentanil"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300, 400, 500))



for(a in 10:12){
  curve(10^(ko2002Q0EXPD)[a] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a] * 10^(ko2002Q0EXPD)[a] * x) - 1)), 
        xlim = c(5, 1200), ylim = c(0, 500), log = "x", from = .01, to = 1200, col = "orange", cex.lab = 1.25, cex.main = 1.25,
        axes = FALSE, xlab = "FR", ylab = "Consumption", main = paste("Monkey", unique(ko2002Data$monkey)[a-9] ), lwd = 1.25)
  curve(10^(ko2002Q0EXPD)[a+3] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a+3] * 10^(ko2002Q0EXPD)[a+3] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^(ko2002Q0EXPD)[a+3*2] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a+3*2] * 10^(ko2002Q0EXPD)[a+3*2] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^ko2002Q0SND[a] * exp(- 10^ko2002AlphaSND[a] * 10^ko2002Q0SND[a] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(10^ko2002Q0SND[a+3] * exp(- 10^ko2002AlphaSND[a+3] * 10^ko2002Q0SND[a+3] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(10^ko2002Q0SND[a+3*2] * exp(- 10^ko2002AlphaSND[a+3*2] * 10^ko2002Q0SND[a+3*2] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  points(subset(ko2002Data, monkey == unique(ko2002Data$monkey)[a-9] & 
                  drug == unique(ko2002Data$drug)[2], select = c(x, y)), 
         pch =rep(c(23, 22, 21), each = 5)) 
  axis(1, at = tickFunc(), tcl = -.25, labels = NA)
  axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
  axis(2, at = c(0, 100, 200, 300, 400, 500))
}

curve(10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[4]) * 10^(ko2002k * 
     (exp(-(10^(coef(ko2002EXPDIndGNLS)[10] + coef(ko2002EXPDIndGNLS)[13])) * 10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[4]) 
          * x) - 1)), log = "x",
      xlim = c(5, 1200), ylim = c(0,500), col = "orange", ylab = "Mean Consumption", lwd = 1.5, 
      main = "Fentanyl", xlab = "FR", cex.main = 2, cex.lab = 1.25, axes = FALSE, from = .01, to = 1200)
curve(10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[5]) * 10^(ko2002k * 
                                                                           (exp(-(10^(coef(ko2002EXPDIndGNLS)[10] + coef(ko2002EXPDIndGNLS)[14])) * 10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[5]) 
                                                                                * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[6]) * 10^(ko2002k * 
                                                                           (exp(-(10^(coef(ko2002EXPDIndGNLS)[10] + coef(ko2002EXPDIndGNLS)[15])) * 10^(coef(ko2002EXPDIndGNLS)[1] + coef(ko2002EXPDIndGNLS)[6]) 
                                                                                * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)

curve(10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[4]) * 
        exp(-(10^(coef(ko2002SNDIndGNLS)[10] + coef(ko2002SNDIndGNLS)[13])*
                10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[4])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[5]) * 
        exp(-(10^(coef(ko2002SNDIndGNLS)[10] + coef(ko2002SNDIndGNLS)[14])*
                10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[5])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[6]) * 
        exp(-(10^(coef(ko2002SNDIndGNLS)[10] + coef(ko2002SNDIndGNLS)[15])*
                10^(coef(ko2002SNDIndGNLS)[1] + coef(ko2002SNDIndGNLS)[6])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)

points(x = rep(ko2002xVals, 3), 
       y =aggregate(y ~ x + dose, subset(ko2002Data, drug == "fentanyl"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 5), cex = 1.5)
aggregate(y ~ x + dose, subset(ko2002Data, drug == "fentanyl"),mean, na.rm = TRUE)$y
arrows(x0 = rep(ko2002xVals), x1 = rep(ko2002xVals), 
       y0= aggregate(y ~ x + dose, subset(ko2002Data, drug == "fentanyl"),mean, na.rm = TRUE)$y - 
         aggregate(y ~ x + dose, subset(ko2002Data, drug == "fentanyl"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(ko2002Data, drug == "fentanyl"),mean, na.rm = TRUE)$y + 
         aggregate(y ~ x + dose, subset(ko2002Data, drug == "fentanyl"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300, 400, 500))


for(a in 19:21){
  curve(10^(ko2002Q0EXPD)[a] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a] * 10^(ko2002Q0EXPD)[a] * x) - 1)), 
        xlim = c(5, 1000), ylim = c(0, 500), log = "x", from = .01, to = 1200, col = "orange", axes = FALSE, cex.main = 1.25,
        xlab = "FR", ylab = "Consumption", main = paste("Monkey", unique(ko2002Data$monkey)[a-18]), lwd = 1.25, cex.lab = 1.25)
  curve(10^(ko2002Q0EXPD)[a+3] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a+3] * 10^(ko2002Q0EXPD)[a+3] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^(ko2002Q0EXPD)[a+3*2] * 10^(ko2002k * (exp(-10^(ko2002AlphaEXPD)[a+3*2] * 10^(ko2002Q0EXPD)[a+3*2] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^ko2002Q0SND[a] * exp(- 10^ko2002AlphaSND[a] * 10^ko2002Q0SND[a] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(10^ko2002Q0SND[a+3] * exp(- 10^ko2002AlphaSND[a+3] * 10^ko2002Q0SND[a+3] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(10^ko2002Q0SND[a+3*2] * exp(- 10^ko2002AlphaSND[a+3*2] * 10^ko2002Q0SND[a+3*2] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  points(subset(ko2002Data, monkey == unique(ko2002Data$monkey)[a-18] & 
                  drug == unique(ko2002Data$drug)[3], select = c(x, y)), 
         pch =rep(c(23, 22, 21), each = 5)) 
  axis(1, at = tickFunc(), tcl = -.25, labels = NA)
  axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
  axis(2, at = c(0, 100, 200, 300, 400, 500))
}

dev.off()


# Ko 2002 R2/MAE ----
ko2002Outputs <- ko2002Outputs |> 
  mutate(
    r2Exp = pmap(list(Exp, data), rSquarer),
    MAEExp = pmap(list(Exp, data), MAEmaker),
    AICcExp = pmap(list(Exp, data), AICcmaker),
    r2Free = pmap(list(FreeK, data), rSquarer),
    MAEFree = pmap(list(FreeK, data), MAEmaker),
    AICcFree = pmap(list(FreeK, data), AICcmaker),
    r2SND = pmap(list(SND, data), rSquarer),
    MAESND = pmap(list(SND, data), MAEmaker),
    AICcSND = pmap(list(SND, data), AICcmaker)
  )

rbind(ko2002Outputs |> unnest(c(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND)) |>
        ungroup() |>
        select(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND) |> data.frame() |>
        sapply(range, na.rm = TRUE) |> round(4),
      median = ko2002Outputs |> unnest(c(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND)) |>
        ungroup() |>
        select(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND) |> data.frame() |>
        sapply(median, na.rm = TRUE) |> round(4)
)


# Winger et al. (2002) ----
winger2002Data <- read.csv("winger2002Organized.csv")
winger2002Data$combDrug <- as.factor(paste(winger2002Data$drug, winger2002Data$dose))
winger2002Data$dose <- as.factor(winger2002Data$dose)
winger2002Data$aloneDose <- rep(c("high", "mid", "low"), each = 18)
winger2002Data$xy <- winger2002Data$x * winger2002Data$y
winger2002Data$k <- beezdemand::GetK(winger2002Data)

winger2002Exp <- function(ind){
    est<- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = data.frame(ind),
                        iter = c(10, 10),
                        start_lower = c(q0 = 1, alpha = -5),
                        start_upper = c(q0 = 3, alpha = 0))
    return(est)
}

winger2002FreeK <- function(ind){
    est<- nls_multstart(y ~ 10^(q0) * 10^(ks * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = data.frame(ind),
                        iter = c(7, 7, 7),
                        start_lower = c(q0 = 1, ks = 1, alpha = -5),
                        start_upper = c(q0 = 3, ks = 3, alpha = 0),
                        algorithm = "LM",
                        lower = c(-Inf, 0, -Inf),
                        upper = c(Inf, 10, Inf))
    return(est)
}

winger2002SND <- function(ind){
    est<- nls_multstart(y ~ 10^q0 * exp(-10^alpha *10^q0 * x), data = data.frame(ind),
                        iter = c(10, 10),
                        start_lower = c(q0 = 1, alpha = -5),
                        start_upper = c(q0 = 3, alpha = 0))
    return(est)
}

winger2002Fits <- winger2002Data |> 
    group_by(monkey, combDrug) |> 
    nest() |> 
    mutate(Exp = map(data, winger2002Exp), 
           FreeK = map(data, winger2002FreeK),
           SND = map(data, winger2002SND))

saveRDS(winger2002Fits, "winger2002fits.rds")
winger2002Fits <- readRDS("winger2002fits.rds")

winger2002Outputs <- winger2002Fits |> mutate(ExpCoef = map(Exp, coef),
                                              FreeKCoef = map(FreeK, coef),
                                              SNDCoef = map(SND, coef))

winger2002Outputs |> unnest(c(FreeKCoef)) |> pull(FreeKCoef)

winger2002AlphaEXPD <- winger2002Outputs |> unnest(c(ExpCoef)) |> filter(ExpCoef < 0) |> pull(ExpCoef)
winger2002Q0EXPD <- winger2002Outputs |> unnest(c(ExpCoef)) |> filter(ExpCoef > 0) |> pull(ExpCoef)
winger2002AlphaFreeK <- winger2002Outputs |> unnest(c(FreeKCoef)) |> filter(FreeKCoef < 0) |> pull(FreeKCoef)
Winger2002FreeValues <- winger2002Outputs |> unnest(c(FreeKCoef)) |> filter(FreeKCoef > 0) |> pull(FreeKCoef)
winger2002Q0FreeK <- Winger2002FreeValues[seq(1,53, by = 2)]
winger2002KFreeK <- Winger2002FreeValues[seq(2,54, by = 2)]
winger2002AlphaSND <- winger2002Outputs |> unnest(c(SNDCoef)) |> filter(SNDCoef < 0) |> pull(SNDCoef)
winger2002Q0SND <- winger2002Outputs |> unnest(c(SNDCoef)) |> filter(SNDCoef > 0) |> pull(SNDCoef)



# Winger 2002 Pairs ----
svg("winger2002pairs.svg", 7.5, 7.5)
pairs.panels(data.frame(
    Q0EXPD = winger2002Q0EXPD,
    Q0Free = winger2002Q0FreeK,
    SNDQ0 = winger2002Q0SND,
    AlphaEXPD = winger2002AlphaEXPD,
    AlphaFree = winger2002AlphaFreeK,
    SNDAlpha = winger2002AlphaSND), 
    labels = c(expression(log[10](italic(Q)[0*EXPD])), 
               expression(log[10](italic(Q)[0*FreeK])), 
               expression(log[10](italic(Q)[0*SND])), 
                 expression(log[10](alpha[EXPD])), 
                 expression(log[10](alpha[FreeK])), 
               expression(log[10](alpha[SND]))),
    col = "black", hist.col = "lightblue1",
    ellipses = FALSE, lm = TRUE, main = "Winger et al. (2002) SND and EXPD Correlations", scale =FALSE)
dev.off()



# Winger 2002 GNLS----

winger2002SNDNLS <- nls_multstart(y ~ (10^q0) * exp(-(10^alpha) * (10^q0) * x),
                                  data = winger2002Data, iter = c(10, 10),
                                  start_lower = c(q0 = 1, alpha = -5),
                                  start_upper = c(q0 = 3, alpha = 1))
winger2002SNDGNLS <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                          data = winger2002Data,
                          params = list(q0 ~ drug*aloneDose,  
                                        alpha ~ drug*aloneDose
                          ),
                          start = list(fixed = c(coef(winger2002SNDNLS)[1], rep(0, 8), coef(winger2002SNDNLS)[2], rep(0, 8))),
                          verbose = 2,
                          control = list(msMaxIter = 50000,
                                         niterEM = 5000,
                                         maxIter = 5000,
                                         pnlsTol = .01,
                                         tolerance = .01,
                                         apVar = T,
                                         minScale = .0000001,
                                         opt = "optim"), na.action = na.omit)

winger2002SNDIndGNLS <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                             data = winger2002Data,
                             params = list(q0 ~ combDrug,  
                                           alpha ~ combDrug
                             ),
                             start = list(fixed = c(coef(winger2002SNDNLS)[1], rep(0, 8), coef(winger2002SNDNLS)[2], rep(0, 8))),
                             verbose = 2,
                             control = list(msMaxIter = 50000,
                                            niterEM = 5000,
                                            maxIter = 5000,
                                            pnlsTol = .01,
                                            tolerance = .01,
                                            apVar = T,
                                            minScale = .0000001,
                                            opt = "optim"), na.action = na.omit)

anova(winger2002SNDGNLS, winger2002SNDIndGNLS)
anova(winger2002SNDGNLS)
coef(winger2002SNDGNLS)

emmeans(winger2002SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(winger2002SNDGNLS, param = "q0", specs =  pairwise~aloneDose|drug)

winger2002EXPDNLS <- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                                   data = winger2002Data, iter = c(10, 10),
                                   start_lower = c(q0 = 1, alpha = -5),
                                   start_upper = c(q0 = 3, alpha = 1))
winger2002EXPDGNLS <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                           data = winger2002Data,
                           params = list(q0 ~ drug*aloneDose,  
                                         alpha ~ drug*aloneDose
                           ),
                           start = list(fixed = c(coef(winger2002EXPDNLS)[1], rep(0, 8), coef(winger2002EXPDNLS)[2], rep(0, 8))),
                           verbose = 2,
                           control = list(msMaxIter = 50000,
                                          niterEM = 5000,
                                          maxIter = 5000,
                                          pnlsTol = .01,
                                          tolerance = .01,
                                          apVar = T,
                                          minScale = .0000001,
                                          opt = "optim"), na.action = na.omit)

winger2002EXPDIndGNLS <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                              data = winger2002Data,
                              params = list(q0 ~ combDrug,  
                                            alpha ~ combDrug
                              ),
                              start = list(fixed = c(coef(winger2002EXPDNLS)[1], rep(0, 8), coef(winger2002EXPDNLS)[2], rep(0, 8))),
                              verbose = 2,
                              control = list(msMaxIter = 50000,
                                             niterEM = 5000,
                                             maxIter = 5000,
                                             pnlsTol = .01,
                                             tolerance = .01,
                                             apVar = T,
                                             minScale = .0000001,
                                             opt = "optim"), na.action = na.omit)

summary(winger2002EXPDGNLS)
anova(winger2002EXPDGNLS)
coef(winger2002EXPDGNLS)

emmeans(winger2002EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(winger2002EXPDGNLS, param = "q0", specs = pairwise~aloneDose|drug, data = winger2002Data)

anova(winger2002EXPDGNLS, winger2002SNDGNLS)

1 - sum((predict(winger2002SNDGNLS) - na.omit(winger2002Data)$y)^2)/sum((mean(winger2002Data$y, na.rm = TRUE) - na.omit(winger2002Data)$y)^2)
1 - sum((predict(winger2002EXPDGNLS) - na.omit(winger2002Data)$y)^2)/sum((mean(winger2002Data$y, na.rm = TRUE) - na.omit(winger2002Data)$y)^2)

anova(winger2002EXPDGNLS, type = "marginal")
anova(winger2002SNDGNLS, type = "marginal")
emmeans(winger2002EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(winger2002SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)


# Winger 2002 Figures ----
svg("winger2002figures.svg", width = 9, height = 7.25)

winger2002k <- GetK(winger2002Data)
winger2002xVals <- c(1, 10, 32, 100, 320, 1000) 
layout(matrix(c(rep(1, 4), 2:13), byrow = TRUE, ncol = 4), heights = c(.06, rep(.3125, 3)))
par(mar = rep(0,4))
plot.new()
text(.5, .5, "Winger et al. (2002)", cex = 3)
legend(-.025, 1.1, c("Low Dose", "Medium Dose", "High Dose", "EXPD", "SND", NA), ncol = 2, xpd = NA,
       bty = "n", pch = c(21, 22, 23, NA, NA, NA), lty = c(NA, NA, NA, 1, 5, NA), col = c(1, 1, 1, "orange", "blue", NA), pt.cex = 1.5,
       xjust = 0)

par(mar = c(5, 4, 2, .5), pty = "m")

curve(10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[4]) * 10^(winger2002k * 
   (exp(-(10^(coef(winger2002EXPDIndGNLS)[10] + coef(winger2002EXPDIndGNLS)[13])) * 10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[4]) 
        * x) - 1)),  log = "x",
   xlim = c(.5, 1200), ylim = c(0,500), ylab = "Mean Consumption", lwd = 1.5, 
   main = "Ketamine", xlab = "FR", cex.main = 2, axes = FALSE, cex.lab = 1.25,
      col = "orange", from = .01, to = 1200)
curve(10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[5]) * 10^(winger2002k * 
                       (exp(-(10^(coef(winger2002EXPDIndGNLS)[10] + coef(winger2002EXPDIndGNLS)[14])) * 10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[5]) 
                            * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[6]) * 10^(winger2002k * 
                       (exp(-(10^(coef(winger2002EXPDIndGNLS)[10] + coef(winger2002EXPDIndGNLS)[15])) * 10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[6]) 
                        * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[4]) * 
        exp(-(10^(coef(winger2002SNDIndGNLS)[10] + coef(winger2002SNDIndGNLS)[13])*
                10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[4])*x)), 
      add = TRUE, from = .01, to = 1200, lty = 5, lwd = 1.5, col = "blue")
curve(10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[5]) * 
        exp(-(10^(coef(winger2002SNDIndGNLS)[10] + coef(winger2002SNDIndGNLS)[14])*
                10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[5])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[6]) * 
        exp(-(10^(coef(winger2002SNDIndGNLS)[10] + coef(winger2002SNDIndGNLS)[15])*
                10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[6])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)

points(x = rep(winger2002xVals, 3), 
       y =aggregate(y ~ x + dose, subset(winger2002Data, drug == "ketamine"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 6), cex = 1.5)
aggregate(y ~ x + dose, subset(winger2002Data, drug == "ketamine"),mean, na.rm = TRUE)$y
arrows(x0 = rep(winger2002xVals), x1 = rep(winger2002xVals), 
       y0= aggregate(y ~ x + dose, subset(winger2002Data, drug == "ketamine"),mean, na.rm = TRUE)$y - 
         aggregate(y ~ x + dose, subset(winger2002Data, drug == "ketamine"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(winger2002Data, drug == "ketamine"),mean, na.rm = TRUE)$y + 
         aggregate(y ~ x + dose, subset(winger2002Data, drug == "ketamine"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300, 400, 500))

for(a in 1:3){
curve(10^(winger2002Q0EXPD)[a] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a] * 10^(winger2002Q0EXPD)[a] * x) - 1)), 
      xlim = c(.5, 1200), ylim = c(0, 500), log = "x", from = .01, to = 1200, col = "orange", axes = FALSE, cex.lab = 1.25, cex.main = 1.25,
      xlab = "FR", ylab = "Consumption", main = paste("Monkey", unique(winger2002Data$monkey)[a]))
curve(10^(winger2002Q0EXPD)[a+3] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a+3] * 10^(winger2002Q0EXPD)[a+3] * x) - 1)), 
      from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
curve(10^(winger2002Q0EXPD)[a+3*2] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a+3*2] * 10^(winger2002Q0EXPD)[a+3*2] * x) - 1)), 
      from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
curve(10^winger2002Q0SND[a] * exp(- 10^winger2002AlphaSND[a] * 10^winger2002Q0SND[a] * x), 
      from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
curve(10^winger2002Q0SND[a+3] * exp(- 10^winger2002AlphaSND[a+3] * 10^winger2002Q0SND[a+3] * x), 
      from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
curve(10^winger2002Q0SND[a+3*2] * exp(- 10^winger2002AlphaSND[a+3*2] * 10^winger2002Q0SND[a+3*2] * x), 
      from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
points(subset(winger2002Data, monkey == unique(winger2002Data$monkey)[a] & 
                drug == unique(winger2002Data$drug)[1], select = c(x, y)), 
       pch =rep(c(23, 22, 21), each = 6)) 
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300, 400, 500))
}


curve(10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[7]) * 10^(winger2002k * 
                             (exp(-(10^(coef(winger2002EXPDIndGNLS)[10] + coef(winger2002EXPDIndGNLS)[16])) * 10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[7]) 
                      * x) - 1)),
      col = "orange", log = "x",
      xlim = c(.5, 1200), ylim = c(0,500), from = .01, to = 1200, cex.lab = 1.25,
      ylab = "Mean Consumption", xlab = "FR", lwd = 1.5, main = "Phencyclidine", cex.main = 2, axes = FALSE)
curve(10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[8]) * 10^(winger2002k * 
                 (exp(-(10^(coef(winger2002EXPDIndGNLS)[10] + coef(winger2002EXPDIndGNLS)[17])) * 10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[8]) 
                      * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[9]) * 10^(winger2002k * 
         (exp(-(10^(coef(winger2002EXPDIndGNLS)[10] + coef(winger2002EXPDIndGNLS)[18])) * 10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[9]) 
                                  * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[7]) * 
        exp(-(10^(coef(winger2002SNDIndGNLS)[10] + coef(winger2002SNDIndGNLS)[16])*
                10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[7])*x)), add  =TRUE,from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[8]) * 
        exp(-(10^(coef(winger2002SNDIndGNLS)[10] + coef(winger2002SNDIndGNLS)[17])*
                10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[8])*x)), add  =TRUE,from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[9]) * 
        exp(-(10^(coef(winger2002SNDIndGNLS)[10] + coef(winger2002SNDIndGNLS)[18])*
                10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[9])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)

points(x = rep(winger2002xVals, 3), 
       y =aggregate(y ~ x + dose, subset(winger2002Data, drug == "phencyclidine"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 6), cex = 1.5)
arrows(x0 = rep(winger2002xVals), x1 = rep(winger2002xVals), 
       y0= aggregate(y ~ x + dose, subset(winger2002Data, drug == "phencyclidine"),mean, na.rm = TRUE)$y - 
         aggregate(y ~ x + dose, subset(winger2002Data, drug == "phencyclidine"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(winger2002Data, drug == "phencyclidine"),mean, na.rm = TRUE)$y + 
         aggregate(y ~ x + dose, subset(winger2002Data, drug == "phencyclidine"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300, 400, 500))

for(a in 10:12){
  curve(10^(winger2002Q0EXPD)[a] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a] * 10^(winger2002Q0EXPD)[a] * x) - 1)), 
        xlim = c(.5, 1200), ylim = c(0, 500), log = "x", from = .01, to = 1200, col = "orange", cex.lab = 1.25,
        axes = FALSE, xlab = "FR", lwd = 1.25, cex.main = 1.25,
        ylab = "Consumption", main = paste("Monkey", unique(winger2002Data$monkey)[a-9] ))
  curve(10^(winger2002Q0EXPD)[a+3] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a+3] * 10^(winger2002Q0EXPD)[a+3] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^(winger2002Q0EXPD)[a+3*2] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a+3*2] * 10^(winger2002Q0EXPD)[a+3*2] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^winger2002Q0SND[a] * exp(- 10^winger2002AlphaSND[a] * 10^winger2002Q0SND[a] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
  curve(10^winger2002Q0SND[a+3] * exp(- 10^winger2002AlphaSND[a+3] * 10^winger2002Q0SND[a+3] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
  curve(10^winger2002Q0SND[a+3*2] * exp(- 10^winger2002AlphaSND[a+3*2] * 10^winger2002Q0SND[a+3*2] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
  points(subset(winger2002Data, monkey == unique(winger2002Data$monkey)[a-9] & 
                  drug == unique(winger2002Data$drug)[2], select = c(x, y)), 
         pch =rep(c(23, 22, 21), each = 6)) 
  axis(1, at = tickFunc(), tcl = -.25, labels = NA)
  axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
  axis(2, at = c(0, 100, 200, 300, 400, 500))
}


curve(10^(coef(winger2002EXPDIndGNLS)[1]) * 10^(winger2002k * 
          (exp(-(10^coef(winger2002EXPDIndGNLS)[10]) * 10^(coef(winger2002EXPDIndGNLS)[1]) * x) - 1)), log = "x",
      xlim = c(.5, 1200), ylim = c(0,500), from = .01, to = 1200, cex.lab = 1.25,
      col = "orange", ylab = "Mean Consumption", xlab = "FR", lwd = 1.5, main = "Dizocilpine", cex.main = 2, axes = FALSE)
curve(10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[2]) * 10^(winger2002k * 
                   (exp(-(10^(coef(winger2002EXPDIndGNLS)[10] + coef(winger2002EXPDIndGNLS)[11])) * 10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[2]) 
                        * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[3]) * 10^(winger2002k * 
                       (exp(-(10^(coef(winger2002EXPDIndGNLS)[10] + coef(winger2002EXPDIndGNLS)[12])) * 10^(coef(winger2002EXPDIndGNLS)[1] + coef(winger2002EXPDIndGNLS)[3]) 
                                * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)
curve(10^coef(winger2002SNDIndGNLS)[1] * exp(-(10^coef(winger2002SNDIndGNLS)[10])*10^coef(winger2002SNDIndGNLS)[1]*x), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[2]) * 
        exp(-(10^(coef(winger2002SNDIndGNLS)[10] + coef(winger2002SNDIndGNLS)[11])*
                10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[2])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[3]) * 
        exp(-(10^(coef(winger2002SNDIndGNLS)[10] + coef(winger2002SNDIndGNLS)[12])*
                10^(coef(winger2002SNDIndGNLS)[1] + coef(winger2002SNDIndGNLS)[3])*x)), add = TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)
points(x = rep(winger2002xVals, 3), 
       y =aggregate(y ~ x + dose, subset(winger2002Data, drug == "dizocilpine"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 6), cex = 1.5)

arrows(x0 = rep(winger2002xVals), x1 = rep(winger2002xVals), 
       y0= aggregate(y ~ x + dose, subset(winger2002Data, drug == "dizocilpine"),mean, na.rm = TRUE)$y - 
         aggregate(y ~ x + dose, subset(winger2002Data, drug == "dizocilpine"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(winger2002Data, drug == "dizocilpine"),mean, na.rm = TRUE)$y + 
         aggregate(y ~ x + dose, subset(winger2002Data, drug == "dizocilpine"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300, 400, 500))

for(a in 19:21){
  curve(10^(winger2002Q0EXPD)[a] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a] * 10^(winger2002Q0EXPD)[a] * x) - 1)), 
        xlim = c(.5, 1000), ylim = c(0, 500), lwd = 1.25, cex.lab = 1.25, cex.main = 1.25,
        log = "x", from = .01, to = 1200, col = "orange", axes = FALSE,
        xlab = "FR", ylab = "Consumption", main = paste("Monkey", unique(winger2002Data$monkey)[a-18]))
  curve(10^(winger2002Q0EXPD)[a+3] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a+3] * 10^(winger2002Q0EXPD)[a+3] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^(winger2002Q0EXPD)[a+3*2] * 10^(winger2002k * (exp(-10^(winger2002AlphaEXPD)[a+3*2] * 10^(winger2002Q0EXPD)[a+3*2] * x) - 1)), 
        from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
  curve(10^winger2002Q0SND[a] * exp(- 10^winger2002AlphaSND[a] * 10^winger2002Q0SND[a] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
  curve(10^winger2002Q0SND[a+3] * exp(- 10^winger2002AlphaSND[a+3] * 10^winger2002Q0SND[a+3] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
  curve(10^winger2002Q0SND[a+3*2] * exp(- 10^winger2002AlphaSND[a+3*2] * 10^winger2002Q0SND[a+3*2] * x), 
        from = .01, to = 1200, col = "blue", add = TRUE, lwd = 1.25, lty = 5)
  points(subset(winger2002Data, monkey == unique(winger2002Data$monkey)[a-18] & 
                  drug == unique(winger2002Data$drug)[3], select = c(x, y)), 
         pch =rep(c(23, 22, 21), each = 6)) 
  axis(1, at = tickFunc(), tcl = -.25, labels = NA)
  axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
  axis(2, at = c(0, 100, 200, 300, 400, 500))
}

dev.off()

# Winger 2002 R2/MAE ----

winger2002Outputs <- winger2002Outputs |> 
    mutate(
        r2Exp = pmap(list(Exp, data), rSquarer),
        MAEExp = pmap(list(Exp, data), MAEmaker),
        AICcExp = pmap(list(Exp, data), AICcmaker),
        r2Free = pmap(list(FreeK, data), rSquarer),
        MAEFree = pmap(list(FreeK, data), MAEmaker),
        AICcFree = pmap(list(FreeK, data), AICcmaker),
        r2SND = pmap(list(SND, data), rSquarer),
        MAESND = pmap(list(SND, data), MAEmaker),
        AICcSND = pmap(list(SND, data), AICcmaker)
        )

rbind(winger2002Outputs |> unnest(c(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND)) |>
    ungroup() |>
    select(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND) |> data.frame() |>
    sapply(range, na.rm = TRUE) |> round(4),
    median = winger2002Outputs |> unnest(c(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND)) |>
    ungroup() |>
    select(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND) |> data.frame() |>
    sapply(median, na.rm = TRUE) |> round(4)
)



# Winger et al. 2006 ----
winger2006Data <- read.csv("winger2006Organized.csv")
winger2006Data$combDrug <- as.factor(paste(winger2006Data$drug, winger2006Data$dose))
winger2006Data$drug <- as.factor(winger2006Data$drug)
winger2006Data$dose <- as.factor(winger2006Data$dose)
winger2006Data$aloneDose <- as.factor(rep(c("low", "mid", "high"), each = 18))
winger2006Data$xy <- winger2006Data$x * winger2006Data$y
winger2006Data$monkey <- as.factor(winger2006Data$monkey)
winger2006Data$k <- GetK(winger2006Data)

winger2006Exp <- function(ind){
    est<- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = data.frame(ind),
                        iter = c(10, 10),
                        start_lower = c(q0 = 1, alpha = -5),
                        start_upper = c(q0 = 3, alpha = 0))
    return(est)
}

winger2006FreeK <- function(ind){
    est<- nls_multstart(y ~ 10^(q0) * 10^(ks * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = data.frame(ind),
                        iter = c(7, 7, 7),
                        start_lower = c(q0 = 1, ks = 1, alpha = -5),
                        start_upper = c(q0 = 3, ks = 3, alpha = 0),
                        algorithm = "LM",
                        lower = c(-Inf, 0, -Inf),
                        upper = c(Inf, 10, Inf))
    return(est)
}

winger2006SND <- function(ind){
    est<- nls_multstart(y ~ 10^q0 * exp(-10^alpha *10^q0 * x), data = data.frame(ind),
                        iter = c(10, 10),
                        start_lower = c(q0 = 1, alpha = -5),
                        start_upper = c(q0 = 3, alpha = 0))
    return(est)
}


winger2006Fits <- winger2006Data |> 
    group_by(monkey, combDrug) |> 
    nest() |> 
    mutate(Exp = map(data, winger2006Exp), 
           FreeK = map(data, winger2006FreeK),
           SND = map(data, winger2006SND))

saveRDS(winger2006Fits, "winger2006fits.rds")
winger2006Fits <- readRDS("winger2006fits.rds")

winger2006Outputs <- winger2006Fits |> mutate(ExpCoef = map(Exp, coef),
                                              FreeKCoef = map(FreeK, coef),
                                              SNDCoef = map(SND, coef))


winger2006AlphaEXPD <- winger2006Outputs |> unnest(c(ExpCoef)) |> filter(ExpCoef < 0) |> pull(ExpCoef)
winger2006Q0EXPD <- winger2006Outputs |> unnest(c(ExpCoef)) |> filter(ExpCoef > 0) |> pull(ExpCoef)
winger2006AlphaFreeK <- winger2006Outputs |> unnest(c(FreeKCoef)) |> filter(FreeKCoef < 0) |> pull(FreeKCoef)
Winger2006FreeValues <- winger2006Outputs |> unnest(c(FreeKCoef)) |> filter(FreeKCoef > 0) |> pull(FreeKCoef)
winger2006Q0FreeK <- Winger2006FreeValues[seq(1,53, by = 2)]
winger2006KFreeK <- Winger2006FreeValues[seq(2,54, by = 2)]
winger2006AlphaSND <- winger2006Outputs |> unnest(c(SNDCoef)) |> filter(SNDCoef < 0) |> pull(SNDCoef)
winger2006Q0SND <- winger2006Outputs |> unnest(c(SNDCoef)) |> filter(SNDCoef > 0) |> pull(SNDCoef)

# Winger 2006 Pairs ----
svg("winger2006pairs.svg", 7.5, 7.5)
pairs.panels(data.frame(
    Q0EXPD = winger2006Q0EXPD,
    Q0EXPD = winger2006Q0FreeK,
    SNDQ0 = winger2006Q0SND,
    AlphaEXPD = winger2006AlphaEXPD,
    AlphaEXPD = winger2006AlphaFreeK,
    SNDAlpha = winger2006AlphaSND), labels = c(expression(log[10](italic(Q)[0*EXPD])),
                                               expression(log[10](italic(Q)[0*FreeK])), 
                                               expression(log[10](italic(Q)[0*SND])), 
                                           expression(log[10](alpha[EXPD])), 
                                           expression(log[10](alpha[FreeK])), 
                                           expression(log[10](alpha[SND]))),
    col = "black", hist.col = "lightblue1",
    ellipses = FALSE, lm = TRUE, main = "Winger et al. (2006) SND and EXPD Correlations", scale =FALSE)
dev.off()


# Winger 2006 GNLS----
winger2006SNDNLS <- nls_multstart(y ~ (10^q0) * exp(-(10^alpha) * (10^q0) * x),
                                  data = winger2006Data, iter = c(10, 10),
                                  start_lower = c(q0 = 1, alpha = -5),
                                  start_upper = c(q0 = 3, alpha = 1))
winger2006SNDGNLS <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                          data = winger2006Data,
                          params = list(q0 ~ drug*aloneDose,  
                                        alpha ~ drug*aloneDose
                          ),
                          start = list(fixed = c(coef(winger2006SNDNLS)[1], rep(0, 8), coef(winger2006SNDNLS)[2], rep(0, 8))),
                          verbose = 2,
                          control = list(msMaxIter = 50000,
                                         niterEM = 5000,
                                         maxIter = 5000,
                                         pnlsTol = .01,
                                         tolerance = .21,
                                         apVar = T,
                                         minScale = .0000001,
                                         opt = "optim"), na.action = na.omit)

winger2006SNDIndGNLS <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                             data = winger2006Data,
                             params = list(q0 ~ combDrug,  
                                           alpha ~ combDrug
                             ),
                             start = list(fixed = c(coef(winger2006SNDNLS)[1], rep(0, 8), coef(winger2006SNDNLS)[2], rep(0, 8))),
                             verbose = 2,
                             control = list(msMaxIter = 50000,
                                            niterEM = 5000,
                                            maxIter = 5000,
                                            pnlsTol = .01,
                                            tolerance = .01,
                                            apVar = T,
                                            minScale = .0000001,
                                            opt = "optim"), na.action = na.omit)



anova(winger2006SNDGNLS, winger2006SNDIndGNLS)
anova(winger2006SNDGNLS)
coef(winger2006SNDGNLS)

emmeans(winger2006SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(winger2006SNDGNLS, param = "q0", specs =  pairwise~aloneDose|drug)

winger2006EXPDNLS <- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                                   data = winger2006Data, iter = c(10, 10),
                                   start_lower = c(q0 = 1, alpha = -5),
                                   start_upper = c(q0 = 3, alpha = 1))

winger2006EXPDGNLS <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                           data = winger2006Data,
                           params = list(q0 ~ drug*aloneDose,  
                                         alpha ~ drug*aloneDose
                           ),
                           start = list(fixed = c(coef(winger2006EXPDNLS)[1], rep(0, 8), coef(winger2006EXPDNLS)[2], rep(0, 8))),
                           verbose = 2,
                           control = list(msMaxIter = 50000,
                                          niterEM = 5000,
                                          maxIter = 5000,
                                          pnlsTol = .01,
                                          tolerance = .03,
                                          apVar = T,
                                          minScale = .0000001,
                                          opt = "optim"), na.action = na.omit)

winger2006EXPDIndGNLS <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                              data = winger2006Data,
                              params = list(q0 ~ combDrug,  
                                            alpha ~ combDrug
                              ),
                              start = list(fixed = c(coef(winger2006EXPDNLS)[1], rep(0, 8), coef(winger2006EXPDNLS)[2], rep(0, 8))),
                              verbose = 2,
                              control = list(msMaxIter = 50000,
                                             niterEM = 5000,
                                             maxIter = 5000,
                                             pnlsTol = .01,
                                             tolerance = .01,
                                             apVar = T,
                                             minScale = .0000001,
                                             opt = "optim"), na.action = na.omit)

summary(winger2006EXPDGNLS)
anova(winger2006EXPDGNLS)
coef(winger2006EXPDGNLS)

emmeans(winger2006EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(winger2006EXPDGNLS, param = "q0", specs = pairwise~aloneDose|drug, data = winger2006Data)

anova(winger2006EXPDGNLS, winger2006SNDGNLS)

1 - sum((predict(winger2006SNDGNLS) - na.omit(winger2006Data)$y)^2)/sum((mean(winger2006Data$y, na.rm = TRUE) - na.omit(winger2006Data)$y)^2)
1 - sum((predict(winger2006EXPDGNLS) - na.omit(winger2006Data)$y)^2)/sum((mean(winger2006Data$y, na.rm = TRUE) - na.omit(winger2006Data)$y)^2)


anova(winger2006EXPDGNLS, type = "marginal")
anova(winger2006SNDGNLS, type = "marginal")
emmeans(winger2006EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)
emmeans(winger2006SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)


# Winger 2006 Figures ----
svg("winger2006figures.svg", width = 9, height = 7.25)

layout(matrix(c(rep(1, 4), 2:13), byrow = TRUE, ncol = 4), heights = c(.06, rep(.3125, 3)))
par(mar = rep(0,4))
plot.new()
text(.5, .5, "Winger et al. (2006)", cex = 3)
legend(-.025, 1.1, c("Low Dose", "Medium Dose", "High Dose", "EXPD", "SND", NA), ncol = 2, xpd = NA,
       bty = "n", pch = c(21, 22, 23, NA, NA, NA), lty = c(NA, NA, NA, 1, 5, NA), col = c(1, 1, 1, "orange", "blue", NA), pt.cex = 1.5,
       xjust = 0)
par(mar = c(5, 4, 2, .5), pty = "m")
winger2006k <- GetK(winger2006Data)
winger2006xVals <- c(10,32,100,320,1000, 1780)
curve(10^(coef(winger2006EXPDIndGNLS)[1]) * 10^(winger2006k * 
         (exp(-(10^(coef(winger2006EXPDIndGNLS)[10])) * 10^(coef(winger2006EXPDIndGNLS)[1]) 
                                          * x) - 1)), log = "x", cex.lab = 1.25,
      xlim = c(5, 2000), ylim = c(0, 300), col = "orange", ylab = "Mean Consumption", lwd = 1.5, 
      main = "Cocaine", xlab = "FR", cex.main = 2, axes = FALSE, from = .01, to = 2000)
curve(10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[2]) * 10^(winger2006k * 
         (exp(-(10^(coef(winger2006EXPDIndGNLS)[10] + coef(winger2006EXPDIndGNLS)[11])) * 
                  10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[2]) 
              * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 2000)
curve(10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[3]) * 10^(winger2006k * 
         (exp(-(10^(coef(winger2006EXPDIndGNLS)[10] + coef(winger2006EXPDIndGNLS)[12])) *
                  10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[3]) 
              * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 2000)
curve(10^(coef(winger2006SNDIndGNLS)[1]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10])*
                    10^(coef(winger2006SNDIndGNLS)[1])*x)), add  =TRUE, from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[2]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10] + coef(winger2006SNDIndGNLS)[11])*
                    10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[2])*x)), add  =TRUE, from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[3]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10] + coef(winger2006SNDIndGNLS)[12])*
                    10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[3])*x)), add  =TRUE, from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)
points(x = rep(winger2006xVals, 3), 
       y =aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaine"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 6), cex = 1.5)
aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaine"),mean, na.rm = TRUE)$y
arrows(x0 = rep(winger2006xVals), x1 = rep(winger2006xVals), 
       y0= aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaine"),mean, na.rm = TRUE)$y - 
           aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaine"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaine"),mean, na.rm = TRUE)$y + 
           aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaine"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(jumps = 7), tcl = -.25, labels = NA)
axis(1, at = 10^(-1:4), labels = as.character(c(.1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300, 400, 500))

for(a in 1:3){
    curve(10^(winger2006Q0EXPD)[a] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a] * 10^(winger2006Q0EXPD)[a] * x) - 1)), 
          xlim = c(5, 2000), ylim = c(0, 300), log = "x", from = .01, to = 2000, col = "orange", axes = FALSE, cex.main = 1.25, cex.lab = 1.25,
          xlab = "FR", ylab = "Consumption", main = paste("Monkey", unique(winger2006Data$monkey)[a]), lwd = 1.25)
    curve(10^(winger2006Q0EXPD)[a+3] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a+3] * 10^(winger2006Q0EXPD)[a+3] * x) - 1)), 
          from = .01, to = 2000, col = "orange", add = TRUE, lwd = 1.25)
    curve(10^(winger2006Q0EXPD)[a+3*2] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a+3*2] * 10^(winger2006Q0EXPD)[a+3*2] * x) - 1)), 
          from = .01, to = 2000, col = "orange", add = TRUE, lwd = 1.25)
    curve(10^winger2006Q0SND[a] * exp(- 10^winger2006AlphaSND[a] * 10^winger2006Q0SND[a] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    curve(10^winger2006Q0SND[a+3] * exp(- 10^winger2006AlphaSND[a+3] * 10^winger2006Q0SND[a+3] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    curve(10^winger2006Q0SND[a+3*2] * exp(- 10^winger2006AlphaSND[a+3*2] * 10^winger2006Q0SND[a+3*2] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(subset(winger2006Data, monkey == unique(winger2006Data$monkey)[a] & 
                      drug == unique(winger2006Data$drug)[1], select = c(x, y)), 
           pch =rep(c(21, 22, 23), each = 6)) 
    axis(1, at = tickFunc(jumps = 7), tcl = -.25, labels = NA)
    axis(1, at = 10^(-1:4), labels = as.character(c(.1, 1, 10, 100, 1000, 2000)))
    axis(2, at = c(0, 100, 200, 300, 400, 500))
    
}

curve(10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[7]) * 10^(winger2006k * 
         (exp(-(10^(coef(winger2006EXPDIndGNLS)[10] + coef(winger2006EXPDIndGNLS)[16])) * 10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[7]) 
              * x) - 1)), log = "x", cex.lab = 1.25,
      xlim = c(5, 2000), ylim = c(0, 300), from = .01, to = 2000,
      col = "orange", ylab = "Mean Consumption", xlab = "FR", lwd = 1.5, main = "Remifentanil", cex.main = 2, axes = FALSE)
curve(10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[8]) * 10^(winger2006k * 
         (exp(-(10^(coef(winger2006EXPDIndGNLS)[10] + coef(winger2006EXPDIndGNLS)[17])) * 10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[8]) 
              * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 2000)
curve(10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[9]) * 10^(winger2006k * 
         (exp(-(10^(coef(winger2006EXPDIndGNLS)[10] + coef(winger2006EXPDIndGNLS)[18])) * 10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[9]) 
              * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 2000)

curve(10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[7]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10] + coef(winger2006SNDIndGNLS)[16])*
                    10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[7])*x)), add  =TRUE,from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[8]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10] + coef(winger2006SNDIndGNLS)[17])*
                    10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[8])*x)), add  =TRUE,from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[9]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10] + coef(winger2006SNDIndGNLS)[18])*
                    10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[9])*x)), add  =TRUE, from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)

points(x = rep(winger2006xVals, 3), 
       y =aggregate(y ~ x + dose, subset(winger2006Data, drug == "remifentanil"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 6), cex = 1.5)
arrows(x0 = rep(winger2006xVals), x1 = rep(winger2006xVals), 
       y0= aggregate(y ~ x + dose, subset(winger2006Data, drug == "remifentanil"),mean, na.rm = TRUE)$y - 
           aggregate(y ~ x + dose, subset(winger2006Data, drug == "remifentanil"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(winger2006Data, drug == "remifentanil"),mean, na.rm = TRUE)$y + 
           aggregate(y ~ x + dose, subset(winger2006Data, drug == "remifentanil"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(jumps = 7), tcl = -.25, labels = NA)
axis(1, at = 10^(-1:4), labels = as.character(c(.1, 1, 10, 100, 1000, 2000)))
axis(2, at = c(0, 100, 200, 300, 400, 500))

for(a in 10:12){
    curve(10^(winger2006Q0EXPD)[a] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a] * 10^(winger2006Q0EXPD)[a] * x) - 1)), 
          xlim = c(5, 2000), ylim = c(0, 300), log = "x", from = .01, to = 2000, col = "orange", cex.main = 1.25, cex.lab = 1.25,
          axes = FALSE, xlab = "FR", ylab = "Consumption", main = paste("Monkey", unique(winger2006Data$monkey)[a-9] ), lwd = 1.25)
    curve(10^(winger2006Q0EXPD)[a+3] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a+3] * 10^(winger2006Q0EXPD)[a+3] * x) - 1)), 
          from = .01, to = 2000, col = "orange", add = TRUE, lwd = 1.25)
    curve(10^(winger2006Q0EXPD)[a+3*2] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a+3*2] * 10^(winger2006Q0EXPD)[a+3*2] * x) - 1)), 
          from = .01, to = 2000, col = "orange", add = TRUE, lwd = 1.25)
    curve(10^winger2006Q0SND[a] * exp(- 10^winger2006AlphaSND[a] * 10^winger2006Q0SND[a] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    curve(10^winger2006Q0SND[a+3] * exp(- 10^winger2006AlphaSND[a+3] * 10^winger2006Q0SND[a+3] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    curve(10^winger2006Q0SND[a+3*2] * exp(- 10^winger2006AlphaSND[a+3*2] * 10^winger2006Q0SND[a+3*2] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(subset(winger2006Data, monkey == unique(winger2006Data$monkey)[a-9] & 
                      drug == unique(winger2006Data$drug)[2], select = c(x, y)), 
           pch =rep(c(21, 22, 23), each = 6)) 
    axis(1, at = tickFunc(jumps = 7), tcl = -.25, labels = NA)
    axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
    axis(2, at = c(0, 100, 200, 300, 400, 500))
}


curve(10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[4]) * 10^(winger2006k * 
     (exp(-(10^(coef(winger2006EXPDIndGNLS)[10] + coef(winger2006EXPDIndGNLS)[13])) *
              10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[4]) 
                      * x) - 1)), log = "x", cex.lab = 1.25,
     xlim = c(5, 2000), ylim = c(0, 300), from = .01, to = 2000,
     col = "orange", ylab = "Mean Consumption", xlab = "FR", lwd = 1.5, 
     main = NA, cex.main = 2, axes = FALSE)
title(main = "Cocaine +\nRemifentanil", xpd = NA, cex.main = 2)
curve(10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[5]) * 10^(winger2006k * 
     (exp(-(10^(coef(winger2006EXPDIndGNLS)[10] + coef(winger2006EXPDIndGNLS)[14])) * 
              10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[5]) 
          * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 2000)
curve(10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[6]) * 10^(winger2006k * 
         (exp(-(10^(coef(winger2006EXPDIndGNLS)[10] + coef(winger2006EXPDIndGNLS)[15])) * 
                  10^(coef(winger2006EXPDIndGNLS)[1] + coef(winger2006EXPDIndGNLS)[6]) 
          * x) - 1)), add  =TRUE, 
  col = "orange", lwd = 1.5, from = .01, to = 2000)
curve(10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[4]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10] + coef(winger2006SNDIndGNLS)[13])*
                    10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[4])*x)), add  =TRUE,from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[5]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10] + coef(winger2006SNDIndGNLS)[14])*
                    10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[5])*x)), add  =TRUE,from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)
curve(10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[6]) * 
          exp(-(10^(coef(winger2006SNDIndGNLS)[10] + coef(winger2006SNDIndGNLS)[15])*
                    10^(coef(winger2006SNDIndGNLS)[1] + coef(winger2006SNDIndGNLS)[6])*x)), add  =TRUE, from = .01, to = 2000,
      col = "blue", lwd = 1.5, lty = 5)

points(x = rep(winger2006xVals, 3), 
       y =aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaremi"),mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 6), cex = 1.5)
arrows(x0 = rep(winger2006xVals), x1 = rep(winger2006xVals), 
       y0= aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaremi"),mean, na.rm = TRUE)$y - 
           aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaremi"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       y1= aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaremi"),mean, na.rm = TRUE)$y + 
           aggregate(y ~ x + dose, subset(winger2006Data, drug == "cocaremi"),\(x) sd(x, na.rm = TRUE)/sqrt(3))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(jumps = 8), tcl = -.25, labels = NA)
axis(1, at = 10^(-1:4), labels = as.character(c(.1, 1, 10, 100, 1000, 2000)))
axis(2, at = c(0, 100, 200, 300, 400, 500))


for(a in 19:21){
    curve(10^(winger2006Q0EXPD)[a] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a] * 10^(winger2006Q0EXPD)[a] * x) - 1)), 
          xlim = c(5, 2000), ylim = c(0, 300), log = "x", from = .01, to = 2000, col = "orange", axes = FALSE, cex.main = 1.25, cex.lab = 1.25,
          xlab = "FR", ylab = "Consumption", main = paste("Monkey", unique(winger2006Data$monkey)[a-18]), lwd = 1.25)
    curve(10^(winger2006Q0EXPD)[a+3] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a+3] * 10^(winger2006Q0EXPD)[a+3] * x) - 1)), 
          from = .01, to = 2000, col = "orange", add = TRUE, lwd = 1.25)
    curve(10^(winger2006Q0EXPD)[a+3*2] * 10^(winger2006k * (exp(-10^(winger2006AlphaEXPD)[a+3*2] * 10^(winger2006Q0EXPD)[a+3*2] * x) - 1)), 
          from = .01, to = 2000, col = "orange", add = TRUE, lwd = 1.25)
    curve(10^winger2006Q0SND[a] * exp(- 10^winger2006AlphaSND[a] * 10^winger2006Q0SND[a] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    curve(10^winger2006Q0SND[a+3] * exp(- 10^winger2006AlphaSND[a+3] * 10^winger2006Q0SND[a+3] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    curve(10^winger2006Q0SND[a+3*2] * exp(- 10^winger2006AlphaSND[a+3*2] * 10^winger2006Q0SND[a+3*2] * x), 
          from = .01, to = 2000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(subset(winger2006Data, monkey == unique(winger2006Data$monkey)[a-18] & 
                      drug == unique(winger2006Data$drug)[3], select = c(x, y)), 
           pch =rep(c(21, 22, 23), each = 6)) 
    axis(1, at = tickFunc(jumps = 8), tcl = -.25, labels = NA)
    axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
    axis(2, at = c(0, 100, 200, 300, 400, 500))
}


dev.off()


# Winger 2006 R2/MAE ----
winger2006Outputs <- winger2006Outputs |> 
  mutate(
    r2Exp = pmap(list(Exp, data), rSquarer),
    MAEExp = pmap(list(Exp, data), MAEmaker),
    AICcExp = pmap(list(Exp, data), AICcmaker),
    r2Free = pmap(list(FreeK, data), rSquarer),
    MAEFree = pmap(list(FreeK, data), MAEmaker),
    AICcFree = pmap(list(FreeK, data), AICcmaker),
    r2SND = pmap(list(SND, data), rSquarer),
    MAESND = pmap(list(SND, data), MAEmaker),
    AICcSND = pmap(list(SND, data), AICcmaker)
  )

rbind(winger2006Outputs |> unnest(c(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND)) |>
        ungroup() |>
        select(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND) |> data.frame() |>
        sapply(range, na.rm = TRUE) |> round(4),
      median = winger2006Outputs |> unnest(c(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND)) |>
        ungroup() |>
        select(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND) |> data.frame() |>
        sapply(median, na.rm = TRUE) |> round(4)
)


# Koffarnus & Woods. (2013) ----
koffarnus2013Data <- fread("koffarnus2013data.csv")
koffarnus2013Data <- melt(koffarnus2013Data, id.vars = c("reinforcer", "Sess", "FR"), variable.name = "rat", value.name = "response")

koffarnus2013Data <- koffarnus2013Data |> 
    transmute(id = rat, x = FR, y = response/FR, session = Sess, reinforcer = reinforcer) |> 
    mutate(id = as.factor(id)) |> filter(session %in% 1:5)

koffarnus2013Data$xy <- koffarnus2013Data$x * koffarnus2013Data$y
koffarnus2013Data$k <- GetK(koffarnus2013Data)

koffarnus2013Exp <- function(ind){
    est<- tryCatch(nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = data.frame(ind),
                        iter = c(10, 10),
                        start_lower = c(q0 = 1, alpha = -5),
                        start_upper = c(q0 = 3, alpha = 0)), error = function(e) NA)
    return(est)
}

koffarnus2013FreeK <- function(ind){
    est<- tryCatch(nls_multstart(y ~ 10^(q0) * 10^(ks * (exp(-10^(alpha) * 10^(q0) * x) - 1)), data = data.frame(ind),
                        iter = c(7, 7, 7),
                        start_lower = c(q0 = 1, ks = 1, alpha = -5),
                        start_upper = c(q0 = 3, ks = 3, alpha = 0),
                        algorithm = "LM",
                        lower = c(-Inf, 0, -Inf),
                        upper = c(Inf, 10, Inf)), error = function(e) NA)
    return(est)
}

koffarnus2013SND <- function(ind){
    est<- tryCatch(nls_multstart(y ~ 10^q0 * exp(-10^alpha *10^q0 * x), data = data.frame(ind),
                        iter = c(10, 10),
                        start_lower = c(q0 = 1, alpha = -5),
                        start_upper = c(q0 = 3, alpha = 0)), error = function(e) NA)
    return(est)
}

koffarnus2013Fits <- koffarnus2013Data |> 
    group_by(id, reinforcer) |> 
    nest() |> 
    mutate(
        Exp = map(data, koffarnus2013Exp),
           FreeK = map(data, koffarnus2013FreeK),
           SND = map(data, koffarnus2013SND))


koffarnus2013Outputs <- koffarnus2013Fits |> mutate(ExpCoef = map(Exp, tryCoef),
                                              FreeKCoef = map(FreeK, tryCoef),
                                              SNDCoef = map(SND, tryCoef))

saveRDS(koffarnus2013Outputs, "koffarnus2013Outputs.rds")
koffarnus2013Outputs <- readRDS("koffarnus2013Outputs.rds")



koffarnus2013Q0EXPD <- na.omit(unlist(koffarnus2013Outputs$ExpCoef)[which(names(unlist(koffarnus2013Outputs$ExpCoef)) == "q0")])
koffarnus2013AlphaEXPD <- na.omit(unlist(koffarnus2013Outputs$ExpCoef)[which(names(unlist(koffarnus2013Outputs$ExpCoef)) == "alpha")])

koffarnus2013Q0SND <- na.omit(unlist(koffarnus2013Outputs$SNDCoef)[which(names(unlist(koffarnus2013Outputs$SNDCoef)) == "q0")])
koffarnus2013AlphaSND <- na.omit(unlist(koffarnus2013Outputs$SNDCoef)[which(names(unlist(koffarnus2013Outputs$SNDCoef)) == "alpha")])


koffarnus2013FreeValues <- koffarnus2013Outputs |> unnest(c(FreeKCoef)) |> pull(FreeKCoef)

koffarnus2013AlphaFreeK <- na.omit(unlist(koffarnus2013FreeValues)[which(names(unlist(koffarnus2013FreeValues)) == "alpha")])
koffarnus2013Q0FreeK <- na.omit(unlist(koffarnus2013FreeValues)[which(names(unlist(koffarnus2013FreeValues)) == "q0")])
koffarnus2013KFreeK <- na.omit(unlist(koffarnus2013FreeValues)[which(names(unlist(koffarnus2013FreeValues)) == "ks")])


# Koffarnus 2013 Pairs----
svg("koffarnus2013pairs.svg", 7.5, 7.5)
pairs.panels(data.frame(
    Q0EXPD = koffarnus2013Q0EXPD,
    Q0Free = koffarnus2013Q0FreeK,
    SNDQ0 = koffarnus2013Q0SND,
    AlphaEXPD = koffarnus2013AlphaEXPD,
    AlphaEXPD = koffarnus2013AlphaFreeK,
    SNDAlpha = koffarnus2013AlphaSND), labels = c(expression(log[10](italic(Q)[0*EXPD])),
                                                  expression(log[10](italic(Q)[0*FreeK])), 
                                                  expression(log[10](italic(Q)[0*SND])), 
                                               expression(log[10](alpha[EXPD])), 
                                               expression(log[10](alpha[FreeK])), 
                                               expression(log[10](alpha[SND]))),
    col = "black", hist.col = "lightblue1",
    ellipses = FALSE, lm = TRUE, main = "Koffarnus & Woods (2013) SND and EXPD Correlations", scale =FALSE)
dev.off()

# Koffarnus & Woods. 2013 GNLS ----

koffarnus2013SNDNLS <- nls_multstart(y ~ (10^q0) * exp(-(10^alpha) * (10^q0) * x),
                                  data = koffarnus2013Data, iter = c(10, 10),
                                  start_lower = c(q0 = 1, alpha = -5),
                                  start_upper = c(q0 = 3, alpha = 1))
koffarnus2013SNDGNLS <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                          data = koffarnus2013Data,
                          params = list(q0 ~ reinforcer,  
                                        alpha ~ reinforcer
                          ),
                          start = list(fixed = c(coef(koffarnus2013SNDNLS)[1], 0,
                                                 coef(koffarnus2013SNDNLS)[2], 0)),
                          verbose = 2,
                          control = list(msMaxIter = 50000,
                                         niterEM = 5000,
                                         maxIter = 5000,
                                         pnlsTol = .01,
                                         tolerance = .01,
                                         apVar = T,
                                         minScale = .0000001,
                                         opt = "optim"), na.action = na.omit)

summary(koffarnus2013SNDGNLS)
anova(koffarnus2013SNDGNLS)
coef(koffarnus2013SNDGNLS)

emmeans(koffarnus2013SNDGNLS, param = "alpha", specs = pairwise~reinforcer)
emmeans(koffarnus2013SNDGNLS, param = "q0", specs =  pairwise~reinforcer)

koffarnus2013EXPDNLS <- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                                   data = koffarnus2013Data, iter = c(10, 10),
                                   start_lower = c(q0 = 1, alpha = -5),
                                   start_upper = c(q0 = 3, alpha = 1))
koffarnus2013EXPDGNLS <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                           data = koffarnus2013Data,
                           params = list(q0 ~ reinforcer,  
                                         alpha ~ reinforcer
                           ),
                           start = list(fixed = c(coef(koffarnus2013EXPDNLS)[1], 0,
                                                  coef(koffarnus2013EXPDNLS)[2], 0)),
                           verbose = 2,
                           control = list(msMaxIter = 50000,
                                          niterEM = 5000,
                                          maxIter = 5000,
                                          pnlsTol = .01,
                                          tolerance = 1,
                                          apVar = T,
                                          minScale = .0000001,
                                          opt = "optim"), na.action = na.omit)


summary(koffarnus2013EXPDGNLS)
anova(koffarnus2013EXPDGNLS)
coef(koffarnus2013EXPDGNLS)

emmeans(koffarnus2013EXPDGNLS, param = "alpha", specs = pairwise~reinforcer)
emmeans(koffarnus2013EXPDGNLS, param = "q0", specs = pairwise~reinforcer, data = koffarnus2013Data)

anova(koffarnus2013EXPDGNLS, koffarnus2013SNDGNLS)

1 - sum((predict(koffarnus2013SNDGNLS) - na.omit(koffarnus2013Data)$y)^2)/sum((mean(koffarnus2013Data$y, na.rm = TRUE) - na.omit(koffarnus2013Data)$y)^2)
1 - sum((predict(koffarnus2013EXPDGNLS) - na.omit(koffarnus2013Data)$y)^2)/sum((mean(koffarnus2013Data$y, na.rm = TRUE) - na.omit(koffarnus2013Data)$y)^2)

emmeans(koffarnus2013EXPDGNLS, param = "alpha", specs = pairwise~reinforcer)
emmeans(koffarnus2013SNDGNLS, param = "alpha", specs = pairwise~reinforcer)

# Koffarnus & Woods 2013 Figures ----
svg("koffarnus2013figure.svg", 11, 11)
koffarnus2013k <- GetK(koffarnus2013Data)
koffarnus2013FoodxVals <- c(1, 3, 10, 32, 100) 
koffarnus2013CocainexVals <- c(1, 3, 10, 18, 32) 

layout(matrix(c(rep(1, 5), 2:26), byrow = TRUE, ncol = 5), heights = c(.05, rep((1-.05)/5, 5)))
par(mar = rep(0,4))
plot.new()
text(.5, .5, "Koffarnus & Woods (2013)", cex = 3)
legend(-.025, 1, c("Cocaine", "Food", "EXPD", "SND"), ncol = 2, xpd = NA,
       bty = "n", pch = c(21, 22, NA, NA), lty = c(NA, NA, 1, 5), col = c(1, 1, "orange", "blue"), pt.cex = 1.5)
par(mar = c(5, 4, 2, .5), pty = "m")

curve(10^(coef(koffarnus2013EXPDGNLS)[1] + coef(koffarnus2013EXPDGNLS)[2]) * 10^(koffarnus2013k * 
         (exp(-(10^(coef(koffarnus2013EXPDGNLS)[3] + coef(koffarnus2013EXPDGNLS)[4])) *
                  10^(coef(koffarnus2013EXPDGNLS)[1] + coef(koffarnus2013EXPDGNLS)[2]) 
              * x) - 1)), log = "x",
      xlim = c(1, 100), ylim = c(0,300), col = "orange", ylab = "Mean Consumption", lwd = 1.5, cex.lab = 1.25,
      main = "Cocaine vs Food", xlab = "FR", cex.main = 2, axes = FALSE, from = .01, to = 1200)
curve(10^(coef(koffarnus2013EXPDGNLS)[1]) * 10^(koffarnus2013k * 
         (exp(-(10^coef(koffarnus2013EXPDGNLS)[3]) * 10^(coef(koffarnus2013EXPDGNLS)[1]) 
              * x) - 1)), add  =TRUE, 
      col = "orange", lwd = 1.5, from = .01, to = 1200)

curve(10^(coef(koffarnus2013SNDGNLS)[1] + coef(koffarnus2013SNDGNLS)[2]) * 
          exp(-(10^(coef(koffarnus2013SNDGNLS)[3] + coef(koffarnus2013SNDGNLS)[4])*
                    10^(coef(koffarnus2013SNDGNLS)[1] + coef(koffarnus2013SNDGNLS)[2])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)

curve(10^(coef(koffarnus2013SNDGNLS)[1]) * 
          exp(-(10^(coef(koffarnus2013SNDGNLS)[3])*
                    10^(coef(koffarnus2013SNDGNLS)[1])*x)), add  =TRUE, from = .01, to = 1200,
      col = "blue", lwd = 1.5, lty = 5)

points(x = c(koffarnus2013CocainexVals, koffarnus2013FoodxVals), 
       y =aggregate(y ~ x + reinforcer, koffarnus2013Data, mean, na.rm = TRUE)$y,
       pch = rep(c(21, 22, 23), each = 6), cex = 1.5)
arrows(x0 = c(koffarnus2013CocainexVals, koffarnus2013FoodxVals), 
       x1 = c(koffarnus2013CocainexVals, koffarnus2013FoodxVals), 
       y0= aggregate(y ~ x + reinforcer, koffarnus2013Data, mean, na.rm = TRUE)$y - 
           aggregate(y ~ x + reinforcer, koffarnus2013Data, \(x) sd(x, na.rm = TRUE)/sqrt(length(x)))$y,
       y1= aggregate(y ~ x + reinforcer, koffarnus2013Data, mean, na.rm = TRUE)$y + 
           aggregate(y ~ x + reinforcer, koffarnus2013Data, \(x) sd(x, na.rm = TRUE)/sqrt(length(x)))$y,
       code = 3, length = .025, angle = 90, lty =1)
axis(1, at = tickFunc(), tcl = -.25, labels = NA)
axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
axis(2, at = c(0, 100, 200, 300))

for(a in unique(koffarnus2013Outputs$id)){
    EXPDModDq0 <- data.table(id = koffarnus2013Outputs$id, est = (unlist(koffarnus2013Outputs$ExpCoef)[which(names(unlist(koffarnus2013Outputs$ExpCoef)) == "q0")]))
    EXPDModDAlpha <- data.table(id = koffarnus2013Outputs$id, est = (unlist(koffarnus2013Outputs$ExpCoef)[which(names(unlist(koffarnus2013Outputs$ExpCoef)) == "alpha")]))
    koffSNDq0 <- data.table(id = koffarnus2013Outputs$id, est = (unlist(koffarnus2013Outputs$SNDCoef)[which(names(unlist(koffarnus2013Outputs$SNDCoef)) == "q0")]))
    koffSNDAlpha <- data.table(id = koffarnus2013Outputs$id, est = (unlist(koffarnus2013Outputs$SNDCoef)[which(names(unlist(koffarnus2013Outputs$SNDCoef)) == "alpha")]))
    curve(10^(EXPDModDq0[id == a][1, est]) * 10^(koffarnus2013k *
        (exp(-10^(EXPDModDAlpha[id == a][1, est]) * 10^(EXPDModDq0[id == a][1, est]) * x) - 1)),
          xlim = c(1, 100), ylim = c(0, 300), log = "x", from = .01, to = 120, col = "orange", axes = FALSE,
          xlab = "FR", ylab = "Consumption", main = paste("Rat", a), lwd = 1.25, cex.lab = 1.25)
    curve(10^koffSNDq0[id == a][1, est] * exp(- 10^koffSNDAlpha[id == a][1, est] * 10^koffSNDq0[id == a][1, est] * x), 
          from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    curve(10^(EXPDModDq0)[id == a][2, est] * 10^(koffarnus2013k *
          (exp(-10^(EXPDModDAlpha[id == a][2, est]) * 10^(EXPDModDq0[id == a][2, est]) * x) - 1)), 
          from = .01, to = 1200, col = "orange", add = TRUE, lwd = 1.25)
    curve(10^koffSNDq0[id == a][2, est] * exp(- 10^koffSNDAlpha[id == a][2, est] * 10^koffSNDq0[id == a][2, est] * x), 
          from = .01, to = 1200, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
    points(subset(koffarnus2013Data, id == a, select = c(x, y)), 
           pch =rep(c(21, 22), each = 5)) 
    axis(1, at = tickFunc(), tcl = -.25, labels = NA)
    axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, NA)))
    axis(2, at = c(0, 100, 200, 300, 400, 500))
}

dev.off()

# Koffarnus 2013 R2/MAE----
koffarnus2013Outputs <- koffarnus2013Outputs |> 
  mutate(
    r2Exp = pmap(list(Exp, data), rSquarer),
    MAEExp = pmap(list(Exp, data), MAEmaker),
    AICcExp = pmap(list(Exp, data), AICcmaker),
    r2Free = pmap(list(FreeK, data), rSquarer),
    MAEFree = pmap(list(FreeK, data), MAEmaker),
    AICcFree = pmap(list(FreeK, data), AICcmaker),
    r2SND = pmap(list(SND, data), rSquarer),
    MAESND = pmap(list(SND, data), MAEmaker),
    AICcSND = pmap(list(SND, data), AICcmaker)
  )

rbind(koffarnus2013Outputs |> unnest(c(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND)) |>
        ungroup() |>
        select(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND) |> data.frame() |>
        sapply(range, na.rm = TRUE) |> round(4),
      median = koffarnus2013Outputs |> unnest(c(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND)) |>
        ungroup() |>
        select(r2Exp,r2Free, r2SND, MAEExp, MAEFree, MAESND, AICcExp, AICcFree, AICcSND) |> data.frame() |>
        sapply(median, na.rm = TRUE) |> round(4)
)


# Free and Fixed k Conversions ----

squareWindowX <- function(a){
  rx <- diff(range(a$x, na.rm = TRUE))
  ry <- diff(range(a$y, na.rm = TRUE))
  if(rx > ry){
    return(range(a$x, na.rm = TRUE))
  } else {
    return(mean(range(a$x, na.rm = TRUE)) + 
             c(-diff(range(a$y, na.rm = TRUE))/2, diff(range(a$y, na.rm = TRUE))/2))
  }
}
squareWindowY <- function(a){
  rx <- diff(range(a$x, na.rm = TRUE))
  ry <- diff(range(a$y, na.rm = TRUE))
  if(rx < ry){
    return(range(a$y, na.rm = TRUE))
  } else {
    return(mean(range(a$y, na.rm = TRUE)) + 
             c(-diff(range(a$x, na.rm = TRUE))/2, diff(range(a$x, na.rm = TRUE))/2))
  }
}

mySeq <- function(a, out = 9){
  b <- round(a, 1)
  pretty(b, n = 8)
}

tenOrLow <- function(x){
  fcase(x == 10, 3,
        x < exp(1)/log(10), 2,
        x < 10 & x >= exp(1)/log(10), 1)
}

svg("conversions.svg", 9, 9.75)

layout(matrix(c(1,1,1,1,
                2,3,4,5,
                6,6,6,6,
                7,8,9,10,
                11,11,11,11,
                12,13,14,15,
                16,16,16,16,
                17,18,19,20,
                21,21,21,21), ncol = 4, byrow = TRUE), 
       heights = c(.035, .84/4, .035, .84/4, .035, .84/4, .035, .84/4, .02))

.035*4 + .02
.74/4
## Ko2002----
par(mar = c(0,0,0,0), pty = "m")
plot.new()
text(.5, .33,"Ko et al. (2002)", cex = 1.75, xpd = NA)

par(mar = c(2,4,.5,1), pty = "s")
ko2002Labels <- rep(1:3, each = 9)
ko2002Unconverted <- data.frame(x = ko2002AlphaFreeK,
                                y = na.omit(ko2002AlphaSND), labels = ko2002Labels)
ko2002SNDConverted <- data.frame(x = log10(10^((ko2002AlphaFreeK))* EXPDToSND((ko2002KFreeK))),
                                 y = na.omit(ko2002AlphaSND), labels = ko2002Labels)

ko2002EXPDConverted <- data.frame(x = log10(10^ko2002AlphaFreeK/-VGAM::lambertW(-1/log(10^na.omit(ko2002KFreeK)))),
                                  y = na.omit(ko2002AlphaSND), labels = ko2002Labels)

ko2002EXPD2Converted <- data.frame(x = log10(10^ko2002AlphaFreeK/-VGAM::lambertW(-1/log(10^na.omit(ko2002KFreeK)))),
                                   y = log10(10^ko2002AlphaEXPD/-VGAM::lambertW(-1/log(10^na.omit(unique(ko2002Data$k))))), labels = ko2002Labels)

plot.new()
plot.window(squareWindowX(ko2002Unconverted),
            squareWindowY(ko2002Unconverted)) 

axis(1, mySeq(squareWindowX(ko2002Unconverted), 8), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(ko2002Unconverted), 8), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()

mtext(expression("SND vs EXPD Free"~italic(k)), cex = .8)
title(ylab = expression(log[10](alpha["SND"])), xlab = expression("Free"~italic(k)~log[10](alpha["EXPD"])), line = 2.1, xpd = NA, cex = 2)
points(ko2002Unconverted[,1:2], col = c("black", "blue", "orange2")[ko2002Unconverted$labels],
       pch = c(1, 4, 5)[tenOrLow(ko2002KFreeK)])
for(a in unique(ko2002Labels)){
  lines(range(subset(ko2002Unconverted, labels == a)$x), 
        range(predict(lm(y ~ x, subset(ko2002Unconverted, labels == a)))),
        col = c("black", "blue", "orange2")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(ko2002Unconverted, labels == 1, select = c(x,y)))[1,2],3))*", Alfentanil"),
                  bquote(italic(r)==.(rd(cor(subset(ko2002Unconverted, labels == 2, select = c(x,y)))[1,2],3))*", Remifentanil"),
                  bquote(italic(r)==.(rd(cor(subset(ko2002Unconverted, labels == 3, select = c(x,y)))[1,2],3))*", Fentanyl")),
       bty = "n", text.col = c("black", "blue", "orange2")[1:3], inset = c(-.075, -.025)
)

plot.new()
plot.window(squareWindowX(ko2002SNDConverted),
            squareWindowY(ko2002SNDConverted)) 

axis(1, mySeq(squareWindowX(ko2002SNDConverted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(ko2002SNDConverted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("EXPD-to-SND"), cex = .8)
title(ylab = expression(log[10](alpha["SND"])), xlab = expression("Free"~italic(k)~log[10](alpha["EXPD-to-SND"])), line = 2.1, xpd = NA, cex = 2)
points(ko2002SNDConverted[,1:2], col = c("black", "blue", "orange2")[ko2002SNDConverted$labels],
       pch = c(1, 4, 5)[tenOrLow(ko2002KFreeK)])
for(a in unique(ko2002Labels)){
  lines(range(subset(ko2002SNDConverted, labels == a)$x), 
        range(predict(lm(y ~ x, subset(ko2002SNDConverted, labels == a)))),
        col = c("black", "blue", "orange2")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(ko2002SNDConverted, labels == 1, select = c(x,y)))[1,2],3))*", Alfentanil"),
                  bquote(italic(r)==.(rd(cor(subset(ko2002SNDConverted, labels == 2, select = c(x,y)))[1,2],3))*", Remifentanil"),
                  bquote(italic(r)==.(rd(cor(subset(ko2002SNDConverted, labels == 3, select = c(x,y)))[1,2],3))*", Fentanyl")),
       bty = "n", text.col = c("black", "blue", "orange2")[1:3], inset = c(-.075, -.025)
)


ko2002kUnconverted <- data.frame(x = ko2002AlphaFreeK,
                                 y = na.omit(ko2002AlphaEXPD), labels = ko2002Labels)

plot.new()
plot.window(squareWindowX(ko2002kUnconverted),
            squareWindowY(ko2002kUnconverted)) 

axis(1, mySeq(squareWindowX(ko2002kUnconverted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(ko2002kUnconverted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("EXPD Fixed vs Free"~italic(k)), cex = .8)
title(ylab = expression("Fixed"~italic(k)~log[10](alpha["EXPD"])), xlab = expression("Free"~italic(k)~log[10](alpha["EXPD"])), line = 2.1, xpd = NA, cex = 2)
points(ko2002kUnconverted[,1:2], col = c("black", "purple", "yellow4")[ko2002kUnconverted$labels],
       pch = c(1, 4, 5)[tenOrLow(ko2002KFreeK)])
for(a in unique(ko2002Labels)){
  lines(range(subset(ko2002kUnconverted, labels == a)$x, na.rm = TRUE), 
        range(predict(lm(y ~ x, subset(ko2002kUnconverted, labels == a))), na.rm = TRUE),
        col = c("black", "purple", "yellow4")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(ko2002kUnconverted, labels == 1, select = c(x,y)))[1,2],3))*", Alfentanil"),
                  bquote(italic(r)==.(rd(cor(subset(ko2002kUnconverted, labels == 2, select = c(x,y)))[1,2],3))*", Remifentanil"),
                  bquote(italic(r)==.(rd(cor(subset(ko2002kUnconverted, labels == 3, select = c(x,y)))[1,2],3))*", Fentanyl")),
       bty = "n", text.col = c("black", "purple", "yellow4")[1:3], inset = c(-.075, -.025)
)


plot.new()
plot.window(squareWindowX(ko2002EXPD2Converted),
            squareWindowY(ko2002EXPD2Converted)) 

axis(1, mySeq(squareWindowX(ko2002EXPD2Converted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(ko2002EXPD2Converted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression(italic(W)(italic(k))), cex = .8)
title(ylab = expression("Fixed"~italic(k)~log[10](alpha[italic(W)(italic(k))])),
      xlab = expression("Free"~italic(k)~log[10](alpha[italic(W)(italic(k))])), line = 2.1, xpd = NA, cex = 2)
points(ko2002EXPD2Converted[,1:2], col = c("black", "purple", "yellow4")[ko2002EXPD2Converted$labels],
       pch = c(1, 4, 5)[tenOrLow(ko2002KFreeK)])
for(a in unique(ko2002Labels)){
  lines(range(subset(ko2002EXPD2Converted, labels == a)$x, na.rm = TRUE), 
        range(predict(lm(y ~ x, subset(ko2002EXPD2Converted, labels == a))), na.rm = TRUE),
        col = c("black", "purple", "yellow4")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(ko2002EXPD2Converted, labels == 1, select = c(x,y)), use = "c")[1,2],3))*", Alfentanil"),
                  bquote(italic(r)==.(rd(cor(subset(ko2002EXPD2Converted, labels == 2, select = c(x,y)), use = "c")[1,2],3))*", Remifentanil"),
                  bquote(italic(r)==.(rd(cor(subset(ko2002EXPD2Converted, labels == 3, select = c(x,y)), use = "c")[1,2],3))*", Fentanyl")),
       bty = "n", text.col = c("black", "purple", "yellow4")[1:3], inset = c(-.075, -.025)
)


## Winger2002 ----
par(mar = c(0,0,0,0), pty = "m")
plot.new()
text(.5, .33,"Winger et al. (2002)", cex = 1.75, xpd = NA)

par(mar = c(2,4,.5,1), pty = "s")


winger2002Labels <- rep(1:3, each = 9)
winger2002Unconverted <- data.frame(x = winger2002AlphaFreeK,
                                    y = na.omit(winger2002AlphaSND), labels = winger2002Labels)
winger2002SNDConverted <- data.frame(x = log10(10^((winger2002AlphaFreeK))* EXPDToSND((winger2002KFreeK))),
                                     y = na.omit(winger2002AlphaSND), labels = winger2002Labels)

winger2002EXPDConverted <- data.frame(x = log10(10^winger2002AlphaFreeK/-VGAM::lambertW(-1/log(10^na.omit(winger2002KFreeK)))),
                                      y = na.omit(winger2002AlphaSND), labels = winger2002Labels)

winger2002kUnconverted <- data.frame(x = winger2002AlphaFreeK,
                                     y = na.omit(winger2002AlphaEXPD), labels = winger2002Labels)

winger2002EXPD2Converted <- data.frame(x = log10(10^winger2002AlphaFreeK/-VGAM::lambertW(-1/log(10^na.omit(winger2002KFreeK)))),
                                       y = log10(10^winger2002AlphaEXPD/-VGAM::lambertW(-1/log(10^na.omit(unique(winger2002Data$k))))), labels = winger2002Labels)



plot.new()
plot.window(squareWindowX(winger2002Unconverted),
            squareWindowY(winger2002Unconverted)) 

axis(1, mySeq(squareWindowX(winger2002Unconverted), 6), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(winger2002Unconverted), 8), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("SND vs EXPD Free"~italic(k)), cex = .8)
title(ylab = expression(log[10](alpha["SND"])),
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD"])), line = 2.1, xpd = NA, cex = 2)

points(winger2002Unconverted[,1:2], col = c("black", "blue", "orange2")[winger2002Unconverted$labels],
       pch = c(1, 4, 5)[tenOrLow(winger2002KFreeK)])
for(a in unique(winger2002Labels)){
  lines(range(subset(winger2002Unconverted, labels == a)$x), 
        range(predict(lm(y ~ x, subset(winger2002Unconverted, labels == a)))),
        col = c("black", "blue", "orange2")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(winger2002Unconverted, labels == 1, select = c(x,y)))[1,2],3))*", Ketamine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2002Unconverted, labels == 2, select = c(x,y)))[1,2],3))*", Phencyclidine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2002Unconverted, labels == 3, select = c(x,y)))[1,2],3))*", Dizocilpine")),
       bty = "n", text.col = c("black", "blue", "orange2")[1:3], inset = c(-.075, -.025)
)



plot.new()
plot.window(squareWindowX(winger2002SNDConverted),
            squareWindowY(winger2002SNDConverted)) 

axis(1, mySeq(squareWindowX(winger2002SNDConverted),8), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(winger2002SNDConverted),6), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("EXPD-to-SND"), cex = .8)
title(ylab = expression(log[10](alpha["SND"])),
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD-to-SND"])), line = 2.1, xpd = NA, cex = 2)

points(winger2002SNDConverted[,1:2], col = c("black", "blue", "orange2")[winger2002SNDConverted$labels],
       pch = c(1, 4, 5)[tenOrLow(winger2002KFreeK)])
for(a in unique(winger2002Labels)){
  lines(range(subset(winger2002SNDConverted, labels == a)$x), 
        range(predict(lm(y ~ x, subset(winger2002SNDConverted, labels == a)))),
        col = c("black", "blue", "orange2")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(winger2002SNDConverted, labels == 1, select = c(x,y)))[1,2],3))*", Ketamine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2002SNDConverted, labels == 2, select = c(x,y)))[1,2],3))*", Phencyclidine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2002SNDConverted, labels == 3, select = c(x,y)))[1,2],3))*", Dizocilpine")),
       bty = "n", text.col = c("black", "blue", "orange2")[1:3], inset = c(-.075, -.025)
)



winger2002kUnconverted <- data.frame(x = winger2002AlphaFreeK,
                                     y = na.omit(winger2002AlphaEXPD), labels = winger2002Labels)

plot.new()
plot.window(squareWindowX(winger2002kUnconverted),
            squareWindowY(winger2002kUnconverted)) 

axis(1, mySeq(squareWindowX(winger2002kUnconverted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(winger2002kUnconverted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("EXPD Fixed vs Free"~italic(k)), cex = .8)
title(ylab = expression("Fixed"~italic(k)~log[10](alpha["EXPD"])), 
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD"])), line = 2.1, xpd = NA, cex = 2)

points(winger2002kUnconverted[,1:2], col = c("black", "purple", "yellow4")[winger2002EXPDConverted$labels],
       pch = c(1, 4, 5)[tenOrLow(winger2002KFreeK)])
for(a in unique(winger2002Labels)){
  lines(range(subset(winger2002kUnconverted, labels == a)$x, na.rm = TRUE), 
        range(predict(lm(y ~ x, subset(winger2002kUnconverted, labels == a))), na.rm = TRUE),
        col = c("black", "purple", "yellow4")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(winger2002kUnconverted, labels == 1, select = c(x,y)))[1,2],3))*", Ketamine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2002kUnconverted, labels == 2, select = c(x,y)))[1,2],3))*", Phencyclidine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2002kUnconverted, labels == 3, select = c(x,y)))[1,2],3))*", Dizocilpine")),
       bty = "n", text.col = c("black", "purple", "yellow4")[1:3], inset = c(-.075, -.025)
)


plot.new()
plot.window(squareWindowX(winger2002EXPD2Converted),
            squareWindowY(winger2002EXPD2Converted)) 

axis(1, mySeq(squareWindowX(winger2002EXPD2Converted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(winger2002EXPD2Converted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression(italic(W)(italic(k))), cex = .8)
title(ylab = expression("Fixed"~italic(k)~log[10](alpha[italic(W)(italic(k))])), 
      xlab = expression("Free"~italic(k)~log[10](alpha[italic(W)(italic(k))])), line = 2.1, xpd = NA, cex = 2)

points(winger2002EXPD2Converted[,1:2], col = c("black", "purple", "yellow4")[winger2002EXPD2Converted$labels],
       pch = c(1, 4, 5)[tenOrLow(winger2002KFreeK)])
for(a in unique(winger2002Labels)){
  lines(range(subset(winger2002EXPD2Converted, labels == a)$x, na.rm = TRUE), 
        range(predict(lm(y ~ x, subset(winger2002EXPD2Converted, labels == a))), na.rm = TRUE),
        col = c("black", "purple", "yellow4")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(winger2002EXPD2Converted, labels == 1, select = c(x,y)), use = "c")[1,2],3))*", Ketamine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2002EXPD2Converted, labels == 2, select = c(x,y)), use = "c")[1,2],3))*", Phencyclidine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2002EXPD2Converted, labels == 3, select = c(x,y)), use = "c")[1,2],3))*", Dizocilpine")),
       bty = "n", text.col = c("black", "purple", "yellow4")[1:3], inset = c(-.075, -.025)
)




## Winger2006 ----
par(mar = c(0,0,0,0), pty = "m")
plot.new()
text(.5, .33,"Winger et al. (2006)", cex = 1.75, xpd = NA)

par(mar = c(2,4,.5,1), pty = "s")

winger2006Labels <- rep(1:3, each = 9)
winger2006Unconverted <- data.frame(x = winger2006AlphaFreeK,
                                    y = na.omit(winger2006AlphaSND), labels = winger2006Labels)
winger2006SNDConverted <- data.frame(x = log10(10^((winger2006AlphaFreeK))* EXPDToSND((winger2006KFreeK))),
                                     y = na.omit(winger2006AlphaSND), labels = winger2006Labels)

winger2006EXPDConverted <- data.frame(x = log10(10^winger2006AlphaFreeK/-VGAM::lambertW(-1/log(10^na.omit(winger2006KFreeK)))),
                                      y = na.omit(winger2006AlphaSND), labels = winger2006Labels)

winger2006EXPD2Converted <- data.frame(x = log10(10^winger2006AlphaFreeK/-VGAM::lambertW(-1/log(10^na.omit(winger2006KFreeK)))),
                                       y = log10(10^winger2006AlphaEXPD/-VGAM::lambertW(-1/log(10^na.omit(unique(winger2006Data$k))))), labels = winger2006Labels)



plot.new()
plot.window(squareWindowX(winger2006Unconverted),
            squareWindowY(winger2006Unconverted)) 

axis(1, mySeq(squareWindowX(winger2006Unconverted), 6), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(winger2006Unconverted), 8), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("SND vs EXPD Free"~italic(k)), cex = .8)
title(ylab = expression(log[10](alpha["SND"])), 
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD"])), line = 2.1, xpd = NA, cex = 2)

points(winger2006Unconverted[,1:2], col = c("black", "blue", "orange2")[winger2006Unconverted$labels],
       pch = c(1, 4, 5)[tenOrLow(winger2006KFreeK)])
for(a in unique(winger2006Labels)){
  lines(range(subset(winger2006Unconverted, labels == a)$x), 
        range(predict(lm(y ~ x, subset(winger2006Unconverted, labels == a)))),
        col = c("black", "blue", "orange2")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(winger2006Unconverted, labels == 1, select = c(x,y)))[1,2],3))*", Cocaine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2006Unconverted, labels == 2, select = c(x,y)))[1,2],3))*", Remifentanil"),
                  bquote(italic(r)==.(rd(cor(subset(winger2006Unconverted, labels == 3, select = c(x,y)))[1,2],3))*", Combined")),
       bty = "n", text.col = c("black", "blue", "orange2")[1:3], inset = c(-.075, -.025)
)



plot.new()
plot.window(squareWindowX(winger2006SNDConverted),
            squareWindowY(winger2006SNDConverted)) 

axis(1, mySeq(squareWindowX(winger2006SNDConverted),8), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(winger2006SNDConverted),6), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("EXPD-to-SND"), cex = .8)
title(ylab = expression(log[10](alpha["SND"])),
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD-to-SND"])), line = 2.1, xpd = NA, cex = 2)

points(winger2006SNDConverted[,1:2], col = c("black", "blue", "orange2")[winger2006SNDConverted$labels],
       pch = c(1, 4, 5)[tenOrLow(winger2006KFreeK)])
for(a in unique(winger2006Labels)){
  lines(range(subset(winger2006SNDConverted, labels == a)$x), 
        range(predict(lm(y ~ x, subset(winger2006SNDConverted, labels == a)))),
        col = c("black", "blue", "orange2")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(winger2006SNDConverted, labels == 1, select = c(x,y)))[1,2],3))*", Cocaine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2006SNDConverted, labels == 2, select = c(x,y)))[1,2],3))*", Remifentanil"),
                  bquote(italic(r)==.(rd(cor(subset(winger2006SNDConverted, labels == 3, select = c(x,y)))[1,2],3))*", Combined")),
       bty = "n", text.col = c("black", "blue", "orange2")[1:3], inset = c(-.075, -.025)
)





winger2006kUnconverted <- data.frame(x = winger2006AlphaFreeK,
                                     y = na.omit(winger2006AlphaEXPD), labels = winger2006Labels)

plot.new()
plot.window(squareWindowX(winger2006kUnconverted),
            squareWindowY(winger2006kUnconverted)) 

axis(1, mySeq(squareWindowX(winger2006kUnconverted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(winger2006kUnconverted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("EXPD Fixed vs Free"~italic(k)), cex = .8)
title(ylab = expression("Fixed"~italic(k)~log[10](alpha["EXPD"])),
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD"])), line = 2.1, xpd = NA, cex = 2)

points(winger2006kUnconverted[,1:2], col = c("black", "purple", "yellow4")[winger2006kUnconverted$labels],
       pch = c(1, 4, 5)[tenOrLow(winger2006KFreeK)])
for(a in unique(winger2006Labels)){
  lines(range(subset(winger2006kUnconverted, labels == a)$x, na.rm = TRUE), 
        range(predict(lm(y ~ x, subset(winger2006kUnconverted, labels == a))), na.rm = TRUE),
        col = c("black", "purple", "yellow4")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(winger2006kUnconverted, labels == 1, select = c(x,y)))[1,2],3))*", Cocaine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2006kUnconverted, labels == 2, select = c(x,y)))[1,2],3))*", Remifentanil"),
                  bquote(italic(r)==.(rd(cor(subset(winger2006kUnconverted, labels == 3, select = c(x,y)))[1,2],3))*", Combined")),
       bty = "n", text.col = c("black", "purple", "yellow4")[1:3], inset = c(-.075, -.025)
)


plot.new()
plot.window(squareWindowX(winger2006EXPD2Converted),
            squareWindowY(winger2006EXPD2Converted)) 

axis(1, mySeq(squareWindowX(winger2006EXPD2Converted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(winger2006EXPD2Converted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression(italic(W)(italic(k))), cex = .8)
title(ylab = expression("Fixed"~italic(k)~log[10](alpha[italic(W)(italic(k))])), 
      xlab = expression("Free"~italic(k)~log[10](alpha[italic(W)(italic(k))])), line = 2.1, xpd = NA, cex = 2)

points(winger2006EXPD2Converted[,1:2], col = c("black", "purple", "yellow4")[winger2006EXPD2Converted$labels],
       pch = c(1, 4, 5)[tenOrLow(winger2006KFreeK)])
for(a in unique(winger2006Labels)){
  lines(range(subset(winger2006EXPD2Converted, labels == a)$x, na.rm = TRUE), 
        range(predict(lm(y ~ x, subset(winger2006EXPD2Converted, labels == a))), na.rm = TRUE),
        col = c("black", "purple", "yellow4")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(winger2006EXPD2Converted, labels == 1, select = c(x,y)), use = "c")[1,2],3))*", Cocaine"),
                  bquote(italic(r)==.(rd(cor(subset(winger2006EXPD2Converted, labels == 2, select = c(x,y)), use = "c")[1,2],3))*", Remifentanil"),
                  bquote(italic(r)==.(rd(cor(subset(winger2006EXPD2Converted, labels == 3, select = c(x,y)), use = "c")[1,2],3))*", Combined")),
       bty = "n", text.col = c("black", "purple", "yellow4")[1:3], inset = c(-.075, -.025)
)




## Koffarnus2013 ----
par(mar = c(0,0,0,0), pty = "m")
plot.new()
text(.5, .33,"Koffarnus & Woods (2013)", cex = 1.75, xpd = NA)

par(mar = c(2,4,.5,1), pty = "s")


koffarnus2013Labels <- as.numeric(as.factor(subset(koffarnus2013Outputs, !is.na(Exp))$reinforcer))

koffarnus2013Unconverted <- data.frame(x = na.omit(koffarnus2013AlphaFreeK),
                                       y = na.omit(koffarnus2013AlphaSND), labels = koffarnus2013Labels)
koffarnus2013SNDConverted <- data.frame(x = log10(10^(na.omit(koffarnus2013AlphaFreeK))* EXPDToSND(na.omit(koffarnus2013KFreeK))),
                                        y = na.omit(koffarnus2013AlphaSND), labels = koffarnus2013Labels)

koffarnus2013EXPDConverted <- data.frame(x = log10(10^na.omit(koffarnus2013AlphaFreeK)/-VGAM::lambertW(-1/log(10^na.omit(koffarnus2013KFreeK)))),
                                         y = na.omit(koffarnus2013AlphaSND), labels = koffarnus2013Labels)

koffarnus2013EXPD2Converted <- data.frame(x = log10(10^na.omit(koffarnus2013AlphaFreeK)/-VGAM::lambertW(-1/log(10^na.omit(koffarnus2013KFreeK)))),
                                          y = log10(10^na.omit(koffarnus2013AlphaEXPD)/-VGAM::lambertW(-1/log(10^na.omit(unique(koffarnus2013Data$k))))), labels = koffarnus2013Labels)



plot.new()
plot.window(squareWindowX(koffarnus2013Unconverted),
            squareWindowY(koffarnus2013Unconverted)) 

lapply(koffarnus2013Unconverted, range)

axis(1, mySeq(squareWindowX(koffarnus2013Unconverted), 6), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(koffarnus2013Unconverted), 8), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("SND vs EXPD Free"~italic(k)), cex = .8)
title(ylab = expression(log[10](alpha["SND"])), 
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD"])), line = 2.1, xpd = NA, cex = 2)

points(koffarnus2013Unconverted[,1:2], col = c("black", "blue", "orange2")[koffarnus2013Unconverted$labels],
       pch = c(1, 4, 5)[tenOrLow(koffarnus2013KFreeK)])
for(a in unique(koffarnus2013Labels)){
  lines(range(subset(koffarnus2013Unconverted, labels == a)$x), 
        range(predict(lm(y ~ x, subset(koffarnus2013Unconverted, labels == a)))),
        col = c("black", "blue", "orange2")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(koffarnus2013Unconverted, labels == 1, select = c(x,y)))[1,2],3))*", Cocaine"),
                  bquote(italic(r)==.(rd(cor(subset(koffarnus2013Unconverted, labels == 2, select = c(x,y)))[1,2],3))*", Food")),
       bty = "n", text.col = c("black", "blue", "orange2")[1:3], inset = c(-.075, -.025)
)


koffarnus2013Outputs
plot.new()
plot.window(squareWindowX(koffarnus2013SNDConverted),
            squareWindowY(koffarnus2013SNDConverted)) 

axis(1, mySeq(squareWindowX(koffarnus2013SNDConverted),8), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(koffarnus2013SNDConverted),6), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("EXPD-to-SND"), cex = .8)
title(ylab = expression(log[10](alpha["SND"])), 
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD-to-SND"])), line = 2.1, xpd = NA, cex = 2)


points(koffarnus2013SNDConverted[,1:2], col = c("black", "blue", "orange2")[koffarnus2013SNDConverted$labels],
       pch = c(1, 4, 5)[tenOrLow(koffarnus2013KFreeK)])
for(a in unique(koffarnus2013Labels)){
  lines(range(subset(koffarnus2013SNDConverted, labels == a)$x), 
        range(predict(lm(y ~ x, subset(koffarnus2013SNDConverted, labels == a)))),
        col = c("black", "blue", "orange2")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(koffarnus2013SNDConverted, labels == 1, select = c(x,y)))[1,2],3))*", Cocaine"),
                  bquote(italic(r)==.(rd(cor(subset(koffarnus2013SNDConverted, labels == 2, select = c(x,y)))[1,2],3))*", Food")),
       bty = "n", text.col = c("black", "blue", "orange2")[1:3], inset = c(-.075, -.025)
)



koffarnus2013kUnconverted <- data.frame(x = na.omit(koffarnus2013AlphaFreeK),
                                        y = na.omit(koffarnus2013AlphaEXPD), labels = koffarnus2013Labels)

plot.new()
plot.window(squareWindowX(koffarnus2013kUnconverted),
            squareWindowY(koffarnus2013kUnconverted)) 

axis(1, mySeq(squareWindowX(koffarnus2013kUnconverted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(koffarnus2013kUnconverted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression("EXPD Fixed vs Free"~italic(k)), cex = .8)
title(ylab = expression("Fixed"~italic(k)~log[10](alpha["EXPD"])), 
      xlab = expression("Free"~italic(k)~log[10](alpha["EXPD"])), line = 2.1, xpd = NA, cex = 2)

points(koffarnus2013kUnconverted[,1:2], col = c("black", "purple", "yellow4")[koffarnus2013kUnconverted$labels],
       pch = c(1, 4, 5)[tenOrLow(koffarnus2013KFreeK)])
for(a in unique(koffarnus2013Labels)){
  lines(range(subset(koffarnus2013kUnconverted, labels == a)$x, na.rm = TRUE), 
        range(predict(lm(y ~ x, subset(koffarnus2013kUnconverted, labels == a))), na.rm = TRUE),
        col = c("black", "purple", "yellow4")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(koffarnus2013kUnconverted, labels == 1, select = c(x,y)))[1,2],3))*", Cocaine"),
                  bquote(italic(r)==.(rd(cor(subset(koffarnus2013kUnconverted, labels == 2, select = c(x,y)))[1,2],3))*", Food")),
       bty = "n", text.col = c("black", "purple", "yellow4")[1:2], inset = c(-.075, -.025)
)

plot.new()
plot.window(squareWindowX(koffarnus2013EXPD2Converted),
            squareWindowY(koffarnus2013EXPD2Converted)) 

axis(1, mySeq(squareWindowX(koffarnus2013EXPD2Converted)), cex.axis = .8, padj = -.5)
axis(2, mySeq(squareWindowY(koffarnus2013EXPD2Converted)), cex.axis = .8, padj = .5)
abline(0,1, col = "grey")
box()
mtext(expression(italic(W)(italic(k))), cex = .8)
title(ylab = expression("Fixed"~italic(k)~log[10](alpha[italic(W)(italic(k))])), 
      xlab = expression("Free"~italic(k)~log[10](alpha[italic(W)(italic(k))])), line = 2.1, xpd = NA, cex = 2)

points(koffarnus2013EXPD2Converted[,1:2], col = c("black", "purple", "yellow4")[koffarnus2013EXPD2Converted$labels],
       pch = c(1, 4, 5)[tenOrLow(koffarnus2013KFreeK)])
for(a in unique(koffarnus2013Labels)){
  lines(range(subset(koffarnus2013EXPD2Converted, labels == a)$x, na.rm = TRUE), 
        range(predict(lm(y ~ x, subset(koffarnus2013EXPD2Converted, labels == a))), na.rm = TRUE),
        col = c("black", "purple", "yellow4")[a])
}
legend("topleft", 
       legend = c(bquote(italic(r)==.(rd(cor(subset(koffarnus2013EXPD2Converted, labels == 1, select = c(x,y)), use = "c")[1,2],3))*", Cocaine"),
                  bquote(italic(r)==.(rd(cor(subset(koffarnus2013EXPD2Converted, labels == 2, select = c(x,y)), use = "c")[1,2],3))*", Food")),
       bty = "n", text.col = c("black", "purple", "yellow4")[1:3], inset = c(-.075, -.025)
)

dev.off()




# EV and Pmax Comparisons ----
sum(exp(1)/log(10) <= c(ko2002KFreeK, winger2002KFreeK, winger2006KFreeK, koffarnus2013KFreeK)) / length(c(ko2002KFreeK, winger2002KFreeK, winger2006KFreeK, koffarnus2013KFreeK))


convertedEVs <- rbind(
  data.frame(SNDEV = 1/(10^ ko2002AlphaSND), 
             CONEV = 1/(10^ko2002AlphaFreeK * 
                          EXPDToSND(ko2002KFreeK))
  ),
  
  data.frame(SNDEV = 1/(10^ winger2002AlphaSND), 
             CONEV = 1/(10^winger2002AlphaFreeK * 
                          EXPDToSND(winger2002KFreeK))
  ),
  
  data.frame(SNDEV = 1/(10^ winger2006AlphaSND), 
             CONEV = 1/(10^winger2006AlphaFreeK * 
                          EXPDToSND(winger2006KFreeK))
  ),
  
  
  data.frame(SNDEV = 1/(10^ koffarnus2013AlphaSND), 
             CONEV = 1/(10^koffarnus2013AlphaFreeK * EXPDToSND(koffarnus2013KFreeK))
  )
)


EVOrganized <- data.frame(EV = (na.omit(c(convertedEVs$SNDEV, convertedEVs$CONEV))),
                          study = rep(c(
                            rep("ako2002", nrow(ko2002Outputs)),
                            rep("bwinger2002", nrow(winger2002Outputs)),
                            rep("cwinger2006", nrow(winger2006Outputs)),
                            rep("dkoffarnus2013", length(na.omit(koffarnus2013AlphaSND)))
                          ),2),
                          model = c(rep("aSND", nrow(na.omit(convertedEVs))),
                                    rep("bCon", nrow(na.omit(convertedEVs))))
)


nrow(ko2002Outputs)
nrow(winger2002Outputs)
nrow(winger2006Outputs)
length(na.omit(koffarnus2013AlphaSND))

convertedPmaxs <- rbind(
  data.frame(SNDPmax = 1/(10^ ko2002AlphaSND * 10^ko2002Q0SND), 
             CONPmax = 1/(10^ko2002AlphaFreeK * 10^ko2002Q0FreeK * 
                            EXPDToSND(ko2002KFreeK))
  ),
  
  data.frame(SNDPmax = 1/(10^ winger2002AlphaSND * 10^winger2002Q0SND), 
             CONPmax = 1/(10^winger2002AlphaFreeK * 10^winger2002Q0FreeK *
                            EXPDToSND(winger2002KFreeK))
  ),
  
  data.frame(SNDPmax = 1/(10^ winger2006AlphaSND * 10^winger2006Q0SND), 
             CONPmax = 1/(10^winger2006AlphaFreeK * 10^winger2006Q0FreeK *
                            EXPDToSND(winger2006KFreeK))
  ),
  
  
  data.frame(SNDPmax = 1/(10^ koffarnus2013AlphaSND * 10^koffarnus2013Q0SND), 
             CONPmax = 1/(10^koffarnus2013AlphaFreeK * 10^koffarnus2013Q0FreeK *
                            EXPDToSND(koffarnus2013KFreeK))
  )
)

PmaxOrganized <- data.frame(Pmax = (na.omit(c(convertedPmaxs$SNDPmax, convertedPmaxs$CONPmax))),
                            study = rep(c(
                              rep("ako2002", nrow(ko2002Outputs)),
                              rep("bwinger2002", nrow(winger2002Outputs)),
                              rep("cwinger2006", nrow(winger2006Outputs)),
                              rep("dkoffarnus2013", length(na.omit(koffarnus2013AlphaSND)))
                            ),2),
                            model = c(rep("aSND", nrow(na.omit(convertedPmaxs))),
                                      rep("bCon", nrow(na.omit(convertedPmaxs)))))

convertedOmaxs <- rbind(
  data.frame(SNDOmax = 1/(10^ ko2002AlphaSND * exp(1)), 
             CONOmax = 1/(10^ko2002AlphaFreeK * exp(1) * 
                            EXPDToSND(ko2002KFreeK))
  ),
  
  data.frame(SNDOmax = 1/(10^ winger2002AlphaSND * exp(1)), 
             CONOmax = 1/(10^winger2002AlphaFreeK * exp(1) *
                            EXPDToSND(winger2002KFreeK))
  ),
  
  data.frame(SNDOmax = 1/(10^ winger2006AlphaSND * exp(1)), 
             CONOmax = 1/(10^winger2006AlphaFreeK * exp(1) *
                            EXPDToSND(winger2006KFreeK))
  ),
  
  
  data.frame(SNDOmax = 1/(10^na.omit(koffarnus2013AlphaSND) * exp(1)), 
             CONOmax = 1/(10^na.omit(koffarnus2013AlphaFreeK) * exp(1) *
                            EXPDToSND(koffarnus2013KFreeK))
  )
)

OmaxOrganized <- data.frame(Omax = (na.omit(c(convertedOmaxs$SNDOmax, convertedOmaxs$CONOmax))),
                            study = rep(c(
                              rep("ako2002", nrow(ko2002Outputs)),
                              rep("bwinger2002", nrow(winger2002Outputs)),
                              rep("cwinger2006", nrow(winger2006Outputs)),
                              rep("dkoffarnus2013", length(na.omit(koffarnus2013AlphaSND)))
                            ),2),
                            model = c(rep("aSND", nrow(na.omit(convertedOmaxs))),
                                      rep("bCon", nrow(na.omit(convertedOmaxs)))))


for(a in unique(EVOrganized$study)){
  print(a)
  print(round(c(unlist((t.test(
    log10(subset(EVOrganized, study == a & model == "aSND")$EV),
    log10(subset(EVOrganized, study == a & model == "bCon")$EV))[c("p.value")])),
    corr.test(log10(subset(EVOrganized, study == a & model == "aSND")$EV),
              log10(subset(EVOrganized, study == a & model == "bCon")$EV))$r), 4)
  )
  print("aSND")
  print(round(median(subset(EVOrganized, study == a & model == "aSND")$EV, na.rm = TRUE), 2))
  print(round(quantile(subset(EVOrganized, study == a & model == "aSND")$EV, na.rm = TRUE, c(.25, .75)), 2))
  print("bCon")
  print(round(median(subset(EVOrganized, study == a & model == "bCon")$EV, na.rm = TRUE), 2))
  print(round(quantile(subset(EVOrganized, study == a & model == "bCon")$EV, na.rm = TRUE,c(.25, .75)), 2))
}

for(a in unique(PmaxOrganized$study)){
  print(a)
  print(round(c(unlist((t.test(
    log10(subset(PmaxOrganized, study == a & model == "aSND")$Pmax),
    log10(subset(PmaxOrganized, study == a & model == "bCon")$Pmax))[c("p.value")])),
    corr.test(log10(subset(PmaxOrganized, study == a & model == "aSND")$Pmax),
              log10(subset(PmaxOrganized, study == a & model == "bCon")$Pmax))$r), 4)
  )
  print("aSND")
  print(round(median(subset(PmaxOrganized, study == a & model == "aSND")$Pmax, na.rm = TRUE), 2))
  print(round(quantile(subset(PmaxOrganized, study == a & model == "aSND")$Pmax, na.rm = TRUE, c(.25, .75)), 2))
  print("bCon")
  print(round(median(subset(PmaxOrganized, study == a & model == "bCon")$Pmax, na.rm = TRUE), 2))
  print(round(quantile(subset(PmaxOrganized, study == a & model == "bCon")$Pmax, na.rm = TRUE,c(.25, .75)), 2))
}


for(a in unique(OmaxOrganized$study)){
  print(a)
  print(round(c(unlist((t.test(
    log10(subset(OmaxOrganized, study == a & model == "aSND")$Omax),
    log10(subset(OmaxOrganized, study == a & model == "bCon")$Omax))[c("p.value")])),
    corr.test(log10(subset(OmaxOrganized, study == a & model == "aSND")$Omax),
              log10(subset(OmaxOrganized, study == a & model == "bCon")$Omax))$r), 4)
  )
  print("aSND")
  print(round(median(subset(OmaxOrganized, study == a & model == "aSND")$Omax, na.rm = TRUE), 2))
  print(round(quantile(subset(OmaxOrganized, study == a & model == "aSND")$Omax, na.rm = TRUE, c(.25, .75)), 2))
  print("bCon")
  print(round(median(subset(OmaxOrganized, study == a & model == "bCon")$Omax, na.rm = TRUE), 2))
  print(round(quantile(subset(OmaxOrganized, study == a & model == "bCon")$Omax, na.rm = TRUE,c(.25, .75)), 2))
}



# Logarithmic Modeling Demonstration ----
winger2006logExp <- function(ind){
  est<- nls_multstart(log(y) ~ log(10^q0) + (k * (exp(-(10^alpha) * (10^q0) * x) - 1)), data = data.frame(ind),
                      iter = c(10, 10),
                      start_lower = c(q0 = 1, alpha = -9),
                      start_upper = c(q0 = 3, alpha = -2))
  return(est)
}


winger2006logSND <- function(ind){
  est<- nls_multstart(log(y) ~ log(10^q0) - (10^alpha) * (10^q0) * x, data = data.frame(ind),
                      iter = c(10, 10),
                      start_lower = c(q0 = 1, alpha = -9),
                      start_upper = c(q0 = 3, alpha = -2))
  return(est)
}



GetlnK <- function(dat, mnrange = TRUE) {
  if (mnrange) {
    dat1 <- aggregate(y ~ x, dat, mean)
    (log(max(dat1$y[dat1$y > 0], na.rm = TRUE)) - log(min(dat1$y[dat1$y > 
                                                                   0], na.rm = TRUE))) + 0.5
  }
  else {
    (log(max(dat$y[dat$y > 0], na.rm = TRUE)) - log(min(dat$y[dat$y > 
                                                                0], na.rm = TRUE))) + 0.5
  }
}
GetlnK(subset(winger2006Data, y != 0))

winger2006logFits <- winger2006Data |> 
  mutate(y = ifelse(y == 0, NA, y), k = GetlnK(subset(winger2006Data, y != 0))) |> 
  group_by(monkey, combDrug) |> 
  nest() |> 
  mutate(logExp = map(data, winger2006logExp),
         logSND = map(data, winger2006logSND))

saveRDS(winger2006logFits, "winger2006logfits.rds")
winger2006logFits <- readRDS("winger2006logfits.rds")


winger2006logFits2 <- winger2006Data |> 
  mutate(y = ifelse(y == 0, NA, y), k = GetlnK(subset(winger2006Data, y != 0))*5) |> 
  group_by(monkey, combDrug) |> 
  nest() |> 
  mutate(logExp = map(data, winger2006logExp),
         logSND = map(data, winger2006logSND))

saveRDS(winger2006logFits2, "winger2006logfits2.rds")
winger2006logFits2 <- readRDS("winger2006logfits2.rds")

winger2006logOutputs <- winger2006logFits |> mutate(logExpCoef = map(logExp, coef),
                                                    logSNDCoef = map(logSND, coef))


winger2006logOutputs2 <- winger2006logFits2 |> mutate(logExpCoef = map(logExp, coef),
                                                      logSNDCoef = map(logSND, coef))


winger2006logAlphaEXPD <- winger2006logOutputs |> unnest(c(logExpCoef)) |> filter(logExpCoef < 0) |> pull(logExpCoef)
winger2006logQ0EXPD <- winger2006logOutputs |> unnest(c(logExpCoef)) |> filter(logExpCoef > 0) |> pull(logExpCoef)
winger2006logAlphaSND <- winger2006logOutputs |> unnest(c(logSNDCoef)) |> filter(logSNDCoef < 0) |> pull(logSNDCoef)
winger2006logQ0SND <- winger2006logOutputs |> unnest(c(logSNDCoef)) |> filter(logSNDCoef > 0) |> pull(logSNDCoef)

winger2006logAlphaEXPD2 <- winger2006logOutputs2 |> unnest(c(logExpCoef)) |> filter(logExpCoef < 0) |> pull(logExpCoef)
winger2006logQ0EXPD2 <- winger2006logOutputs2 |> unnest(c(logExpCoef)) |> filter(logExpCoef > 0) |> pull(logExpCoef)
winger2006logAlphaSND2 <- winger2006logOutputs2 |> unnest(c(logSNDCoef)) |> filter(logSNDCoef < 0) |> pull(logSNDCoef)
winger2006logQ0SND2 <- winger2006logOutputs2 |> unnest(c(logSNDCoef)) |> filter(logSNDCoef > 0) |> pull(logSNDCoef)



svg("logwinger.svg", width = 7.5, height = 5.75)
par(mfrow = c(2,3), mar = c(4.75, 3.75, 3.75, 1.5), pty = "m")

logWingerK1 <- winger2006logOutputs$data[[1]]$k[1]
logWingerK1
plot(winger2006logQ0SND, winger2006logQ0EXPD,
     main = bquote(italic(Q)[0*EXPL]~vs~italic(Q)[0*LSND]~"when"~italic(k)==.(logWingerK1)),
     xlab = expression(log[10](italic(Q)[0*LSND])),
     ylab = expression(log[10](italic(Q)[0*EXPL])), pch = 22, mgp = c(2.5,1,0))
legend("topleft", 
       legend = bquote(italic(r)==.(rd(cor.test(winger2006logQ0EXPD, winger2006logQ0SND)$estimate, 3))
       ), bty = "n"
)

plot(winger2006logAlphaSND, winger2006logAlphaEXPD,
     main = bquote(alpha[EXPL]~vs~alpha[LSND]~when~italic(k)==.(logWingerK1)),
     xlab = expression(log[10](alpha[LSND])),
     ylab = expression(log[10](alpha[EXPL])), mgp = c(2.5,1,0))
legend("topleft", 
       legend = bquote(italic(r)==.(rd(cor.test(winger2006logAlphaEXPD, winger2006logAlphaSND)$estimate, 3))
       ), bty = "n"
)

for(a in 1){
  curve(log((10^winger2006logQ0EXPD)[a]) + (logWingerK1 * (exp(-(10^winger2006logAlphaEXPD)[a] * (10^winger2006logQ0EXPD)[a] * x) - 1)), 
        xlim = c(5, 2000), ylim = c(-2, 6), log = "x", from = .01, to = 10000, col = "orange", axes = FALSE,
        xlab = "FR", ylab = expression(ln("Consumption")), main = paste0("Cocaine, Monkey ", unique(winger2006Data$monkey)[a], ", 'Typical' Span"), 
        lwd = 1.25, mgp = c(2.5,1,0))
  curve(log((10^winger2006logQ0EXPD)[a+3]) + (logWingerK1 * (exp(-(10^winger2006logAlphaEXPD)[a+3] * (10^winger2006logQ0EXPD)[a+3] * x) - 1)), 
        from = .01, to = 10000, col = "orange", add = TRUE, lwd = 1.25)
  curve(log((10^winger2006logQ0EXPD)[a+3*2]) + (logWingerK1 * (exp(-(10^winger2006logAlphaEXPD)[a+3*2] * (10^winger2006logQ0EXPD)[a+3*2] * x) - 1)), 
        from = .01, to = 10000, col = "orange", add = TRUE, lwd = 1.25)
  curve(log(10^winger2006logQ0SND[a]) - 10^winger2006logAlphaSND[a] * 10^winger2006logQ0SND[a] * x, 
        from = .01, to = 10000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(log(10^winger2006logQ0SND[a+3]) - 10^winger2006logAlphaSND[a+3] * 10^winger2006logQ0SND[a+3] * x, 
        from = .01, to = 10000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(log(10^winger2006logQ0SND[a+3*2]) - 10^winger2006logAlphaSND[a+3*2] * 10^winger2006logQ0SND[a+3*2] * x, 
        from = .01, to = 10000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  points(rep(winger2006xVals, 3), log(subset(winger2006Data, monkey == unique(winger2006Data$monkey)[a] & 
                                               drug == unique(winger2006Data$drug)[1])$y), 
         pch =rep(c(21, 22, 23), each = 6)) 
  axis(1, at = tickFunc(jumps = 9), tcl = -.25, labels = NA)
  axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, "10k")))
  axis(2)
}

logWingerK2 <- winger2006logOutputs2$data[[1]]$k[1]
logWingerK2
plot(winger2006logQ0SND2, winger2006logQ0EXPD2,
     main = bquote(italic(Q)[0*EXPL]~vs~italic(Q)[0*LSND]~"when"~italic(k)==.(logWingerK2)),
     xlab = expression(log[10](italic(Q)[0*LSND])),
     ylab = expression(log[10](italic(Q)[0*EXPL])), pch = 22, mgp = c(2.5,1,0))
legend("topleft", 
       legend = bquote(italic(r)==.(rd(cor.test(winger2006logQ0EXPD2, winger2006logQ0SND2)$estimate, 3))
       ), bty = "n"
)
plot(winger2006logAlphaSND2, winger2006logAlphaEXPD2,
     main = bquote(alpha[EXPL]~vs~alpha[LSND]~"when"~italic(k)==.(logWingerK2)),
     xlab = expression(log[10](alpha[LSND])),
     ylab = expression(log[10](alpha[EXPL])), mgp = c(2.5,1,0))
legend("topleft", 
       legend = bquote(italic(r)==.(rd(cor.test(winger2006logAlphaEXPD2, winger2006logAlphaSND2)$estimate, 3))
       ), bty = "n"
)


for(a in 1){
  curve(log((10^winger2006logQ0EXPD2)[a]) + (logWingerK2 * (exp(-(10^winger2006logAlphaEXPD2)[a] * (10^winger2006logQ0EXPD2)[a] * x) - 1)), 
        xlim = c(5, 2000), ylim = c(-2, 6), log = "x", from = .01, to = 10000, col = "orange", axes = FALSE,
        xlab = "FR", ylab = expression(ln("Consumption")), main = paste0("Cocaine, Monkey ", unique(winger2006Data$monkey)[a], ", Span * 5"), 
        lwd = 1.25, mgp = c(2.5,1,0))
  curve(log((10^winger2006logQ0EXPD2)[a+3]) + (logWingerK2 * (exp(-(10^winger2006logAlphaEXPD2)[a+3] * (10^winger2006logQ0EXPD2)[a+3] * x) - 1)), 
        from = .01, to = 10000, col = "orange", add = TRUE, lwd = 1.25)
  curve(log((10^winger2006logQ0EXPD2)[a+3*2]) + (logWingerK2 * (exp(-(10^winger2006logAlphaEXPD2)[a+3*2] * (10^winger2006logQ0EXPD2)[a+3*2] * x) - 1)), 
        from = .01, to = 10000, col = "orange", add = TRUE, lwd = 1.25)
  curve(log(10^winger2006logQ0SND2[a]) - 10^winger2006logAlphaSND2[a] * 10^winger2006logQ0SND2[a] * x, 
        from = .01, to = 10000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(log(10^winger2006logQ0SND2[a+3]) - 10^winger2006logAlphaSND[a+3] * 10^winger2006logQ0SND2[a+3] * x, 
        from = .01, to = 10000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  curve(log(10^winger2006logQ0SND2[a+3*2]) - 10^winger2006logAlphaSND2[a+3*2] * 10^winger2006logQ0SND2[a+3*2] * x, 
        from = .01, to = 10000, col = "blue", add = TRUE, lty = 5, lwd = 1.25)
  points(rep(winger2006xVals, 3), log(subset(winger2006Data, monkey == unique(winger2006Data$monkey)[a] & 
                                               drug == unique(winger2006Data$drug)[1])$y), 
         pch =rep(c(21, 22, 23), each = 6)) 
  axis(1, at = tickFunc(jumps = 9), tcl = -.25, labels = NA)
  axis(1, at = 10^(-3:4), labels = as.character(c(.001, .01, .1, 1, 10, 100, 1000, "10k")))
  axis(2)
}



dev.off()

# Alpha Comparisons and Additional Supplemental ----

anova(ko2002EXPDGNLS, type = "marginal")
format(data.frame(emmeans(ko2002EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)$contrasts),
        digits = 1, nsmall = 3)
anova(ko2002SNDGNLS, type = "marginal")
format(data.frame(emmeans(ko2002SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)$contrasts),
        digits = 1, nsmall = 3)

anova(winger2002EXPDGNLS, type = "marginal")
format(data.frame(emmeans(winger2002EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)$contrasts),
        digits = 1, nsmall = 3)
anova(winger2002SNDGNLS, type = "marginal")
format(data.frame(emmeans(winger2002SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)$contrasts),
        digits = 1, nsmall = 3)

anova(winger2006EXPDGNLS, type = "marginal")
format(data.frame(emmeans(winger2006EXPDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)$contrasts),
        digits = 1, nsmall = 3)
anova(winger2006SNDGNLS, type = "marginal")
format(data.frame(emmeans(winger2006SNDGNLS, param = "alpha", specs = pairwise~aloneDose|drug)$contrasts),
        digits = 1, nsmall = 3)

## Ko 2002 EXPL vs EXPD vs SND ----


# HS EXPL model
ko2002EXPLNLSTestAlf <- nls_multstart(log10(y) ~ log10(10^q0) + k * (exp(-(10^alpha) * (10^q0) * x)-1),
                                    data = subset(ko2002Data, drug == "alfentanil"), iter = c(10, 10),
                                    start_lower = c(q0 = 1, alpha = -5),
                                    start_upper = c(q0 = 3, alpha = 1))

ko2002EXPLGNLSTestAlf1 <- gnls(log10(y) ~ log10(10^q0) + k * (exp(-(10^alpha) * (10^q0) * x)-1),
                             data = subset(ko2002Data, drug == "alfentanil"),
                             params = list(q0 ~ aloneDose,  
                                           alpha ~ 1
                             ),
                             start = list(fixed = c(coef(ko2002EXPLNLSTestAlf)[1], rep(0, 2), coef(ko2002EXPLNLSTestAlf)[2])),
                             verbose = 2,
                             control = list(msMaxIter = 50000,
                                            niterEM = 5000,
                                            maxIter = 5000,
                                            pnlsTol = .01,
                                            tolerance = .01,
                                            apVar = T,
                                            minScale = .0000001,
                                            opt = "optim"), na.action = na.omit)

ko2002EXPLGNLSTestAlfVar <- gnls(log10(y) ~ log10(10^q0) + k * (exp(-(10^alpha) * (10^q0) * x)-1),
                               data = subset(ko2002Data, drug == "alfentanil"),
                               params = list(q0 ~ aloneDose,  
                                             alpha ~ aloneDose
                               ),
                               start = list(fixed = c(coef(ko2002EXPLNLSTestAlf)[1], rep(0, 2), coef(ko2002EXPLNLSTestAlf)[2], rep(0, 2))),
                               verbose = 2,
                               control = list(msMaxIter = 50000,
                                              niterEM = 5000,
                                              maxIter = 5000,
                                              pnlsTol = .01,
                                              tolerance = .01,
                                              apVar = T,
                                              minScale = .0000001,
                                              opt = "optim"), na.action = na.omit)

# Koffarnus EXPD model
ko2002EXPDNLSTestAlf <- nls_multstart(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                                      data = subset(ko2002Data, drug == "alfentanil"), iter = c(10, 10),
                                      start_lower = c(q0 = 1, alpha = -5),
                                      start_upper = c(q0 = 3, alpha = 1))

ko2002EXPDGNLSTestAlf1 <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                               data = subset(ko2002Data, drug == "alfentanil"),
                               params = list(q0 ~ aloneDose,  
                                             alpha ~ 1
                               ),
                               start = list(fixed = c(coef(ko2002EXPDNLSTestAlf)[1], rep(0, 2), coef(ko2002EXPDNLSTestAlf)[2])),
                               verbose = 2,
                               control = list(msMaxIter = 50000,
                                              niterEM = 5000,
                                              maxIter = 5000,
                                              pnlsTol = .01,
                                              tolerance = .01,
                                              apVar = T,
                                              minScale = .0000001,
                                              opt = "optim"), na.action = na.omit)

ko2002EXPDGNLSTestAlfVar <- gnls(y ~ 10^(q0) * 10^(k * (exp(-10^(alpha) * 10^(q0) * x) - 1)),
                                 data = subset(ko2002Data, drug == "alfentanil"),
                                 params = list(q0 ~ aloneDose,  
                                               alpha ~ aloneDose
                                 ),
                                 start = list(fixed = c(coef(ko2002EXPDNLSTestAlf)[1], rep(0, 2), coef(ko2002EXPDNLSTestAlf)[2], rep(0, 2))),
                                 verbose = 2,
                                 control = list(msMaxIter = 50000,
                                                niterEM = 5000,
                                                maxIter = 5000,
                                                pnlsTol = .01,
                                                tolerance = .01,
                                                apVar = T,
                                                minScale = .0000001,
                                                opt = "optim"), na.action = na.omit)


# SND model
ko2002SNDNLSTestAlf <- nls_multstart(y ~ (10^q0) * exp(-(10^alpha) * (10^q0) * x),
                                     data = subset(ko2002Data, drug == "alfentanil"), iter = c(10, 10),
                                     start_lower = c(q0 = 1, alpha = -5),
                                     start_upper = c(q0 = 3, alpha = 1))

# One alpha estimated
ko2002SNDGNLSTestAlf1 <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                              data = subset(ko2002Data, drug == "alfentanil"),
                              params = list(q0 ~ aloneDose,  
                                            alpha ~ 1
                              ),
                              start = list(fixed = c(coef(ko2002SNDNLSTestAlf)[1], rep(0, 2), coef(ko2002SNDNLSTestAlf)[2])),
                              verbose = 2,
                              control = list(msMaxIter = 50000,
                                             niterEM = 5000,
                                             maxIter = 5000,
                                             pnlsTol = .01,
                                             tolerance = .01,
                                             apVar = T,
                                             minScale = .0000001,
                                             opt = "optim"), na.action = na.omit)
# Multiple alphas estimated
ko2002SNDGNLSTestAlfVar <- gnls(y ~ (10^q0) * exp(-(10^alpha) *(10^q0) * x),
                                data = subset(ko2002Data, drug == "alfentanil"),
                                params = list(q0 ~ aloneDose,  
                                              alpha ~ aloneDose
                                ),
                                start = list(fixed = c(coef(ko2002SNDNLSTestAlf)[1], rep(0, 2), coef(ko2002SNDNLSTestAlf)[2], rep(0, 2))),
                                verbose = 2,
                                control = list(msMaxIter = 50000,
                                               niterEM = 5000,
                                               maxIter = 5000,
                                               pnlsTol = .01,
                                               tolerance = .01,
                                               apVar = T,
                                               minScale = .0000001,
                                               opt = "optim"), na.action = na.omit)



anova(ko2002EXPLGNLSTestAlf1, ko2002EXPLGNLSTestAlfVar)
anova(ko2002EXPDGNLSTestAlf1, ko2002EXPDGNLSTestAlfVar)
anova(ko2002SNDGNLSTestAlf1, ko2002SNDGNLSTestAlfVar)





