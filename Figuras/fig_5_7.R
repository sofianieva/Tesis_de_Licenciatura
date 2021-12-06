library('SynchWave')
library("forecast")
library("fGarch")
library('Metrics')


t <- seq(0,10, (1/100))
t2 <- seq(-4,14, (1/100))

A1 <- function(t) (1+0.1*cos(t))*(atan(1132/87 - 200*t/87))/2 + 2
A2 <- function(t) (3.5*(t<=7.5)+2*(t>7.5))

phi1 <- function(t) (t+0.1*sin(t))
phi2 <- function(t) (3.4*t-0.02*(abs(t)^(2.3)))

s21 <- function(t) (A1(t)*cos(2*pi*phi1(t)))
s22 <- function(t) (A2(t)*cos(2*pi*phi2(t)))
s2 <- function(t) s21(t) + s22(t)

T1 <- function(t) (8*(1/(1+(t/5)^2)+exp(-t/10)))
T2 <- function(t) (2*t+10*exp(-(t-4)^2/6))

Xarma1f <- function(t) arima.sim(list(order=c(1,0,1), ar=c(-.5), ma=c(-.4)), n=length(t), rand.gen = rt, df=4)
sigma <- function(t) (1+.1*cos(pi*t))
X1f <- function(t) (2*sigma(t)*Xarma1f(t))
Xarma2f <- function(t) arima.sim(list(order=c(1,0,1), ar=c(.2), ma=c(-.51)), n=length(t), rand.gen = rt, df=4)
X2f <- function(t) sigma(t)*(4*Xarma1f(t)*(t<=5)+Xarma2f(t)*(t>5)) 
l <- list(omega = c(1), alpha = c(0.2), beta=c(0.2, 0.3))
spec = garchSpec(model = l)
Xgarchf <- function(t) garchSim(spec, n = length(t))
X3f <- function(t) 2*as.numeric(Xgarchf(t))

trend <- T1
Xf <- X2f
sigma0 <- 1

set.seed(6)
SX <- Xf(t2)
wavelet <- "hhhat"
# "bump", "hhhat" y "morlet" funcionan

Y0 <- function(t) s2(t) + trend(t) + sigma0*SX
Y0t <- Y0(t2)
dt <- t[2]-t[1]
nv <- 32
opt <- list(type = wavelet)

# Ridge extraction automatico

# Synchrosqueezed wavelet transform
sstfit <- synsq_cwt_fw(t, Y0t, nv, opt)
# Ridge extraction
lambda <- 1e+04
nw <- 16
imtfit <- curve_ext_multi(sstfit$Tx, log2(sstfit$fs), 2, lambda, nw)
# Reconstruction
curvefit <- curve_ext_recon(sstfit$Tx, sstfit$fs, imtfit$Cs, opt, nw)

#temp <- synsq_filter_pass(sstfit$Wx, sstfit$asc, 0.1, 40);
cwtfit <- cwt_fw(Y0t, opt$type, nv, dt, opt)
Y0_sinT <- cwt_iw(cwtfit$Wx, opt$type, opt)
T_hat <- Y0t - Y0_sinT  
r_hat = Y0t - curvefit[,1] - curvefit[,2] - T_hat

# FIGURA 5.7
par(mfcol=c(5,2),  mai=c(0.25,0.3,0.1,0.2))

plot(t, Y0t[401:1401], type="l", xaxt = "n")
plot(t, curvefit[,2][401:1401], type="l", col="red", xaxt = "n")
lines(t, s21(t), lty=1)
plot(t, curvefit[,1][401:1401], type="l", col="red", xaxt = "n")
lines(t, s22(t), lty=1)
plot(t, trend(t), type="l", xaxt = "n")
lines(t, T_hat[401:1401], col="red", lty=1)
plot(t, SX[401:1401], type="l")
lines(t, r_hat[401:1401], col="red", lty=1)

par(mfrow=c(1,1), mai=c(2,2,1.5,1.5))
image.plot(list(x=t2, y=sstfit$fs, z=t(abs(sstfit$Tx))),
           xlab="Tiempo", ylab="Frecuencia", main="Representación tiempo-frecuencia aplicando SST",
           col=designer.colors(64, c("grey100", "grey50", "grey25", "grey0")), ylim=c(0.5, 5), xlim=c(0, 10))
