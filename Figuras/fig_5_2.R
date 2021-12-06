library('SynchWave')

# Definimos las los componentes y el dominio
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

f <- function(t) s2(t) + T1(t)

# Definimos los parámetros de Synchrosqueezing
set.seed(7)
ft <- f(t2)
dt <- t[2]-t[1]
nv <- 32
wavelet <- "hhhat"
opt <- list(type = wavelet)

# Señal extendida
# Synchrosqueezed wavelet transform
sstfit <- synsq_cwt_fw(t2, ft, nv, opt)

# Extraccion automática de curvas
lambda <- 1e+04
nw <- 16
nc <- 2
imtfit <- curve_ext_multi(sstfit$Tx, log2(sstfit$fs), nc, lambda, nw)

# Reconstruction componentes estacionales
curvefit <- curve_ext_recon(sstfit$Tx, sstfit$fs, imtfit$Cs, opt, nw)
s21_hat <- curvefit[,2][401:1401]
s22_hat <- curvefit[,1][401:1401]

# Reconstruction tendencia
cwtfit <- cwt_fw(ft, opt$type, nv, dt, opt)
#opt$gamma <- est_riskshrink_thresh(cwtfit$Wx, nv)
ft_sinT1 <- cwt_iw(cwtfit$Wx, opt$type, opt)
T1_hat <- ft - ft_sinT1  
T1_hat <- T1_hat[401:1401]


# Señal no extendida
# Synchrosqueezed wavelet transform
sstfit2 <- synsq_cwt_fw(t, f(t), nv, opt)

# Extraccion automática de curvas
lambda <- 1e+04
nw <- 16
nc <- 2
imtfit2 <- curve_ext_multi(sstfi2t$Tx, log2(sstfit2$fs), nc, lambda, nw)

# Reconstruction componentes estacionales
curvefit <- curve_ext_recon(sstfit2$Tx, sstfit2$fs, imtfit2$Cs, opt, nw)
s21_hat <- curvefit[,2]
s22_hat <- curvefit[,1]

# Reconstruction tendencia
cwtfit2 <- cwt_fw(f(t), opt$type, nv, dt, opt)
ft_sinT12 <- cwt_iw(cwtfit2$Wx, opt$type, opt)
T1_hat2 <- f(t) - ft_sinT12  


# FIGURA 5.2
par(mfrow=c(1,2), mai=c(2.2,0.5,2.2,0.2))
plot(t, T1(t), type="l", ylim=c(4,17), xlab="")
lines(t, T1_hat2, col="red", lty=1)
plot(t, T1(t), type="l", ylim=c(4,17), xlab="")
lines(t, T1_hat, col="red", lty=1)