library('SynchWave')

t <- seq(0,10*pi, (1/30))

s1 <- function(t) sin(t)
s2 <- function(t) sin(t*5)


f <- function(t) s1(t) + s2(t)
ft <- f(t)
dt <- t[2]-t[1]
nv <- 32
opt <- list(type = "mhat")

#plot(t, ft, type="l", xlab='', ylab='')

# transformada wavelet continua 
cwtfit <- cwt_fw(ft, opt$type, nv, dt, opt)
par(mfrow=c(3,2), oma = c(2, 2, 2, 2), 
    mar = c(2, 2, 2, 2), 
    mgp = c(2, 1, 0),   
    xpd = NA)

thresh = 0.6
indices = cwtfit$asc[which(cwtfit$asc <= thresh)]

R1 = Re(cwtfit$Wx)
image.plot(list(x=t, y=cwtfit$asc, z=t(R1)), log="y",
           xlab="", ylab="", main="",
           col=designer.colors(64, c("grey90", "grey65", "grey35", "grey10")),
           ylim=c(10, 0.04), xlim=c(3, 27))
lines(t, rep(thresh, length(t)), type='l', col='red')

fcwt1 <- cwt_iw(R1/2, opt$type, opt);
plot(t, fcwt1, type="l", xlab='', ylab='')


R2 = Re(cwtfit$Wx)
R2[1:length(indices), 1:943] <- 0.0
image.plot(list(x=t, y=cwtfit$asc, z=t(R2)), log="y",
           xlab="", ylab="", main="",
           col=designer.colors(64, c("grey90", "grey65", "grey35", "grey10")),
           ylim=c(10, 0.04), xlim=c(3, 27))
lines(t, rep(thresh, length(t)), type='l', col='red')

fcwt2 <- cwt_iw(R2/2, opt$type, opt);
plot(t, fcwt2, type="l", xlab='', ylab='', xlim=c(3, 27))
lines(t, s1(t), lty=3)


indices_h = cwtfit$asc[which(cwtfit$asc > thresh)]
R3 = Re(cwtfit$Wx)
R3[length(indices):320, 1:943] <- 0.0
image.plot(list(x=t, y=cwtfit$asc, z=t(R3)), log="y",
           xlab="", ylab="", main="",
           col=designer.colors(64, c("grey90", "grey65", "grey35", "grey10")),
           ylim=c(10, 0.04), xlim=c(3, 27))
lines(t, rep(thresh, length(t)), type='l', col='red')

fcwt3 <- cwt_iw(R3/2, opt$type, opt);
plot(t, fcwt3, type="l", xlab='', ylab='', xlim=c(3, 27))
lines(t, s2(t), lty=3)
