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


par(mfcol=c(5,2), mai=c(0.25,0.3,0.1,0.2))

# columna izquierda FIGURA 5.8
trend <- T1
Xf <- X2f
sigma0 <- 1

set.seed(86)
TX <- Xf(t)
Y0 <- function(t) s2(t) + trend(t) + sigma0*TX
TY0t <- Y0(t)

#TBATS en Y0
fit <- tbats(TY0t, seasonal.periods=c(100, 100/pi), use.box.cox=TRUE, use.trend=FALSE) 
components <- tbats.components(fit)
r_hat_TBATBS <- TY0t - components[,'level'] - components[,'season1'] - components[,'season2']

plot(t, TY0t, type="l", xaxt = "n")
plot(t,as.double(components[,'season2']), type = "l", ylab="s11",col="red", xaxt = "n")
lines(t, s21(t), lty=1)
plot(t,as.double(components[,'season1']), type = "l", ylab="s12",col="red", xaxt = "n", ylim = c(-4,4))
lines(t, s22(t), lty=1)
plot(t,as.double(components[,'level']), type = "l", ylab="tendencia",col="red", xaxt = "n")
lines(t, trend(t), lty=1)
plot(t, TX, type="l")
lines(t, r_hat_TBATBS, col="red", lty=1)


# columna derecha FIGURA 5.8
trend <- T2
Xf <- X3f
sigma0 <- 0.5

set.seed(95)
TX <- Xf(t)
Y0 <- function(t) s2(t) + trend(t) + sigma0*TX
TY0t <- Y0(t)

#TBATS en Y0
fit <- tbats(TY0t, seasonal.periods=c(100, 100/pi), use.box.cox=TRUE, use.trend=FALSE) 
components <- tbats.components(fit)
r_hat_TBATBS <- TY0t - components[,'level'] - components[,'season1'] - components[,'season2']

plot(t, TY0t, type="l", xaxt = "n")
plot(t,as.double(components[,'season2']), type = "l", ylab="s11",col="red", xaxt = "n")
lines(t, s21(t), lty=1)
plot(t,as.double(components[,'season1']), type = "l", ylab="s12",col="red", xaxt = "n", ylim = c(-4,4))
lines(t, s22(t), lty=1)
plot(t,as.double(components[,'level']), type = "l", ylab="tendencia",col="red", xaxt = "n")
lines(t, trend(t), lty=1)
plot(t, TX, type="l")
lines(t, r_hat_TBATBS, col="red", lty=1)