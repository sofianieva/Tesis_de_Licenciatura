# Definamos las siguientes funciones para modelar la estacionalidad
# Primero hay que definir el dominio
t <- seq(0,10, (1/100))
# estos valores fueron sacados de Chen pag 18 y 19

s11 <- function(t) (2.5*cos(2*pi*t))
s12 <- function(t) (3*cos(2*pi^2*t))
s1 <- function(t) s11(t) + s12(t)

A1 <- function(t) (1+0.1*cos(t))*rev((atan((1:length((t)))/43.5-10)))/2 + 2
#A1 <- function(t) (1+0.1*cos(t))*(atan(1132/87 - 200*t/87))/2 + 2
A2 <- function(t) (3.5*(t<=7.5)+2*(t>7.5))

phi1 <- function(t) (t+0.1*sin(t))
phi2 <- function(t) (3.4*t-0.02*(t^(2.3)))
#para graficar voy a necesitar las derivadas
dphi1 <- function(t) (1+0.1*cos(t))
dphi2 <- function(t) (3.4-0.046*t^(1.3))

s21 <- function(t) (A1(t)*cos(2*pi*phi1(t)))
s22 <- function(t) (A2(t)*cos(2*pi*phi2(t)))
s2 <- function(t) s21(t) + s22(t)

# Ahora consideremos las siguientes dos funciones de tendencia

T1 <- function(t) (8*(1/(1+(t/5)^2)+exp(-t/10)))
T2 <- function(t) (2*t+10*exp(-(t-4)^2/6))

#GRAFICOS
par(mfrow=c(6,2), mai=c(0.3,0.55,0.1,0.3))
plot(t,s11(t), type = "l", ylim=c(-4, 4))
plot(t,s12(t), type = "l", ylim=c(-4, 4))
plot(t,A1(t), type = "l", ylim=c(0, 4))
plot(t,A2(t), type = "l", ylim=c(0, 4))
plot(t,dphi1(t), type = "l", ylim=c(0, 4))
plot(t,dphi2(t), type = "l", ylim=c(0, 4))
plot(t,s21(t), type = "l", ylim=c(-4, 4))
plot(t,s22(t), type = "l", ylim=c(-4, 4))
plot(t,T1(t), type = "l", ylim=c(0,20))
plot(t,T2(t), type = "l", ylim=c(0,20))
plot(t,s2(t)+T1(t), type = "l", ylim=c(-5,25))
plot(t,s1(t)+T2(t), type = "l", ylim=c(-5,25))
par(mfrow = c(1, 1))