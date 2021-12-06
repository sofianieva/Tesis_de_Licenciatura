library('SynchWave')

t <- seq(0,10*pi, (1/30))

s1 <- function(t) sin(t)
s2 <- function(t) sin(t*5)

f <- function(t) s1(t) + s2(t)
ft <- f(t)
dt <- t[2]-t[1]
nv <- 32
opt <- list(type = "morlet")

# transformada synchrosqueezing
sstfit <- synsq_cwt_fw(t, ft, nv, opt)

# Set plot layout
par(mai=c(0.5,1,0.5,1))
layout(mat = matrix(c(1, 2), 
                    nrow = 2, 
                    ncol = 1),
       heights = c(1, 2),    # Heights of the two rows
       widths = c(3))     # Widths of the two columns

plot(t, ft, type="l", xlab='', ylab='')

image.plot(list(x=t, y=sstfit$asc, z=t(abs(sstfit$Tx))), log="y",
           xlab="", ylab="", main="",
           col=designer.colors(64, c("grey0", "grey50", "grey75", "grey100")),
           ylim=c(0.1, 5), xlim=c(8, 22))

# La imagen presente en la tesis fue posteriormente retocada en paint 
# Los otros dos gráficos son los de la figura 3.13