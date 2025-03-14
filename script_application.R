
library(unitARMA)

data("humidity_max")

dados = humidity_max

h=7 # forecasting horizon

# separate treaning and test data
n = nrow(dados)-h
y = dados[1:n, "hr", drop=T]
y_new = dados[(n+1):(n+h), "hr", drop=T]
x = dados[1:n, c("radiacionGlobalInst","ffInst"), drop=F]
x_new = dados[(n+1):(n+h), c("radiacionGlobalInst","ffInst"), drop=F]


# fitting models

fit_ULARMA <- ULARMA(ts(y), ar = 1, ma = 1, link = "logit", h = 7, resid = 3,
                 X = x, X_hat = x_new)

fit_KARMA <- KARMA(ts(y), ar = 1, ma = 1, link = "logit", h = 7, resid = 3,
                     X = x, X_hat = x_new)

fit_BARMA <- BARMA(ts(y), ar = 1, ma = 1, link = "logit", h = 7, resid = 3,
                     X = x, X_hat = x_new)


# fitted graphics
a = cbind(y, fit_ULARMA$fitted)
plot.ts(a, plot.type = "s", col=1:2, main="ULARMA fitted")

b = cbind(y, fit_KARMA$fitted)
plot.ts(b, plot.type = "s", col=1:2, main="KARMA fitted")

c = cbind(y, fit_BARMA$fitted)
plot.ts(c, plot.type = "s", col=1:2, main="BARMA fitted")


# forecasting
a = cbind( c(y, rep(NA,h)) , c(rep(NA,n), fit_ULARMA$forecast) )
plot.ts( tail(a,100), plot.type = "s", col=1:2, main="ULARMA fitted")

b = cbind( c(y, rep(NA,h)) , c(rep(NA,n), fit_KARMA$forecast) )
plot.ts( tail(b,100), plot.type = "s", col=1:2, main="KARMA fitted")

c = cbind( c(y, rep(NA,h)) , c(rep(NA,n), fit_BARMA$forecast) )
plot.ts( tail(c,100), plot.type = "s", col=1:2, main="BARMA fitted")



