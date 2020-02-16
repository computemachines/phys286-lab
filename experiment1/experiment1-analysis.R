library(tidyverse)

indium_timeseries <- read_csv("raw-data/Halflife CRTP.csv",
                              col_names = c("time", "geiger counts"),
                              col_types = cols_only("time"=col_number(), "geiger counts"=col_number()), skip = 2)

indium_background <- read_csv("raw-data/Halflife CRTP.csv",
                           col_names = c(NA, NA, NA, NA, "time", "geiger counts"),
                           col_types = cols_only("time"=col_number(), "geiger counts"=col_integer()), skip = 2)

ib <- indium_background[indium_background %>% complete.cases(),]
ib[,1] <- ib[,1]/30

write.csv(it, "clean-data/indium-halflife")
write.csv(it, "clean-data/indium-background")

names(ib)[1] <- "samples"
plot(ib, main="Background Geiger Counts per 30s Sample")
fit <- lm(ib$`geiger counts` ~ ib$samples)
mean(ib$`geiger counts`)
var(ib$`geiger counts`)
(background_stderr <- var(ib$`geiger counts`)/max(ib$samples))

it <- indium_timeseries
it[,1] <- it[,1]/30
names(it)[1] <- "samples"

theta.0 <- mean(ib$`geiger counts`)
model.0 <- lm(log(`geiger counts` - theta.0) ~ samples, data=it)
coef(model.0)
alpha.0 <- exp(coef(model.0)[1])
tau.0 <- -1/coef(model.0)[2]

start <- list(alpha=alpha.0, theta=theta.0, tau=tau.0)
model <- nls(`geiger counts` ~ alpha * exp(-samples/tau) + theta, data=it, start = start)
model


plot(it, main="Indium Source Geiger Counts per 30s sample over time",
     ylim=c(0, max(it$`geiger counts`)))
lines(it$samples, predict(model, list(time=it$samples)), lwd=3)
model
options(scipen=999)
summary(model)
confint(model)
