## Ref: http://blog.liang2.tw/2013-R-Statistics/ossf.html?full#story
library(ggplot2)

bus.gamma <- function(t) dgamma(t, 15, 1)
curve(bus.gamma, from = 0, to = 50)

ggplot(data.frame(bus.time = c(0, 50)), aes(bus.time)) +
    stat_function(fun=bus.gamma) +
    labs(x="Bus Arrival Time", y="") + theme_bw()


## Multiple distributions on same plot
bus.gamma.1 <- function(t) dgamma(t, 7.5, 0.5)
bus.gamma.2 <- function(t) dgamma(t, 15, 1)
bus.gamma.3 <- function(t) dgamma(t, 30, 2)
curve(bus.gamma.3, from=0, to=50, col='green')
curve(bus.gamma.2, from=0, to=50, add=TRUE, col='blue')
curve(bus.gamma.1, from=0, to=50, add=TRUE, col='red')

## Using ggplot2
plot.fn1 <- data.frame(x=c(0, 50), Functions=factor(1))
plot.fn2 <- data.frame(x=c(0, 50), Functions=factor(2))
plot.fn3 <- data.frame(x=c(0, 50), Functions=factor(3))

g <- ggplot(NULL, aes(x=x, colour=Functions)) +
    stat_function(data=plot.fn1, fun=dgamma,
                  args=list(shape=7.5, rate=0.5)) +
    stat_function(data=plot.fn2, fun=dgamma,
                  args=list(shape=15, rate=1)) +
    stat_function(data=plot.fn3, fun=dgamma,
                  args=list(shape=30, rate=2))

g + labs(x="Bus Arrival Time", color="Parameters", y="") +
    scale_color_brewer(palette='Set1', labels=c(
        "a=7.5, b=0.5", "a=15, b=1", "a=30, b=2"
    )) + theme_bw(20) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))


# Number of simulations for each girl
EXP_REP <- 5000

early_arrival_p <- function(girl_n, th) {
    # If any t_i in 1 ... girl_n such that t_i < th,
    # it means at least one bus has arrived in th minutes.
    # We mark this success bus arrival as TRUE event.
    # By repeating EXP_REP times, the rate having TRUE event
    # should approximate to the real possibility.
    exp.sim <- replicate(
        EXP_REP,
        any(rgamma(girl_n, 15, 1) <= th)
    )
    mean(exp.sim)  # return the rate having TRUE event
}

least.prob <- 0.8

# First we guess that N ranges from 1 to 30
which.max(
    sapply(1:30, function(n) early_arrival_p(n, 10)) > least.prob
)

# Theoretically
ceiling(log(1 - least.prob) / log(1 - pgamma(10, 15, 1)))
