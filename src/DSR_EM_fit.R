#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(readr)
library(mixtools)
# Arguments
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file.default one columne, contain DSR value"),
  make_option(c("-n", "--num"), dest = "num", default = "2",
              help="number of mixed normal distributions. default 2"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "test.pdf",
              help = "[opt] output file name.")
)

parser <- OptionParser(usage = "mapdrawer [options]", option_list = option_list)
arguments <- parse_args(parser, positional_arguments=c(0,Inf))
num <- as.numeric(arguments$options$num)
infile <- arguments$options$infile
dt_all <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)

mixmdla_2 <- normalmixEM(log(dt_all[[1]]+10,10),k=num,maxit = 30) ## may need to use mean.constr parameter to constrain EM fitting analysis
mu = c()
sigma = c()
lambd <- c()
for (i in seq(1,num)) {
  mu = c(mu,mixmdla_2$mu[i])
  sigma = c(sigma,mixmdla_2$sigma[i])
  lambd <- c(lambd,mixmdla_2$lambda[i])
}

turn_label_log <- function(x) {paste0("10e-", 6-x)}
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
p1 <- data.frame(x = mixmdla_2$x) %>% 
  ggplot() +
  geom_histogram(aes(x, ..density..), bins=50, colour = "black", fill = "grey90") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdla_2$mu[1], mixmdla_2$sigma[1], lam = mixmdla_2$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdla_2$mu[2], mixmdla_2$sigma[2], lam = mixmdla_2$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdla_2$mu[3], mixmdla_2$sigma[3], lam = mixmdla_2$lambda[3]),
                colour = "green", lwd = 1.5)+
  ylab("Density") + xlab("Variant density (per bp)") +
  cowplot::theme_cowplot()
ggsave(outfile,p1, height = 5, width =6, limitsize = FALSE)


