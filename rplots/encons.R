#!/usr/local/bin/Rint

# See: https://thesquareplanet.com/blog/interactive-r-scripts/
library(grDevices)
library(utils)
library(stats)
library(graphics)
library(datasets)
library(methods)
library(base)

# I make a path conveniently, so that I can run this script from
# the rplots folder inside my simulation dir or directly from the simulation dir.
if ( basename(getwd()) == "rplots" ) {
  pluto_data <- paste(dirname(getwd()),"/out/energy_cons.dat", sep="")
} else {
  pluto_data <- paste(getwd(),"/out/energy_cons.dat", sep="")
}

pl <- read.table(pluto_data)
# pl_vars <- c(1, 2, 3)
Eerr_rel <- (pl$Etot-(pl$Etot[1] + pl$E_adv_in + pl$E_tc_in + pl$E_res_in))/pl$Etot

# set device driver for interactive plotting, I hope!
X11()

# plot
par(mfrow=c(3,2)) # all plots on one page

plot(pl$t, pl$Etot, type="o")
grid()
plot(pl$t, pl$E_tc_in, type="o")
grid()
plot(pl$t, pl$E_res_in, type="o")
grid()
plot(pl$t, pl$E_adv_in, type="o")
grid()
plot(pl$t, Eerr_rel, type="o")
grid()