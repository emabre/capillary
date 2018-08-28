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

# set device driver for interactive plotting, I hope!
X11()

# plot
par(mfrow=c(3,1)) # all plots on one page

plot(pl$t, pl$dt, type="o")
grid()
plot(pl$t, type="o")
grid()
plot(pl$dt, type="o")
grid()