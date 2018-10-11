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
  curr_tab <- paste(dirname(getwd()),"/current_table.dat", sep="")
} else {
  curr_tab <- paste(getwd(),"/current_table.dat", sep="")
}

curr <- read.table(curr_tab)

# set device driver for interactive plotting, I hope!
X11()

# plot
plot(curr$V1, curr$V2, type="o")
grid()
