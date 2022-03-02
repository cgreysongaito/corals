rm(list=ls())

library(QPot)
library(zeallot)
library(tidyverse)
library(viridis)

setwd("/home/chrisgg/Documents/Guelph/McCannLab/corals/")

## Set up isoclines
resiso_data_05 <- read_csv("data/resiso_data_05.csv")
resiso_data_07 <- read_csv("data/resiso_data_07.csv")

#Model Set Up
coral_model <- "x * ( r * (1 - c * N) * (1 - y - x) - m - a * y * ( N / (n + N) )) "
macro_model <- "y * ( a * x * ( N / (n + N) ) - g / (y + (1 - y - x)) + l * (1 - y - x))"

#General Functions
param_model <- function(param) {
  c(graze, nut) %<-% param
  model_parms <- c(a = 2.3, g = graze, l = 0.7, r = 1.8, m = 0.15, n = 0.5, N = nut, c = 0.25)
  coral_model_param <- Model2String(coral_model, parms = model_parms, supress.print = TRUE)
  macro_model_param <- Model2String(macro_model, parms = model_parms, supress.print = TRUE)
  parms_eqn <- list(coral_model_param, macro_model_param)
  return(parms_eqn)
}

# Set bounds and step numbers for QPotential calculation
bounds.x <- c(-0.05, 1)
bounds.y <- c(-0.05, 1)
step.number.x <- 1000
step.number.y <- 1000

param_g04_nut05 <- list(graze = 0.4, nut = 0.5)


model.state <- c(x = 0.2, y = 0.4)
model.sigma <- 0.005
model.time <- 1000
# we used 12500 in the figures
model.deltat <- 0.025
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat,
                 x.rhs = param_model(param_g04_nut05)[[1]], y.rhs = param_model(param_g04_nut05)[[2]], sigma = model.sigma)

TSPlot(ts.ex1, deltat=model.deltat, ylim=c(0.0, 1.0))
TSPlot(ts.ex1, deltat=model.deltat, dim = 2)
TSDensity(ts.ex1, dim=1)
TSDensity(ts.ex1, dim=2)

#Quasi-potential for default parameters (in Carling code)

eq1.x <- 0.9047619047619048
eq1.y <- 0.0
eq2.x <- 0.1801318204669627
eq2.y <- 0.41882288305961357

eq1_qp_g04_nut05 <- QPotential(x.rhs = param_model(param_g04_nut05)[[1]],
                                      x.start = eq1.x,
                                      x.bound = bounds.x,
                                      x.num.steps = step.number.x,
                                      y.rhs = param_model(param_g04_nut05)[[2]],
                                      y.start = eq1.y,
                                      y.bound = bounds.y,
                                      y.num.steps = step.number.y)

eq2_qp_g04_nut05 <- QPotential(x.rhs = param_model(param_g04_nut05)[[1]],
                               x.start = eq2.x,
                               x.bound = bounds.x,
                               x.num.steps = step.number.x,
                               y.rhs = param_model(param_g04_nut05)[[2]],
                               y.start = eq2.y,
                               y.bound = bounds.y,
                               y.num.steps = step.number.y)

qp_g04_nut05_global <- QPGlobal(local.surfaces = list(eq1_qp_g04_nut05, eq2_qp_g04_nut05),
                       unstable.eq.x = c(0.4291005457635031), unstable.eq.y = c(0.27492353780416257),
                       x.bound = bounds.x, y.bound = bounds.y)

QPContour(surface = eq1_qp_g04_nut05, dens = c(1000, 1000), x.bound = bounds.x,
          y.bound = bounds.y, c.parm = 5) 
QPContour(surface = eq2_qp_g04_nut05, dens = c(1000, 1000), x.bound = bounds.x,
          y.bound = bounds.y, c.parm = 5) 


#Quasi-potential for Nt=0.45

eq1.x <- 0.9061032863849767
eq1.y <- 0.0

param_nut045 <- list(graze = 0.4, nut = 0.45)
eq1_qp_nut045 <- QPotential(x.rhs = param_model(param_nut045)[[1]],
                               x.start = eq1.x,
                               x.bound = bounds.x,
                               x.num.steps = step.number.x,
                               y.rhs = param_model(param_nut045)[[2]],
                               y.start = eq1.y,
                               y.bound = bounds.y,
                               y.num.steps = step.number.y)


QPContour(surface = eq1_qp_nut045, dens = c(1000, 1000), x.bound = bounds.x,
          y.bound = bounds.y, c.parm = 5) 

#Quasi-potential for default parameters (in Carling code)

eq1.x <- 0.898989898989899
eq1.y <- 0.0
eq2.x <- 0.052093574708453154
eq2.y <- 0.4449201849700071

param_nut07 <- list(graze = 0.4, nut = 0.7)

eq1_qp_g04_nut07 <- QPotential(x.rhs = param_model(param_nut07)[[1]],
                               x.start = eq1.x,
                               x.bound = bounds.x,
                               x.num.steps = step.number.x,
                               y.rhs = param_model(param_nut07)[[2]],
                               y.start = eq1.y,
                               y.bound = bounds.y,
                               y.num.steps = step.number.y)

eq2_qp_g04_nut07 <- QPotential(x.rhs = param_model(param_nut07)[[1]],
                               x.start = eq2.x,
                               x.bound = bounds.x,
                               x.num.steps = step.number.x,
                               y.rhs = param_model(param_nut07)[[2]],
                               y.start = eq2.y,
                               y.bound = bounds.y,
                               y.num.steps = step.number.y)

qp_g04_nut07_global <- QPGlobal(local.surfaces = list(eq1_qp_g04_nut07, eq2_qp_g04_nut07),
                                unstable.eq.x = c(0.0, 0.5819530444904251), unstable.eq.y = c(0.0, 0.16655650787678872),
                                x.bound = bounds.x, y.bound = bounds.y)

QPContour(surface = eq1_qp_g04_nut07, dens = c(1000, 1000), x.bound = bounds.x,
          y.bound = bounds.y, c.parm = 5) 
QPContour(surface = eq2_qp_g04_nut07, dens = c(1000, 1000), x.bound = bounds.x,
          y.bound = bounds.y, c.parm = 5) 

