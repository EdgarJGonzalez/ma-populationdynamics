# Population dynamics of Mimosa acantholoba
# open r terminal in desired work directory or substitute dir by directory path
dir <- getwd()
setwd(dir) # setting work directory
# load requiered libraries
library(AICcmodavg)
library(betareg)
library(EnvStats)
library(gamm4)
library(ggplot2)
library(lme4)
library(mgcv)
library(reshape2)
library(plotly)
library(gridExtra)
library(plot3D)
# read data
{
  s.db <- read.csv("./Data/survival.csv") # survival with estimated h1
  s1.db <- read.csv("./Data/survival-1.csv") # survival without estimated h1
  g.db <- read.csv("./Data/growth.csv")
  f1.db <- read.csv("./Data/reproduction-probability.csv")
  f2.db <- read.csv("./Data/reproduction-fruits.csv")
  f3.db <- read.csv("./Data/reproduction-seeds.csv")
  f4.db <- read.csv("./Data/establishment.csv")
  f5.db <- read.csv("./Data/understory.csv")
 }
# definitions
{
 Age.mature <- 100
 m <- 100 # small-matrix dimension
 M <- 200000 # big-matrix dimension
 min.lh1 <- min(s.db$ln.h1, na.rm = TRUE)-0.0001
 min.lh2 <- min(s.db$ln.h1, na.rm = TRUE)-0.0001
 max.lh1 <- max(s.db$ln.h1, na.rm = TRUE)+0.0001
 max.lh2 <- max(s.db$ln.h2, na.rm = TRUE)+0.0001
 E.pred <- seq(min(min.lh1, min.lh2), max(max.lh1, max.lh2), length.out = M+1)
 h <- E.pred[2]-E.pred[1]
 X.pred <- (E.pred[2:(M+1)]+E.pred[1:M])/2
 t1.pred <- rep(0, length.out = M)
 sM <- M/m
 e.pred <- seq(min(X.pred), max(X.pred), length.out = m+1)
 x.pred <- (e.pred[2:(m+1)]+e.pred[1:m])/2
 f4.pred <- Age.pred <- 1:Age.mature
 nplot <- nlevels(s.db$PLOT)
}
# modelling vital rates
{
 glmm.s <- glmer(sup~ln.h1*Age+(0+ln.h1| Census)+(0+Age+ln.h1|plot), s.db, binomial) # survival: s
 lmm.g <- lmer(ln.h2~ln.h1*Age+(1+ln.h1| Census)+(1+ln.h1|plot), g.db) # growth: g
 glmm.f1 <- glmer(Rep~Age*ln.h1 + (1| plot), f1.db, binomial) # reproduction probability: f1
 gam.f2 <- gamm4(TotFrut~t2(ln.h1, Age, k = 4), random = ~ (0 + Age + ln.h1|plot), family = negbin(1.179556), data = f2.db) # fruit number per individual: f2
 glmm.f3 <- glmer(N.seed~ln.h1+Age + (1| plot), f3.db, poisson) # seed number per fruit: f3
 den.f5 <- density(x = f5.db$h1, n = m, na.rm = TRUE, from = min(f5.db$h1, na.rm = TRUE), to = max(exp(max.lh1), exp(max.lh2)))
 den.f5$y <- den.f5$y/sum(den.f5$y)
 # gam.f4
 {
	und.ages <- f4.db$Age <- as.factor(as.character(f4.db$Age)) # get ages with register in understory data
	stb.n <- as.data.frame(table(f4.db$Age)) # get number of recruits per year
	names(stb.n) <- c("Age", "SNB") # SNB = recruits
	und.ages <- as.numeric(levels(und.ages))
	new.stb.n <- data.frame(Age = und.ages[order(und.ages)], SNB = rep(NA, length(und.ages))) # create new dataframe
	for (i in 1:nrow(stb.n)) { # order dataframe
		if (stb.n$Age[i] == new.stb.n$Age[which(new.stb.n$Age == stb.n$Age[i])]) {
			new.stb.n$SNB[which(new.stb.n$Age == stb.n$Age[i])] <- stb.n$SNB[i]
		}
	}
	new.stb.n$SNB[which(is.na(new.stb.n$SNB) == TRUE)] <- 0 # setting 0 to ages with register but no recruits
	stb.n <- new.stb.n
	stb.n <- transform(stb.n, Age = as.numeric(as.character(Age)))
	n.ages <- stb.n$Age
	f4.a <- rep(NA, length(n.ages)) # create new dataframe for establishment probability
	for (a in n.ages) {
		stb.a <- subset(s.db, Age == a) # create new data for each year
		n.str.a <- as.data.frame(cbind(x.pred, hist(log(as.numeric(as.character(stb.a$h1))), breaks = e.pred, plot = FALSE)$counts))
		names(n.str.a) <- c("ln.h1", "n")
		n.vec.a <- n.str.a$ln.h1
		pred.f1.a <- predict(glmm.f1, newdata = data.frame(ln.h1 = x.pred, Age = rep(a-1, length(x.pred))), type = "response", re.form = NA)
		pred.f2.a <- predict(gam.f2$gam, newdata = data.frame(ln.h1 = x.pred, Age = rep(a-1, length(x.pred))), type = "response", re.form = NA)
		pred.f3.a <- predict(glmm.f3, newdata = data.frame(ln.h1 = x.pred, Age = rep(a-1, length(x.pred))), type = "response", re.form = NA)
		pred.f.a <- pred.f1.a*pred.f2.a*pred.f3.a # predicting and multiplying reproduction functions
		new.n.a <- sum(pred.f.a) # seed number produced in t-1
		stb.n.a <- stb.n$SNB[which(stb.n$Age == a)] # observed recruits number
		f4.a[which(stb.n$Age == a)] <- stb.n.a/new.n.a # establishment probability from t-1 to t
		}
	f4.a[which(f4.a == 0)] <- 0.00000001 # change 0 to perform beta regression
	x <- n.ages
	gam.f4 <- gam(f4.a~s(x, k = 3), family = betar(link = "logit"))
}
# discretizing vital rates functions
{
 # g array
 GGG <- matrix(NA, ncol = Age.mature, nrow = M)
 for (a in 1:Age.mature)
  GGG[, a] <- predict(lmm.g, newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), type = "response", re.form = NA)
 MUG <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {
  init.g <- 1
  for (i in 1:m) {
     MUG[i, a] <- mean(GGG[init.g:(init.g+sM-1), a])
   init.g <- init.g + sM
  }
 }
 sdg <- sigma(lmm.g)
 G <- array(NA, dim = c(Age.mature, m, m))
 for (a in 1:Age.mature) {
  for (i in 1:m) {
   G[a, 1, i] <- pnorm(e.pred[2], MUG[i, a], sdg)
   for (j in 2:(m-1))
    G[a, j, i] <- pnorm(e.pred[j+1], MUG[i, a], sdg) - pnorm(e.pred[j], MUG[i, a], sdg)
   G[a, m, i] <- 1-pnorm(e.pred[m], MUG[i, a], sdg) # avoiding probability eviction
  }
 } 
 # s matrix
 SSS <- matrix(NA, ncol = Age.mature, nrow = M)
 for (a in 1:Age.mature)
  SSS[ ,a] <- predict(glmm.s, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
 MUS <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:(Age.mature)) {
  init.s <- 1
  for (i in 1:m) {
     MUS[i, a] <- mean(SSS[init.s:(init.s+sM-1), a])
     init.s <- init.s + sM
  }
 }
 S <- MUS
 # f1 matrix
 FF1 <- matrix (NA, nrow = M, ncol = Age.mature)
 for (a in 1:Age.mature)
  FF1[ ,a] <- predict(glmm.f1, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
 MUF1 <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {
    init.f1 <- 1
    for (i in 1:m) {
     MUF1[i, a] <- mean(FF1[init.f1:(init.f1+sM-1), a])
     init.f1 <- init.f1 + sM
  }
 }
 F1 <- MUF1
 # f2 matrix
 FF2 <- matrix (NA, nrow = M, ncol = Age.mature)
 for (a in 1:Age.mature)
  FF2[ ,a] <- predict(gam.f2$gam, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
 MUF2 <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {
    init.f2 <- 1
    for (i in 1:m) {
     MUF2[i, a] <- mean(FF2[init.f2:(init.f2+sM-1), a])
     init.f2 <- init.f2 + sM
  }
 }
 F2 <- MUF2
 # f3 matrix
 FF3 <- matrix (NA, nrow = M, ncol = Age.mature)
 for (a in 1:Age.mature)
  FF3[ ,a] <- predict(glmm.f3, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
 MUF3 <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {
    init.f3 <- 1
    for (i in 1:m) {
     MUF3[i, a] <- mean(FF3[init.f3:(init.f3+sM-1), a])
     init.f3 <- init.f3 + sM
  }
 }
 F3 <- MUF3
 # f4 vector
 F4 <- predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response")
 # f5 vector
 F5 <- den.f5$y
}

# creating kernel
{
 k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m)) 
 s.i.a <- S
 f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
 for (a in 1:Age.mature)
  f.i.a[, a] <- F1[, a]*F2[, a]*F3[, a]*F4[a]
 for (a in 1:Age.mature)
  for (i in 1:m) {
   s.i.j.a[a, i, ] <- s.i.a[i, a]
   g.i.j.a[a, , ] <- G[a, , ]
   p.i.j.a[a, i, ] <- p.i <- s.i.j.a[a, i, ]*G[a, i, ]
   f.i.j.a[a, i, ] <- f.i.a[i, a]
   F5.i.j.a[a, i, ] <- F5
   f.i <- f.i.a[i, a]*F5
   f.i.j.a[a, i, ] <- f.i.j.a[a, i, ]*F5
   k.i.j.a[a, i, ] <- p.i+f.i
  }
}
# total predicted lambda
{
	lam.list <- c()
	n.0 <- which(s1.db$Age == 0) # get trees in first year (row number)
	n.0.h <- log(s1.db[n.0,]$h2) # get hight ... since this dataframe has a new estimated h1, h2 is the observed first height for the first year
	n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts # count number of trees in each size class
	init.n.a.v <- n.0.v # renameF vector
	size.v.a.NM <- list(NA) # create list for size structure change, NM = no migration
	n.list <- c(sum(init.n.a.v)) # create population size vector & setting first value
	for (a in 1:Age.mature) {
		 n.a.v <- k.i.j.a[a,,]%*%init.n.a.v # get kernel for first value (year)
		 lam.a <- sum(n.a.v)/sum(init.n.a.v) # get lambda for first year
		 size.v.a.NM[[a]] <- init.n.a.v 
		 init.n.a.v <- n.a.v # rename initial vector to start again
		 n.list <- c(n.list, sum(n.a.v)) # add next population size value to vector
		 lam.list <- c(lam.list, lam.a) # add next lambda value to vector
		 }
		tot.lam.pred <- lam.list # lam.list is the transitory lambda vector
		lam.list.NM <- lam.list # set lambda list with no migration
} 
# total observed lambda
{
   	   ob.lam.list <- c()				# create observed lambda vector
 	   N.0 <- subset(s1.db, Age == 0)		# get row number of trees in first year
 	   N.0 <- subset(N.0, sup == 1)			# filtering to obtain only alive trees
 	   N.init <- nrow(N.0)				# get number of alive trees in first year
 	   N.list <- c(N.init)				# create a new list of observed population sizes and setting first value
 	   for (i in 1:Age.mature) {			
 		   ns.i <- subset(s1.db, Age == i)	# create a set of trees in year a, starting from a = 1
 		   ns.i <- subset(ns.i, sup == 1)	# filtering to get only alive trees
 		   N.ns.i <- nrow(ns.i)			# counting trees
 		   if (N.ns.i != 0 & N.init != 0) {	# if it is not an empty vector:
 			   lambda <- N.ns.i/N.init	# get lambda
 			 } else {			# if it is:
 				 lambda <- NA		# set lambda to NA (avoiding to get 0 instead of NA)
 			 }
 		   N.init <- N.ns.i			# set a new initial population for next lambda
 		   N.list <- c(N.list, N.init)		# add new population size to vector
 		   ob.lam.list <- c(ob.lam.list, lambda)# aff new lambda to observed lambda vector
 	    }
     ns.mature <- subset(s1.db, Age == Age.mature) 	# get mature forest subset (avoid population overestimating due to individuals repeat at Age.mature)
     ns.mature <- subset(ns.mature, sup == 1)		# get alive trees only
     ns.mature <- transform(ns.mature, Census = as.factor(as.character(Census))) # transform census year as factor
     ns.mature <- droplevels(ns.mature) # drop useless levels
     N.mat.list <- c()  # create new vector of population sizes
     for (i in levels(ns.mature$Census)) {		# for each year in mature forest:
     	ns.mat.yr <- subset(ns.mature, Census == i)	# get a subset of population in census year = i
     	ns.mat.N <- nrow(ns.mat.yr)			# count trees
     	N.mat.list <- c(N.mat.list, ns.mat.N)		# add population size to vector
     }
     lam.m.list <- c()					# create a new vector of lambdas of mature forest
     init.n.mat <- N.mat.list[1]			# set first population size
     for (j in 2:length(N.mat.list)) {			# for every j in population size list:
     	lam.j <- N.mat.list[j]/init.n.mat		# get lambda
     	lam.m.list <- c(lam.m.list,lam.j)		# add lambda to list
     	init.n.mat <- N.mat.list[j]			# reset initial population size
     }
     lam.mature <- geoMean(lam.m.list)			# apply geometric mean to get mature forest lambda
     ob.lam.list <- c(ob.lam.list, lam.mature)		# add observed lambda in mature forest to the whole succession list
     ob.lam.df <- as.data.frame(matrix(ncol = 2, nrow = Age.mature+1)) # create database for observed lambdas
     names(ob.lam.df) <- c("Age", "ob.lambda")		# naming database
     ob.lam.df$Age <- 1:(Age.mature+1)			# set age range
     ob.lam.df$ob.lambda[1:(Age.mature+1)] <- ob.lam.list # set observed lambda value
     lambda.df <- as.data.frame(list(lambda = lam.list, Age = 1:Age.mature)) # create lambda data frame for plotting
}
# predicted lambda by plot
{   
     load("./Data/k-p-list.RData")						# read list of kernels calculated by plot: poner link del código?
     s1.db <- transform(s1.db, PLOT = as.factor(as.character(plot)))		# transform plot variable to factor
     nplot <- nlevels(s1.db$PLOT)						# get number of plots
     est.lam.df <- as.data.frame(matrix(nrow = Age.mature, ncol = nplot, NA))	# create new database for estimated lambda values
     names(est.lam.df) <- 1:nplot						# plot names
     for (i in 1:nplot) {								# for each plot:
       	 Mim.s.p <- subset(s1.db, PLOT == i & !is.na(ln.h2))			# get a subset of plot i omitting NA's
         if (Mim.s.p$PLOT[1] != 5) {						# if plot number is 5 (mature forest):
          if (all(!is.na(Mim.s.p$Age))) {					# if it is not an empty list:
             lam.list <- c()							# create new empty list for lambda values
             min.a <- min(Mim.s.p$Age, na.rm = TRUE)				# set maximum height for that plot
             max.a <- max(Mim.s.p$Age, na.rm = TRUE)				# set minimum height for that plot
             n.1.i <- which(Mim.s.p$Age == min.a)				# get row number of trees
             n.1.h.i <- log(Mim.s.p[n.1.i,]$h2)[order(log(Mim.s.p[n.1.i,]$h2))] # get height of trees 
             init.n.a.v.i <- n.1.v.i <- hist(n.1.h.i, breaks = e.pred, plot = FALSE)$counts # get size structure vector and set initial vector
             for (a in ((min.a+1):max.a)) {					# for every year:
                  n.a.v.i <- k.p.list[[as.integer(i)]][(a-min.a),,]%*%init.n.a.v.i # multiply size vector by plot kernel
                  lam.a.i <- sum(n.a.v.i)/sum(init.n.a.v.i)			# get lambda 
                  init.n.a.v.i <- n.a.v.i					# reset initial vector
                  lam.list <- c(lam.list, lam.a.i)				# add lambda to list
              }
             est.lam.df[(min.a+1):max.a, as.numeric(i)] <- lam.list		# fill database
        }
       }
     }
     est.lam.list <- list(NA)							# create new list
     for (i in 1:nplot) {							# for each plot
      if (i != 5) {								# if nor plot 5 (not mature forest):
        plot.l.i <- as.data.frame(matrix(nrow = length(which(!is.na(est.lam.df[,i]) == TRUE)), ncol = 2)) # create new database for that plot
        names(plot.l.i) <- c("Age", "lambda")					# name columns
        plot.l.i$Age <- which(!is.na(est.lam.df[,i]))				# set ages
        plot.l.i$lambda <- est.lam.df[,i][which(!is.na(est.lam.df[,i]) == TRUE)]# set lambda values
        est.lam.list[[i]] <- plot.l.i 						# set plot number
        } else {								# if it is plot 5:
        plot.l.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))			# create new database
        names(plot.l.5) <- c("Age", "lambda")					# name columns
        plot.l.5$Age[1] <- Age.mature						# set age
        plot.l.5$lambda[1] <- lam.mature					# set lambda value
        est.lam.list[[5]] <- plot.l.5						# set plot
      }
     }
}
# observed lambda by plot
{   
           nplot <- nlevels(s1.db$PLOT)						# get number of plots
           s1.db <- transform(s1.db, PLOT = as.factor(as.character(plot)))	# transoform plots as factor
     	   plot.lam.list <- as.data.frame(matrix(nrow = Age.mature, ncol = length(levels(s1.db$PLOT)), NA)) # create new database for observed lambda list by plot
 	   names(plot.lam.list) <- 1:nplot					# name plot columns
 	   for (i in levels(s1.db$PLOT)) {					# for every plot:
  	 	 Mim.s.p <- subset(s1.db, PLOT == i)				# get a subset of the whole database
		 if (Mim.s.p$PLOT[1] != 5) {					# if plot number is not 5: !!! (tengo que ver por qué, está vacío)		
		  if (all(!is.na(Mim.s.p$Age))) {				# if subset has all ages registered:
		     lam.list <- c()						# create a new lambda list
		     min.a <- min(Mim.s.p$Age, na.rm = TRUE)			# get maximum age
		     max.a <- max(Mim.s.p$Age, na.rm = TRUE)			# get minimum age
		     n.1.i <- which(Mim.s.p$Age == min.a)			# get rows containing trees in year a
		     n.1.h.i <- log(Mim.s.p[n.1.i,]$h2)[order(log(Mim.s.p[n.1.i,]$h2))] # get tree hight and order
		     init.n.a.v.i <- n.1.v.i <- hist(n.1.h.i, breaks = e.pred, plot = FALSE)$counts # count trees in each size class and set a population structure vector
		     for (a in ((min.a+1):max.a)) {				# for each age in plot subset:
			  n.a.db <- droplevels(subset(Mim.s.p, Age == a))	# drop levels
			  n.a.v.i <- nlevels(n.a.db$id)				# count trees
			  lam.a.i <- n.a.v.i/sum(init.n.a.v.i)			# get lambda
			  init.n.a.v.i <- n.a.v.i				# reset initial n
			  lam.list <- c(lam.list, lam.a.i)			# add lambda to list
		      }
		     plot.lam.list[(min.a+1):max.a, as.numeric(i)] <- lam.list	# set lambda values and plot
		}
	       }
	   }
	   nplot <- nlevels(s1.db$PLOT)	
 	   obs.lam.plot <- list(NA)						# create lambda list
 	   for (i in 1:nplot) {							# for each plot:
      		if (i != 5) {							# if it is not 5:
        		plot.l.i <- as.data.frame(matrix(nrow = length(which(!is.na(plot.lam.list[,i]) == TRUE)), ncol = 2)) # create new database: cols = 2 & rows = number of not NA in list
        		names(plot.l.i) <- c("Age", "lambda")			# set names
       			plot.l.i$Age <- which(!is.na(plot.lam.list[,i]))	# set valid ages
        		plot.l.i$lambda <- plot.lam.list[,i][which(!is.na(plot.lam.list[,i]) == TRUE)] # set lambda values
        		obs.lam.plot[[i]] <- plot.l.i 				# add lambda matrix to list
			}							
      		else if (i == 5) {						# if plot is 5:
        		plot.l.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))	# create new database for mature forest
        		names(plot.l.5) <- c("Age", "lambda")			# name it
        		plot.l.5$Age[1] <- Age.mature				# set age as mature
        		plot.l.5$lambda[1] <- lam.mature			# set lambda as calculated for mature forest above
        		obs.lam.plot[[5]] <- plot.l.5				# set plot number
      			}
 	   	       }
 }
# estimating inmigration
{
  lambda.obs <- obs.lam.plot
  lambda.est <- est.lam.list
  distance.l <- function(c) {
  	  dist.l <- 0
  	  lam.est <- c()
	  for (p in 1:nplot)
		if (p != 5) {
		    S.p <- subset(s1.db, plot == p & !is.na(ln.h2))
		    S.p <- droplevels(S.p)
		    min.a <- min(S.p$Age)
		    max.a <- max(S.p$Age)
		    n.1 <- which(S.p$Age == min.a)
		    n.1.h <- log(S.p[n.1,]$h2)
		    n.1.v <- hist(n.1.h, breaks = e.pred, plot = FALSE)$counts
		    init.n.a.v <- n.1.v+c*F5
		    for (a in (min.a+1):max.a) {
			n.a.v <- k.p.list[[p]][(a-min.a),,]%*%init.n.a.v
			lam.a <- (sum(n.a.v)+c)/sum(init.n.a.v)
			init.n.a.v <- n.a.v+c*F5
			lam.est <- c(lam.est, lam.a)
		    }
		    lambda.obs.plot <- rep(NA, nrow(lambda.obs[[p]]))
		    for (j in (1:nrow(lambda.obs[[p]])))
			dist.l <- dist.l + (lam.est[j]-lambda.obs[[p]]$lambda[j])^2
		dist.l <- sqrt(dist.l)
		cat(paste0("c = ", c, ", dist = ", dist.l, "\n"))
	        }
  	  return(dist.l)
 } 
 # Try different estimators
 {
    opt.brent <- optim(0, distance.l, method = "Brent", lower = 0, upper = 100)
    opt.brent$par # 10.79303, dist = 1.702837
    opt.lbf <- optim(0, distance.l, method = "L-BFGS-B")
    c <- opt.lbf$par # 10.79292, dist = 1.702837
    opt.bfgs <- optim(0, distance.l, method = "BFGS")
    opt.bfgs$par # 10.64919, dist = 1.702837
 }
}

# plots
# total, no migration
{
     lambda.df.NM <- as.data.frame(list(lambda = lam.list.NM, Age = 1:Age.mature))
     lambda.edad <- qplot(x = Age, y = log(lambda), data = lambda.df.NM, color = "red", geom = "line", xlab = "Edad sucesional (años)", ylab = expression(italic(r))) +
     	scale_fill_discrete(guide=FALSE) +
     	geom_point(data = ob.lam.df, mapping = aes(x = Age, y = log(ob.lambda)))
     	lambda.edad + 
     	theme_minimal() +
     	theme(legend.position = "none")
     lambda.edad
     #ggsave("lambda-edad.png", lambda.edad, device = "png", width = 9, height = 7, units = "in", dpi = 180*2)
}
# by plot, no migration
 {
 	 lambda.df.2.NM <- transform(lambda.df.NM)
 	 p.lam.p.2 <- as.data.frame(matrix(nrow = 0, ncol = 2))
 	 names(p.lam.p.2) <- c("Age", "lambda")
 	 for (i in 1:nplot) {
 		 p.lam.p.i <-transform(obs.lam.plot[[i]])
 		 p.lam.p.2  <- rbind(p.lam.p.2, p.lam.p.i)
 	 }
 	 p.lam.ob <- ggplot(p.lam.p.2, aes(x = Age, y = lambda, alpha = 1/2), show.legend = FALSE) + 
 	 theme_minimal() +
 	 geom_line(size = 1) + 
     	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) + 
 	 geom_line(data = lambda.df.2.NM, aes(x = Age, y = lambda), col = "red", size = 2, alpha = 1/3, show.legend = FALSE) +
 	 labs(x = expression(paste("Abandonment time ", italic(t), " (years)")), y = expression(lambda))+
 	 scale_alpha(guide = "none") +
	 scale_y_continuous(limits = c(0, 5.8))
 	 p.lam.ob
 	 #ggsave("no-migration.png", p.lam.ob, device = "png", width = 18, height = 12, units = "in", dpi = 180*2)
 }
# by plot, with migration
{
	obs.lam.p <- as.data.frame(matrix(nrow = 0, ncol = 3))
	names(obs.lam.p) <- c("Age", "lambda", "plot")
	for (i in 1:nplot) {
		obs.lam.p.i <-transform(obs.lam.plot[[i]], plot = as.factor(i))
		obs.lam.p  <- rbind(obs.lam.p, obs.lam.p.i)
		}
 # calculate new lambda list considering migration
 {
     	lam.list <- c()
     	n.0 <- which(s1.db$Age == 0)
     	n.0.h <- log(s1.db[n.0,]$h2)
     	n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
     	init.n.a.v <- n.0.v
	size.v.a.WM <- list(NA) # size structure vector by year, WM = with migration
     	n.list <- c(sum(init.n.a.v))
     	for (a in 1:Age.mature) {
       		n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
		size.v.a.WM[[a]] <- init.n.a.v
       		lam.a <- (sum(n.a.v)+c)/sum(init.n.a.v) # add c individuals
       		init.n.a.v <- n.a.v+c*F5		# add size structure
       		n.list <- c(n.list, sum(n.a.v))
       		lam.list <- c(lam.list, lam.a)
       		}
      	tot.lam.pred <- lam.list
	lambda.df <- as.data.frame(list(lambda = lam.list, Age = 1:Age.mature))
   	lambda.df.2 <- transform(lambda.df, plot = as.factor("total"))
   } 

	plot.l.graf.c <- ggplot(lambda.df.2, aes(x = Age, y = lambda)) +
				theme_minimal() +
				xlab(expression(paste("Abandonment time ", italic(t), " (years)"))) +
				ylab(expression(lambda)) +
				theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) + 
				geom_line(size = 2, alpha = 1/3, col = "red") + 
				geom_line(data = obs.lam.p, aes(x = Age, y = lambda), alpha = 1/2, size = 1) +
				scale_alpha(guide = "none") +
				scale_y_continuous(limits = c(0, 5.8))
	plot.l.graf.c
	#ggsave("with-migration.png", p.lam.ob, device = "png", width = 18, height = 12, units = "in", dpi = 180*2)
}
# size structure change over time
{
	pop.size.NM <- c()
	pop.size.WM <- c()
	size.v.a.NM.s <- list(NA)
	size.v.a.WM.s <- list(NA)
	for (a in 1:Age.mature) {
		pop.size.NM <- c(pop.size.NM, sum(size.v.a.NM[[a]]))		# get population size 
		size.v.a.NM.s[[a]] <- size.v.a.NM[[a]]/sum(size.v.a.NM[[a]])		# standardized sized vector, NO MIGRATION
		pop.size.WM <- c(pop.size.WM, sum(size.v.a.WM[[a]]))
		size.v.a.WM.s[[a]] <- size.v.a.WM[[a]]/sum(size.v.a.WM[[a]])	
		}
	plot(1:Age.mature, pop.size.NM, main = "Population vector, no migration") # plot population size
	plot(1:Age.mature, pop.size.WM, main = "Population vector, with migration")
	#
	dat.NM <- as.data.frame(matrix(nrow = 10*Age.mature, ncol = 3)) # create new df for size structure change, NO MIGRATION
	names(dat.NM) <- c("Age", "Prob.j", "size")
	count <- 1:10
	for (i in count) {
		dat.NM[((i-1)*Age.mature+1):(i*Age.mature),1] <- rep((count[i]*10)-9, Age.mature)
		dat.NM[((i-1)*Age.mature+1):(i*Age.mature),2] <- size.v.a.NM.s[[(count[i]*10)-9]]
		dat.NM[((i-1)*Age.mature+1):(i*Age.mature),3] <- exp(x.pred)
		}
	dat.NM[901:1000,1] <- rep(Age.mature,Age.mature)
	dat.NM[901:1000,2] <- size.v.a.NM.s[[Age.mature]]
	dat.NM$Age <- as.factor(as.character(dat.NM$Age))
	#
	dat.WM <- as.data.frame(matrix(nrow = 10*Age.mature, ncol = 3)) # create new df for size structure change, WITH MIGRATION
	names(dat.WM) <- c("Age", "Prob.j", "size")
	count <- 1:10
	for (i in count) {
		dat.WM[((i-1)*Age.mature+1):(i*Age.mature),1] <- rep((count[i]*10)-9, Age.mature)
		dat.WM[((i-1)*Age.mature+1):(i*Age.mature),2] <- size.v.a.WM.s[[(count[i]*10)-9]]
		dat.WM[((i-1)*Age.mature+1):(i*Age.mature),3] <- exp(x.pred)
		}
	dat.WM[901:1000,1] <- rep(Age.mature,Age.mature)
	dat.WM[901:1000,2] <- size.v.a.WM.s[[Age.mature]]
	#
	size.str.mat.NM <- matrix(unlist(size.v.a.NM.s), ncol = 100, byrow = TRUE)
	size.str.mat.WM <- matrix(unlist(size.v.a.WM.s), ncol = 100, byrow = TRUE)
	hist3D(exp(x.pred), (1:Age.mature), size.str.mat.NM, col = "grey", border = "black", xlab = "Size", ylab = "Age", zlab = "Probability", main = "Size structure change without migration")	
	hist3D(exp(x.pred), (1:Age.mature), size.str.mat.WM, col = "grey", border = "black", xlab = "Size", ylab = "Age", zlab = "Probability", main = "Size structure change with migration")


# esto ya no está limpio

# plots F5
  f5.p.df <- data.frame(f5.db$Age, f5.db$h1)
  f5.p.gg <- ggplot(f5.p.df, aes(x=f5.db.h1)) + 
    theme_minimal() +
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
    xlab(expression(paste("Recruits height (m)"))) +
    ylab(expression("Probability"))
   #ggsave("art-f5.png", f5.p.gg, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)

	#plots F4
	
 	#png("f4Plot.png", width = 900, height = 600, units = "px", pointsize = 20)
 	f4.pred <- seq(0, Age.mature)
 	plot(x, f4.a, xlim = c(0, Age.mature), xlab = "Edad sucesional (años)", ylab = "Probabilidad de establecimiento")
 	#lines(f4.pred, predict(glm.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "green")
    #lines(f4.pred, predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "blue")
    lines(f4.pred, predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "red")
 	#dev.off()
}

