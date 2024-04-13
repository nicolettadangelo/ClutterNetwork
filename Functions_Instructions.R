
#########################################
########## Loading libraries ############
#########################################

library(MASS)
library(spatstat)
library(segmented)


#########################################
##### Creating working directory ########
#########################################

source("rates_localdetect.R")
source("class_net.R")
source("nncleanEngine2.R")
source("dknn2.R")
source("calculate_statistics.R")
source("class_net_Engine.R")


##### Instructions for using the classification method ######


## 1. The method will be applied to a "X" object of class lpp

class(X)[1] == "lpp"


## 2. The linear network must be connected

is.connected(domain(X)) == FALSE

X <- lpp(as.ppp(X), subnet) 
# where "subnet" is the largest connected component in the list
# connected(domain(X), what = "components")

is.connected(domain(X)) == TRUE


## 3. The linear network must have a non-sparse representation

domain(X)[["sparse"]]==TRUE

X <- lpp(as.ppp(X), as.linnet(domain(X)[["lines"]]))

domain(X)[["sparse"]]==FALSE


## 4. Applying the classification method

classification <- class_net(X, n_feat = 1, verbose = F)


## 5. Plotting the classified point patterns

marks(X) <- as.factor(classification$id)

plot(X, col = "grey", cols = c("#56B4E9", "#E69F00"), 
       chars = c(1, 1), main = "")  
# The points of the 'feature' mark are those that belong to the clusters 
# according to the classification method
