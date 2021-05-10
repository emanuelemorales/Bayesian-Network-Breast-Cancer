#import libraries
library(gRbase)
library(Rgraphviz)
library(igraph)
library(bnlearn)
library(dplyr)
library(heatmaply)

#Import dataset and pre-processing:______________________________________________
#df_cleaned is the dataframe with the variables categorized.
df_cleaned = read.csv("data_cleaned.csv", sep = ";", encoding = "UTF-8 BOM")
colnames(df_cleaned)[1] <- "CENTRO"

#Import data about dietary patterns
df_diet = read.csv("pesi.csv", encoding = "UTF-8 BOM")
colnames(df_diet)[1] <- "X1"

#union of the two datasets
df = cbind(df_cleaned, df_diet)
head(df)

#Cancer = 1/ no cancer = 0
df$V2[df$V2 == 2] = 0
summary(df)

#Scale data
df$X1 = scale(df$X1)
df$X2 = scale(df$X2)
df$X3 = scale(df$X3)
df$X4 = scale(df$X4)

#categorize the diets
df[,c("X1", "X2", "X3", "X4")] = as.data.frame(sapply(df[,c("X1", "X2", "X3", "X4")], function (x) {
  qx <- quantile(x)
  cut(x, qx, include.lowest = TRUE,
      labels = 1:4)
}))

#Factorize the categorical variables
names <- c(1:ncol(df))
df[,names] <- lapply(df[,names] , factor)
str(df)

#Diet_tot
diet = c("AGE", "EDU", "GIN4", "ANTR0", "CHILD", "SMOKE", "ALCOHOL", "FIS4", "X2","X3","X4", "V2")
diet.df = df[,diet]

#whitelist and blacklist based on information of the logistic analysis
wl = data.frame(from =c("AGE", "AGE", "AGE","EDU","EDU","EDU","GIN4","GIN4","GIN4","GIN4", "ANTR0", "ANTR0", "CHILD","CHILD", "SMOKE",
                        "SMOKE","SMOKE","SMOKE","ALCOHOL","ALCOHOL","ALCOHOL","ALCOHOL","FIS4","FIS4","FIS4", "X2", "X3", "X4"),
                
                to = c("X2", "X3", "X4","V2","X2","X3", "V2", "X2", "X3", "X4","X2","X4", "V2","X2", "V2", "X2", "X3", "X4", 
                       "V2", "X2", "X3", "X4", "X2", "X3", "X4", "V2", "V2", "V2"))

bl = data.frame(from = c("V2"), to = c("EDU"))

#APPLICATION OF LEARNING PROCEDURES______________________________________________________________________
#TRY BOOTSTRAP IN PC MODEL
set.seed(1)
str.diff_pc = boot.strength(diet.df, R = 100, algorithm = "pc.stable", 
                            algorithm.args = list(whitelist = wl, blacklist = bl))

avg.diff_pc = averaged.network(str.diff_pc, threshold = 0.5)
strength.plot(avg.diff_pc, str.diff_pc, shape = "ellipse", highlight = list(arcs = wl, lwd = 0.5))

pc_mat = amat(avg.diff_pc)
pc_mat

#TRY BOOTSTRAP IN HC MODEL_______________________________________________________________________
set.seed(1)
str.diff_hc = boot.strength(diet.df, R = 100, algorithm = "hc",  
                            algorithm.args = list(score="bic", whitelist = wl, blacklist = bl))
avg.diff_hc = averaged.network(str.diff_hc, threshold = 0.5)
strength.plot(avg.diff_hc, str.diff_hc, shape = "ellipse", highlight = list(arcs = wl, lwd = 0.5))

head(str.diff_hc)
hc_mat = amat(avg.diff_hc)
hc_mat

#TRY BOOTSTRAP IN GS MODEL______________________________________________________________________
set.seed(1)
str.diff_gs = boot.strength(diet.df, R = 100, algorithm = "gs",  
                            algorithm.args = list(whitelist = wl, blacklist = bl))
avg.diff_gs = averaged.network(str.diff_gs, threshold = 0.5)
strength.plot(avg.diff_gs, str.diff_gs, shape = "ellipse", highlight = list(arcs = wl, lwd = 0.5))

gs_mat = amat(avg.diff_gs)
gs_mat

#TRY BOOTSTRAP IN TABU GREEDY MODEL_____________________________________________________________
set.seed(1)
str.diff_tabu = boot.strength(diet.df, R = 100, algorithm = "tabu",  
                            algorithm.args = list(whitelist = wl, blacklist = bl))
attr(str.diff_tabu, "threshold")
strength.plot(avg.diff_tabu, str.diff_tabu, shape = "ellipse", highlight = list(arcs = wl, lwd = 0.5))

tabu_mat = amat(avg.diff_tabu)
tabu_mat

#TRY BOOTSTRAP IN MAX-MIN HILL CLIMBING MODEL_______________________________________________________
set.seed(1)
str.diff_rsmax2 = boot.strength(diet.df, R = 100, algorithm = "rsmax2",  
                              algorithm.args = list(whitelist = wl, blacklist = bl))

avg.diff_rsmax2 = averaged.network(str.diff_rsmax2, threshold = 0.5)
strength.plot(avg.diff_rsmax2, str.diff_rsmax2, shape = "ellipse", highlight = list(arcs = wl, lwd = 0.5))

rsmax2_mat = amat(avg.diff_rsmax2)
rsmax2_mat

#Union of the matrices
mat_sum = hc_mat + gs_mat + pc_mat + tabu_mat + rsmax2_mat
mat_sum

#create a bayesian network based on the results of the structure learning algorithms matrix:
bn_from_mat = function(adj_mat, tshl){
  adj = adj_mat
  adj[which(mat_sum < tshl)] = 0
  adj[which(mat_sum >= tshl)] = 1
  
  model = empty.graph(colnames(adj))
  amat(model) = adj
  
  return(model)
}
plt_bn = function(model){
  gR = graphviz.plot(model, layout = "dot", shape = "rectangle", 
                     highlight = list(nodes = c("V2"), col = c("tomato"), fill = "orange"))
  node.attrs = nodeRenderInfo(gR)
  node.attrs$textCol[c("X3","X2","X4")] ="tomato"
  node.attrs$shape[c("X3","X2","X4")] ="ellipse"
  nodeRenderInfo(gR) = node.attrs
  
  return(renderGraph((gR)))
}

#Base model:
var = colnames(rsmax2_mat)
e = empty.graph(var)
arcs(e) = wl

#model_2: threshold = 2
model_2 = bn_from_mat(mat_sum, 2)
model_2 = cextend(model_2)

#model_3: threshold = 3
model_3 = bn_from_mat(mat_sum, 3)
model_3 = cextend(model_3)

#model_5 = threshold = 5
model_5 = bn_from_mat(mat_sum, 5)
model_5 = cextend(model_5)

par(mfrow = c(2,2))
graphviz.compare(e, model_5, model_3, model_2, shape = "ellipse", 
                 main = c("knowledge based", "knowledge based + threshold = 5", "knowledge based + threshold = 3", "knowledge based + threshold 2"))


#We can use the this bayesian model as an expert model and predict the values of one or more
#variables for new individuals, based on the values of some other variables.
#We can measure the predictive accuracy by using the cross-validation:
#here it is applied a k-fold on the mode_5 fixed structure to cross-validate the parameters
#learned.

#Comparison of the CV loss function of the four models:
cv.bic_base = bn.cv(diet.df, bn = e, runs = 10, 
                    algorithm.args = list(score = "bic"))

cv.bic_2 = bn.cv(diet.df, bn = model_2, runs = 10, 
                 algorithm.args = list(score = "bic"))

cv.bic_3 = bn.cv(diet.df, bn = model_3, runs = 10, 
                 algorithm.args = list(score = "bic"))

cv.bic_5 = bn.cv(diet.df, bn = model_5, runs = 10, 
                 algorithm.args = list(score = "bic"))


plot(cv.bic_base, cv.bic_5, cv.bic_3, cv.bic_2, xlab = c("BIC KB", "BIC t=5", "BIC t = 3", "BIC t = 2"))


#Chose the model with thresshold = 2 to prune arcs
fit = bn.fit(cextend(model_2), diet.df)


par(mfrow = c(1,1))
gR = graphviz.plot(fit, layout = "dot", shape = "rectangle", 
                   highlight = list(nodes = c("V2"), col = c("tomato"), fill = "orange"))

node.attrs = nodeRenderInfo(gR)
node.attrs$textCol[c("X3","X2","X4")] ="tomato"
node.attrs$shape[c("X3","X2","X4")] ="ellipse"
nodeRenderInfo(gR) = node.attrs
renderGraph((gR))

### MODEL STRATIFICATION ###
#Looking at the network structure observing the path exploiting the causality of the model.
#Now we can stratify by confounding variables:
#Values of the variables:
age_values = unique(df$AGE)
edu_values = unique(df$EDU)
smoke_values = unique(df$SMOKE)
alcohol_values = unique(df$ALCOHOL)
bmi_values = unique(df$ANTR0)
child_values = unique(df$CHILD)
gin_values = unique(df$GIN4)
fis_values = unique(df$FIS4)


#X2 PLOT:_____________________________________________________________________________________________________

par(mfrow = c(1,2))
#Probability of cancer by varying X2:
probs_X2 = c()
for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X2 == value))
  probs_X2[i] = prob
}

#Stratification by alcohol:
X2_ALCOHOL = list()
for (j in 1:length(alcohol_values)){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & ALCOHOL == alcohol_values[j]))
    probs[i] = prob
  }
  X2_ALCOHOL[[j]] = probs
}

#PLOT X2/ALCOHOL:
g_range = range(0.1, 0.6)


plot(probs_X2, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n") 
axis(1, at = 1:4) 
#abline(h = prob_cancer, col = "grey") +
lines(X2_ALCOHOL[[1]], type = "o", col = "green",lwd = 1) 
lines(X2_ALCOHOL[[2]], type = "o", col = "red", pch = 19, lwd = 1) 
lines(X2_ALCOHOL[[3]], type = "o", col = "orange", pch = 19, lwd = 1) 
lines(X2_ALCOHOL[[4]], type = "o", col = "blue", pch = 19, lwd = 1) 
title(main="Prob. of cancer given X2 and Alcohol Status", col.main="red") 
title(xlab="Quartiles of diet X2") 
title(ylab="Probability of cancer") 
legend(1, g_range[2], c("0 unit/week ","1-5 unit/week", "5-10 unit/week", ">10 unit/week", "no strat."), cex=0.8, 
       col=c("blue","orange", "red", "green", "black"), pch = 21, lty = 1, title = "Stratified by Alcohol Status:", horiz = FALSE)



#Stratification by child:
X2_CHILD = list()
for (j in 1:length(child_values)){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & CHILD == child_values[j]))
    probs[i] = prob
  }
  
  X2_CHILD[[j]] = probs
}

#PLOT X2/CHILD:
g_range = range(0.1, 0.6)
plot(probs_X2, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X2_CHILD[[1]], type = "o", col = "green",lwd = 1)
lines(X2_CHILD[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X2_CHILD[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
title(main="Probabilities of cancer given X2 and number of children", col.main="red")
title(xlab="Quartiles of diet X2")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("0","1", "2+", "no strat."), cex=0.8, 
       col=c("green","red", "orange", "black"), pch = 21, lty = 1, title = "Stratified by num. of children:", horiz = FALSE)

#X2: Probability with confounding BMI:
X2_BMI = list()
for (j in 1:length(bmi_values)){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X2 == value & ANTR0 == bmi_values[i]))
    probs[i] = prob
  }
  
  X2_BMI[[j]] = probs
}

#PLOT X2/BMI:
g_range = range(0.2, 0.6)
plot(probs_X2, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X2_BMI[[1]], type = "o", col = "green",lwd = 1)
lines(X2_BMI[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X2_BMI[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
lines(X2_BMI[[4]], type = "o", col = "blue", pch = 19, lwd = 1)
title(main="Prob. of cancer given X2 and BMI", col.main="red")
title(xlab="Quartiles of diet X2")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("underweight","normal weight", "overweight", "obesity", "no strat."), cex=0.8, 
       col=c("blue","orange", "green", "red", "black"), pch = 21, lty = 1, title = "Stratified by BMI:", horiz = FALSE)


#X3 PLOT_________________________________________________________________________________________________________________________________
#Probability of cancer by varying X3
par(mfrow = c(1,2))

probs_X3 = c()

for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X3 == value))
  probs_X3[i] = prob
}
probs_X3

#Stratification by EDU:
X3_EDU = list()
for (j in 1:length(edu_values)){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & EDU == edu_values[j]))
    probs[i] = prob
  }
  X3_EDU[[j]] = probs
}

#PLOT X3/EDU:
g_range = range(0.0, 0.6)
plot(probs_X3, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X3_EDU[[1]], type = "o", col = "green",lwd = 1)
lines(X3_EDU[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X3_EDU[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
title(main="Prob. of cancer given X3 and Edu.", col.main="red")
title(xlab="Quartiles of diet X3")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("elementary","middle/high", "university", "no strat."), cex=0.8, 
       col=c("green","red", "orange", "black"), pch = 21, lty = 1, title = "Stratified by Education :", horiz = FALSE)


#X3: Probability with confounding GIN4:
X3_GIN = list()
for (j in 1:length(gin_values)){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X3 == value & GIN4 == toString(j)))
    probs[i] = prob
  }
  X3_GIN[[j]] = probs
}

#PLOT X3/GIN4:
g_range = range(0.0, 0.6)
plot(probs_X3, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X3_GIN[[1]], type = "o", col = "green",lwd = 1)
lines(X3_GIN[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X3_GIN[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
title(main="Prob. of cancer given X3 and Menop. status", col.main="red")
title(xlab="Quartiles of diet X3")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("pre-menopause","peri-menopause", "post-menopause", "no strat."), cex=0.8, 
       col=c("green","red", "orange", "black"), pch = 21, lty = 1, title = "Stratified by menopausal status:", horiz = FALSE)

#X4 PLOT_________________________________________________________________________________________________________________________________
#Probability of cancer by varying X4
par(mfrow = c(1,2))
probs_X4 = c()
for (i in 1:4){
  value = toString(i)
  prob = cpquery(fit, (V2 == "1") , (X4 == value))
  probs_X4[i] = prob
}

#PLOT X4/SMOKE
str(df$SMOKE)
X4_SMOKE = list()
for (j in 1:length(smoke_values)){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X4 == value & SMOKE == smoke_values[j]))
    probs[i] = prob
  }
  
  X4_SMOKE[[j]] = probs
}

#plot X4/SMOKE
g_range = range(0.1, 0.5)
plot(probs_X4, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X4_SMOKE[[1]], type = "o", col = "green",lwd = 1)
lines(X4_SMOKE[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X4_SMOKE[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
title(main="Prob. of cancer given X4 and SMOKE", col.main="red")
title(xlab="Quartiles of diet X4")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("non smoker","smoker for <=25y", "smoker for >25y", "No strat."), cex=0.8, 
       col=c("green","red", "orange", "black"), pch = 21, lty = 1, title = "Stratified by smoking status", horiz = FALSE)


#X4: Probability with confounding FIS4:
X4_FIS4 = list()
for (j in 1:length(fis_values)){
  probs = c()
  for (i in 1:4){
    value = toString(i)
    prob = cpquery(fit, (V2 == "1") , (X4 == value & FIS4 == toString(j)))
    probs[i] = prob
  }
  
  X4_FIS4[[j]] = probs
}

g_range = range(0.1, 0.5)
plot(probs_X4, type = "o", col = "black", lty = 4, ylim = g_range, ann = FALSE, x = 1:4, xaxt = "n")
axis(1, at = 1:4)
#abline(h = prob_cancer, col = "grey")
lines(X4_FIS4[[1]], type = "o", col = "green",lwd = 1)
lines(X4_FIS4[[2]], type = "o", col = "red", pch = 19, lwd = 1)
lines(X4_FIS4[[3]], type = "o", col = "orange", pch = 19, lwd = 1)
lines(X4_FIS4[[4]], type = "o", col = "blue", pch = 19, lwd = 1)
title(main="Prob. of cancer given X4 and Physical Activity", col.main="red")
title(xlab="Quartiles of diet X4")
title(ylab="Probability of cancer")
legend(1, g_range[2], c("<2H","2-4 H", "5-7 H", ">7H", "No strat."), cex=0.8, 
       col=c("red","green", "blue","orange", "black"), pch = 21, lty = 1, title = "Stratified by Physical Activity (19-25 y.o.)", horiz = FALSE)

