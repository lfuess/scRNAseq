##this is the code to run statistics looking for significant differences in cell type populations across populations##

##briefly we run separate GLMs for each cell type, looking at the effect of population on cell type abundance##
##we then use Tukey's post-hoc tests to look at pairwise comparisons across population types##

#install.packages("multcomp")
library(multcomp)

##HCs##
HC_data = read.csv("HCs.csv")
HC_data$Population = as.factor(HC_data$Population)

lrfit = glm( cbind(HCs, Others) ~ Population, 
             data = HC_data, family = binomial)
summary(lrfit)
summary(glht(lrfit, mcp(Population="Tukey")))

##Neutrophils##
Neut_data = read.csv("Neutrophils.csv")
Neut_data$Population = as.factor(Neut_data$Population)
lrfit = glm( cbind(Neutrophils, Others) ~ Population, 
             data = Neut_data, family = binomial)
summary(lrfit)
summary(glht(lrfit, mcp(Population="Tukey")))

##APCs##
APC_data = read.csv("APCs.csv")
APC_data$Population = as.factor(APC_data$Population)
lrfit = glm( cbind(APCs, Others) ~ Population, 
             data = APC_data, family = binomial)
summary(lrfit)
summary(glht(lrfit, mcp(Population="Tukey")))

##Fibs- NS##
Fib_data = read.csv("Fibroblasts.csv")
Fib_data$Population = as.factor(Fib_data$Population)
lrfit = glm( cbind(Fibroblasts, Others) ~ Population, 
             data = Fib_data, family = binomial)
summary(lrfit)
summary(glht(lrfit, mcp(Population="Tukey")))

##Bcells##
B_data = read.csv("Bcells.csv")
B_data$Population = as.factor(B_data$Population)
lrfit = glm( cbind(B.cells, Others) ~ Population, 
             data = B_data, family = binomial)
summary(lrfit)
summary(glht(lrfit, mcp(Population="Tukey")))

##RBCs##
RBC_data = read.csv("RBCs.csv")
RBC_data$Population = as.factor(RBC_data$Population)
lrfit = glm( cbind(RBCs, Others) ~ Population, 
             data = RBC_data, family = binomial)
summary(lrfit)
summary(glht(lrfit, mcp(Population="Tukey")))

##Platelets##
Plat_data = read.csv("Platelets.csv")
Plat_data$Population = as.factor(Plat_data$Population)
lrfit = glm( cbind(Platelets, Others) ~ Population, 
             data = Plat_data, family = binomial)
summary(lrfit)
summary(glht(lrfit, mcp(Population="Tukey")))

##NKs##
NK_data = read.csv("NKs.csv")
NK_data$Population = as.factor(NK_data$Population)
lrfit = glm( cbind(NKs, Others) ~ Population, 
             data = NK_data, family = binomial)
summary(lrfit)
summary(glht(lrfit, mcp(Population="Tukey")))
