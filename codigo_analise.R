# clean workspace
rm(list=ls())

library(confintr) # package for confidence intervals
library(ggplot2)
library(car) # <--- Install if needed with install.packages("car")
library(lmtest)
library(multcomp)
library(knitr)
library(stringr)


#===================
####
# Load data
####

data <- read.table("C:/Users/Ariani/Downloads/ec03_ufmg-main/dados_final.csv",
                   header = TRUE)

summary(data)


# Aggregate data (algorithm means by instance group)
aggdata <- with(data,
                aggregate(x   = Resultado,
                          by  = list(ConfiguraC'C#o, DimensC#o),
                          FUN = mean))

# Rename columns
names(aggdata) <- c("ConfiguraC'C#o", 
                    "DimensC#o",
                    "Y")

# Coerce categorical variables to factors
for (i in 1:2){
  aggdata[, i] <- as.factor(aggdata[, i])
}

summary(aggdata)

####
## Exploratory data analysis: plot observations by Configuration and Dependency.
####

# png(filename = "../figs/algo_lineplot.png",
#     width = 1000, height = 400, 
#     bg = "transparent")
p <- ggplot(aggdata, aes(x = DimensC#o, 
                         y = Y, 
                         group = ConfiguraC'C#o, 
                         colour = ConfiguraC'C#o))
p + geom_line(linetype=2) + geom_point(size=5)
# dev.off()

#Boxplot
ggplot(aggdata,
       aes(x = ConfiguraC'C#o,
           y = Y,
           fill = ConfiguraC'C#o))+
  geom_boxplot() + geom_point() +
  ggtitle("DimensC#o ConfiguraC'C#o", "(original data + boxplots)") + theme(legend.position = "name")

####
## Statistical modeling
####

# First model
# Response is result of the overall mean plus configuration effect (treatment)
# and dependency (block) effect

model <- aov(Y~ConfiguraC'C#o+DimensC#o,
             data = aggdata)

summary(model)
summary.lm(model)$r.squared

#Verificando a independC*ncia dos dados
inde <- dwtest(model)

#Verificando normalidade dos resC-duo
shap = shapiro.test(model$residuals)

qqPlot(model$residuals,
       pch = 16,
       lwd = 3,
       cex = 2,
       las = 1)

#Verificando a premissa de homoscedasticidade
#teste de Fligner-Killen e grC!fico dos resC-duos pelos valores ajustados.
fligner <- fligner.test(Y~ConfiguraC'C#o, data = aggdata)
plot(x = model$fitted.values,
     y = model$residuals)

#Estimativa do poder do teste
n <- 64
sigma = 3.257e+14
power_result <- power.anova.test(groups = model$rank , between.var = 3.257e+14 ,
                                 within.var = (sigma^2), sig.level = 0.05, power = NULL, n = n)
# A potC*ncia do teste C) dada por power_result$power.
power_value <- power_result$power

# Graphical test of assumptions
# png(filename = "../figs/algo_res1.png",
#     width = 1000, height = 500, 
#     bg = "transparent")
par(mfrow = c(2, 2))
plot(model, pch = 20, las = 1)
# dev.off()

################################################################################
# Try log-transformed data

model2 <- aov(log(Y)~ConfiguraC'C#o+DimensC#o,
              data = aggdata)
summary(model2)
summary.lm(model2)$r.squared

#Verificando a independC*ncia dos dados
inde2 <- dwtest(model2)

#Verificando normalidade dos resC-duo
shap2 = shapiro.test(model2$residuals)

qqPlot(model2$residuals,
       pch = 16,
       lwd = 3,
       cex = 2,
       las = 1)

#Verificando a premissa de homoscedasticidade
#teste de Fligner-Killen e grC!fico dos resC-duos pelos valores ajustados.
fligner <- fligner.test(Y~ConfiguraC'C#o, data = aggdata)
plot(x = model2$fitted.values,
     y = model2$residuals)

#Estimativa do poder do teste
n <- 64
sigma =  669
power_result2 <- power.anova.test(groups = model2$rank , between.var = 669.5 ,
                                 within.var = (sigma^2), sig.level = 0.05, power = NULL, n = n)
# A potC*ncia do teste C) dada por power_result$power.
power_value2 <- power_result2$power

# Graphical test of assumptions
# png(filename = "../figs/algo_res2.png",
#     width = 1000, height = 500, 
#     bg = "transparent")
par(mfrow = c(2, 2))
plot(model2, pch = 20, las = 1)
# dev.off()

################################################################################

# model diagnostic
par(mfrow = c(1, 1))
qqPlot(model$residuals, pch = 20, las = 1)

# model2 diagnostic
par(mfrow = c(1, 1))
qqPlot(model2$residuals, pch = 20, las = 1)

####
## Blocking efficiency 
####

#For model
mydf        <- as.data.frame(summary(model)[[1]])
MSblocks    <- mydf["DimensC#o","Mean Sq"]
MSe         <- mydf["Residuals","Mean Sq"]
a           <- length(unique(aggdata$ConfiguraC'C#o))
b           <- length(unique(aggdata$DimensC#o))
((b - 1) * MSblocks + b * (a - 1) * MSe) / ((a * b - 1) * MSe)

#For model2

mydf        <- as.data.frame(summary(model2)[[1]])
MSblocks    <- mydf["DimensC#o","Mean Sq"]
MSe         <- mydf["Residuals","Mean Sq"]
a           <- length(unique(aggdata$Algorithm))
b           <- length(unique(aggdata$Instance_Group))
((b - 1) * MSblocks + b * (a - 1) * MSe) / ((a * b - 1) * MSe)

####
## Post-hoc multiple comparisons
####

# using model 
duntest     <- glht(model,
                    linfct = mcp(ConfiguraC'C#o = "Dunnett"))

summary(duntest)

duntestCI   <- confint(duntest)
# png(filename = "../figs/algo_mcp.png",
#     width = 1000, height = 500, 
#     bg = "transparent")
par(mar = c(5, 8, 4, 2), las = 1)
plot(duntestCI,
     xlab = "Mean difference (log scale)")
# dev.off()

# using model 2 - log-transformed data
duntest2     <- glht(model2,
                    linfct = mcp(ConfiguraC'C#o = "Dunnett"))

summary(duntest2)

duntestCI2   <- confint(duntest2)
# png(filename = "../figs/algo_mcp.png",
#     width = 1000, height = 500, 
#     bg = "transparent")
par(mar = c(5, 8, 4, 2), las = 1)
plot(duntestCI2,
     xlab = "Mean difference (log scale)")
# dev.off()




