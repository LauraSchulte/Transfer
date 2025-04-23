###Transfer experiment
#authors: Laura Schulte, Pia Oswald, Eva Rousselle, Max Mühlenhaupt, Manuela Schmidt, Barbara A. Caspers
#Created: 15.10.2024
#last modified: 23.04.2025




###############################################################
#Is the recapture rate evenly distributed among the different locations and types of habitat?


####2024 
#transfer locations

# Observed recaptures per location
observed <- c(S1 = 11, S2 = 8, P1 = 23, P3 = 6)

# Numbers of captures from day of transfer as basis 
basis <- c(S1 = 67, S2 = 45, P1 = 76, P3 = 46)

# total number of recaptures 
n_observed <- sum(observed)

# expected values based on relative abundacne from first capture
expected <- (basis / sum(basis)) * n_observed

# expected and observed values 
cat("Beobachtete Werte:\n")
print(observed)
cat("\nErwartete Werte:\n")
print(round(expected, 2))

# calculate chi2 test
chi2_wert <- sum((observed - expected)^2 / expected)

# df
df <- length(observed) - 1

# calculate p-value
p_wert <- pchisq(chi2_wert, df = df, lower.tail = FALSE)

# results
cat("\nChi²-Wert:", round(chi2_wert, 3), "\n")
cat("Freiheitsgrade:", df, "\n")
cat("p-Wert:", round(p_wert, 4), "\n")

chisq.test(x = observed, p = basis / sum(basis))



####2023 
#transfer locations


# Observed recaptures per location
observed <- c(S1 = 20, S2 = 8, P1 = 15, P3 = 0)

# Numbers of captures from day of transfer as basis 
basis <- c(S1 = 49, S2 = 42, P1 = 78, P3 = 52)

# total number of recaptures 
n_observed <- sum(observed)

# expected values based on relative abundacne from first capture
expected <- (basis / sum(basis)) * n_observed

# expected and observed values 
cat("Beobachtete Werte:\n")
print(observed)
cat("\nErwartete Werte:\n")
print(round(expected, 2))

# calculate chi2 test
chi2_wert <- sum((observed - expected)^2 / expected)

# df
df <- length(observed) - 1

# calculate p-value
p_wert <- pchisq(chi2_wert, df = df, lower.tail = FALSE)

# results
cat("\nChi²-Wert:", round(chi2_wert, 3), "\n")
cat("Freiheitsgrade:", df, "\n")
cat("p-Wert:", round(p_wert, 4), "\n")

chisq.test(x = observed, p = basis / sum(basis))


###############################################################
#Is there a difference in the size of pond and stream larvae and the gill size at the beginning of the experiment?

####body size beginning of the exp####
#import dataset
library(readxl)
data <- read_excel("body_size_increase_23_24.xlsx", na="NA",  
                   col_types = c("text", "text", "text","text","text","numeric", "numeric", "numeric", "numeric", "numeric")) 
head(data)

library(dplyr)
library(car)

#inspect raw data 
hist(data$initial_body_size)
boxplot(data$initial_body_size ~ data$original_habitat)

#test for normal distribution & homogeneity of variance
shapiro.test(data$initial_body_size)

var.test(data$initial_body_size ~ data$original_habitat)


#transform data to normal distribution
# first, create objects with the transformations, e.g. logarithm etc.
library(bestNormalize)
(arcsinh_initial_body_size <- arcsinh_x(data$initial_body_size))
(boxcox_initial_body_size <- boxcox(data$initial_body_size))
(centerscale_initial_body_size <- center_scale(data$initial_body_size))
(orderNorm_initial_body_size <- orderNorm(data$initial_body_size))
(yeojohnson_initial_body_size <- yeojohnson(data$initial_body_size))
(sqrt_initial_body_size <- sqrt(data$initial_body_size))
(log_initial_body_size <- log(data$initial_body_size))

#histograms first impression
par(mfrow = c(2,4)) # display all histograms in one window
hist(data$initial_body_size)
MASS::truehist(arcsinh_initial_body_size$x.t, main = "Arcsinh transformation", nbins = 12) # x.t stands for the transformed variable
MASS::truehist(boxcox_initial_body_size$x.t, main = "Box Cox transformation", nbins = 12)
MASS::truehist(centerscale_initial_body_size$x.t, main = "center_scale transformation", nbins = 12)
MASS::truehist(orderNorm_initial_body_size$x.t, main = "orderNorm transformation", nbins = 12)
MASS::truehist(yeojohnson_initial_body_size$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(sqrt_initial_body_size, main = "squareroot transformation", nbins = 12)
MASS::truehist(log_initial_body_size, main = "log transformation", nbins = 12)

# let R suggest best transformation
bn.initial_body_size<-bestNormalize(data$initial_body_size, out_of_sample = FALSE) 
bn.initial_body_size

# create an object with automatically the best transformation using the following expression (e.g.) any.name <- orderNorm$x.t  # (i.e result of best normalization)
(yeojohnson_initial_body_size <- yeojohnson(data$initial_body_size))
yj.initial_body_size <- yeojohnson_initial_body_size$x.t

# add the new object as column to your data frame
initial_body_size.transformed <- cbind(data, yj.initial_body_size)
hist(initial_body_size.transformed$yj.initial_body_size)

# test for normal distribution and homogeneity of variance of the transformed variable
shapiro.test(initial_body_size.transformed$yj.initial_body_size)

var.test(initial_body_size.transformed$yj.initial_body_size ~ initial_body_size.transformed$original_habitat)


qqnorm(data$initial_body_size)
qqline(data$initial_body_size)


library(Matrix)
library(lme4)
library(lmerTest)
library(performance)
library(see)


modelperformance1<-lm(yj.initial_body_size ~ original_habitat + year, data=initial_body_size.transformed) 
modelperformance2<-lmer(yj.initial_body_size ~ original_habitat + (1|year), data=initial_body_size.transformed) 


AIC(modelperformance1,modelperformance2)

compare_performance(modelperformance1,modelperformance2, rank=T)


summary(modelperformance1)

# check model
check_model(modelperformance1)

############################################
#### gill size beginning of the exp ####

#import dataset
library(readxl)
data <- read_excel("gill_size_change_23_24.xlsx", na="NA",  
                   col_types = c("text", "text", "text","text","text", "numeric", "numeric", "numeric", "numeric", "numeric")) 
head(data)

library(dplyr)
library(car)

#inspect raw data 
hist(data$initial_gill_size)
boxplot(data$initial_gill_size ~ data$original_habitat)

#test for normal distribution & homogeneity of variance
shapiro.test(data$initial_gill_size)

var.test(data$initial_gill_size ~ data$original_habitat)


#transform data to normal distribution
# first, create objects with the transformations, e.g. logarithm etc.
library(bestNormalize)
(arcsinh_initial_gill_size <- arcsinh_x(data$initial_gill_size))
(boxcox_initial_gill_size <- boxcox(data$initial_gill_size))
(centerscale_initial_gill_size <- center_scale(data$initial_gill_size))
(orderNorm_initial_gill_size <- orderNorm(data$initial_gill_size))
(yeojohnson_initial_gill_size <- yeojohnson(data$initial_gill_size))
(sqrt_initial_gill_size <- sqrt(data$initial_gill_size))
(log_initial_gill_size <- log(data$initial_gill_size))

#histograms first impression
par(mfrow = c(2,4)) # display all histograms in one window
hist(data$initial_gill_size)
MASS::truehist(arcsinh_initial_gill_size$x.t, main = "Arcsinh transformation", nbins = 12) # x.t stands for the transformed variable
MASS::truehist(boxcox_initial_gill_size$x.t, main = "Box Cox transformation", nbins = 12)
MASS::truehist(centerscale_initial_gill_size$x.t, main = "center_scale transformation", nbins = 12)
MASS::truehist(orderNorm_initial_gill_size$x.t, main = "orderNorm transformation", nbins = 12)
MASS::truehist(yeojohnson_initial_gill_size$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(sqrt_initial_gill_size, main = "squareroot transformation", nbins = 12)
MASS::truehist(log_initial_gill_size, main = "log transformation", nbins = 12)

# let R suggest best transformation
bn.initial_gill_size<-bestNormalize(data$initial_gill_size, out_of_sample = FALSE) 
bn.initial_gill_size

# create an object with automatically the best transformation using the following expression (e.g.) any.name <- orderNorm$x.t  # (i.e result of best normalization)
(oderNorm_initial_gill_size <- orderNorm (data$initial_gill_size))
on.initial_gill_size <- orderNorm_initial_gill_size$x.t

# add the new object as column to your data frame
initial_gill_size.transformed <- cbind(data, on.initial_gill_size)
hist(initial_gill_size.transformed$on.initial_gill_size)

# test for normal distribution and homogeneity of variance of the transformed variable
shapiro.test(initial_gill_size.transformed$on.initial_gill_size)
#p-value = 1 > 0.05 Normal distribution given 
var.test(initial_gill_size.transformed$on.initial_gill_size ~ initial_gill_size.transformed$original_habitat)
#p-value = 0.4393 > 0.05 homogeneity of variance is given 

qqnorm(data$initial_gill_size)
qqline(data$initial_gill_size)


library(Matrix)
library(lme4)
library(lmerTest)
library(performance)
library(see)


modelperformance1<-lm(on.initial_gill_size ~ original_habitat + year, data=initial_gill_size.transformed) 
modelperformance2<-lmer(on.initial_gill_size ~ original_habitat + (1|year), data=initial_gill_size.transformed) 


AIC(modelperformance1,modelperformance2)

compare_performance(modelperformance1,modelperformance2, rank=T)


summary(modelperformance1)


# check model
check_model(modelperformance1)





####################################
#Is there a correlation between gill size and total body length across larvae from 2023 and 2024 and both habitat types?


library(readxl)
data <- read_excel("gill_body_correlation_23_24.xlsx", na="NA",  
                   col_types = c("text", "numeric", "numeric")) 
head(data)

#test for normal distribution 
shapiro.test(data$initial_body_size)


corr <- cor.test(data$initial_body_size, data$initial_gill_size, method = 'pearson')
corr

library("ggpubr")
ggscatter(data, x = "initial_body_size", y = "initial_gill_size", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", color = "#228B22", 
          xlab = "Initial body size (cm)", ylab = "Initial gill size (cm)") 







##########
# Is there a difference in the daily growth rate between pond and stream larvae from 2023 and 2024?

#import dataset
library(readxl)
data <- read_excel("body_size_increase_23_24.xlsx", na="NA",  
                   col_types = c("text","text","text","text", "text", "numeric","numeric","numeric","numeric","numeric"))

head(data)

library(dplyr)
library(car)

#set categories
as.factor(data$original_habitat)
as.factor(data$transfer_habitat)
as.factor(data$Group)
as.factor(data$Match)
as.numeric(data$nb_days_btw)
as.numeric(data$weekly)
as.numeric(data$daily)
as.numeric(data$initial_body_size)
as.numeric(data$year)

#inspect raw data 
hist(data$daily)
boxplot(data$daily ~ data$Group)
boxplot(data$daily ~ data$Match)


#test for normal distribution & homogeneity of variance
shapiro.test(data$daily)

var.test(data$daily ~ data$transfer_habitat)

var.test(data$daily ~ data$Match)


#transform data
# first, create objects with the transformations, e.g. logarithm etc.
library(bestNormalize)
(arcsinh_daily <- arcsinh_x(data$daily))
(boxcox_daily <- boxcox(data$daily))
(centerscale_daily <- center_scale(data$daily))
(orderNorm_daily <- orderNorm(data$daily))
(yeojohnson_daily <- yeojohnson(data$daily))
(sqrt_daily <- sqrt(data$daily))
(log_daily <- log(data$daily))

# then have a look at the histograms of the different transformations for first impression
par(mfrow = c(2,4)) # display all histograms in one window
hist(data$daily)
MASS::truehist(arcsinh_daily$x.t, main = "Arcsinh transformation", nbins = 12) # x.t stands for the transformed variable
MASS::truehist(boxcox_daily$x.t, main = "Box Cox transformation", nbins = 12)
MASS::truehist(centerscale_daily$x.t, main = "center_scale transformation", nbins = 12)
MASS::truehist(orderNorm_daily$x.t, main = "orderNorm transformation", nbins = 12)
MASS::truehist(yeojohnson_daily$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(sqrt_daily, main = "squareroot transformation", nbins = 12)
MASS::truehist(log_daily, main = "log transformation", nbins = 12)

# let R suggest best transformation
bn.daily<-bestNormalize(data$daily, out_of_sample = FALSE) 
bn.daily

# create an object with automatically the best transformation 
(orderNorm_daily <- orderNorm(data$daily))
on.daily <- orderNorm_daily$x.t

# add the new object as column to your data frame
daily.transformed <- cbind(data, on.daily)
hist(daily.transformed$on.daily)

# test for normal distribution and homogeneity of variance of the transformed variable
shapiro.test(daily.transformed$on.daily)
 
var.test(daily.transformed$on.daily ~ daily.transformed$transfer_habitat)


# run model and investigate model output
library(Matrix)
library(lme4)
library(lmerTest)
library(performance)
library(see)

#habitat
moddaily1<-lmer(on.daily ~ Match + initial_body_size + (1|year), data=daily.transformed)
moddaily2<-lmer(on.daily ~ Match + (1|year), data=daily.transformed)
moddaily3<-lm(on.daily ~ Match, data=daily.transformed)


AIC(moddaily1,moddaily2,moddaily3)

compare_performance(moddaily1,moddaily2,moddaily3, rank=T)

summary(moddaily2)
check_model(moddaily2)


daily<-
  ggplot(data=daily.transformed, aes(Match,daily), fill=Match)+
  geom_boxplot(aes(fill=Match), outlier.shape = NA)+
  geom_jitter(aes(fill=Match), shape=21)+
  #facet_wrap(~year)+
  scale_fill_manual(values=c("#8B864E","#228B22"))+
  stat_summary(fun=mean, shape=8, show.legend=FALSE) +
  # annotate("text", x = 1.5, y = 25, size =10, label = "*")+
  scale_x_discrete(labels=c("Match","Mismatch"))+
  theme_classic(base_size=18, base_family="Arial")+
  labs(x="Treatment group", y="Daily growth rate in total size (cm)")+
  #scale_y_continuous(breaks=seq(0,25,5), limits = c(0, 0.6))+
  theme(axis.text.x=element_text(family="Arial", size=18, color="black"), 
        axis.text.y=element_text(family="Arial", size=18, color="black"),
        legend.position="none")
dev.off()

daily

###################OR###################################

moddaily4<-lmer(on.daily ~ Group + initial_body_size + (1|year), data=daily.transformed)
moddaily5<-lmer(on.daily ~ Group + (1|year), data=daily.transformed)


AIC(moddaily4,moddaily5)

compare_performance(moddaily4,moddaily5, rank=T)


summary(moddaily5)

#posthoc test for groups
library(emmeans)
emmeans(moddaily5, list(pairwise ~ Group), adjust = "tukey")

#check model
check_model(moddaily5)


#posthoc test for groups
library(emmeans)
emmeans(moddaily5, list(pairwise ~ transfer_habitat), adjust = "tukey")


# plot as boxplot
library(plyr)
library(ggpubr)
count(data, "Group")


daily<-
  ggplot(data=daily.transformed, aes(Group,daily), fill=group)+
  geom_boxplot(aes(fill=Group), outlier.shape = NA)+
  geom_jitter(aes(fill=Group), shape=21)+
  #facet_wrap(~year)+
  scale_fill_manual(values=c("#cdc673", "#228B22", "#8B864E","#66CD00"))+
  stat_summary(fun=mean, shape=8, show.legend=FALSE) +
  # annotate("text", x = 1.5, y = 25, size =10, label = "*")+
  scale_x_discrete(labels=c("P/P","P/St", "St/P", "St/St"))+
  theme_classic(base_size=18, base_family="Arial")+
  labs(x="Treatment group", y="Daily growth rate (cm)")+
  #scale_y_continuous(breaks=seq(0,25,5), limits = c(0, 0.6))+
  theme(axis.text.x=element_text(family="Arial", size=18, color="black"), 
        axis.text.y=element_text(family="Arial", size=18, color="black"),
        legend.position="none")
dev.off()

daily






##########Is there a difference in the daily growth rate of larvae from ponds and streams in 2019?

#import dataset
library(readxl)
data <- read_excel("Body_size_increase_2019.xlsx", na="NA",  
                   col_types = c("text", "text","text","text", "text", "numeric","numeric","numeric","numeric","numeric"))

head(data)

library(dplyr)
library(car)

#set categories
as.factor(data$original_habitat)
as.factor(data$transfer_habitat)
as.factor(data$Group)
as.factor(data$Match)
as.numeric(data$nb_days_btw)
as.numeric(data$daily)
as.numeric(data$initial_body_size)


#inspect raw data 
hist(data$daily)
boxplot(data$daily ~ data$Group)
boxplot(data$daily ~ data$Match)

#test for normal distribution & homogeneity of variance
shapiro.test(data$daily)

var.test(data$daily ~ data$transfer_habitat)

var.test(data$daily ~ data$Match)


#transform data
# first, create objects with the transformations, e.g. logarithm etc.
library(bestNormalize)
(arcsinh_daily <- arcsinh_x(data$daily))
(boxcox_daily <- boxcox(data$daily))
(centerscale_daily <- center_scale(data$daily))
(orderNorm_daily <- orderNorm(data$daily))
(yeojohnson_daily <- yeojohnson(data$daily))
(sqrt_daily <- sqrt(data$daily))
(log_daily <- log(data$daily))

# then have a look at the histograms of the different transformations for first impression
par(mfrow = c(2,4)) # display all histograms in one window
hist(data$daily)
MASS::truehist(arcsinh_daily$x.t, main = "Arcsinh transformation", nbins = 12) # x.t stands for the transformed variable
MASS::truehist(boxcox_daily$x.t, main = "Box Cox transformation", nbins = 12)
MASS::truehist(centerscale_daily$x.t, main = "center_scale transformation", nbins = 12)
MASS::truehist(orderNorm_daily$x.t, main = "orderNorm transformation", nbins = 12)
MASS::truehist(yeojohnson_daily$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(sqrt_daily, main = "squareroot transformation", nbins = 12)
MASS::truehist(log_daily, main = "log transformation", nbins = 12)

# let R suggest best transformation
bn.daily<-bestNormalize(data$daily, out_of_sample = FALSE) 
bn.daily

# create an object with automatically the best transformation 
(orderNorm_daily <- orderNorm(data$daily))
on.daily <- orderNorm_daily$x.t

# add the new object as column to your data frame
daily.transformed <- cbind(data, on.daily)
hist(daily.transformed$on.daily)

# test for normal distribution and homogeneity of variance of the transformed variable
shapiro.test(daily.transformed$on.daily)

var.test(daily.transformed$on.daily ~ daily.transformed$transfer_habitat)


# run model and investigate model output
library(Matrix)
library(lme4)
library(lmerTest)
library(performance)
library(see)

moddaily10<-lm(on.daily ~ Group + initial_body_size, data=daily.transformed)
moddaily20<-lm(on.daily ~ Group, data=daily.transformed)


AIC(moddaily10,moddaily20)

compare_performance(moddaily10,moddaily20, rank=T)

summary(moddaily10)


#check model
check_model(moddaily10)


# plot as boxplot
library(plyr)
library(ggpubr)
count(data, "Group")


daily<-
  ggplot(data=daily.transformed, aes(Group,daily), fill=group)+
  geom_boxplot(aes(fill=Group), outlier.shape = NA)+
  geom_jitter(aes(fill=Group), shape=21)+
  #facet_wrap(~year)+
  scale_fill_manual(values=c("#cdc673", "#228B22", "#8B864E","#66CD00"))+
  stat_summary(fun=mean, shape=8, show.legend=FALSE) +
  # annotate("text", x = 1.5, y = 25, size =10, label = "*")+
  scale_x_discrete(labels=c("P/P","P/St", "St/P", "St/St"))+
  theme_classic(base_size=18, base_family="Arial")+
  labs(x="Treatment group", y="Daily growth rate (cm)")+
  #scale_y_continuous(breaks=seq(0,25,5), limits = c(0, 0.6))+
  theme(axis.text.x=element_text(family="Arial", size=18, color="black"), 
        axis.text.y=element_text(family="Arial", size=18, color="black"),
        legend.position="none")
dev.off()

daily






##########
# Is there a difference in the gill size change between pond and stream larvae from 2023 and 2024?

#import dataset
library(readxl)
data <- read_excel("gill_size_change_23_24.xlsx", na="NA",  
                   col_types = c("text","text","text","text", "text", "numeric","numeric","numeric","numeric","numeric"))

head(data)

library(dplyr)
library(car)

#set categories
as.factor(data$original_habitat)
as.factor(data$transfer_habitat)
as.factor(data$Group)
as.factor(data$Match)
as.numeric(data$nb_days_btw)
as.numeric(data$weekly)
as.numeric(data$daily)
as.numeric(data$initial_gill_size)
as.numeric(data$year)

#inspect raw data 
hist(data$daily)
boxplot(data$daily ~ data$Group)
boxplot(data$daily ~ data$Match)


#test for normal distribution & homogeneity of variance
shapiro.test(data$daily)

var.test(data$daily ~ data$transfer_habitat)

var.test(data$daily ~ data$Match)


#transform data
# first, create objects with the transformations, e.g. logarithm etc.
library(bestNormalize)
(arcsinh_daily <- arcsinh_x(data$daily))
(boxcox_daily <- boxcox(data$daily))
(centerscale_daily <- center_scale(data$daily))
(orderNorm_daily <- orderNorm(data$daily))
(yeojohnson_daily <- yeojohnson(data$daily))
(sqrt_daily <- sqrt(data$daily))
(log_daily <- log(data$daily))

# then have a look at the histograms of the different transformations for first impression
par(mfrow = c(2,4)) # display all histograms in one window
hist(data$daily)
MASS::truehist(arcsinh_daily$x.t, main = "Arcsinh transformation", nbins = 12) # x.t stands for the transformed variable
MASS::truehist(boxcox_daily$x.t, main = "Box Cox transformation", nbins = 12)
MASS::truehist(centerscale_daily$x.t, main = "center_scale transformation", nbins = 12)
MASS::truehist(orderNorm_daily$x.t, main = "orderNorm transformation", nbins = 12)
MASS::truehist(yeojohnson_daily$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(sqrt_daily, main = "squareroot transformation", nbins = 12)
MASS::truehist(log_daily, main = "log transformation", nbins = 12)

# let R suggest best transformation
bn.daily<-bestNormalize(data$daily, out_of_sample = FALSE) 
bn.daily

# create an object with automatically the best transformation 
(orderNorm_daily <- orderNorm(data$daily))
on.daily <- orderNorm_daily$x.t

# add the new object as column to your data frame
daily.transformed <- cbind(data, on.daily)
hist(daily.transformed$on.daily)

# test for normal distribution and homogeneity of variance of the transformed variable
shapiro.test(daily.transformed$on.daily)

var.test(daily.transformed$on.daily ~ daily.transformed$transfer_habitat)

# run model and investigate model output
library(Matrix)
library(lme4)
library(lmerTest)
library(performance)
library(see)

#habitat
moddaily1<-lmer(on.daily ~ Match + initial_gill_size + (1|year), data=daily.transformed)
moddaily2<-lmer(on.daily ~ Match + (1|year), data=daily.transformed)
moddaily3<-lm(on.daily ~ Match, data=daily.transformed)


AIC(moddaily1,moddaily2,moddaily3)

compare_performance(moddaily1,moddaily2,moddaily3, rank=T)
 

summary(moddaily1)
check_model(moddaily1)


daily<-
  ggplot(data=daily.transformed, aes(Match,daily), fill=Match)+
  geom_boxplot(aes(fill=Match), outlier.shape = NA)+
  geom_jitter(aes(fill=Match), shape=21)+
  #facet_wrap(~year)+
  scale_fill_manual(values=c("#8B864E","#228B22"))+
  stat_summary(fun=mean, shape=8, show.legend=FALSE) +
  # annotate("text", x = 1.5, y = 25, size =10, label = "*")+
  scale_x_discrete(labels=c("Match","Mismatch"))+
  theme_classic(base_size=18, base_family="Arial")+
  labs(x="Treatment group", y="Daily growth rate in total size (cm)")+
  #scale_y_continuous(breaks=seq(0,25,5), limits = c(0, 0.6))+
  theme(axis.text.x=element_text(family="Arial", size=18, color="black"), 
        axis.text.y=element_text(family="Arial", size=18, color="black"),
        legend.position="none")
dev.off()

daily



###################OR###################################



moddaily4<-lmer(on.daily ~ Group + initial_gill_size + (1|year), data=daily.transformed)
moddaily5<-lmer(on.daily ~ Group + (1|year), data=daily.transformed)


AIC(moddaily4,moddaily5)

compare_performance(moddaily4,moddaily5, rank=T)


summary(moddaily4)

#posthoc test for groups
library(emmeans)
emmeans(moddaily5, list(pairwise ~ Group), adjust = "tukey")

#check model
check_model(moddaily4)


#posthoc test for groups
library(emmeans)
emmeans(moddaily5, list(pairwise ~ transfer_habitat), adjust = "tukey")


# plot as boxplot
library(plyr)
library(ggpubr)
count(data, "Group")


daily<-
  ggplot(data=daily.transformed, aes(Group,daily), fill=group)+
  geom_boxplot(aes(fill=Group), outlier.shape = NA)+
  geom_jitter(aes(fill=Group), shape=21)+
  #facet_wrap(~year)+
  scale_fill_manual(values=c("#cdc673", "#228B22", "#8B864E","#66CD00"))+
  stat_summary(fun=mean, shape=8, show.legend=FALSE) +
  # annotate("text", x = 1.5, y = 25, size =10, label = "*")+
  scale_x_discrete(labels=c("P/P","P/St", "St/P", "St/St"))+
  theme_classic(base_size=18, base_family="Arial")+
  labs(x="Treatment group", y="Daily change in gill size (cm)")+
  #scale_y_continuous(breaks=seq(0,25,5), limits = c(0, 0.6))+
  theme(axis.text.x=element_text(family="Arial", size=18, color="black"), 
        axis.text.y=element_text(family="Arial", size=18, color="black"),
        legend.position="none")
dev.off()

daily



#############
# Can pond larvae be recorded over a longer time than stream larvae? 

#import dataset
library(readxl)
setwd("D:/Uni_Backup_20250308/fire_salamander/12_transfer_experiment/2024/Statistics")
data <- read_excel("body_size_increase_23_24.xlsx", na="NA",  
                   col_types = c("text", "text", "text","text","text","numeric", "numeric", "numeric", "numeric", "numeric")) 
head(data)

#inspect raw data 
hist(data$nb_days_btw)
boxplot(data$nb_days_btw ~ data$transfer_habitat)
boxplot(data$nb_days_btw ~ data$Group)


#test for normal distribution & homogeneity of variance
shapiro.test(data$nb_days_btw)

var.test(data$nb_days_btw ~ data$transfer_habitat)


#Mann Whitney U test to test differences between habitat types
install.packages("psych")
library(psych)
wilcox.test(nb_days_btw~transfer_habitat, data = data, 
            exact = FALSE, 
            correct = FALSE, 
            conf.int = FALSE)



#Kruskal-Wallis test to test differences between Groups
kruskal.test(data$nb_days_btw~data$Group)


library(rstatix)
dunn_test(nb_days_btw~Group, data=data, 
          p.adjust.method = "bonferroni")


library("ggpubr")

data %>%
  ggplot(aes(x=nb_days_btw, 
             y=daily, 
             color=Group))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Time between first and last capture (days)", y="Daily growth rate (cm)")+
  scale_color_manual(values=c("#8B864E", "#228B22", "#cdc673", "#66cd00")) +
  theme(axis.text.x=element_text(family="Arial", size = 24, color="black"), 
        axis.text.y=element_text(family="Arial", size = 24, color="black"),
        axis.title=element_text(size = 24),
        legend.text=element_text(size = 24),
        legend.title=element_text(size = 24),
        legend.key.size = unit(1,'cm'))
ggsave("scatterplot_Group.png")



#transform data
# first, create objects with the transformations, e.g. logarithm etc.
library(bestNormalize)
(arcsinh_nb_days_btw <- arcsinh_x(data$nb_days_btw))
(boxcox_nb_days_btw <- boxcox(data$nb_days_btw))
(centerscale_nb_days_btw <- center_scale(data$nb_days_btw))
(orderNorm_nb_days_btw <- orderNorm(data$nb_days_btw))
(yeojohnson_nb_days_btw <- yeojohnson(data$nb_days_btw))
(sqrt_nb_days_btw <- sqrt(data$nb_days_btw))
(log_nb_days_btw <- log(data$nb_days_btw))

# then have a look at the histograms of the different transformations for first impression
par(mfrow = c(2,4)) # display all histograms in one window
hist(data$nb_days_btw)
MASS::truehist(arcsinh_nb_days_btw$x.t, main = "Arcsinh transformation", nbins = 12) # x.t stands for the transformed variable
MASS::truehist(boxcox_nb_days_btw$x.t, main = "Box Cox transformation", nbins = 12)
MASS::truehist(centerscale_nb_days_btw$x.t, main = "center_scale transformation", nbins = 12)
MASS::truehist(orderNorm_nb_days_btw$x.t, main = "orderNorm transformation", nbins = 12)
MASS::truehist(yeojohnson_nb_days_btw$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(sqrt_nb_days_btw, main = "squareroot transformation", nbins = 12)
MASS::truehist(log_nb_days_btw, main = "log transformation", nbins = 12)

# let R suggest best transformation
bn.nb_days_btw<-bestNormalize(data$nb_days_btw, out_of_sample = FALSE) 
bn.nb_days_btw

# create an object with automatically the best transformation 
(orderNorm_nb_days_btw <- orderNorm(data$nb_days_btw))
on.nb_days_btw <- orderNorm_nb_days_btw$x.t

# add the new object as column to your data frame
nb_days_btw.transformed <- cbind(data, on.nb_days_btw)
hist(nb_days_btw.transformed$on.nb_days_btw)

# test for normal distribution and homogeneity of variance of the transformed variable
shapiro.test(nb_days_btw.transformed$on.nb_days_btw)
#p-value = 1 > 0.05 Normal distribution given 
var.test(nb_days_btw.transformed$on.nb_days_btw ~ nb_days_btw.transformed$transfer_habitat)
#p-value = 0.4288 > 0.05 homogeneity of variance is given 

# run model and investigate model output
library(Matrix)
library(lme4)
library(lmerTest)
library(performance)
library(see)

#number of days between first and last capture vs habitat of transfer
modnb_days_btw1<-lm(on.nb_days_btw ~ Group, data=nb_days_btw.transformed)

summary(modnb_days_btw1)
check_model(modnb_days_btw1)

library(emmeans)
emmeans(modnb_days_btw1, list(pairwise ~ Group), adjust = "tukey")

ANCOVA(modnb_days_btw1) 





#############
# Do the abiotic parameters differ between the habitat types? 

#import dataset
library(readxl)
data <- read_excel("Abiotic_factors_20250319.xlsx", na="NA",  
                   col_types = c("text", "text", "date", "numeric", "numeric", "numeric", "numeric", "text")) 
head(data)

as.factor(data$Year)

#########Watertemp
#inspect raw data 
hist(data$Watertemp)
boxplot(data$Watertemp ~ data$Habitat_type)


#test for normal distribution & homogeneity of variance
shapiro.test(data$Watertemp)

var.test(data$Watertemp ~ data$Habitat_type)


library(Matrix)
library(lme4)
library(lmerTest)
library(performance)
library(see)


modelwatertemp1<-lmer(Watertemp ~ Habitat_type + (1|Year), data=data) 
modelwatertemp2<-lm(Watertemp ~ Habitat_type + Year, data=data) 

AIC(modelwatertemp1,modelwatertemp2)

compare_performance(modelwatertemp1,modelwatertemp2, rank=T)


summary(modelwatertemp2)


# check model
check_model(modelwatertemp2)


library("ggpubr")
ggscatter(data, x = "Session", y = "Watertemp",
         add = "reg.line", conf.int = TRUE,
         cor.coef = TRUE, cor.method = "pearson", color = "#228B22",
         xlab = "Date", ylab = "Watertemperature (°C)") +
  facet_wrap(~Year)


data %>%
  ggplot(aes(x=Session, 
             y=Watertemp, 
             color=Habitat_type))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Number of capture session", y="Water temperature (°C)")+
  scale_color_manual(values=c("#8B864E", "#228B22")) +
    theme(axis.text.x=element_text(family="Arial", size = 24, color="black"), 
        axis.text.y=element_text(family="Arial", size = 24, color="black"),
        axis.title=element_text(size = 24),
        legend.text=element_text(size = 24),
        legend.title=element_text(size = 24),
        legend.key.size = unit(1,'cm'))
ggsave("scatterplot_Group.png")





#########O2
#inspect raw data 
hist(data$O2)
boxplot(data$O2 ~ data$Habitat_type)


#test for normal distribution & homogeneity of variance
shapiro.test(data$O2)

var.test(data$O2 ~ data$Habitat_type)


#transform data
# first, create objects with the transformations, e.g. logarithm etc.
library(bestNormalize)
(arcsinh_O2 <- arcsinh_x(data$O2))
(boxcox_O2 <- boxcox(data$O2))
(centerscale_O2 <- center_scale(data$O2))
(orderNorm_O2 <- orderNorm(data$O2))
(yeojohnson_O2 <- yeojohnson(data$O2))
(sqrt_O2 <- sqrt(data$O2))
(log_O2 <- log(data$O2))

# then have a look at the histograms of the different transformations for first impression
par(mfrow = c(2,4)) # display all histograms in one window
hist(data$O2)
MASS::truehist(arcsinh_O2$x.t, main = "Arcsinh transformation", nbins = 12) # x.t stands for the transformed variable
MASS::truehist(boxcox_O2$x.t, main = "Box Cox transformation", nbins = 12)
MASS::truehist(centerscale_O2$x.t, main = "center_scale transformation", nbins = 12)
MASS::truehist(orderNorm_O2$x.t, main = "orderNorm transformation", nbins = 12)
MASS::truehist(yeojohnson_O2$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(sqrt_O2, main = "squareroot transformation", nbins = 12)
MASS::truehist(log_O2, main = "log transformation", nbins = 12)

# let R suggest best transformation
bn.O2<-bestNormalize(data$O2, out_of_sample = FALSE) 
bn.O2

# create an object with automatically the best transformation 
(orderNorm_O2 <- orderNorm(data$O2))
on.O2 <- orderNorm_O2$x.t

# add the new object as column to your data frame
O2.transformed <- cbind(data, on.O2)
hist(O2.transformed$on.O2)

# test for normal distribution and homogeneity of variance of the transformed variable
shapiro.test(O2.transformed$on.O2)

var.test(O2.transformed$on.O2 ~ O2.transformed$Habitat_type)


library(Matrix)
library(lme4)
library(lmerTest)
library(performance)
library(see)

modelO21<-lmer(on.O2 ~ Habitat_type + (1|Year), data=O2.transformed) 
modelO22<-lm(on.O2 ~ Habitat_type + Year, data=O2.transformed) 

AIC(modelO21,modelO22)

compare_performance(modelO21,modelO22, rank=T)


summary(modelO22)



data %>%
  ggplot(aes(x=Session, 
             y=O2, 
             color=Habitat_type))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Number of capture session", y="Dissolved oxygen (%)")+
  scale_color_manual(values=c("#8B864E", "#228B22")) +
  theme(axis.text.x=element_text(family="Arial", size = 24, color="black"), 
        axis.text.y=element_text(family="Arial", size = 24, color="black"),
        axis.title=element_text(size = 24),
        legend.text=element_text(size = 24),
        legend.title=element_text(size = 24),
        legend.key.size = unit(1,'cm'))
ggsave("scatterplot_Group.png")






######################################################

#Plot total body length over sessions per group for both years


# Paket laden
library(ggplot2)

# at which sessions do we have measuremnts per group? and what are the total sizes 
data <- data.frame(
  Group = c(rep("P/P", 9), rep("P/St", 8), rep("St/P", 8), rep("St/St", 8)),
  Ereignis = c(1, 2, 3, 4, 5, 6, 7, 8, 9,   # P/P
               1, 2, 3, 4, 5, 6, 7, 8, # P/St
               1, 3, 4, 5, 6, 7, 8, 9, # St/P
               1, 2, 3, 4, 5, 6, 7, 8), #St/St     
  AbsoluteGroesse = c(3.98, 4.03, 4.43, 3.9, 4.3, 4.94, 4.7, 4.9, 5.2,
                      3.57, 3.56, 3.73, 3.6, 3.98, 3.93, 3.93, 4.85,
                      3.54, 3.63, 3.35, 3.63, 4, 4.56, 4.48, 4.7,
                      3.49, 3.13, 3.90, 3.73, 3.94, 3.95, 4.8, 4.7)
)


print(data)

#Plot
ggplot(data, aes(x = Ereignis, y = AbsoluteGroesse, color = Group, group = Group)) +
  geom_line(size = 1) +     # connect lines 
  geom_point(size = 3) +    # shows single measurements
  scale_x_continuous(breaks = 1:9) +  # x achsis with sessions 1-9 
  scale_color_manual(values = c("P/P" = "#8B864E", "P/St" = "#228B22", "St/P" = "#cdc673", "St/St" = "#66cd00")) +
  labs(
    x = "Session",
    y = "Total body length (cm)"
  ) +
  theme_minimal()

