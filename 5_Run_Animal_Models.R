#
# Animal models for meiosis traits
#
#
#

library(GenABEL)
library(magrittr)
library(plyr)
library(ggplot2)
library(asreml)
library(reshape2)

load("data/GRM.RData")

source("r/ASReml4.EstEffects.R")

attr(grminv, which = "INVERSE") <- TRUE

recsumm <- read.table("results/3_20200721_Individual_Recombination_Rates.txt", header = T, stringsAsFactors = F)

recsumm$RRID <- factor(recsumm$RRID)
recsumm$RRID.SEX <- factor(recsumm$RRID.SEX)

recsumm <- subset(recsumm, !is.na(RRID.SEX) & !is.na(RRID.Fhat3))

recsumm.m <- subset(recsumm, RRID.SEX == "M") %>% droplevels
recsumm.f <- subset(recsumm, RRID.SEX == "F") %>% droplevels

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Plot the Data                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(recsumm, aes(TotalRecombCount)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~RRID.SEX, ncol = 2) +
  ggtitle("Total Crossover Count")


ggplot(recsumm, aes(RRID.Age.2, TotalRecombCount, col = RRID.SEX)) + 
  geom_point(alpha = 0.2) + 
  stat_smooth() +
  labs(x = "Age", y = "Autosomal Crossover Count", col = "Sex") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, 16, 2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. How is the heritability?                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

anmod.acc.all <- asreml(fixed = TotalRecombCount ~ RRID.SEX + RRID.Fhat3,
                        random = ~ vm(RRID, grminv) + ide(RRID),
                        data = recsumm,
                        residual = ~idv(units))

asreml4pin(anmod.acc.all)
summary.asreml(anmod.acc.all, coef = T)$coef.fixed # Fixed effects
wald.asreml(anmod.acc.all)


anmod.acc.males <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                          random = ~ vm(RRID, grminv) + ide(RRID),
                          data = recsumm.m,
                          residual = ~idv(units))


asreml4pin(anmod.acc.males)
summary.asreml(anmod.acc.males, coef = T)$coef.fixed # Fixed effects
wald.asreml(anmod.acc.males)

anmod.acc.females <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                            random = ~ vm(RRID, grminv) + ide(RRID),
                            data = recsumm.f,
                            residual = ~idv(units))

asreml4pin(anmod.acc.females)
summary.asreml(anmod.acc.females, coef = T)$coef.fixed # Fixed effects
wald.asreml(anmod.acc.females)

##~~ Test significance of heritability


anmod.acc.all.temp <- asreml(fixed = TotalRecombCount ~ RRID.SEX + RRID.Fhat3,
                        random = ~ RRID,
                        data = recsumm,
                        residual = ~idv(units))

anmod.acc.males.temp <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                          random = ~ RRID,
                          data = recsumm.m,
                          residual = ~idv(units))

anmod.acc.females.temp <- asreml(fixed = TotalRecombCount ~ RRID.Fhat3,
                            random = ~ RRID,
                            data = recsumm.f,
                            residual = ~idv(units))

2*(anmod.acc.all$loglik - anmod.acc.all.temp$loglik)
2*(anmod.acc.males$loglik - anmod.acc.males.temp$loglik)
2*(anmod.acc.females$loglik - anmod.acc.females.temp$loglik)

pchisq(2*(anmod.acc.all$loglik - anmod.acc.all.temp$loglik), 1, lower.tail = F)
pchisq(2*(anmod.acc.males$loglik - anmod.acc.males.temp$loglik), 1, lower.tail = F)
pchisq(2*(anmod.acc.females$loglik - anmod.acc.females.temp$loglik), 1, lower.tail = F)

#~~ Get the model information

vartab <- rbind(cbind(Model = "Both Sexes", vpredict(anmod.acc.all, VP ~ V1+V2+V3)) ,
                cbind(Model = "Males"     ,vpredict(anmod.acc.males, VP ~ V1+V2+V3)) ,
                cbind(Model = "Females"   , vpredict(anmod.acc.females, VP ~ V1+V2+V3)))

names(vartab) <- c("Model", "Vp", "Vp.se")

vartab$Vobs = c(var(recsumm$TotalRecombCount, na.rm = T),
                var(recsumm.m$TotalRecombCount, na.rm = T),
                var(recsumm.f$TotalRecombCount, na.rm = T))
vartab$Meanobs = c(mean(recsumm$TotalRecombCount, na.rm = T),
                   mean(recsumm.m$TotalRecombCount, na.rm = T),
                   mean(recsumm.f$TotalRecombCount, na.rm = T))

vartab$N <- c(length(recsumm$TotalRecombCount),
              length(recsumm.m$TotalRecombCount),
              length(recsumm.f$TotalRecombCount))

vartab$Nids <- c(recsumm$RRID %>% unique %>% length,
                 recsumm.m$RRID %>% unique %>% length,
                 recsumm.f$RRID %>% unique %>% length)

vartab$Nxovers <- c(sum(recsumm$TotalRecombCount),
                    sum(recsumm.m$TotalRecombCount),
                    sum(recsumm.f$TotalRecombCount))




#~~ Format the fixed effects

fixef.tab <- rbind(cbind(Model = "Both Sexes", summary.asreml(anmod.acc.all, coef = T)$coef.fixed,
                         variable = row.names(summary.asreml(anmod.acc.all, coef = T)$coef.fixed)),
               cbind(Model = "Males"     , summary.asreml(anmod.acc.males, coef = T)$coef.fixed,
                     variable = row.names(summary.asreml(anmod.acc.males, coef = T)$coef.fixed)),
               cbind(Model = "Females"   , summary.asreml(anmod.acc.females, coef = T)$coef.fixed,
                     variable = row.names(summary.asreml(anmod.acc.females, coef = T)$coef.fixed))) %>% data.frame
for(i in 2:4) fixef.tab[,i] <- as.numeric(as.character(fixef.tab[,i]))

# #~~ Format the Random Effects

ranef.tab <- rbind(cbind(Model = "Both Sexes", Trait = "TotalRecombCount", asreml4pin(anmod.acc.all)),
                   cbind(Model = "Males"     , Trait = "TotalRecombCount", asreml4pin(anmod.acc.males)),
                   cbind(Model = "Females"   , Trait = "TotalRecombCount", asreml4pin(anmod.acc.females)))


head(ranef.tab)
ranef.tab$Trait <- NULL
ranef.tab$variable <- rep(c("h2", "pe2", "e2"), times = 3)
ranef.tab$Effect <- NULL
ranef.tab <- ranef.tab[,c(1, 4, 2, 3)]

temp1 <- dcast(ranef.tab, Model ~ variable, value.var = c("Estimate"))
temp2 <- dcast(ranef.tab, Model ~ variable, value.var = c("SE"))

names(temp2)[2:4] <- paste0(names(temp2)[2:4], ".se")

ranef.tab <- join(temp1, temp2)
ranef.tab <- subset(ranef.tab, select = c(Model, h2, h2.se, pe2, pe2.se, e2, e2.se))

vartab <- join(vartab, ranef.tab)

#~~ Wald tests

wald.res <- rbind(cbind(Model = "Both Sexes", wald.asreml(anmod.acc.all)),
                  cbind(Model = "Males", wald.asreml(anmod.acc.males)),
                  cbind(Model = "Females", wald.asreml(anmod.acc.females))) %>% data.frame

for(i in 3:5) wald.res[,i] <- as.numeric(as.character(wald.res[,i]))
wald.res$variable <- c("(Intercept)", "RRID.SEX", "RRID.Fhat3", "residual", 
                       "(Intercept)", "RRID.Fhat3", "residual",
                       "(Intercept)", "RRID.Fhat3", "residual")
wald.res <- subset(wald.res, variable != "residual")


write.table(vartab, "results/5_Animal_Model_Results.txt", row.names = F, sep = "\t", quote = F)
write.table(vartab, "results/5_Animal_Model_Results_LaTeX.txt", row.names = F, sep = " & ", quote = F)

write.table(fixef.tab, "results/5_Animal_Model_Results_Fixef.txt", row.names = F, sep = "\t", quote = F)
write.table(wald.res, "results/5_Animal_Model_Results_Wald.txt", row.names = F, sep = "\t", quote = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Bivariate Animal Models                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Create female and male specific columns to allow the model to run:

recsumm$Male.RR   <- ifelse(recsumm$RRID.SEX == "M"  , recsumm$TotalRecombCount, NA)
recsumm$Female.RR <- ifelse(recsumm$RRID.SEX == "F", recsumm$TotalRecombCount, NA)

#~~ Run the raw model.

init_test <- c(0.5, 1, 1)
names(init_test) <- c("U", "U", "U")

grm.summ.bivar.raw <- asreml(fixed  = cbind(Female.RR, Male.RR) ~ trait + trait:RRID.Fhat3, 
                             random = ~ corgh(trait, init = init_test):vm(RRID, grminv),
                             residual   = ~ units:us(trait),
                             data = recsumm,
                             maxit = 20, workspace = 500e+6, pworkspace = 500e+6)

summary(grm.summ.bivar.raw)



init_test <- c(0, 1, 1)
names(init_test) <- c("F", "U", "U")

grm.summ.bivar.ra0 <- asreml(fixed  = cbind(Female.RR, Male.RR) ~ trait + trait:RRID.Fhat3, 
                        random = ~ corgh(trait, init = init_test):vm(RRID, grminv),
                        residual   = ~ units:us(trait),
                        data = recsumm,
                        maxit = 20)

summary(grm.summ.bivar.ra0)

2*(grm.summ.bivar.raw$loglik - grm.summ.bivar.ra0$loglik)
pchisq(2*(grm.summ.bivar.raw$loglik - grm.summ.bivar.ra0$loglik), 1, lower.tail = F)


init_test <- c(0.999, 1, 1)
names(init_test) <- c("F", "U", "U")

grm.summ.bivar.ra1 <- asreml(fixed  = cbind(Female.RR, Male.RR) ~ trait + trait:RRID.Fhat3, 
                             random = ~ corgh(trait, init = init_test):vm(RRID, grminv),
                             residual   = ~ units:us(trait),
                             data = recsumm,
                             maxit = 20, workspace = 500e+6, pworkspace = 500e+6)

summary(grm.summ.bivar.ra1)

2*(grm.summ.bivar.raw$loglik - grm.summ.bivar.ra1$loglik)
pchisq(2*(grm.summ.bivar.raw$loglik - grm.summ.bivar.ra1$loglik), 1, lower.tail = F)

