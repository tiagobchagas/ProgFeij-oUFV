#install.packages("emmeans")
library(emmeans)
library(dplyr)

setwd("/Users/tiagobchagas/Desktop/Qualificação")
list.files()
carioca = read.csv("conjuntacarioca.csv", header = T, sep = ";")

carioca <- carioca %>% 
  mutate_at(vars(Trials, Blocks, Rep, Treat, Families), as.factor)

carioca2019 = carioca %>% filter(Trials == "LS19")
str(carioca2019)
carioca2020 = carioca %>% filter(Trials == "LT20")
carioca2021 = carioca %>% filter(Trials == "LT21")
carioca2022 = carioca %>% filter(Trials == "LT22")
#install.packages("agricolae")

library(agricolae)
#install.packages("MASS")
library(MASS)
#install.packages("nlme")
library(nlme)
#X<-carioca # Leitura dos dados X



library(lmerTest)
m1f <- lm(PROD ~ Families + Rep + Blocks, data = carioca2019)
anova(m1f)
m1 <- lmer(PROD ~ Families + Rep + (1 | Blocks), data = carioca2019)
anova(m1)

#summary(glht(fit.ibd2, linfct = mcp(Treat = c(1, 0, 0, 0, 0, -1))))

print(m1)
emm1 <- emmeans(m1, specs = ~Families)
#v <- c("emmean", "SE")
#bind_cols( as.data.frame(emm0)[, v] %>% rename_all(paste0, "_m0"))          
#summary(fit.ibd2)

#str(emm0)
# prod2019 = print(emm0)

#system("say Chama, Acunha!")
#library(notifier)
#notify(  title = "Terminou macho",  msg = c("Chama","Acunha"))

m2 <- lmer(PROD ~ Families + Rep + (1 | Blocks), data = carioca2020)
m2f <- lm(PROD ~ Families + Rep +  Blocks, data = carioca2020)
m2r <- lmer(PROD ~ (1| Families) + Rep + (1 | Blocks), data = carioca2020)
lrt
m2r
anova(m2r)
plot(m2f)
hist(m2f$residuals)
anova(m2f)
levels(carioca2020$Blocks)
m2 <- lmer(PROD ~ Families + (1 | Blocks), data = carioca2020)
emm2 <- emmeans(m2, specs = ~Families)
prod2020 = print(emm2)



m3 <- lmer(AG ~ Families + Rep + (1 |Blocks), data = carioca2020)
emm3 <- emmeans(m3, specs = ~Families)
ag2020 = print(emm3)


m4 <- lmer(PROD ~ Families + (1 | Blocks), data = carioca2021)
emm4 <- emmeans(m4, specs = ~Families)
prod2021 = print(emm4)


m5 <- lmer(AG ~ Families + (1 | Blocks), data = carioca2021)
emm5 <- emmeans(m5, specs = ~Families)
ag2021 = print(emm5)

m6 <- lmer(Arq ~ Families + (1 | Blocks), data = carioca2021)
emm6 <- emmeans(m6, specs = ~Families)
arq2021 = print(emm6)

m7 <- lmer(PROD ~ Families + (1 | Blocks), data = carioca2022)
emm7 <- emmeans(m7, specs = ~Families)
prod2022 = print(emm7)

m8 <- lmer(AG ~ Families + (1 | Blocks), data = carioca2022)
emm8 <- emmeans(m8, specs = ~Families)
ag2022 = print(emm8)

m9 <- lmer(Arq ~ Families + (1 | Blocks), data = carioca2022)
emm9 <- emmeans(m9, specs = ~Families)
arq2022 = print(emm9)
system("say Chama, Acunha!")
#renomear 
colnames(prod2019)[2] = "prod2019"
colnames(prod2020)[2] = "prod2020"
colnames(ag2020)[2] = "ag2020"
colnames(prod2021)[2] = "prod2021"
colnames(ag2021)[2] = "ag2021"
colnames(arq2021)[2] = "arq2021"
colnames(prod2022)[2] = "prod2022"
colnames(ag2022)[2] = "ag2022"
colnames(arq2022)[2] = "arq2022"
#filtrar as famílias anteriormente selecionadas nas safras anteriores
MeanProd2019= select(filter(carioca19, Families %in% c(target)), c(Families,prod2019))
MeanProd2020= select(filter(prod2020, Families %in% c(target)), c(Families,prod2020))
MeanAg2020= select(filter(ag2020, Families %in% c(target)), c(Families,ag2020))
MeanProd2021= select(filter(prod2021, Families %in% c(target)), c(Families,prod2021))
MeanAg2021= select(filter(ag2021, Families %in% c(target)), c(Families,ag2021))
MeanArq2021= select(filter(arq2021, Families %in% c(target)), c(Families,arq2021))
MeanProd2022= select(filter(prod2022, Families %in% c(target)), c(Families,prod2022))
MeanAg2022= select(filter(ag2022, Families %in% c(target)), c(Families,ag2022))
MeanArq2022= select(filter(arq2022, Families %in% c(target)), c(Families,arq2022))


#juntas os arquivos
carioca19=MeanProd2019
carioca20=bind_cols(MeanProd2020,MeanAg2020)
carioca21=bind_cols(MeanProd2021,MeanAg2021,MeanArq2021)
carioca22=bind_cols(MeanProd2022,MeanAg2022,MeanArq2022)



cariocameans=bind_rows(list(carioca19,carioca20, carioca21, carioca22), .id = "id")

target <- c(prod2022[,1])

#filter(carioca19, Families %in% target)  # equivalently, dat %>% filter(name %in% target)
#filter(df,date=='Sunday'| date=='Monday')
#target == carioca19["Families"]
#prod2019= select(filter(carioca19, Families %in% c(target)), c(Families,prod2019))
## Merging Shapefile informations
#setwd("/Users/tiagobchagas/Desktop/HTP CREA 2022/")
#dir()
#results <- read_excel("Data prova Eduardo.xlsx")
#Index <- cbind(Red["Red"],Green["Green"], Blue["Blue"], NGRDI["NGRDI"])
#df = merge(x =cariocameans, y = carioca["Treat"], by = "Families")
#df

#data_all <- merge(carioca,          # Properly merge data
#                cariocameans,
#               by.x = "Treat",
#              by.y = "Families*")
#data_all                       

#results <- bind_cols(carioca,cariocameans)
#head(results)
#tail(results)

## Write Excel Table
#install.packages("openxlsx")
library(openxlsx)
#setwd("/Users/tiagobchagas/Desktop/HTP CREA 2022/")
write.xlsx(cariocameans, "CariocaMeans.xlsx", sheetName = "MEANS", 
           colNames = TRUE, rowNnames = F, append = FALSE)

