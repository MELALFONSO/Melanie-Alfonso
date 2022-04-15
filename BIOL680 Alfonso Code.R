cell.data=read.csv("FISH_cell_intensity_scramble.csv")
View(cell.data)
cell.data[,"HSM_flask"]=factor(cell.data[,"HSM_flask"],levels=1:3, ordered=FALSE)
str(cell.data,give.attr=FALSE)
library(ggplot2)
ggplot(cell.data, aes(x = probe_name, y = avg_intensity, group = HSM_flask, col = HSM_flask)) + 
  geom_point() + stat_summary(fun = mean, geom = "line") + theme_bw() + xlab("gene-specific FISH probe")+ylab("Average cell signal intensity")
with(cell.data, interaction.plot(x.factor = probe_name, 
                                trace.factor = HSM_flask, 
                                response = avg_intensity))
options(contrasts = c("contr.treatment", "contr.poly"))
library(lmerTest)
fit.cell.data = lmer(avg_intensity ~ probe_name + (probe_name| HSM_flask), data = cell.data)
anova(fit.cell.data)
#fixed effect of probe_name is significant
summary(fit.cell.data)
fixef(fit.cell.data)
#purely fixed effects model below
fit.cell.data.aov <- aov(avg_intensity ~ probe_name * HSM_flask, data = cell.data)
summary(fit.cell.data.aov)
#check assumptions
plot(fit.cell.data,which=3,caption="",col="black")
qqnorm(resid(fit.cell.data))
qqline(resid(fit.cell.data))
shapiro.test(residuals(fit.cell.data))
library(car)
leveneTest(residuals(fit.cell.data) ~ cell.data$HSM_flask)
