library(data.table)
library(performance)
library(marginaleffects)
library(plyr)
library(lme4)
library(ggplot2)
library(readxl)
library(car)

data <- read_excel("Procrustes_residuals_default_vs_custom_chat.xlsx")
setDT(data)

data[ , method := factor(sub(" default vs.*", "", x = Comparison)) ] 
data$isolation_source <- factor(data$isolation_source)
data$sampleID <- as.numeric( factor(data$sample) )

table(data$method)
table(data$isolation_source)
table(data$sample)

mod1 <- lmer(residual ~ method * isolation_source + (1|sampleID), data = data)
mod2 <- lmer(residual ~ method + isolation_source + (1|sampleID), data = data)

compare_performance(mod1, mod2)
# mod1 AICc = -9.0, R2 = 0.611   --> better model
# mod2 AICc = -2.0, R2 = 0.426

test_performance(mod1, mod2)  # p < 0.001 --> models are not equivalent

car::Anova(mod1, test.statistic = "F")
# -> Interaction is significant


predss <- predictions(
    model = mod1,
    newdata = datagrid(
        isolation_source = unique,
        method = unique),
    re.form = NA)

print(names(predss))
print(head(predss))
preds <- as.data.table(predss)
print(names(preds))
#preds <- as.data.table(predss[ c(11:12, 2:7) ])

ggplot(data = preds, aes(x = method, y = estimate, color = isolation_source)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_line(aes(group = isolation_source), position = position_dodge(width = 0.2)) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.2)) +
    geom_point(size = 2, position = position_dodge(width = 0.2)) +
    theme_classic() +
    labs(x = "Method", y = "Residual")


## Pairwise comparisons
avg_comparisons(model = mod1, variables = list(method = "pairwise"), re.form = NA)


avg_comparisons(model = mod1, variables = list(method = "pairwise"), re.form = NA)
