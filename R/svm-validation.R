# https://blog.devgenius.io/ensemble-models-using-caret-in-r-d54e4e646968
# https://rviews.rstudio.com/2019/06/19/a-gentle-intro-to-tidymodels/

library(tidyverse)
library(caret)
library(patchwork)
# library(caretEnsemble)
library(showtext)
library(ggokabeito)
library(viridis)
# library(xgboost)
library(ggforce)

set.seed(42)

font_add("Fuzzy Bubbles", regular = "fonts/ABeeZee-Regular.ttf")
showtext_auto()
theme_clean <- function() {
  theme_minimal(base_family = "Fuzzy Bubbles") +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 16, family = "Fuzzy Bubbles"),
          plot.background = element_rect(fill = "white", color = NA),
          axis.text = element_text(size = 16),
          axis.title = element_text(face = "bold", size = 20),
          strip.text = element_text(face = "bold", size = rel(0.8), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.text = element_text(size = 16))
}

z1 <- readRDS("data/z_4500_clean_extra") %>% 
  drop_na()
z1_extra <- readRDS("data/z_4500_clean_extra_extra_unnormalised") %>% 
  mutate(cor_logit = gtools::logit(cor, min = -1, max = 1) %>%
           scale())
names(z1_extra) <- names(z1)
z1_extra <- z1_extra |> 
  mutate(objid = as.character(objid),
         objid1 = as.character(objid1),
         specobjid = as.character(specobjid),
         specobjid1 = as.character(specobjid1))
z1 <- z1 |> 
  mutate(objid = as.character(objid),
         objid1 = as.character(objid1),
         specobjid = as.character(specobjid),
         specobjid1 = as.character(specobjid1))

validation <- z1_extra #|> 
 # filter(subclass != "CV", subclass1 != "CV")

svm1 <- readRDS("Locus-elsevier1/models/svm-sigma008-010C18-20")
svm1
plot(svm1)
svm1pred<-predict(svm1, validation)
svm1values<-data.frame(obs=validation$cor_logit, 
                       pred=svm1pred, 
                       res=svm1pred-validation$cor_logit)
yardstick::metrics(data = svm1values, truth = obs, estimate = pred)
defaultSummary(svm1values)

plot1 <- ggplot(svm1values, aes(y=obs,
                       x=pred, 
                       colour=abs(res))) +
  geom_point(alpha=0.9, show.legend = F, size = 0.3) + 
#  geom_smooth(se=FALSE,colour="red", linetype="dashed", size=0.5)+ 
  geom_abline(slope=1, linetype="dashed") +
  scale_color_gradient(low = "#7BA0B4", high = "#0A2D46") + 
  labs(y = "Observered Logit Correlation Value",
       x = "Predicted Logit Correlation Value") +
  theme_clean() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

plot2 <- ggplot(svm1values, aes(y = res,
                       x = pred)) +
  geom_point(alpha=0.9, show.legend = F, col = "#44728C", size = 0.3) + 
#  geom_smooth(se=FALSE,colour="red", linetype="dashed", size=0.5)+ 
#  geom_abline(slope=1, linetype="dashed") +
  labs(y = "Residuals",
       x = "Predicted Logit Correlation Value") +
  theme_clean() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))


svm1values |> 
  mutate(extreme = res>2) |> 
  ggplot(aes(res)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white", bins = 50) +
  geom_circle(aes(x0 = 3.42, y0 = 0.02, r = 0.08), 
              inherit.aes = FALSE, colour = "red") +
  geom_density(lwd = 0.5, colour = 4,
               fill = 4, alpha = 0.25) +
  labs(x = "Residuals from SVR Model Fit") +
  theme_clean() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())
# ggsave(filename = "Locus-elsevier1/figures/histogram.png", width = 515, height = 355, units = "px")

qqplot(svm1values$pred, 
       svm1values$obs,
       main = "",
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles")
abline(a = 0, b = 1, col = "darkblue", lwd = 2)

plot3 <- svm1values %>% 
  ggplot(aes(sample = res)) + 
  geom_qq(size = 0.3) +
#  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis_d() +
  theme_clean()

M <- 0.3
rating <- function(gr, rr, ir, gt, rt, it) {
  delta.CS <- (gt - rt) - (gr - rr)
  delta.CL <- (rt - it) - (rr - ir)
  RS <- 1 - abs(delta.CS / M)
  RL <- 1 - abs(delta.CL / M)
  RS * RL
}
z1 <- z1 %>% 
  mutate(oisin = rating(g, r, i, g1, r1, i1))
cor(z1$oisin, z1$cor_logit)

layout <- "
AAAA#CCCC
BBBB#CCCC
"
plot1 + plot2 + plot3 + 
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ") ")
