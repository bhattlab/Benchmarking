library(ggplot2)
library(here)
library(RColorBrewer)

results <- read.table(here("test/maaslin_test.txt"), sep="\t", header=TRUE)
results <- mutate(results, signif=ifelse(qval < 0.05, "TRUE", "FALSE"))
limit <- max(abs(results$coef)) * c(-1, 1)


ggplot(results, aes(x=comparison, y=feature, fill=coef)) + 
  geom_tile() + 
  theme_bw() +
  scale_x_discrete(position = "top") +
  facet_grid(~metadata, scales="free_x") + 
  scale_fill_distiller(type="div", palette="RdBu", limit=limit) +
  geom_point(aes(alpha=signif), color = "white", size = 5, show.legend = F, shape=8) + 
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0)) + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) + 
  theme(axis.text.y = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.title.y = element_blank()) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(
    strip.text.x = element_text(
      size = 18, color = "black", face = "bold"
    ),
    strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid"
    ),
    strip.placement = 'outside'
  ) 
  

heatmap(results, Colv = comparison, Rowv = coef, scale="column")


results2 <- results %>% select(comparison, feature, coef)
as.matrix(results2)


##### TESTING AGGREGATION #####
donor <- c("D1", "D1", "D1", "D1", "D1", "D1", "D2", "D2", "D2", "D2", "D2", "D2", "D3", "D3","D3","D3","D3","D3")
replicate <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
condition <- c("NF", "NF", "NF", "OF", "OF", "OF","NF", "NF", "NF", "OF", "OF", "OF","NF", "NF", "NF", "OF", "OF", "OF")
diversity <- c(1,2,3,4,5,6,7,8,9,2,4,6,7,8,9,1,9,2)
testdf <- data.frame(donor, replicate, condition, diversity)

testdf %>%
  group_by(donor, condition) %>%
  summarise_at(vars(diversity), list(name = mean))

