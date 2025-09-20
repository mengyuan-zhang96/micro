##########
library(agricolae)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggpubr)
library(ranger)
library(ggplot2)
design <- read.delim("design.txt",row.names = 1)
DX <- read.delim("DX.txt")
a <- DX$sample
design <- design[rownames(design) %in% a, ]
colnames(design)
####The significance of the impact of environmental factors on output####
data <- design[,c(3,6,8,9,10,12)]
data$Slope <- as.factor(data$Slope)
#
library(car)
lm_model <- lm(Production ~ Slope + Soil + Adder + Temperature + Altitude, data = design)
# VIF
vif(lm_model)
# ）
rf_model <- ranger(Production ~ ., data = data, importance = 'impurity')
rf_model <- ranger(Production ~ ., data = data, importance = 'permutation', num.trees = 500, max.depth = 10)
rf_model
# Get and print the importance of the variable
var_importance <- importance(rf_model)
print(var_importance)
importance_df <- data.frame(Variable = names(var_importance),
                            Importance = var_importance)

# Sort the variables by their importance
importance_df <- importance_df[order(-importance_df$Importance),]

# Visualize the importance of variables
P1 <- ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_segment(aes(xend = Variable, yend = 0), color = "#4CAF50", size = 1) +
  geom_point(color = "#4CAF50", size = 4) +
  coord_flip() +
  labs(title = 'Variable Importance',
       x = 'Variables',
       y = 'Importance') +
  theme_minimal(base_size = 15) +
  theme(
    axis.text = element_text(size = 12, family = "Arial", color = "black"),  #
    axis.title = element_text(size = 14, family = "Arial", face = "bold", color = "black"),  #
    plot.title = element_text(size = 18, family = "Arial", face = "bold", hjust = 0.5, color = "black"),  #
    plot.margin = margin(10, 10, 10, 10),  #
    panel.grid.major = element_blank(),  #
    panel.grid.minor = element_blank(),  #
    axis.line = element_line(size = 0.8, color = "black")  #
  ) +
  scale_y_continuous(expand = c(0, 0))  #
P1
####Mixed-effect model####
library(lme4)
colnames(data)
model <- lmer(Production ~ Slope + Soil + Adder + Temperature + Altitude + (1 | Name), data = design)
model_interaction <- lmer(Production ~ Slope * Soil * Adder * Temperature * Altitude + (1 | Name), data = design)
summary(model_interaction)
library(emmeans)

#install.packages("emmeans")
library(emmeans)
emm_slope <- emmeans(model, ~ Slope)
emm_soil <- emmeans(model, ~ Soil)
emm_adder <- emmeans(model, ~ Adder)
emm_temp <- emmeans(model, ~ Temperature)
emm_altitude <- emmeans(model, ~ Altitude)


summary(emm_slope)
summary(emm_soil)
summary(emm_adder)
summary(emm_temp)
summary(emm_altitude)
emm_adder_slope <- emmeans(model, ~ Adder | Temperature)
summary(emm_adder_slope)


pairs(emm_adder_slope)



####Difference test####
data_diff_aov <- design[,c(3,7)]
colnames(data_diff_aov) <- c("Abundance","group")
data_diff_aov$group
data_diff_aov$group <- as.factor(data_diff_aov$group)
#data_diff_aov <-data_diff_aov[-c(7,17),]
shapiro.test(data_diff_aov$Abundance)
kruskal.test(Abundance~group,data = data_diff_aov)
library(dunn.test)
dunn.test(data_diff_aov$Abundance, data_diff_aov$group, method = "none")
dunn.test(data_diff_aov$Abundance, data_diff_aov$group,method = "holm")
pairwise.wilcox.test(data_diff_aov$Abundance, data_diff_aov$group, p.adjust.method = "none")
#
summary(aov(Abundance~group,data = data_diff_aov))
model= aov(Abundance~group,data = data_diff_aov)
Tukey_HSD <-TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
Tukey_HSD_tble <- as.data.frame(Tukey_HSD$group)
out <- LSD.test(model, "group",p.adj = "none")
stat <- out$groups
stat
#LSD.test(y, trt, DFerror, MSerror, alpha = 0.05, p.adj=c("none","holm","hommel",
#"hochberg", "bonferroni", "BH", "BY", "fdr"), group=TRUE, main = NULL,console=FALSE)
#","#DBA9A8", "#E43030", "#E99B78", "#FF8831"
labels <- c("a", "ab", "b")
####Take altitude as an example####
data_diff_aov <- design[,c(3,7)]
colnames(data_diff_aov) <- c("Abundance","group")
data_diff_aov$group
data_diff_aov$group <- as.factor(data_diff_aov$group)
labels <- c("a", "ab", "b")
medians <- aggregate(Abundance ~ group, data = data_diff_aov, median)

medians$label <- labels  #
P1 <- ggplot(data_diff_aov, aes(x = group, y = Abundance, fill = group)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = NA,  #
               color = c("#00aeba", "#e7b800", "#5d9a3b"),  #
               size = 0.7) +  # 保持中位数线
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(fill = group), color = "black", shape = 21, show.legend = FALSE) +  #
  theme_test() +
  labs(x = "Altitude", y = "Production") +  #
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank()) +  #
  scale_fill_manual(values = c("#00aeba", "#e7b800", "#5d9a3b")) +  #
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), #
    plot.caption = element_text(size = 12),  #
    axis.text = element_text(size = 12)  #
  ) +
  #
  geom_text(data = medians, aes(x = group, y = Abundance + 200, label = label), size = 5, vjust = 0)

#
print(P1)
####2. alpha####
####....2.1####
####........Difference test验####
data <- read.delim("alpha_diversity_indices-fun.txt")
design <- read.delim("design.txt",row.names = 1)
design$ID <- rownames(design)
data_diff = merge(data,design, by="ID")
data_diff1 <- data_diff[,c(3,15,16,17,18,20)]
data_diff_aov <- data_diff[,c(3,20)]
colnames(data_diff_aov) <- c("Abundance","group")
data_diff_aov$group
data_diff_aov$group <- as.factor(data_diff_aov$group)
library(dplyr)

data_cleaned <- data_diff_aov %>%
  group_by(group) %>%
  filter(Abundance >= quantile(Abundance, 0.25) - 1.5 * IQR(Abundance) &
           Abundance <= quantile(Abundance, 0.75) + 1.5 * IQR(Abundance)) %>%
  ungroup()


head(data_cleaned)
data_diff_aov <- data_cleaned
#data_diff_aov <-data_diff_aov[-c(7,17),]
shapiro.test(data_diff_aov$Abundance)
kruskal.test(Abundance~group,data = data_diff_aov)
library(dunn.test)

dunn.test(data_diff_aov$Abundance, data_diff_aov$group, method = "none")

dunn.test(data_diff_aov$Abundance, data_diff_aov$group,method = "holm")

pairwise.wilcox.test(data_diff_aov$Abundance, data_diff_aov$group, p.adjust.method = "none")

summary(aov(Abundance~group,data = data_diff_aov))
model= aov(Abundance~group,data = data_diff_aov)
Tukey_HSD <-TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
Tukey_HSD_tble <- as.data.frame(Tukey_HSD$group)
out <- LSD.test(model, "group",p.adj = "none")
stat <- out$groups
stat

labels <- c("a", "a", "a", "a")
data_diff1$Altitude1 <- as.factor(data_diff1$Altitude1)
data_diff1$Slope<- as.factor(data_diff1$Slope)
data_long <- data_diff1 %>%
  pivot_longer(cols = c(Altitude1, Temperature, Soil, Slope, Adder),
               names_to = "GroupType",
               values_to = "GroupValue")

fill_colors <- c(
  "1" = "#00aeba", "2" = "#e7b800", "3" = "#5d9a3b",  # Altitude1
  "4500-5000" = "#5d9a3b", "5500+" = "#DBA9A8", "5000-5500" = "#a8452f", "4000-4500" = "#e7b800", "-4000" = "#00aeba",  #
  "HR" = "#00aeba", "SH" = "#e7b800", "SD" = "#5d9a3b", "ZS" = "#a8452f",  # Soil
  "1" = "#00aeba", "2" = "#e7b800", "3" = "#5d9a3b", "4" = "#a8452f",  # Slope
  "B" = "#00aeba", "X" = "#e7b800", "DB" = "#5d9a3b", "N" = "#a8452f",
  "DN" = "#DBA9A8", "D" = "#E43030", "XB" = "#E99B78", "XN" = "#FF8831"  # Adder
)


group_counts <- data_long %>%
  group_by(GroupType) %>%
  summarise(count = n_distinct(GroupValue)) %>%
  ungroup()


data_long <- data_long %>%
  left_join(group_counts, by = "GroupType") %>%
  mutate(width_ratio = count / max(count))


P1 <- ggplot(data_long, aes(x = GroupValue, y = Chao1, color = GroupValue, fill = GroupValue)) +
  geom_boxplot(aes(width = 0.6), alpha = 0.5, color = "black", position = position_dodge(width = 0.8)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  facet_grid(~ GroupType, scales = "free_x", space = "free") +  #
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  #
  labs(title = "Boxplot of Chao1 by Group1 to Group5",
       x = "Group Type", y = "Chao1") +  #
  theme_bw() +  #
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # x
    axis.text.y = element_text(size = 12),  #
    strip.text = element_text(size = 14, face = "bold"),  #
    strip.background = element_rect(fill = "lightgray"),  #
    panel.grid.major = element_blank(),  #
    panel.grid.minor = element_blank(),  #
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  #
    plot.margin = margin(10, 10, 10, 10),  #
    panel.spacing.x = unit(1, "lines")  #
  ) +
  scale_fill_manual(values = fill_colors) +  #
  scale_color_manual(values = fill_colors)  #
ggsave("../01alpha/chao1.pdf", plot =P1, width = 12, height = 4)
####Take altitude as an example####
data_diff_aov <- data_diff[,c(5,15)]
colnames(data_diff_aov) <- c("Abundance","group")
data_diff_aov$group
data_diff_aov$group <- as.factor(data_diff_aov$group)
medians <- aggregate(Abundance ~ group, data = data_diff_aov, median)
P1 <- ggplot(data_diff_aov, aes(x = group, y = Abundance, fill = group)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = NA,
               color = c("#00aeba", "#e7b800","#5d9a3b"),  ),
               size=0.7) +  # 保持中位数线
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(fill = group), color = "black", shape = 21, show.legend = FALSE) +
  theme_test() +
  labs(x = "soil", y = "shannon") +
  #scale_y_continuous(expand = c(0, 0), limits = c(1000, 5000)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank()) +  # 移除方格线
  scale_fill_manual(values =  c("#00aeba", "#e7b800","#5d9a3b")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.caption = element_text(size = 12),
    axis.text = element_text(size = 12)
  )

#geom_text(data = medians, aes(x = group, y = Abundance + 200, label = labels), size = 5, vjust = 0)
P1

##混合效应模型##
install.packages("lmerTest")
library(lme4)  # 用于混合效应模型
library(lmerTest)  # 用于显著性检验
library(ggplot2)  # 用于可视化
# 构建混合效应模型
model <- lmer(Production ~ Shannon + Altitude + Temperature + Soil +
                Slope + Adder + (1 | Name), data = data_diff)
colnames(data_diff)
model <- lmer(Shannon ~ Altitude + Temperature + Soil +
                Slope + Adder + (1 | Name), data = data_diff)
summary(model)
anova(model)
##多样性变化##
data_diff_aov1 <- data_diff[,c(5,13)]
colnames(data_diff_aov1) <- c("Abundance","group")
data_diff_aov1$group <- as.character.factor(data_diff_aov1$group)
model= aov(Abundance~group,data = data_diff_aov1)
summary(model)
#dunnTest(data_diff_aov$Abundance ~ data_diff_aov$group, method = "bh")
Tukey_HSD <-TukeyHSD(model)
Tukey_HSD_tble <- as.data.frame(Tukey_HSD$group)
Tukey_HSD_tble
library(doBy)
out <- LSD.test(model, "group",p.adj = "none")
stat <- out$groups
stat
P3 <- ggplot(data_diff_aov, aes(x = group, y = Abundance, fill = group)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = NA,  #
               color = c("#00aeba", "#e7b800","#5d9a3b","#a8452f"), ),
               size=0.7) +  # 保持中位数线
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(fill = group), color = "black", shape = 21, show.legend = FALSE) +  # 添加数据点
  theme_test() +
  labs(x = "data_R", y = "evenness") + #
  #scale_y_continuous(expand = c(0, 0), limits = c(0.725, 0.745)) +  #
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid = element_blank()) +  #
  scale_fill_manual(values =  c("#00aeba", "#e7b800","#5d9a3b","#a8452f")) +  #
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), #
    plot.caption = element_text(size = 12), #
    axis.text = element_text(size = 12) #
  )
P3
####....2.3beta####
library(vegan)
library(ggplot2)

otu_table <- read.delim("otu_table_norm.txt",row.names = 1)
design<- read.delim("design.txt",row.names = 1)
colSums(otu_table)
otu_relative_abundance <- sweep(otu_table, 2, colSums(otu_table), FUN = "/")
colSums(otu_relative_abundance)
otu_table <-as.data.frame(t(otu_table))
####Take altitude as an example####
metadata = data.frame(Group = design$Altitude1,
                      SampleID =row.names(design),
                      row.names = row.names(design))

otu_table <- otu_table[rownames(otu_table) %in% metadata$SampleID, ]
metadata <- metadata[rownames(otu_table), ]  #
metadata$Group <- as.factor(metadata$Group)

bray_curtis_dm <- vegdist(otu_table, method = "bray")

library(vegan)

permanova_result <- adonis2(bray_curtis_dm ~ Group, data = metadata, permutations = 999)
print(permanova_result)


print(pairwise_Canberra)

anosim_result <- anosim(bray_curtis_dm, metadata$Group, permutations = 999)
print(anosim_result)


pcoa_results <- cmdscale(bray_curtis_dm, eig = TRUE, k = 2)


pcoa_df <- as.data.frame(pcoa_results$points)
colnames(pcoa_df) <- c("PC1", "PC2")
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df$Group <- metadata$Group  #

P1 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  stat_ellipse(aes(fill = Group), level = 0.4,type = "norm", alpha = 0.2, geom = "polygon", show.legend = FALSE) +  #
  #geom_text(aes(label = SampleID), vjust = -0.5, size = 3) +
  xlab(paste("PCoA1 (", round(pcoa_results$eig[1] / sum(pcoa_results$eig) * 100, 2), "%)", sep = "")) +
  ylab(paste("PCoA2 (", round(pcoa_results$eig[2] / sum(pcoa_results$eig) * 100, 2), "%)", sep = "")) +
  ggtitle("PCoA Plot (Bray-Curtis Distance)") +
  theme_test() +
  #xlim(c(-0.1, 0.1)) +  #
  #ylim(c(-0.1, 0.1)) +  #
  scale_color_manual(values = c("#00aeba", "#e7b800", "#5d9a3b")) +  #
  scale_fill_manual(values = c("#00aeba", "#e7b800", "#5d9a3b"))  #
P1

####3.Relative abundance bar chart####
package_list = c("reshape2","ggplot2","vegan")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

package_list = c("digest","ggrepel")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
package_list = c("kassambara/ggpubr")
for(p in package_list){
  q=unlist(strsplit(p,split = "/"))[2]
  if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install_github(p)
    suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
####....3.1 Phylum####

otu_table <- read.delim("P.txt",row.names = 1)
otu_table <- otu_table[, colnames(otu_table) %in% rownames(design)]
colSums(otu_table)
tax_sample = otu_table
otu_table <-as.data.frame(t(otu_table))


colnames(design)

sampFile = data.frame(group = design$Soil,
                      sample =row.names(otu_table),
                      row.names = row.names(otu_table))

idx = rownames(sampFile) %in% colnames(tax_sample) # match design with alpha
sampFile = sampFile[idx,]
tax_sample = tax_sample[,rownames(sampFile)]

colSums(tax_sample)

mean_sort = tax_sample[(order(-rowSums(tax_sample))), ]
mean_sort = as.data.frame(mean_sort)
colSums(mean_sort)

other = colSums(mean_sort[11:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:10, ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[11] = c("Low abundance")

colSums(mean_sort)

merge_tax=mean_sort

mat_t = t(merge_tax)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,c(-1,-3)]


mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno


mean_sort=as.data.frame(mat_mean_final)

mean_sort$tax = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("tax")))

data_all$value <- as.numeric(data_all$value)

custom_colors <- c("#A2D2E7", "#67A8CD","#FFC17F","#CF9F88", "#6FB3A8", "#B3E19B", "#50AA4B","#FF9D9F","#F36569", "#3581B7", "#CDB6DA", "#704BA3","#9A7FBD","#DBA9A8", "#E43030", "#E99B78", "#FF8831")

p = ggplot(data_all, aes(x=variable, y=value, fill=tax)) +
  geom_bar(stat="identity", position="stack", width=0.7) +
  #scale_y_continuous(limits=c(0, 800)) +
  scale_fill_manual(values=custom_colors) +
  xlab("Groups") +
  ylab("Abundance") +
  theme_classic()

p

