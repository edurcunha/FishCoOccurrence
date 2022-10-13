################################################################################
#  Rodrigues et al. 2022
#  Functional dissimilarity determines the co-occurrence of native and non-native species
################################################################################

# co-occurrence analysis 
# hypothesis 1 = functionally similar native and non-native species are less likely to co-occur than expected by chance
# hypothesis 2 = native and non-native species co-occur less when compared to pairs of native species

# Packages
library(cooccur) # estimate the co-occurrence
library(ggplot2) # perform graphs
library(lme4) # perform models
library(lmerTest) #obtain p values
library(MuMIn) # obtain r²
library(car) # Anova

# Open the presence=absence matrix of species by samples
spt.ab <- read.table("spatial_matrix.txt", h=T) #abundance data
spt.in <- read.table("spatial_matrix.txt", h=T) #incidence data
spt.in[,8:85][spt.in[,8:85]>0] <- 1  #transform in Presence/Absence
colnames(spt.ab)
colnames(spt.in)
rownames(spt.in) <- spt.in$sites

# CO-OCCUR PACKAGE 
# Veech, J. A. (2013) A probabilist model for analysing species co-occurrence.
# Griffith et al. (2016) cooccur: Probabilistic species co-occurrence analysis in R
library(cooccur)

# Create a matrix sites x species with 0 and 1 only
m <- spt.in[,8:85]

# Transpose the matrix to species x sites
m.t <- t(m)

# Calculate the co-occurrence
cooccur.spp <- cooccur(mat=m.t, type = "spp_site", thresh = FALSE, 
                       spp_names = TRUE)
summary(cooccur.spp)
plot(cooccur.spp)
hist(cooccur.spp$results$prob_cooccur)
pair.profile(cooccur.spp)
obs.v.exp(cooccur.spp)
write.table(cooccur.spp$results, 'cooccur_package.txt')
#cooccur.res <- cooccur.spp$results 

# Add the status (native and non-native) of species, 
# the group (native x native; native x non-native; non-native x non-native) 
# and the scaled_cooccur (obs_cooccur/exp_cooccur)
cooccur.res <- read.table("cooccur_package.txt", h=T)

hist(cooccur.res$scaled_cooccur)

cooccur.res <- cooccur.res[1:2703,] # select only native x native and native x non-native groups
summary(cooccur.res)
cooccur.res <- cooccur.res[cooccur.res$obs_cooccur != 0,] #remove obs_cooccur = 0
plot(x = cooccur.res$group, y = cooccur.res$scaled_cooccur)

# hypothesis 1 = co-occurrence ~ functional distance
# open the script that estimates the functional distance (Shared_script2)
# generalized linear mixed models (lme4)

# random factor: sp1_name (only native sp.)
# y: log(scaled_cooccur)
# x: dist_mt_gow

# Open the matrix with distance values
dataset <- read.table('cooccur_package.txt', h=T)

# Select native x nonnative
sub <- dataset[dataset$group == "nativenonnative" & dataset$obs_cooccur != 0,]

sub$sp1_name <- as.factor(sub$sp1_name)
sub$sp2_name <- as.factor(sub$sp2_name)
sub$sp1 <- as.factor(sub$sp1)

sub <- na.omit(sub)

mod2.lmer.a <- lmer(log(scaled_cooccur) ~ dist_mt_gow + (1|sp1_name), 
                    data = sub)
summary(mod2.lmer.a)
plot(mod2.lmer.a)
AIC(mod2.lmer.a)
r.squaredGLMM(mod2.lmer.a)

# hypothesis 2 = co-occurrence ~ groups
# generalized linear mixed models (lme4)

# random factor: sp1_name (only native sp.)
# y: log(scaled_cooccur)
# x: group

with(cooccur.res, hist(scaled_cooccur))
with(cooccur.res, hist(log(scaled_cooccur)))
colnames(cooccur.res)

cooccur.res <- na.omit(cooccur.res)
mod1.lmer.a <- lmer(log(scaled_cooccur) ~ group + (1|sp1_name), 
                    data = cooccur.res)
summary(mod1.lmer.a)
plot(mod1.lmer.a)
library(car)
Anova(mod1.lmer.a, type = 'III')
AIC(mod1.lmer.a)
r.squaredGLMM(mod1.lmer.a)

##### Figures #####

# Figure 1 is at Shared_script2

#### Figure 2 ####
#install BiocManager
#install ComplexHeatmap
biocLite("ComplexHeatmap")
library(BiocManager)
library(ComplexHeatmap) #Heatmap function not heatmap function
library(RColorBrewer)
library(circlize)

cooccur_heatmap <- pairwise_cooccur[1:78, 1:53] 
#pairwise_cooccur = species x species matrix with scaled_cooccur values

column_ha_bar <- rowAnnotation("Number \n of sampling \n events" = anno_barplot(c(4,8,107,95,15,40,4,8,19,9,19,6,4,7,7,3,18,1,2,3,52,53,7,61,55,38,98,1,43,16,3,68,10,77,5,130,83,21,16,1,2,1,13,4,40,14,32,9,86,6,69,8,1,52,16,61,1,1,1,5,72,1,14,41,138,9,3,6,1,120,19,11,94,19,40,40,126,41),  
                                                                                height = unit(3, "cm"))) #create a vector with each species count in order based on column
split <- rep(c("Native", "Nonnative"), c(53, 25))
split

col_fun <- colorRamp2(c(0,1,10,20), c("khaki1","salmon1","orangered3", "orangered4"))

col_column <- c("seagreen4","seagreen4","seagreen4","seagreen4","seagreen4",
                "deepskyblue4","deepskyblue4","seagreen4","seagreen4",
                "deepskyblue4","darkorange","darkorange","seagreen4",
                "darkorange", "darkorange",
                "indianred1","indianred1","indianred1","indianred1","indianred1","indianred1","indianred1","indianred1",
                "seagreen4","seagreen4","seagreen4","seagreen4","seagreen4","seagreen4","seagreen4",
                "indianred1","seagreen4","indianred1","indianred1","indianred1","seagreen4",
                "indianred1","seagreen4","indianred1","seagreen4","indianred1","seagreen4",
                "indianred1","indianred1","seagreen4","seagreen4","seagreen4","seagreen4","seagreen4",
                "darkorange","seagreen4","seagreen4","indianred1")

col_row <- c("seagreen4","seagreen4","seagreen4","seagreen4","seagreen4",
             "deepskyblue4","deepskyblue4","seagreen4","seagreen4",
             "deepskyblue4","darkorange","darkorange","seagreen4",
             "darkorange", "darkorange",
             "indianred1","indianred1","indianred1","indianred1","indianred1","indianred1","indianred1","indianred1",
             "seagreen4","seagreen4","seagreen4","seagreen4","seagreen4","seagreen4","seagreen4",
             "indianred1","seagreen4","indianred1","indianred1","indianred1","seagreen4",
             "indianred1","seagreen4","indianred1","seagreen4","indianred1","seagreen4",
             "indianred1","indianred1","seagreen4","seagreen4","seagreen4","seagreen4","seagreen4",
             "darkorange","seagreen4","seagreen4","indianred1",
             "deepskyblue4","indianred1","indianred1","indianred1","cyan3",
             "seagreen4","indianred1","seagreen4","indianred1","indianred1","seagreen4",
             "indianred1","indianred1","indianred1","red",
             "indianred1","indianred1","indianred1","indianred1","seagreen4","darkorange",
             "seagreen4","indianred1","seagreen4","indianred1")

order.row <- c("characiformes","characiformes","characiformes","characiformes","characiformes",
               "perciformes","perciformes","characiformes","characiformes",
               "perciformes","gymnotiformes","gymnotiformes","characiformes",
               "gymnotiformes", "gymnotiformes",
               "siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes",
               "characiformes","characiformes","characiformes","characiformes","characiformes","characiformes","characiformes",
               "siluriformes","characiformes","siluriformes","siluriformes","siluriformes","characiformes",
               "siluriformes","characiformes","siluriformes","characiformes","siluriformes","characiformes",
               "siluriformes","siluriformes","characiformes","characiformes","characiformes","characiformes","characiformes",
               "gymnotiformes","characiformes","characiformes","siluriformes",
               "perciformes","siluriformes","siluriformes","siluriformes","pleuronectiformes",
               "characiformes","siluriformes","characiformes","siluriformes","siluriformes","characiformes",
               "siluriformes","siluriformes","siluriformes","myliobatiformes",
               "siluriformes","siluriformes","siluriformes","siluriformes","characiformes","gymnotiformes",
               "characiformes","siluriformes","characiformes","siluriformes")
order.row <- as.factor(order.row)
order.row <- data.frame(order = order.row, species = rownames(cooccur_heatmap))
order.row <- order.row[order(order.row$order), ]

order.col <- c("characiformes","characiformes","characiformes","characiformes","characiformes",
               "perciformes","perciformes","characiformes","characiformes",
               "perciformes","gymnotiformes","gymnotiformes","characiformes",
               "gymnotiformes", "gymnotiformes",
               "siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes",
               "characiformes","characiformes","characiformes","characiformes","characiformes","characiformes","characiformes",
               "siluriformes","characiformes","siluriformes","siluriformes","siluriformes","characiformes",
               "siluriformes","characiformes","siluriformes","characiformes","siluriformes","characiformes",
               "siluriformes","siluriformes","characiformes","characiformes","characiformes","characiformes","characiformes",
               "gymnotiformes","characiformes","characiformes","siluriformes")
order.col <- as.factor(order.col)
order.col <- data.frame(order = order.col, species = colnames(cooccur_heatmap))
order.col <- order.col[order(order.col$order), ]

ht_cooccur <- Heatmap(cooccur_heatmap, na_col = "grey70",
                      col= col_fun,
                      column_order = order.col$species,
                      row_names_side = "left",
                      row_names_centered = TRUE,
                      row_order = order.row$species,
                      row_split = split,
                      left_annotation = column_ha_bar,
                      cluster_columns = FALSE, cluster_rows = FALSE, 
                      column_names_gp = gpar(fontsize = 10, col = col_column),
                      row_names_gp = gpar(fontsize = 10, col = col_row),
                      width = unit(15, "cm"), height = unit(25, "cm"),
                      heatmap_legend_param = list(title = "Scaled co-occurrence", 
                                                  at = c(0, 1, 10, 20), legend_width = unit(4, "cm"),
                                                  direction = "horizontal"),
                      row_gap = unit(3, "mm"),
                      border_gp = gpar(col = "black", lwd = 1))
ht_cooccur
draw(ht_cooccur, heatmap_legend_side = "bottom")
#save png with width 1000 and aspect ratio on

#seagreen4 = characiformes
#deepskyblue4 = cichliformes
#darkorange = gymnotiformes
#indianred1 = siluriformes
#cyan3 = pleuronectiformes
#red = myliobatiformes  

#### Figure 3 ####
#mod2.lmer.a <- lmer(log(scaled_cooccur) ~ dist_mt_gow2 + (1|sp1_name), 
#                    data = sub)
#OBTAIN FITTED VALUES AND ADD A NEW COLUM IN DATASET
sub$fitted_values <- predict(mod2.lmer.a, type = "response")

#PLOT THE PREDICTED DATA
gg <- ggplot(sub, aes(x=dist_mt_gow2, y=fitted_values)) + 
  geom_point(size = 4, col = rgb(0,0,0,.2)) +
  geom_smooth(color="black", method="glm", se=T) + 
  labs(x = 'Functional distance', y = 'Scaled co-occurrence \n (predicted values)') +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))
gg

ggsave(plot = gg, width = 4, height = 3, dpi = 600, 
       filename = "figure3.png")

#### Figure 4 ####
#OBTAIN FITTED VALUES AND ADD A NEW COLUM IN DATASET
#mod1.lmer.a <- lmer(log(scaled_cooccur) ~ group + (1|sp1_name), 
#                    data = cooccur.res)
cooccur.res$fitted_values <- predict(mod1.lmer.a, type = "response")

library(PupillometryR)
library(ggplot2)
pointcolors <- setNames( c("black","burlywood3","chartreuse4","darkorange2"),levels(cooccur.res$group))
p_fitted <- ggplot(cooccur.res, aes(x = group, y = fitted_values)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2,
                   fill='gray85',color='gray85', width = 0.9) + 
  scale_x_discrete(labels = c("nativenative" = "Native \n x \n Native", 
                              "nativenonnative" = "Native \n x \n Nonnative")) + 
  scale_y_continuous(name="Scaled Co-occurrence \n (predicted values)",
                     limits=c(-0.2,1.3)) + 
  geom_jitter(width=0.1, aes(colour=factor(group)), size=2, 
              alpha=0.5) + 
  scale_color_manual(values=pointcolors) + 
  theme_bw() + theme(axis.text.x = element_text(size = 12), 
                     axis.title.x = element_blank(), #remove the name of the factor
                     axis.text.y = element_text(size = 12),
                     panel.border = element_blank(), #remove the borders
                     axis.line = element_line(colour = "black", size=0.7), #add border just on axis
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(size=12), 
                     legend.position = 'none')

p_fitted
ggsave(plot = p_fitted, width = 4, height = 3, dpi = 600, 
       filename = "figure4.png")

#### Figure S2 ####
### graph number of sites
spt.sites <- read.table("species_sites.txt", h=T) #matrix with the number of sites for each sp.
library(dplyr)
library(ggplot2)
library(stringr)
remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

plt <- ggplot(spt.sites[1:10,]) +
  geom_col(aes(x = reorder(species, site), y = site, fill = site), position = "dodge2", show.legend = TRUE,
           alpha = .9) +
  coord_polar() + 
  # New fill and legend title for number of tracks per region
  scale_fill_gradientn("Number of samples",
                       colours = c( "#6C5B7B","#C06C84","#F67280","#F8B195")) +
  theme_test() +
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 12)) 
plt
ggsave(plot = plt, width = 6, height = 4, dpi = 600, 
       filename = "Fig. S2.png")

#### Figure S3 ####
traits_heatmap <- pairwise_traits[1:53, 54:78]

col_column <- c("deepskyblue4","indianred1","indianred1","indianred1",
                "cyan3","seagreen4","indianred1","seagreen4","indianred1",
                "indianred1","seagreen4","indianred1","indianred1","indianred1",
                "red","indianred1","indianred1","indianred1","indianred1",
                "seagreen4","darkorange","seagreen4","indianred1",
                "seagreen4","indianred1")

col_row <- c("seagreen4","seagreen4","seagreen4","seagreen4","seagreen4",
             "deepskyblue4","deepskyblue4","seagreen4","seagreen4",
             "deepskyblue4","darkorange","darkorange","seagreen4",
             "darkorange","darkorange","indianred1","indianred1","indianred1",
             "indianred1","indianred1","indianred1","indianred1","indianred1",
             "seagreen4","seagreen4","seagreen4","seagreen4","seagreen4",
             "seagreen4","seagreen4","indianred1","seagreen4","indianred1",
             "indianred1","indianred1","seagreen4","indianred1","seagreen4",
             "indianred1","seagreen4","indianred1","seagreen4","indianred1",
             "indianred1","seagreen4","seagreen4","seagreen4","seagreen4",
             "seagreen4","darkorange","seagreen4","seagreen4","indianred1")

order.row.traits <- c("characiformes","characiformes","characiformes","characiformes","characiformes",
                      "perciformes","perciformes","characiformes","characiformes",
                      "perciformes","gymnotiformes","gymnotiformes","characiformes",
                      "gymnotiformes", "gymnotiformes",
                      "siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes","siluriformes",
                      "characiformes","characiformes","characiformes","characiformes","characiformes","characiformes","characiformes",
                      "siluriformes","characiformes","siluriformes","siluriformes","siluriformes","characiformes",
                      "siluriformes","characiformes","siluriformes","characiformes","siluriformes","characiformes",
                      "siluriformes","siluriformes","characiformes","characiformes","characiformes","characiformes","characiformes",
                      "gymnotiformes","characiformes","characiformes","siluriformes")
order.row.traits <- as.factor(order.row.traits)
order.row.traits <- data.frame(order = order.row.traits, species = rownames(traits_heatmap))
order.row.traits <- order.row.traits[order(order.row.traits$order), ]

order.col.traits <- c("perciformes","siluriformes","siluriformes","siluriformes","pleuronectiformes",
                      "characiformes","siluriformes","characiformes","siluriformes","siluriformes","characiformes",
                      "siluriformes","siluriformes","siluriformes","myliobatiformes",
                      "siluriformes","siluriformes","siluriformes","siluriformes","characiformes","gymnotiformes",
                      "characiformes","siluriformes","characiformes","siluriformes")
order.col.traits <- as.factor(order.col.traits)
order.col.traits <- data.frame(order = order.col.traits, species = colnames(traits_heatmap))
order.col.traits <- order.col.traits[order(order.col.traits$order), ]

#sort(colnames(traits_heatmap))
#sort(rownames(traits_heatmap))
#col = col_column
#col = col_row

ht_traits <- Heatmap(traits_heatmap, column_order = order.col.traits$species,
                     row_order = order.row.traits$species,
                     cluster_columns = FALSE, cluster_rows = FALSE, 
                     column_names_gp = gpar(fontsize = 10, col = col_column),
                     row_names_gp = gpar(fontsize = 10, col = col_row),
                     width = unit(15, "cm"), height = unit(20, "cm"),
                     heatmap_legend_param = list(title = "Functional \n distance", legend_height = unit(4, "cm"),
                                                 direction = "horizontal"),
                     border = TRUE,
                     color = colorRampPalette(c("navy","white","firebrick3"))(50),
                     column_title = "Non-native species", 
                     column_title_side = "bottom",
                     row_title = "Native species",
                     row_title_side = "right")
ht_traits
draw(ht_traits, heatmap_legend_side = "bottom")