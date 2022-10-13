################################################################################
#  Rodrigues et al. 2022
#  Functional dissimilarity determines the co-occurrence of native and non-native species
################################################################################

# Functional analysis 
# hypothesis 1 = functionally similar native and non-native species are less likely to co-occur than expected by chance
# hypothesis 2 = native and non-native species co-occur less when compared to pairs of native species

#  FUNCTIONAL DATA: Brosse et al. 2021

# Open the species by traits matrix
spt.traits <- read.table("spatial_traits.txt", h=T) #all spatial traits 
rownames(spt.traits) <- spt.traits$species1
spt.m.traits <- spt.traits[,c(8:17)] 

# Open the sites by species matrix (abundance)
spt.ab <- read.table("spatial_matrix.txt", h=T) #abundance data
rownames(spt.ab) <- spt.ab$sites

#### Estimate the quality of functional space ####

# Load packages
library(cluster)
source('quality_funct_space_fromdist.R')

# Standardizes the variables
scale.mt <- scale(spt.m.traits, center = T, scale = T)
# it is not necessary to standardize if it's going to use Gower's distance

# Creates a dissimilarity matrix with Gower's distance
gow.mt <- daisy(spt.m.traits, metric = "gower")

# Estimates the quality of functional spaces
quality.mt <- quality_funct_space_fromdist(gow.mt, nbdim = 10,
                                           plot = "quality_morpho")
quality.mt$meanSD #provides the mSD of each functional space

# RESULT: the best functional space is the one with 8 dimensions (mSD = 0.000475)

# Provides the functional distance in the best FS (similarity)
# standardized distances between species
dist.sp.mt <- quality.mt$details_funct_space$dist_st$m_8D
dist.sp.mt <- data.frame(as.matrix(dist.sp.mt))
write.table(dist.sp.mt, 'distance_mt_gow.txt')

# coordinates of species in the nbdim multidimensional functional space (PCoA axis)
coord.sp.mt <- quality.mt$details_funct_space$mat_coord[,1:8] #to build de FS

# The functional similarity between species was measured as the 
# standardized distance between each pair of species in the 
# functional space. The best functional space was estimated
# according to Maire et al. 2015. For our data, the best functional
# space presented eight dimensions (mSD = 0.000475).

#### Build functional space ####

# create a vector with the status of species
status <- rep( c("native", "non-native"), c(53,25))
coord.sp.mt <- data.frame(status, coord.sp.mt)

### Figure 1 ###
library(ggplot2)
library(tidyverse)
theme_set(theme_bw(16))

p1 <- coord.sp.mt %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()
p1

hull_data <- coord.sp.mt %>%
  drop_na() %>%
  group_by(status) %>%
  slice(chull(PC1, PC2))

p1 +
  geom_polygon(data = hull_data, aes(fill = status,
                                     colour = status),
               alpha = 0.3, show.legend = FALSE)

axis12 <- ggplot(coord.sp.mt, aes(x = PC1, y = PC2, color = status)) +
  geom_point(size = 3) + theme_test() +
  scale_color_manual(values=c("#046c9aff", "#ffcbaeff")) +
  labs(x = "PCoA 1", y = "PCoA 2") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size=12),
        legend.position = c(0.85,0.85)) +
  scale_y_continuous(limits = c(-0.3, 0.4)) + scale_x_continuous(limits = c(-0.4, 0.6)) +
  geom_polygon(data = hull_data, aes(fill = status, colour = status),
               alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values=c("#046c9aff", "#ffcbaeff"))
axis12

p1 <- coord.sp.mt %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_point()
p1

hull_data <- coord.sp.mt %>%
  drop_na() %>%
  group_by(status) %>%
  slice(chull(PC3, PC4))

p1 +
  geom_polygon(data = hull_data, aes(fill = status,
                                     colour = status),
               alpha = 0.3, show.legend = FALSE)

axis34 <- ggplot(coord.sp.mt, aes(x = PC3, y = PC4, color = status)) +
  geom_point(size = 3) + theme_test() +
  scale_color_manual(values=c("#046c9aff", "#ffcbaeff")) +
  labs(x = "PCoA 3", y = "PCoA 4") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size=12),
        legend.position = c(0.85,0.85)) +
  scale_y_continuous(limits = c(-0.3, 0.4)) + scale_x_continuous(limits = c(-0.4, 0.6)) +
  geom_polygon(data = hull_data, aes(fill = status, colour = status),
               alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values=c("#046c9aff", "#ffcbaeff")) 
axis34

p1 <- coord.sp.mt %>%
  ggplot(aes(x = PC5, y = PC6)) +
  geom_point()
p1

hull_data <- coord.sp.mt %>%
  drop_na() %>%
  group_by(status) %>%
  slice(chull(PC5, PC6))

p1 +
  geom_polygon(data = hull_data, aes(fill = status,
                                     colour = status),
               alpha = 0.3, show.legend = FALSE)

axis56 <- ggplot(coord.sp.mt, aes(x = PC5, y = PC6, color = status)) +
  geom_point(size = 3) + theme_test() +
  scale_color_manual(values=c("#046c9aff", "#ffcbaeff")) +
  labs(x = "PCoA 5", y = "PCoA 6") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size=12),
        legend.position = c(0.85,0.85)) +
  geom_polygon(data = hull_data, aes(fill = status, colour = status),
               alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values=c("#046c9aff", "#ffcbaeff")) +
  scale_y_continuous(limits = c(-0.3, 0.4)) + scale_x_continuous(limits = c(-0.4, 0.6))
axis56

p1 <- coord.sp.mt %>%
  ggplot(aes(x = PC7, y = PC8)) +
  geom_point()
p1

hull_data <- coord.sp.mt %>%
  drop_na() %>%
  group_by(status) %>%
  slice(chull(PC7, PC8))

p1 +
  geom_polygon(data = hull_data, aes(fill = status,
                                     colour = status),
               alpha = 0.3, show.legend = FALSE)

axis78 <- ggplot(coord.sp.mt, aes(x = PC7, y = PC8, color = status)) +
  geom_point(size = 3) + theme_test() +
  scale_color_manual(values=c("#046c9aff", "#ffcbaeff")) +
  labs(x = "PCoA 7", y = "PCoA 8") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size=12),
        legend.position = c(0.85,0.85)) +
  geom_polygon(data = hull_data, aes(fill = status, colour = status),
               alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values=c("#046c9aff", "#ffcbaeff")) +
  scale_y_continuous(limits = c(-0.3, 0.4)) + scale_x_continuous(limits = c(-0.4, 0.6))
axis78

library(ggpubr) #combine multiple plots
funct_space <- ggarrange(axis12, axis34, axis56, axis78,
                         ncol = 2, nrow = 2,
                         common.legend = TRUE, legend = "bottom")
funct_space

ggsave(plot = funct_space, width = 6, height = 5, dpi = 600, 
       filename = "figure1b.pdf")

#PERMANOVA
library(vegan)
#Step 1: create the distance matrix
gow.mt

#Step 2: visualize the groups in a PCoA
pcoa <- cmdscale(gow.mt, eig = TRUE)
pcoa
ordiplot(scores(pcoa)[ ,c(1,2)], type = "t")

#Step 3: PERMANOVA
?adonis2

status <- (rep(c("native","nonnative"), c(53,25)))
spt.m.traits2 <- data.frame(spt.m.traits, status = status)

#checking homogeneity 
resu.permdisp <- betadisper(gow.mt, spt.m.traits2$status)
resu.permdisp
permutest(resu.permdisp)

#permanova
resu.permanova <- adonis2(gow.mt ~ spt.m.traits2$status)
resu.permanova
#there is already a p value
#report degrees of freedom, F-value and p-value
