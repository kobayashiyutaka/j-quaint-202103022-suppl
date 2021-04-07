library(ggplot2)
library(geosphere)
library(philentropy)
library(ggrepel)
library(sf)
library(xtable)
library(gridExtra)
library(scales)

set.seed(123)

########################
## GRAPHIC COMPONENTS ##
########################

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalette_biplot <- alpha(cbPalette[c(7,5,3,4)],0.8)

## Plot themes

plt_theme1 <- theme(panel.background=element_rect(fill='transparent'),
                    panel.border=element_rect(fill='transparent',
                                              color='black'),
                    panel.grid.major=element_line(color=NA,
                                                  size=0.5),
                    panel.grid.minor=element_line(color=NA,
                                                  size=0.25),
                    axis.text=element_text(size=12,
                                           color='black'))

map_theme1 <-
  theme(panel.background = element_rect(fill='transparent'),
        panel.border = element_rect(fill='transparent',
                                    color='black'),
        axis.title = element_blank())

# Reading the world-map shape file
lakes_highres <- rgdal::readOGR('ne_50m_lakes/ne_50m_lakes.shp')
lakes_lowres <- rgdal::readOGR('ne_110m_lakes/ne_110m_lakes.shp')
world_highres <- rgdal::readOGR('ne_50m_land/ne_50m_land.shp')
world_lowres <- rgdal::readOGR('ne_110m_land/ne_110m_land.shp')
is_eurasia <- vector(length=nrow(world_highres))
is_eurasia_lowres <- vector(length=nrow(world_lowres))
is_eurasia_lakes <- vector(length=nrow(lakes_highres))
is_eurasia_lakes_lowres <- vector(length=nrow(lakes_lowres))

# Extracting Eurasia from the world map (high resolution)
for(i in 1:length(is_eurasia)) {
  bbox <- world_highres[i,]@bbox
  if(bbox[1,1] > -20 &&
     bbox[2,1] > -48) {
    is_eurasia[i] <- TRUE
  }
}

# Extracting Eurasia from the world map (low resolution)
for(i in 1:length(is_eurasia_lowres)) {
  bbox <- world_lowres[i,]@bbox
  if(bbox[1,1] > -20 &&
     bbox[2,1] > -48) {
    is_eurasia_lowres[i] <- TRUE
  }
}

# Removing small lakes (high resolution)
cut_area <- quantile(areaPolygon(lakes_highres),probs=c(0.9))
for(i in 1:length(is_eurasia_lakes)) {
  bbox <- lakes_highres[i,]@bbox
  lake_area <- areaPolygon(lakes_highres[i,])
  if(bbox[1,1] > -20 &&
     bbox[2,1] > -48 &&
     lake_area > cut_area) {
    is_eurasia_lakes[i] <- TRUE
  } else {
    is_eurasia_lakes[i] <- FALSE
  }
}

# Removing small lakes (low resolution)
for(i in 1:length(is_eurasia_lakes_lowres)) {
  bbox <- lakes_lowres[i,]@bbox
  lake_area <- areaPolygon(lakes_lowres[i,])
  if(bbox[1,1] > -20 &&
     bbox[2,1] > -48 &&
     lake_area > cut_area) {
    is_eurasia_lakes_lowres[i] <- TRUE
  } else {
    is_eurasia_lakes_lowres[i] <- FALSE
  }
}

eurasia <- world_highres[is_eurasia,]
eurasia_lowres <- world_lowres[is_eurasia_lowres,]
eurasia_lakes <- lakes_highres[is_eurasia_lakes,]
eurasia_lakes_lowres <- lakes_lowres[is_eurasia_lakes_lowres,]

# Converting SPDF (SpatialPolygonDataFrame) to simple features (sf)
eurasia_sf <- st_as_sf(eurasia)
eurasia_lowres_sf <- st_as_sf(eurasia_lowres)
eurasia_lakes_sf <- st_as_sf(eurasia_lakes)
eurasia_lakes_lowres_sf <- st_as_sf(eurasia_lakes_lowres)

##################
## READING DATA ##
##################

## Reading the Paleo-Asia dataset
paleo <- read.csv('PADBMode200801_excluding_Africa.csv',
                  header=T,na.strings='',stringsAsFactors = F)
## Removing incomplete cases
## But actually there are no such cases in the latest dataset, 
## so that the following code just copies paleo to paleo_comp.
paleo_comp <- paleo[complete.cases(paleo),]

## Column numbers for mode information
mode_seq <- 6:29
mode_names <- names(paleo_comp)[mode_seq]

# Computing the number of mode per layer
paleo_comp$mode_num <- apply(paleo_comp[,mode_seq],1,sum)

####################
## REGION SYSTEMS ##
####################

## Defining regions and subregions
region_names <- c('East Africa','Levant',
                  'Arabia','South Caucasus','Zagros',
                  'Central Asia','South Asia','South-East Asia',
                  'East Asia','Russia')

west_asia <- c('Levant','Arabia','South Caucasus','Zagros')

EastAfrica <- c('Ethiopia','Kenya','Tanzania, United Republic of')
Levant <- c('Syrian Arab Republic',
            'Lebanon','Jordan','Israel','Egypt')
Arabia <- c('Yemen','United Arab Emirates')
SouthCaucasus <- c('Armenia','Azerbaijan','Georgia')
Zagros <- c('Iran, Islamic Republic of','Iraq','Turkey')
CentralAsia <- c('Uzbekistan','Tajikistan')
SouthAsia <- c('Afghanistan','India','Pakistan','Sri Lanka')
SouthEastAsia <- c('Indonesia','Malaysia','Philippines',
                   'Thailand','Timor-Leste','Viet Nam')
EastAsia <- c('China','Korea, Republic of')
Russia <- c('Russian Federation')

paleo_comp$Region <- NA
paleo_comp$Region[paleo_comp$Country %in% EastAfrica] <- 'East Africa'
paleo_comp$Region[paleo_comp$Country %in% Levant] <- 'Levant'
paleo_comp$Region[paleo_comp$Country %in% Arabia] <- 'Arabia'
paleo_comp$Region[paleo_comp$Country %in% SouthCaucasus] <- 'South Caucasus'
paleo_comp$Region[paleo_comp$Country %in% Zagros] <- 'Zagros'
paleo_comp$Region[paleo_comp$Country %in% CentralAsia] <- 'Central Asia'
paleo_comp$Region[paleo_comp$Country %in% SouthAsia] <- 'South Asia'
paleo_comp$Region[paleo_comp$Country %in% SouthEastAsia] <- 'South-East Asia'
paleo_comp$Region[paleo_comp$Country %in% EastAsia] <- 'East Asia'
paleo_comp$Region[paleo_comp$Country %in% Russia] <- 'Russia'
paleo_comp <- transform(paleo_comp,
                        Region_f=factor(paleo_comp$Region))

## Five-region system
paleo_comp$Region5 <- ifelse(paleo_comp$Region=='Russia','North Asia',
                      ifelse(paleo_comp$Region=='East Asia','East Asia',
                      ifelse(paleo_comp$Region %in%
                            c('South Asia','South-East Asia'),
                            'South/South-East Asia',
                      ifelse(paleo_comp$Region=='East Africa',
                            'East Africa',
                            'West/Central Asia'))))

paleo_comp$Region5_f <- factor(paleo_comp$Region5,
                               levels=c('East Asia',
                                        'South/South-East Asia',
                                        'North Asia',
                                        'West/Central Asia',
                                        'East Africa'))

## Seven-region system
paleo_comp$Region7 <- ifelse(paleo_comp$Region %in% west_asia,
                             'West Asia',
                             paleo_comp$Region)

paleo_comp$Region7_f <- factor(paleo_comp$Region7)

## Want to include Africa ? 
## If yes, comment out the following
paleo_comp <- paleo_comp[paleo_comp$Region5 != 'East Africa',]

# Mode frequencies (popularities)
mode_freq <- apply(paleo_comp[,mode_seq],2,sum)
mode_freq_MP <- apply(paleo_comp[paleo_comp$Era=='MP',mode_seq],2,sum)
mode_freq_MUP <- apply(paleo_comp[paleo_comp$Era=='MUP',mode_seq],2,sum)
mode_freq_UP <- apply(paleo_comp[paleo_comp$Era=='UP',mode_seq],2,sum)

## PLOT
## Plotting the relative popularities of modes
# in Periods 1-3 separately

# Preparing a data frame
freq_df <- data.frame(mode=rep(mode_names,3),
                       Era=rep(c('MP','MUP','UP'),
                                  each=length(mode_names)),
                       freq=c(mode_freq_MP/sum(paleo_comp$Era=='MP'),
                              mode_freq_MUP/sum(paleo_comp$Era=='MUP'),
                              mode_freq_UP/sum(paleo_comp$Era=='UP'))
)

freq_df$Era <- factor(freq_df$Era,
                          ordered=T,
                          levels=c('MP','MUP','UP'))

# Relative popularities for all modes
plt_popularity <- ggplot(freq_df,
                          aes(x=mode,y=freq)) +
  geom_bar(stat='identity',width=0.8,color="#444444",size=0.4,
           aes(fill=Era),position='dodge') +
  scale_fill_viridis_d(limits=c('MP','MUP','UP'),
                       labels=1:3) +
  scale_x_discrete(limits=mode_names)+
  scale_y_continuous(limits=c(0,1.0),expand=c(0,0)) +
  theme_bw() +
  theme(
    legend.key.size=unit(1,'cm'),
    text=element_text(family='serif',size=25),
    axis.title.x=element_text(margin=margin(t=10,r=0,l=0,b=0,unit='pt')),
    axis.title.y=element_text(margin=margin(t=0,r=10,l=0,b=0,unit='pt')),
    axis.text.x=element_text(size=16),
    plot.margin=unit(c(2,2,2,2),'cm'),
    aspect.ratio=0.5
  )+
  labs(fill='Period',x='Mode',
       y='Frequency') 

## Figure 1 of Nishiaki et al. (A4, Landscape)
## To plot Figure 1, uncomment the following code:
#plt_popularity

## Spatial distributions of all modes

# Plot theme
mode_dist_theme <-
  theme(legend.position='none',
        text=element_text(family='serif',
                          size=14),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid=element_blank()
  ) 

# All periods
plt_mode_distrib <- list()
for(mode in mode_names){
  plt_mode_distrib[[mode]] <- 
    ggplot(eurasia_lowres_sf) + geom_sf(color=NA,
                                 fill='white',
                                 size=0.4) + 
    geom_sf(data=eurasia_lakes_lowres_sf,color=NA,
            fill='skyblue',size=0.4) +
    map_theme1 +
    theme(panel.background=element_rect(fill='skyblue'))+
    geom_point(data=paleo_comp[paleo_comp[[mode]]==1,],
               aes(x=Longitude,y=Latitude),
               col='black',alpha=0.5,size=1) +
    xlim(c(30,150)) + ylim(c(-10,80)) +
    mode_dist_theme +
    labs(title=paste(mode,':',mode_freq[mode]))
}

# Figure S1:4 of Nishiaki et al.
## To plot Figure S1:4, uncomment the following code:
#grid.arrange(grobs=plt_mode_distrib,nrow=4)

plt_mode_distribUP <- list()
for(mode in mode_names){
  plt_mode_distribUP[[mode]] <- 
    ggplot(eurasia_lowres_sf) + geom_sf(color=NA,
                                 fill='white',
                                 size=0.4) +
    geom_sf(data=eurasia_lakes_lowres_sf,color=NA,
            fill='skyblue',size=0.4) +
    map_theme1 +
    theme(panel.background=element_rect(fill='skyblue'))+
    geom_point(data=paleo_comp[paleo_comp[[mode]]==1 &
                                 (paleo_comp$Era=='UP'),],
               aes(x=Longitude,y=Latitude),
               col='black',alpha=0.5,size=1) +
    ylim(c(-10,80)) +
    xlim(c(30,150)) + 
    mode_dist_theme +
    labs(title=paste(mode,':',mode_freq_UP[mode]))
}

## Figure S1:3 of Nishiaki et al.
## To plot Figure S1:3, uncomment the following code:
#grid.arrange(grobs=plt_mode_distribUP,
#             nrow=4)

plt_mode_distribMP <- list()
for(mode in mode_names){
  plt_mode_distribMP[[mode]] <- 
    ggplot(eurasia_lowres_sf) + geom_sf(color=NA,
                                 fill='white',
                                 size=0.4) +
    geom_sf(data=eurasia_lakes_lowres_sf,color=NA,
            fill='skyblue',size=0.4) +
    map_theme1 +
    theme(panel.background=element_rect(fill='skyblue'))+
    geom_point(data=paleo_comp[paleo_comp[[mode]]==1 &
                                 (paleo_comp$Era=='MP'),],
               aes(x=Longitude,y=Latitude),
               col='black',alpha=0.5,size=1) +
    mode_dist_theme +ylim(c(-10,80)) +
    xlim(c(30,150)) + 
    labs(title=paste(mode,':',mode_freq_MP[mode]))
}
## Figure S1:1 of Nishiaki et al.
## To plot Figure S1:1, uncomment the following code:
#grid.arrange(grobs=plt_mode_distribMP,nrow=4)

plt_mode_distribMUP <- list()
for(mode in mode_names){
  plt_mode_distribMUP[[mode]] <- 
    ggplot(eurasia_lowres_sf) + geom_sf(color=NA,
                                 fill='white',
                                 size=0.4) + 
    geom_sf(data=eurasia_lakes_lowres_sf,color=NA,
            fill='skyblue',size=0.4) +
    map_theme1 +
    theme(panel.background=element_rect(fill='skyblue'))+
    geom_point(data=paleo_comp[paleo_comp[[mode]]==1 &
                                 (paleo_comp$Era=='MUP'),],
               aes(x=Longitude,y=Latitude),
               col='black',alpha=0.5,size=1) +
    mode_dist_theme +ylim(c(-10,80)) +
    xlim(c(30,150)) + 
    labs(title=paste(mode,':',mode_freq_MUP[mode]))
}

## Figure S1:2 of Nishiaki et al.
## To plot Figure S1:2, uncomment the following code:
#grid.arrange(grobs=plt_mode_distribMUP,nrow=4)

## unscaled PCA

pca_res <- prcomp(paleo_comp[,mode_seq])
paleo_comp$PC1 <- pca_res$x[,1]
paleo_comp$PC2 <- pca_res$x[,2]

pr_scores <- data.frame(pca_res$x[,1:2])
names(pr_scores) <- c('PC1','PC2')

# Computing raw loadings (not standardized by SD)
pr_loadings <- pca_res$rotation * rep(pca_res$sdev,each=24)

## Make raw loadings 3 times bigger
pr_loadings_x3 <- data.frame(3*pr_loadings[,1:2])

## Standard deviations of modes
sd_modes <- apply(paleo_comp[,mode_names],2,sd)

## Standardized loadings (will be 1.5 times magnified)
pr_sd_loadings <- data.frame(pr_loadings[,1:2] / sd_modes)

pr_loadings_plot <- pr_sd_loadings*1.5
## Want to use raw loadings instead of standardized loadings??
## If yes, uncomment the following line:
#pr_loadings_plot <- pr_loadings_x3

# data frame for drawing arrows
arrow_loadings <- data.frame(matrix(data=0,nrow=24,ncol=2))
names(arrow_loadings) <- c('origin_x','origin_y')
arrow_loadings <- cbind(arrow_loadings,pr_loadings_plot)

# Biplot (only for MP)
pca_ylimits <- c(-1.7,1.11)
pca_xlimits <- c(-1.8,1.31)

pca_vars <- pca_res$sdev * pca_res$sdev 
pca_explained <- pca_vars/sum(pca_vars)

biplot_region_labels <- c('E. Asia','S./S.E. Asia','N. Asia',
                          'W./C. Asia','E. Africa')
biplot_region_labels_full <- c('East Asia',
                               'South/Southeast Asia',
                               'North Asia',
                               'West/Central Asia',
                               'East Africa')

biplot_theme <-  theme(text=element_text(family='serif',
                                         size=16),
                       legend.title=element_text(size=14),
                       plot.margin=margin(c(0.5,0.5,0.5,0.5),
                                          unit='cm'))

plt_PCAbiplotMP <- ggplot(pr_scores[paleo_comp$Era=='MP',],aes(x=PC1,y=PC2)) +
  geom_vline(xintercept=0,col='gray')+
  geom_hline(yintercept=0,col='gray')+
  geom_count(
    aes(color=paleo_comp[paleo_comp$Era=='MP',]$Region5_f))+
  geom_point(data=pr_loadings_plot,size=0.5,color='tomato')+
  scale_size_area(max_size=8)  +
  scale_color_manual('Regions',
                     values=cbPalette_biplot,
                     label=biplot_region_labels)+
  geom_segment(data=arrow_loadings,
               aes(x=origin_x,y=origin_y,
                   xend=PC1,yend=PC2),color='tomato',
               arrow=arrow(length=unit(0.03,'npc'))) +
  scale_y_continuous(limits=pca_ylimits)+
  scale_x_continuous(limits=pca_xlimits) +
  labs(x=paste('PC1 (',
               format(100*round(pca_explained[1],3),nsmall=1),'%)',sep=''),
       y=paste('PC2 (',100*round(pca_explained[2],3),'%)',sep=''))+
  geom_text_repel(data=pr_loadings_plot,aes(label=mode_names),family='serif',
                  size=4,max.overlaps=20)+
  guides(size=guide_legend('Count')) +
  guides(color = guide_legend(override.aes = list(size=4)))+
  plt_theme1 + biplot_theme +
  annotate("text",x=1.2,y=-1.5,label='1',size=8) +
  coord_fixed() 

plt_PCAbiplotUP <- ggplot(pr_scores[paleo_comp$Era=='UP',],aes(x=PC1,y=PC2)) +
  geom_vline(xintercept=0,col='gray')+
  geom_hline(yintercept=0,col='gray')+
  geom_count(
    aes(color=paleo_comp[paleo_comp$Era=='UP',]$Region5_f)) +
  geom_point(data=pr_loadings_plot,size=0.5,color='tomato') + 
  scale_size_area(max_size=8,breaks=seq(5,25,10)) + 
  scale_color_manual('Regions',values=cbPalette_biplot,
                     label=biplot_region_labels)+
  geom_segment(data=arrow_loadings,
               aes(x=origin_x,y=origin_y,
                   xend=PC1,yend=PC2),color='tomato',
               arrow=arrow(length=unit(0.03,'npc'))) +
  scale_y_continuous(limits=pca_ylimits)+
  scale_x_continuous(limits=pca_xlimits)+
  labs(x=paste('PC1 (',
               format(100*round(pca_explained[1],3),nsmall=1),'%)',sep=''),
       y=paste('PC2 (',100*round(pca_explained[2],3),'%)',sep=''))+
  geom_text_repel(data=pr_loadings_plot,aes(label=mode_names),family='serif',
                  size=4,max.overlaps = 20)+
  guides(size=guide_legend('Count')) +
  guides(color = guide_legend(override.aes = list(size=4)))+
  plt_theme1 + biplot_theme +
  annotate("text",x=1.2,y=-1.5,label='3',size=8) +
  coord_fixed() 

plt_PCAbiplotMUP <- ggplot(pr_scores[paleo_comp$Era=='MUP',],aes(x=PC1,y=PC2)) +
  geom_vline(xintercept=0,col='gray')+
  geom_hline(yintercept=0,col='gray')+
  geom_count(
    aes(color=paleo_comp[paleo_comp$Era=='MUP',]$Region5_f)) +
  geom_point(data=pr_loadings_plot,size=0.5,color='tomato') + 
  scale_size_area(max_size=8,breaks=c(2,5,8)) + 
  scale_color_manual('Regions',values=cbPalette_biplot,
                     label=biplot_region_labels)+
  geom_segment(data=arrow_loadings,
               aes(x=origin_x,y=origin_y,
                   xend=PC1,yend=PC2),color='tomato',
               arrow=arrow(length=unit(0.03,'npc'))) +
  scale_y_continuous(limits=pca_ylimits)+
  scale_x_continuous(limits=pca_xlimits)+
  labs(x=paste('PC1 (',
               format(100*round(pca_explained[1],3),nsmall=1),'%)',sep=''),
       y=paste('PC2 (',100*round(pca_explained[2],3),'%)',sep=''))+
  geom_text_repel(data=pr_loadings_plot,aes(label=mode_names),family='serif',
                  size=4,max.overlaps = 20)+
  guides(size=guide_legend('Count'))+
  guides(color = guide_legend(override.aes = list(size=4)))+
  plt_theme1 + biplot_theme +
  annotate("text",x=1.2,y=-1.5,label='2',size=8) +
  coord_fixed() 

## Figure 4 of Nishiaki et al. (A4, Portrait)
## To plot Figure 4, uncomment the following code:
#grid.arrange(grobs=list(plt_PCAbiplotMP,
#                       plt_PCAbiplotMUP,
#                       plt_PCAbiplotUP),nrow=3)

plt_PCAbiplot <- ggplot(pr_scores,aes(x=PC1,y=PC2)) +
  geom_vline(xintercept=0,col='gray')+
  geom_hline(yintercept=0,col='gray')+
  geom_count(
    aes(color=paleo_comp$Region5_f)) +
  geom_point(data=pr_loadings_plot,size=0.5,color='tomato')+
  scale_size_area(max_size=11) + 
  scale_color_manual('Regions',
                     values=cbPalette_biplot,
                     label=biplot_region_labels)+
  geom_segment(data=arrow_loadings,
#               size=0.7,
               aes(x=origin_x,y=origin_y,
                   xend=PC1,yend=PC2),color='tomato',
               arrow=arrow(length=unit(0.03,'npc'))) +
  scale_y_continuous(limits=pca_ylimits)+
  scale_x_continuous(limits=pca_xlimits)+
  labs(x=paste('PC1 (',
               format(100*round(pca_explained[1],3),nsmall=1),'%)',sep=''),
       y=paste('PC2 (',100*round(pca_explained[2],3),'%)',sep=''))+
  geom_text_repel(data=pr_loadings_plot,aes(label=mode_names),family='serif',
#                  fontface='bold',
                  size=4,max.overlaps = 100)+
  guides(size=guide_legend('Count'))+
  plt_theme1 + theme(text=element_text(family='serif',
                                       size=18),
                     legend.key.size=unit(0.8,'cm'),
                     legend.title=element_text(size=18),
                     axis.text=element_text(family='serif',
                                            size=16),
                     axis.title.x=element_text(margin=margin(t=10,r=0,l=0,b=0,unit='pt')),
                     plot.margin=unit(c(2,2,2,2),unit='cm'))+
  guides(color = guide_legend(override.aes = list(size=4)))+
  coord_fixed() 

## Figure 3 of Nishiaki et al. (A4, Portrait)
## To plot Figure 3, uncomment the following code:
#plt_PCAbiplot

## SCALED PCA
scpca_res <- prcomp(paleo_comp[,mode_names],scale.=T)
scpr_scores <- data.frame(scpca_res$x[,1:2])
names(scpr_scores) <- c('PC1','PC2')
scpr_loadings <- scpca_res$rotation * rep(scpca_res$sdev,each=24)
paleo_comp[which(scpr_scores$PC2 < -5),]

## Make loadings 4 times bigger
scpr_loadings_x4 <- data.frame(4*scpr_loadings[,1:2])
# data frame for drawing arrows
scarrow_loadings <- data.frame(matrix(data=0,nrow=24,ncol=2))
names(scarrow_loadings) <- c('origin_x','origin_y')
scarrow_loadings <- cbind(scarrow_loadings,scpr_loadings_x4)

#scpca_ylimits <- c(-1.7,1.1)
#scpca_xlimits <- c(-1.8,1.31)

scpca_vars <- scpca_res$sdev * scpca_res$sdev 
scpca_explained <- scpca_vars/sum(scpca_vars)

plt_scPCAbiplot <- ggplot(scpr_scores,aes(x=PC1,y=PC2)) +
  geom_vline(xintercept=0,col='gray')+
  geom_hline(yintercept=0,col='gray')+
  geom_count(
    aes(color=paleo_comp$Region5_f)) +
  geom_point(data=scpr_loadings_x4,size=0.5,color='tomato') + 
  scale_size_area(max_size=11) + 
  scale_color_manual('Regions',values=cbPalette_biplot,
                     label=biplot_region_labels)+
  geom_segment(data=scarrow_loadings,
               aes(x=origin_x,y=origin_y,
                   xend=PC1,yend=PC2),color='tomato',
               arrow=arrow(length=unit(0.03,'npc'))) +
#  scale_y_continuous(limits=scpca_ylimits)+
#  scale_x_continuous(limits=scpca_xlimits)+
  labs(x=paste('PC1 (',
               format(100*round(scpca_explained[1],3),nsmall=1),'%)',sep=''),
       y=paste('PC2 (',100*round(scpca_explained[2],3),'%)',sep=''))+
  geom_text_repel(data=scpr_loadings_x4,aes(label=mode_names),
                  family='serif',size=4)+
  guides(size=guide_legend('Count'))+
  guides(color = guide_legend(override.aes = list(size=4)))+
  plt_theme1 + theme(text=element_text(family='serif',
                                       size=18),
                     legend.key.size=unit(0.8,'cm'),
                     legend.title=element_text(size=18),
                     axis.text=element_text(family='serif',
                                            size=16),
                     axis.title.x=element_text(margin=margin(t=10,r=0,l=0,b=0,unit='pt')),
                     plot.margin=unit(c(2,2,2,2),unit='cm'))+
  coord_fixed() 

## Figure S2 of Nishiaki et al. (A4, Portrait)
## To plot Figure S2, uncomment the following code:
#plt_scPCAbiplot

## PLOT
## Geographic representation of PC scores
geo_pca <- data.frame(pca_res$x[,c(1,2)])
names(geo_pca) <- c('PC1','PC2')
geo_pca <- cbind(geo_pca,
                 paleo_comp[,c('Era','Latitude','Longitude')])

pc2min <- min(geo_pca[,'PC2'])
pc1min <- min(geo_pca[,'PC1'])

pca_map_xlimits <- c(20,160)
pca_map_ylimits <- c(-20,80)
pc2_limits <- c(-1.6,1.0)
pc1_limits <- c(-1.8,1.3)
pc1_breaks <- seq(-1.6,1.2,0.4)
pc2_breaks <- seq(-1.6,1.0,0.4)

geo_theme <- theme(text=element_text(size=16,family='serif'),
      legend.title=element_blank(),
      legend.position='right',
      legend.key.height = unit(0.1,units='npc'),
      title=element_text(size=16),
      axis.text=element_text(size=10),
      plot.margin=margin(c(0.5,0.5,0.5,0.5),unit='cm'))

plt_geo_pc2UP <- ggplot(eurasia_sf) + geom_sf(color='#999999',
                                              fill='white',
                                           size=0.2) + 
  geom_sf(data=eurasia_lakes_sf,color='#999999',
          fill='white',size=0.2) +
          map_theme1 +
  geom_point(
    data=geo_pca[paleo_comp$Era=='UP',],
    aes(x=Longitude,
        y=Latitude,
    color=geo_pca[paleo_comp$Era=='UP','PC2']),
    size=2, 
    shape=16,alpha=0.8)+
  scale_color_viridis_c(limits=pc2_limits,
                        breaks=pc2_breaks,
                        oob=squish)+
  geo_theme +
#  labs(title='(f) PC2 scores in Period 3') +
  annotate('text',x=90,y=-15,label='3',size=9) +
  scale_x_continuous(limits=pca_map_xlimits)+
  scale_y_continuous(limits=pca_map_ylimits) 

plt_geo_pc2MP <- ggplot(eurasia_sf) + geom_sf(color='#999999',
                                              fill='white',
                                              size=0.2) + 
  geom_sf(data=eurasia_lakes_sf,color='#999999',
          fill='white',size=0.2) +
  map_theme1 +
  geom_point(
    data=geo_pca[paleo_comp$Era=='MP',],
    aes(x=Longitude,
        y=Latitude,
        color=geo_pca[paleo_comp$Era=='MP','PC2']),
    size=2,#stroke=0.05,color="#444444",
    shape=16,alpha=0.8) +
  scale_color_viridis_c(limits=pc2_limits,
                      breaks=pc2_breaks,
                      oob=squish)+
#  labs(title='(d) PC2 scores in Period 1')+
  annotate('text',x=90,y=-15,label='1',size=9) +
  geo_theme +
  scale_x_continuous(limits=pca_map_xlimits) + 
  scale_y_continuous(limits=pca_map_ylimits)

plt_geo_pc2MUP <- ggplot(eurasia_sf) + geom_sf(color='#999999',
                                              fill='white',
                                              size=0.2) + 
  geom_sf(data=eurasia_lakes_sf,color='#999999',
          fill='white',size=0.2) +
  map_theme1 +
  geom_point(
    data=geo_pca[paleo_comp$Era=='MUP',],
    aes(x=Longitude,
        y=Latitude,
        color=geo_pca[paleo_comp$Era=='MUP','PC2']),
    size=2,#stroke=0.05,color="#444444",
    shape=16,alpha=0.8) +
  scale_color_viridis_c(limits=pc2_limits,
                        breaks=pc2_breaks,
                        oob=squish)+
#  labs(title='(e) PC2 scores in Period 2')+
  annotate('text',x=90,y=-15,label='2',size=9) +
  geo_theme +
  scale_x_continuous(limits=pca_map_xlimits)+
  scale_y_continuous(limits=pca_map_ylimits)

# Figure 6 of Nishiaki et al. (A4, Portrait)
# To plot Figure 6, uncomment the following code:
#grid.arrange(grobs=list(plt_geo_pc2MP,
#                        plt_geo_pc2MUP,
#                        plt_geo_pc2UP),
#             nrow=3)

## PC1

plt_geo_pc1UP <- ggplot(eurasia_sf) + geom_sf(color='#999999',
                                              fill='white',
                                              size=0.2) + 
  geom_sf(data=eurasia_lakes_sf,color='#999999',
          fill='white',size=0.2) +
  map_theme1 +
  geom_point(
    data=geo_pca[paleo_comp$Era=='UP',],
    aes(x=Longitude,
        y=Latitude,
        color=geo_pca[paleo_comp$Era=='UP','PC1']),
    size=2,
    shape=16,alpha=0.8)+
  scale_color_viridis_c(limits=pc1_limits,
                        breaks=pc1_breaks,
                        oob=squish)+
  geo_theme +
#  labs(title='(c) PC1 scores in Period 3') +
  annotate("text",x=90,y=-15,label='3',size=9)+
  scale_x_continuous(limits=pca_map_xlimits)+
  scale_y_continuous(limits=pca_map_ylimits)

plt_geo_pc1MP <- ggplot(eurasia_sf) + geom_sf(color='#999999',
                                              fill='white',
                                              size=0.2) +  
  geom_sf(data=eurasia_lakes_sf,color='#999999',
          fill='white',size=0.2) +
  map_theme1 +
  geom_point(
    data=geo_pca[paleo_comp$Era=='MP',],
    aes(x=Longitude,
        y=Latitude,
        color=geo_pca[paleo_comp$Era=='MP','PC1']),
    size=2,#stroke=0.05,color="#444444",
    shape=16,alpha=0.8) +
  scale_color_viridis_c(limits=pc1_limits,
                        breaks=pc1_breaks,
                        oob=squish)+
#  labs(title='(a) PC1 scores in Period 1')+
  annotate("text",x=90,y=-15,label='1',size=9)+
  geo_theme +
  scale_x_continuous(limits=pca_map_xlimits) + 
  scale_y_continuous(limits=pca_map_ylimits)

plt_geo_pc1MUP <- ggplot(eurasia_sf) + geom_sf(color='#999999',
                                               fill='white',
                                               size=0.2) +
  geom_sf(data=eurasia_lakes_sf,color='#999999',
          fill='white',size=0.2) +
  map_theme1 +
  geom_point(
    data=geo_pca[paleo_comp$Era=='MUP',],
    aes(x=Longitude,
        y=Latitude,
        color=geo_pca[paleo_comp$Era=='MUP','PC1']),
    size=2,#stroke=0.05,color="#444444",
    shape=16,alpha=0.8) +
  scale_color_viridis_c(limits=pc1_limits,oob=squish,
                        breaks=pc1_breaks)+
#  labs(title='(b) PC1 scores in Period 2')+
  annotate("text",x=90,y=-15,label='2',size=9) +
  geo_theme +
  scale_x_continuous(limits=pca_map_xlimits)+
  scale_y_continuous(limits=pca_map_ylimits)

## Figure 5 of Nishiaki et al. (A4, Portrait)
## To plot Figure 5, uncomment the following code:
#grid.arrange(grobs=list(plt_geo_pc1MP,
#                        plt_geo_pc1MUP,
#                        plt_geo_pc1UP),nrow=3)

# Table 3 (Regional difference in recorded mode number)
modnum_summary <- aggregate(paleo_comp$mode_num,by=list(paleo_comp$Region),FUN=summary)
layno_region <- table(paleo_comp$Region)
modnum_summary_df <- data.frame(Region=modnum_summary[,1],layno=as.vector(layno_region),round(modnum_summary[,-1],2))[,c(1,2,3,5,6,8)]
modnum_total <- summary(paleo_comp$mode_num)
names(modnum_summary_df) <- c('Region','#Assemblages',
                              'Min','Median','Mean','Max')
modnum_summary_df$Region <- factor(modnum_summary_df$Region,
                                   levels=c(levels(modnum_summary_df$Region),'Total'))
modnum_summary_df <- rbind(modnum_summary_df,
      c('Total',sum(modnum_summary_df[,'#Assemblages']),round(modnum_total,2)[c(1,3,4,6)]))

write(print(xtable(modnum_summary_df,
                   digits=c(0,0,0,0,0,2,0)),include.rownames=F),
      file='table3.tex')

# Table 4 (Spatial structure)
# Distance matrix calculation
# Geographic (crow-fly distance)
d <- distm(paleo_comp[,c('Longitude','Latitude')],fun=distGeo)
# Hamming (Manhattan) distance
d_hamming <- distance(paleo_comp[,mode_seq],method='manhattan') 
within_region <- outer(paleo_comp$Region7,paleo_comp$Region7,FUN='==')


## Distances between layers

dMP <- d[paleo_comp$Era=='MP',paleo_comp$Era=='MP']
dUP <- d[paleo_comp$Era=='UP',paleo_comp$Era=='UP']
dMUP <- d[paleo_comp$Era=='MUP',paleo_comp$Era=='MUP']

dMP_hamming <- d_hamming[paleo_comp$Era=='MP',paleo_comp$Era=='MP']
dUP_hamming <- d_hamming[paleo_comp$Era=='UP',paleo_comp$Era=='UP']
dMUP_hamming <- d_hamming[paleo_comp$Era=='MUP',paleo_comp$Era=='MUP']

within_regionMP <- outer(paleo_comp[paleo_comp$Era=='MP','Region7'],
                           paleo_comp[paleo_comp$Era=='MP','Region7'],
                           FUN='==')
within_regionMUP <- outer(paleo_comp[paleo_comp$Era=='MUP','Region7'],
                           paleo_comp[paleo_comp$Era=='MUP','Region7'],
                           FUN='==')
within_regionUP <- outer(paleo_comp[paleo_comp$Era=='UP','Region7'],
                           paleo_comp[paleo_comp$Era=='UP','Region7'],
                           FUN='==')

dMP_df <- data.frame(dMP=dMP[upper.tri(dMP)]/10^6,
                     dMP_hamming=dMP_hamming[upper.tri(dMP_hamming)],
                     within_region=within_regionMP[upper.tri(within_regionMP)])
dUP_df <- data.frame(dUP=dUP[upper.tri(dUP)]/10^6,
                     dUP_hamming=dUP_hamming[upper.tri(dUP_hamming)],
                     within_region=within_regionUP[upper.tri(within_regionUP)])
dMUP_df <- data.frame(dMUP=dMUP[upper.tri(dMUP)]/10^6,
                     dMUP_hamming=dMUP_hamming[upper.tri(dMUP_hamming)],
                     within_region=within_regionMUP[upper.tri(within_regionMUP)])
d_df <- data.frame(d=d[upper.tri(d)]/10^6,
                   d_hamming=d_hamming[upper.tri(d_hamming)],
                   within_region=within_region[upper.tri(within_region)])

geom_ham_cor <- list(
  MP=cor(dMP_df,method='spearman'),
  UP=cor(dUP_df,method='spearman'),
  MUP=cor(dMUP_df,method='spearman'),
  total=cor(d_df,method='spearman')
)

sum(paleo_comp$Region7=='Russia')/895
## Computation of within-region correlation
round(cor(d_df[d_df$within_region,1:2]),2)

## Computation of correlation excluding inter-regional pairs
between_EastWest <- outer(paleo_comp$Region7=='West Asia',paleo_comp$Region7=='East Asia') |
  outer(paleo_comp$Region7=='East Asia',paleo_comp$Region7=='West Asia',FUN='&')
between_EastRussia <- outer(paleo_comp$Region7=='Russia',paleo_comp$Region7=='East Asia') |
  outer(paleo_comp$Region7=='East Asia',paleo_comp$Region7=='Russia',FUN='&')
between_WestRussia <- outer(paleo_comp$Region7=='West Asia',paleo_comp$Region7=='Russia') |
  outer(paleo_comp$Region7=='Russia',paleo_comp$Region7=='West Asia',FUN='&')

d_df_exc <- data.frame(d=d[upper.tri(d)]/10^6,
                       d_hamming=d_hamming[upper.tri(d_hamming)],
                       between_EastWest=between_EastWest[upper.tri(between_EastWest)],
                       between_EastRussia=between_EastRussia[upper.tri(between_EastRussia)],
                       between_WestRussia=between_WestRussia[upper.tri(between_WestRussia)])

round(cor(d_df_exc[!d_df_exc$between_EastWest,1:2],method='spearman'),2)
round(cor(d_df_exc[!d_df_exc$between_EastRussia,1:2],method='spearman'),2)
round(cor(d_df_exc[!d_df_exc$between_WestRussia,1:2],method='spearman'),2)

## Computation of the distribution of geographic distances
mean(d_df[d_df$within_region,1])
mean(d_df_exc[d_df_exc$between_EastWest,1])
mean(d_df_exc[d_df_exc$between_WestRussia,1])
mean(d_df_exc[d_df_exc$between_EastRussia,1])

### 
d_hamming_summary <- matrix(
c(quantile(dMP_df$dMP_hamming,probs=c(0.25,0.5,0.75,0.95,1)),
quantile(dMUP_df$dMUP_hamming,probs=c(0.25,0.5,0.75,0.95,1)),
quantile(dUP_df$dUP_hamming,probs=c(0.25,0.5,0.75,0.95,1)),
quantile(d_df$d_hamming,probs=c(0.25,0.5,0.75,0.95,1))),
nrow=4,byrow=T)
d_hamming_summary_df <- data.frame(d_hamming_summary)
names(d_hamming_summary_df) <- c('25%','Median','75%','95%','Max')
d_hamming_summary_df$Period <- c('1','2','3','All')
d_hamming_summary_df$Mean <-
  c(mean(dMP_df$dMP_hamming),
    mean(dMUP_df$dMUP_hamming),
    mean(dUP_df$dUP_hamming),
    mean(d_df$d_hamming))
d_hamming_summary_df[,'#Assemblages'] <-
  c(sum(paleo_comp$Era=='MP'),
    sum(paleo_comp$Era=='MUP'),
    sum(paleo_comp$Era=='UP'),
    nrow(paleo_comp))
d_hamming_summary_df <- d_hamming_summary_df[,c(6,8,1,2,7,3:5)]
d_hamming_summary_df$Corr <-
  format(round(c(geom_ham_cor$MP[1,2],
    geom_ham_cor$MUP[1,2],
    geom_ham_cor$UP[1,2],
    geom_ham_cor$total[1,2]),2),digits=2)
write(print(xtable(d_hamming_summary_df,
             digits=c(0,0,0,0,0,2,0,0,0,2)),include.rownames=F),
      file='table4.tex')

## Fig. 2 (histogram of Hamming distance)
hist_theme <-   theme(
                      axis.text=element_text(family='serif',size=16),
                      plot.title=element_text(hjust=0.5),
                      title=element_text(family='serif',size=22),
                      axis.title=element_text(family='serif',size=20),
                      legend.text=element_text(family='serif',size=12),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank(),
                      legend.position='bottom',
                      plot.margin=margin(c(0.5,0.5,0.5,0.5),unit='cm')) 

plt_hamming_hist <- ggplot(d_df,aes(x=d_hamming,fill=within_region)) +
  geom_bar(position=position_dodge2(preserve='single',reverse=T,
                                    padding=0.2),
           color=NA,width=0.8)+
  scale_fill_manual('Regions',limits=c(T,F),labels=c('Same','Different'),
                    values=cbPalette[c(7,3)],guide=F) +
  theme_bw() + hist_theme +
  labs(x='Cultural distance',y='Count',
       title='4')+
  scale_y_continuous(expand=c(0,0),limits=c(0,50000))

plt_hamming_histMP <- ggplot(dMP_df,aes(x=dMP_hamming,fill=within_region)) +
  geom_bar(position=position_dodge2(preserve='single',reverse=T,
                                   padding=0.2),
           color=NA,width=0.8)+
  scale_fill_manual('Regions',limits=c(T,F),labels=c('Same','Different'),
                    values=cbPalette[c(7,3)],guide=F)+
  theme_bw() + hist_theme +
  labs(x='Cultural distance',y='Count',
       title='1') +
  scale_y_continuous(expand=c(0,0),limits=c(0,8000)) 

plt_hamming_histMUP <- ggplot(dMUP_df,
                              aes(x=dMUP_hamming,fill=within_region)) +
  geom_bar(position=position_dodge2(preserve='single',reverse=T,
                                    padding=0.2),
                                   color=NA,width=0.8)+
  scale_fill_manual('Regions',limits=c(T,F),labels=c('Same','Different'),
                    values=cbPalette[c(7,3)],guide=F) +
  theme_bw() + hist_theme +
  labs(x='Cultural distance',y='Count',
       title='2')  +
  scale_y_continuous(expand=c(0,0),limits=c(0,700))

plt_hamming_histUP <- ggplot(dUP_df,
                             aes(x=dUP_hamming,fill=within_region)) +
  geom_bar(position=position_dodge2(preserve='single',reverse=T,
                                    padding=0.2),
           color=NA,width=0.8)+
  scale_fill_manual('Regions',limits=c(T,F),labels=c('Same','Different'),
                    values=cbPalette[c(7,3)],guide=F) +
  theme_bw() + hist_theme +
  labs(x='Cultural distance',y='Count',title='3')+
  scale_y_continuous(expand=c(0,0),limits=c(0,11000))

## Figure 2 of Nishiaki et al. (A4, Landscape)
## To plot Figure 2, uncomment the code below:
#grid.arrange(grobs=list(plt_hamming_histMP,plt_hamming_histMUP,
#                        plt_hamming_histUP,plt_hamming_hist),
#             nrow=2)

# Generate information for Table 1
paleo_count <- paleo_comp[c('SiteID','Era','Region','Country')]
paleo_count$one <- 1
aggregate(paleo_count$one,by=list(paleo_count$Region),FUN='sum')
paleo_MP <- paleo_count[paleo_count$Era=='MP',]
paleo_MUP <- paleo_count[paleo_count$Era=='MUP',]
paleo_UP <- paleo_count[paleo_count$Era=='UP',]
paleo_unknown <- paleo_count[paleo_count$Era=='unknouwn',]

paleo_table1 <- aggregate(paleo_MP$one,by=list(paleo_MP$Region),drop=F,FUN='sum')
names(paleo_table1) <- c('Region','period1')
paleo_table1$period2 <- 0
paleo_table1$period3 <- 0
paleo_table1$unknown <- 0

paleo_MUPagg <- aggregate(paleo_MUP$one,by=list(paleo_MUP$Region),drop=T,FUN='sum')
paleo_table1[paleo_table1$Region %in% paleo_MUPagg$Group.1,'period2'] <- paleo_MUPagg$x
paleo_UPagg <- aggregate(paleo_UP$one,by=list(paleo_UP$Region),drop=F,FUN='sum')
paleo_table1[paleo_table1$Region %in% paleo_UPagg$Group.1,'period3'] <- paleo_UPagg$x
paleo_unknownagg <- aggregate(paleo_unknown$one,by=list(paleo_unknown$Region),drop=F,FUN='sum')
paleo_table1[paleo_table1$Region %in% paleo_unknownagg$Group.1,'unknown'] <- paleo_unknownagg$x
paleo_table1$total <- aggregate(paleo_count$one,by=list(paleo_comp$Region),drop=F,FUN='sum')$x

paleo_count_site <- unique(paleo_count[c('Region','SiteID','one')])
paleo_table1$sites <- aggregate(paleo_count_site$one,by=list(paleo_count_site$Region),FUN='sum')$x
paleo_table1 <- paleo_table1[c(1,7,2:6)]
## To print information for Table 1, uncomment the following line
#print(paleo_table1)



