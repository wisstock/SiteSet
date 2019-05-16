# Copyright © 2019 Institure of Molecular Biology and Genetics of NASU,
# Systems Biology Research Group
# Copyright © 2019 Borys Olifirov


##### INIT #####
require(ggplot2)
require(RColorBrewer)
require(HDInterval)
require(fitdistrplus)

# Mode <- function(x) {  # mode calculation function
#  ux <- unique(x)
#  ux[which.max(tabulate(match(x, ux)))]
#}

setwd("/home/astria/Bio/Ctools/SiteSet/LocPars")
position.df <- read.csv('semisite_df.csv')


gene.list <- levels(position.df$gene)

for (current.gene in gene.list) {
  nam <- paste(current.gene, 'df', sep = '.')
  inner.df <- subset(position.df, gene == current.gene)
  assign(nam, inner.df)
}
rm(inner.df)



##### GART #####
h <- hdi(density(GART.df$location[GART.df$factor == 'YTHDC1']),
         credMass = .25,
         allowSplit = TRUE)  # расчет интервала плотности вероятности для 0.25


ggplot(GART.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +  # low='blue', high='red'
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  labs(title = 'GART') +
  xlab('Позиція у інтроні (bp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність імовірності \n розподілу сайтів')) +
  geom_segment(aes(x = 744 , y = .5,
                   xend = 744, yend = 3.5),
               size = .5,
               colour = 'grey') +
  annotate('rect', xmin = h[[1]],
                   xmax = h[[2]],
           ymin = .5, ymax = 3.5, alpha = .3,
           fill = 'white')


##### NAP1L #####
h <- hdi(density(GART.df$location[GART.df$factor == 'YTHDC1']),
         credMass = .25,
         allowSplit = TRUE)
ggplot(NAP1L.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  labs(title = 'NAP1L') +
  xlab('Позиція у інтроні (bp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів')) +
  geom_vline(xintercept = 384, colour = '#5B588E', size = 1) +
  annotate('rect', xmin = 384-250, xmax = 384+250,
           ymin = .5, ymax = 3.5, alpha = .2,
           fill = 'yellow')


# PCIF1 plot
ggplot(PCIF1.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  labs(title = 'PCIF1') +
  xlab('Позиція у інтроні (bp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів')) +
  geom_vline(xintercept = 270, colour = '#5B588E', size = 1) +
  annotate('rect', xmin = 270-250, xmax = 270+250,
           ymin = .5, ymax = 3.5, alpha = .2,
           fill = 'yellow')


# CSTF3 plot
ggplot(CSTF3.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  labs(title = 'CSTF3') +
  xlab('Позиція у інтроні (bp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів')) +
  geom_vline(xintercept = 2304, colour = '#5B588E', size = 1) +
  annotate('rect', xmin = 2304-250, xmax = 2304+250,
           ymin = .5, ymax = 3.5, alpha = .2,
           fill = 'yellow')


# ZMYM3
ggplot(ZMYM3.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  labs(title = 'ZMYM3') +
  xlab('Позиція у інтроні (bp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів')) +
  geom_vline(xintercept = 65, colour = '#5B588E')
