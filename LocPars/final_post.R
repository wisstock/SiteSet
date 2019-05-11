# Copyright © Borys Olifirov 2019
# Institure of Molecular Biology and Genetics of NASU,
# Systems Biology Research Group

##### INIT #####
require(ggplot2)
require(RColorBrewer)

setwd("/home/astria/Bio/Ctools/SiteSet/LocPars")
position.df <- read.csv('semisite_df.csv')
position.df$location <- position.df$location/1000


gene.list <- levels(position.df$gene)

for (current.gene in gene.list) {
  nam <- paste(current.gene, 'df', sep = '.')
  inner.df <- subset(position.df, gene == current.gene)
  assign(nam, inner.df)
}
rm(inner.df)


##### PLOT #####


# GART plot
ggplot(GART.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  xlab('Позиція у інтроні (kbp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів')) +
  geom_vline(xintercept = 0.744)


# NAP1L plot
ggplot(NAP1L.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  xlab('Позиція у інтроні (kbp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів'))


# PCIF1 plot
ggplot(PCIF1.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  xlab('Позиція у інтроні (kbp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів'))


# CSTF3 plot
ggplot(CSTF3.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  xlab('Позиція у інтроні (kbp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів'))


# ZMYM3
ZMYM3.df$location <-ZMYM3.df$location*1000 
ggplot(ZMYM3.df,
       aes(x = location,
           y = factor(factor))) +
  stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
  scale_fill_gradient(low='blue', high='red') +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono') +
  xlab('Позиція у інтроні (kbp)') + 
  ylab('Фактор') +
  guides(fill = guide_legend(title='Щільність сайтів'))
