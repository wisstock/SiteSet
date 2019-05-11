##### INIT SECTION ####

require(ggplot2)
require(gridExtra)
require(cluster)

setwd("/home/astria/Bio/Ctools/SiteSet/BSite")
position.df <- read.csv('bsite_df.csv')

##### ASSHOLE SECTION #####

gene.name <- 'ZMYM3'

gart.df <- subset(position.df, gene == gene.name & b_factor == 'SRSF3' |
                    gene == gene.name & b_factor == 'YTHDC1' |
                    gene == gene.name & b_factor == 'SRSF10' |
                    gene == gene.name & b_factor == 'PAS',
                  select = c(b_factor, location, score, site))

tapply(gart.df$b_factor, 
       gart.df$b_factor,
       length)

apply(gart.df[2], 2, length(gart.df$b_factor))


          
##### VISUALISATION SECTION #####

gart.df <- subset(position.df, gene == gene.name & b_factor == 'SRSF3' |
                    gene == gene.name & b_factor == 'YTHDC1' |
                    gene == gene.name & b_factor == 'SRSF10' |
                    gene == gene.name & b_factor == 'PAS',
                  select = c(b_factor, location, score, site))


nk <- length(gart.df$score[gart.df$b_factor == 'YTHDC1'])

gart.k <- kmeans(gart.df$location, nk)
gart.df$clusters <- factor(gart.k$cluster)

ggplot(gart.df, aes(x = location, 
                    colour = b_factor, 
                    fill = b_factor)) +
  geom_dotplot(alpha = 0.8) +
  geom_dotplot(data = subset(gart.df,
                             b_factor == 'YTHDC1',
                             select = c('location', 'clusters')),
               fill = '#5B544E',
               colour = '#5B544E') +
  scale_x_continuous(breaks = seq(0, max(gart.df$location), 25),
                     limits = c(0, max(gart.df$location))) +
  labs(fill = 'Factor', 
       colour = 'Factor',
       x = 'Location, bp',
       y = 'Count',
       title = gene.name)
