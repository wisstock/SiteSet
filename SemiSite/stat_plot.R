# Copyright © 2019 Institure of Molecular Biology and Genetics of NASU,
# Systems Biology Research Group
# Copyright © 2019 Borys Olifirov

require(ggplot2)
require(RColorBrewer)
require(HDInterval)


setwd("/home/astria/Bio/Ctools/SiteSet/SemiSite")
position.df <- read.csv('semisite_df.csv')


gene.list <- levels(position.df$gene)

position.stat <- data.frame(gene = character(),
                            comp = character(),
                            ks.p = integer(),
                            crit.p = integer())

for (current.gene in gene.list) {
  nam <- paste(current.gene, 'df', sep = '.')
  inner.df <- subset(position.df, gene == current.gene)
  
  # CDF plot
  gart.srsf3 <- inner.df$location[inner.df$factor == 'SRSF3']
  gart.srsf10 <- inner.df$location[inner.df$factor == 'SRSF10']
  gart.ythdc1 <- inner.df$location[inner.df$factor == 'YTHDC1']
  gart.uni <- unique(sort(inner.df$location))
  
  get_df.ecdf <- function(x, group, level = 0.05) { 
    n <- length(x)
    x.sort <- sort(x)
    y <- (1:n)/n 
    # CI по теореме Дворецкого-Кифера-Вольфовица (ДКВ)
    #  epsilon = sqrt(log(2/level)/(2*n))
    #  L = pmax(y - epsilon, 0)
    #  U = pmin(y + epsilon, 1)
    #  D <- approx.ksD(n)
    #  U3 <- pmin(y + D, 1)
    #  L3 <- pmax(y - D, 0)
    # CI на основе центральной предельной теоремы (ЦПТ) confidence band
    z <- qnorm(1-level/2)
    U = pmin(y + z*sqrt(y*(1-y)/n ),1)
    L = pmax(y - z*sqrt(y*(1-y)/n ),0)
    data.frame(x=x.sort, y, group, z, U, L) 
  }
  
  df.srsf3 <- get_df.ecdf(gart.srsf3,'SRSF3')
  df.all <- rbind(df.srsf3, get_df.ecdf(gart.srsf10, 'SRSF10'))
  df.all <- rbind(df.all, get_df.ecdf(gart.ythdc1, 'YTHDC1'))
  
  ggplot(df.all, aes(x=x, y=y, colour=group)) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin = L, ymax = U, fill = group), alpha = .3) +
    theme_minimal(base_size = 12,
                  base_family = 'ubuntu mono') +
    labs(title = 'GART') +
    xlab('Позиція у інтроні (bp)') + 
    ylab('Накопичена імовірність') +
    guides(fill = guide_legend(title='Фактор'),
           colour = guide_legend(title='Фактор'))
  
  
  
  # GART density bar plot
  h <- hdi(density(inner.df$location[inner.df$factor == 'YTHDC1']),
           credMass = .25,
           allowSplit = TRUE)  # calculate higest dens interval for a 0.25 prob
  
  ggplot(inner.df,
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
  
  # GART stat test
  ks.3.10 <- ks.test(inner.df$location[inner.df$factor == 'SRSF3'],
                     inner.df$location[inner.df$factor == 'SRSF10'],
                     alternative = 'two.sided')
  ks.y.10 <- ks.test(inner.df$location[inner.df$factor == 'YTHDC1'],
                     inner.df$location[inner.df$factor == 'SRSF10'],
                     alternative = 'two.sided')
  ks.y.3 <- ks.test(inner.df$location[inner.df$factor == 'YTHDC1'],
                    inner.df$location[inner.df$factor == 'SRSF3'],
                    alternative = 'two.sided')
  ks.y.3.l <- ks.test(inner.df$location[inner.df$factor == 'YTHDC1'],
                      inner.df$location[inner.df$factor == 'SRSF3'],
                      alternative = 'less')
  ks.y.3.g <- ks.test(inner.df$location[inner.df$factor == 'YTHDC1'],
                      inner.df$location[inner.df$factor == 'SRSF3'],
                      alternative = 'greater')
  
  critical.p <- p.adjust(c(ks.3.10$p.value,
                           ks.y.10$p.value,
                           ks.y.3.l$p.value),
                         "bonferroni")[1]
  
}
