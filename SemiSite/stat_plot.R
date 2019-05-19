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
                            crit.p = integer(),
                            stringsAsFactors=FALSE)

apa.pos <- data.frame(gene = c('GART', 'NAP1L', 'ZMYM3', 'CSTF3', 'PCIF1'),
                      position = c(744, 384, 65, 2304, 279))

test <- c()
# loop over all genes in input data frame
for (current.gene in gene.list) {
  inner.df <- subset(position.df, gene == current.gene)
  inner.apa <- apa.pos$position[apa.pos$gene == current.gene]
  test <- append(test, current.gene)
  
  # subsetting by factors
  inner.srsf3 <- inner.df$location[inner.df$factor == 'SRSF3']
  inner.srsf10 <- inner.df$location[inner.df$factor == 'SRSF10']
  inner.ythdc1 <- inner.df$location[inner.df$factor == 'YTHDC1']
  inner.uni <- unique(sort(inner.df$location))
  
  
  # KS stat test
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
  
  
  get_df.ecdf <- function(x, group, level = 0.05) { 
    
    n <- length(x)
    x.sort <- sort(x)
    y <- (1:n)/n 
    
    # confidence band calculated bu Central limit theorem
    z <- qnorm(1-level/2)
    U = pmin(y + z*sqrt(y*(1-y)/n ),1)
    L = pmax(y - z*sqrt(y*(1-y)/n ),0)
    data.frame(x=x.sort, y, group, z, U, L) 
  }
  
  df.srsf3 <- get_df.ecdf(inner.srsf3,'SRSF3')
  df.all <- rbind(df.srsf3, get_df.ecdf(inner.srsf10, 'SRSF10'))
  df.all <- rbind(df.all, get_df.ecdf(inner.ythdc1, 'YTHDC1'))
  
  
  # calculate higest density interval for a 0.25 prob
  h <- hdi(density(inner.df$location[inner.df$factor == 'YTHDC1']),
           credMass = .25,
           allowSplit = TRUE)
  
  
  # ECDF plot
  ecdf.nam <- paste(current.gene, 'ecdf', sep = '.')  # generating plot name in air
  
  ecdf <- ggplot(df.all, aes(x=x, y=y, colour=group)) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin = L, ymax = U, fill = group), alpha = .2) +
    theme_minimal(base_size = 12,
                  base_family = 'ubuntu mono') +
    labs(title = current.gene) +
    xlab('Позиція у інтроні (bp)') + 
    ylab('Кумулятивна імовірність') +
    guides(fill = guide_legend(title='Фактор'),
           colour = guide_legend(title='Фактор'))
  
  assign(ecdf.nam, ecdf)
  

  # density bar plot
  dens.nam <- paste(current.gene, 'dens', sep = '.')  # generating plot name in air

  dens <- ggplot(inner.df,
         aes(x = location,
             y = factor(factor))) +
    stat_density(aes(fill = stat(density)), geom = "raster", position = "identity") +
    scale_fill_gradient(low='blue', high='red') +  # low='blue', high='red'
    geom_segment(aes(x = 744,     # PAS location
                     xend = 744,
                     y = .5,
                     yend = 3.5),
                 size = .5,
                 colour = 'grey') +
    geom_text(aes(label = 'XYI',
                  x = 744,
                  y = 3.4),
              family = 'ubuntu mono')
    annotate('rect', xmin = h[[1]],              # interval bar
             xmax = h[[2]],
             ymin = .5, ymax = 3.5, alpha = .3,
             fill = 'white') +
    theme_minimal(base_size = 12,
                  base_family = 'ubuntu mono') +
    labs(title = current.gene) +
    xlab('Позиція у інтроні (bp)') + 
    ylab('Фактор') +
    guides(fill = guide_legend(title='Щільність імовірності \n розподілу сайтів'))
      
  
  assign(dens.nam, dens)
  
  
}

