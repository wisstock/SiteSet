# Copyright © 2019 Institure of Molecular Biology and Genetics of NASU,
# Systems Biology Research Group
# Copyright © 2019 Borys Olifirov

require(ggplot2)
require(RColorBrewer)
require(HDInterval)


setwd("/home/astria/Bio/Ctools/SiteSet/SemiSite")


position.df <- read.csv('semisite_df.csv')

gene.list <- levels(position.df$gene)
factor.list <- levels(position.df$factor)
factor.comb <- combn(factor.list, 2)

position.stat <- data.frame(gene = character(),
                            comb = character(),
                            ks.hyp = character(),
                            ks.p = integer(),
                            crit.p = integer(),
                            stringsAsFactors=FALSE)

apa.pos <- data.frame(gene = c('GART', 'NAP1L', 'ZMYM3', 'CSTF3', 'PCIF1'),
                      position = c(744, 384, 65, 2304, 279))


ks_stat <- function(current.gene, input.factor, input.data) {
  # KS stat calculation
  
  alt.hyp <- c('two.sided', 'less', 'greater')
  inner.df <- subset(input.data, gene == current.gene)
  
  ks_loop <- function(two.factors, input.df, hyp.list) {
    target.p <- 0
    for (current.hyp in hyp.list) {  # loop over the alt hyp variance
      
      current.ks <-  ks.test(input.df$location[input.df$factor == two.factors[1]],
                             input.df$location[input.df$factor == two.factors[2]],
                             alternative = current.hyp)
      
      if(target.p < current.ks$p.value) {
        target.p <- current.ks$p.value
        target.hyp <- current.hyp
        next
      }
      return(c(paste(two.factors, collapse='_'), target.hyp, target.p))
    }
  }

  
  ks.res <- apply(input.factor, MARGIN = 2,
                  FUN = ks_loop,
                  input.df = inner.df, hyp.list = alt.hyp)
  
  current.df <- data.frame(current.gene, t(ks.res))
  colnames(current.df) <- c('gene', 'comb', 'ks.hyp', 'ks.p')
  
  return(current.df)
}

ecdf_plot <- function(input.gene, input.data, input.stat) {
  
  inner.df <- subset(input.data, gene == input.gene)
  inner.stat <- subset(input.stat, gene == input.gene)
  
  # subsetting by factors
  inner.srsf3 <- inner.df$location[inner.df$factor == 'SRSF3']
  inner.srsf10 <- inner.df$location[inner.df$factor == 'SRSF10']
  inner.ythdc1 <- inner.df$location[inner.df$factor == 'YTHDC1']
  inner.uni <- unique(sort(inner.df$location))
  
  
  critical.p <- p.adjust(c(inner.stat$ks.p[inner.stat$comb == 'SRSF10_SRSF3'],
                           inner.stat$ks.p[inner.stat$comb == 'SRSF10_YTHDC1'],
                           inner.stat$ks.p[inner.stat$comb == 'SRSF3_YTHDC1']),
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
  rm(df.srsf3)
  
  # ECDF plot
  ecdf.nam <- paste(input.gene, 'ecdf', sep = '.')  # generating plot name in air
  
  ecdf <- ggplot(df.all, aes(x=x, y=y, colour=group)) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin = L, ymax = U, fill = group), alpha = .2) +
    theme_minimal(base_size = 12,
                  base_family = 'ubuntu mono') +
    labs(title = input.gene) +
    xlab('Позиція у інтроні (bp)') + 
    ylab('Кумулятивна імовірність') +
    guides(fill = guide_legend(title='Фактор'),
           colour = guide_legend(title='Фактор'))
  
  return(assign(ecdf.nam, ecdf))
}

# dens_plot <- function() {
  # calculate higest density interval for a 0.25 prob
  h <- hdi(density(inner.df$location[inner.df$factor == 'YTHDC1']),
           credMass = .25,
           allowSplit = TRUE)
  
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


stat.res <- lapply(gene.list, ks_stat,
                   input.factor = factor.comb,
                   input.data = position.df)
position.stat <- as.data.frame(do.call(rbind, stat.res))

rm(factor.comb)
rm(stat.res)

lapply(gene.list, ecdf_plot,
       input.data = position.df,
       input.stat = position.stat)
