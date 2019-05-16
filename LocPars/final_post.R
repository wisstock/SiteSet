# Copyright © 2019 Institure of Molecular Biology and Genetics of NASU,
# Systems Biology Research Group
# Copyright © 2019 Borys Olifirov


##### INIT #####
require(ggplot2)
require(RColorBrewer)
require(HDInterval)
# require(fitdistrplus)

# Mode <- function(x) {  # mode calculation function
#   ux <- unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }

setwd("/home/astria/Bio/Ctools/SiteSet/LocPars")
position.df <- read.csv('semisite_df.csv')


gene.list <- levels(position.df$gene)

position.stat <- data.frame(gene = character(),
                            comp = character(),
                            ks.p = integer(),
                            crit.p = integer())

for (current.gene in gene.list) {
  nam <- paste(current.gene, 'df', sep = '.')
  inner.df <- subset(position.df, gene == current.gene)
  assign(nam, inner.df)
}
rm(inner.df)



##### GART #####
# gart.srsf3.list <- GART.df$location[GART.df$factor == 'SRSF3']  # CDF plot for
# gart.srsf3.list <- sort(gart.srsf3.list)                        # SRSF3 sites
# n <- length(gart.srsf3.list)                                    # locstion distr
# plot(ecdf(gart.srsf3.list))                                     #
# lines(gart.srsf3.list, (1:n)/n, type = 's', col="blue")         #

# CDF plot for GART data
gart.srsf3 <- GART.df$location[GART.df$factor == 'SRSF3']
gart.srsf10 <- GART.df$location[GART.df$factor == 'SRSF10']
gart.ythdc1 <- GART.df$location[GART.df$factor == 'YTHDC1']
gart.uni <- unique(sort(GART.df$location))

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
h <- hdi(density(GART.df$location[GART.df$factor == 'YTHDC1']),
         credMass = .25,
         allowSplit = TRUE)  # calculate higest dens interval for a 0.25 prob

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

# GART stat test
ks.3.10 <- ks.test(GART.df$location[GART.df$factor == 'SRSF3'],
                   GART.df$location[GART.df$factor == 'SRSF10'],
                   alternative = 'two.sided')
ks.y.10 <- ks.test(GART.df$location[GART.df$factor == 'YTHDC1'],
                   GART.df$location[GART.df$factor == 'SRSF10'],
                   alternative = 'two.sided')
ks.y.3 <- ks.test(GART.df$location[GART.df$factor == 'YTHDC1'],
                   GART.df$location[GART.df$factor == 'SRSF3'],
                   alternative = 'two.sided')
ks.y.3.l <- ks.test(GART.df$location[GART.df$factor == 'YTHDC1'],
                  GART.df$location[GART.df$factor == 'SRSF3'],
                  alternative = 'less')
ks.y.3.g <- ks.test(GART.df$location[GART.df$factor == 'YTHDC1'],
                  GART.df$location[GART.df$factor == 'SRSF3'],
                  alternative = 'greater')
 
critical.p <- p.adjust(c(ks.3.10$p.value,
                         ks.y.10$p.value,
                         ks.y.3.l$p.value),
                       "bonferroni")[1]





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
