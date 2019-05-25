# Copyright © 2019 Institure of Molecular Biology and Genetics of NASU,
# Systems Biology Research Group
# Copyright © 2019 Borys Olifirov

require(HDInterval)


setwd("/home/astria/Bio/Ctools/SiteSet/SemiSite")


position.df <- read.csv('semisite_df.csv')

gene.list <- levels(position.df$gene)
factor.list <- levels(position.df$factor)

hdi.stat <- data.frame(lower = character(),
                       upper = character(),
                       factor = character(),
                       gene = character())


hdi.limit <- .2
i <- 1

hdi.loop <- function(input.factor, input.data, input.hdi) {
  h <-lapply(input.factor, function(x, data, mass){
                 hdi(density(data$location[data$factor == x]),
                 credMass = mass,
                 allowSplit = FALSE)
             },
             mass = input.hdi,
             data = input.data)

  h.df <- as.data.frame(do.call(rbind, h))
  h.df$factor <- as.vector(input.factor)
  return(h.df)
}

for(gene in gene.list) {
  current.gene <- gene
  current.df <- subset(position.df, gene == current.gene)
  out.hdi <- hdi.loop(input.factor = factor.list,
                      input.data = current.df,
                      input.hdi = hdi.limit)
  out.hdi$gene <- gene
  hdi.stat <- rbind(hdi.stat, out.hdi)
}

write.csv(hdi.stat, file = 'HDI_stat.csv')
