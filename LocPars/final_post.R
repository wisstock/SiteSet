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

gart.count <- c()
gart.count <- append(length(GART.df$location[GART.df$factor == 'SRSF10']))
gart.count <- append(GART.df$location[GART.df$factor == 'SRSF3'])
gart.count <- append(length(GART.df$location[GART.df$factor == 'YTHDC1']))

cstf3.count <- c()
cstf3.count <- append(length(CSTF3.df$location[CSTF3.df$factor == 'SRSF3']),
                      length(CSTF3.df$location[CSTF3.df$factor == 'SRSF10']),
                      length(CSTF3.df$location[CSTF3.df$factor == 'YTHDC1']))

nap1l.count <- c()
nap1l.count <- append(length(NAP1L.df$location[NAP1L.df$factor == 'SRSF3']),
                      length(NAP1L.df$location[NAP1L.df$factor == 'SRSF10']),
                      length(NAP1L.df$location[NAP1L.df$factor == 'YTHDC1']))

pcif1.count <- c()
pcif1.count <- append(length(PCIF1.df$location[PCIF1.df$factor == 'SRSF3']),
                      length(PCIF1.df$location[PCIF1.df$factor == 'SRSF10']),
                      length(PCIF1.df$location[PCIF1.df$factor == 'YTHDC1']))

zmym3.count <- c()
zmym3.count <- append(length(ZMYM3.df$location[ZMYM3.df$factor == 'SRSF3']),
                      length(ZMYM3.df$location[ZMYM3.df$factor == 'SRSF10']),
                      length(ZMYM3.df$location[ZMYM3.df$factor == 'YTHDC1']))





