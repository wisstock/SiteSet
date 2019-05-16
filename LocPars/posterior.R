
###### INIT #####

require(ggplot2)
require(moments)
require(wesanderson)

setwd("/home/astria/Bio/Ctools/SiteSet/LocPars")
z.loc <- read.csv('semisite_df.csv')
# rel.loc <- read.csv('RelLoc.csv', header = T)
# z.loc <- read.csv('ZLoc.csv', header = T)
# summ <- read.csv('semisite_df.csv', header = T)

# gene.list <- c('NAP1L', 'ZMYM3', 'PCIF1', 'CSTF3', 'GART')


##### GART #####

gart.z <- subset(z.loc, gene == 'GART',
                  select = c('b_factor',
                             'location'))

gart.z <- subset(gart.z, location<=1325+350 & location>=1325-350,
                 select = c('b_factor',
                            'location'))


tapply(gart.z$dep_factor_position,
                  as.factor(gart.z$dep_factor),
                  skewness, na.rm = T)
tapply(gart.z$z_score,
                  as.factor(gart.z$indep_site_position),
                  kurtosis, na.rm = T)


ggplot(gart.z, aes(x = location,
                   fill = b_factor,
                   color = b_factor)) +
  geom_dotplot(alpha = .65) +
  geom_vline(xintercept = 1324, 
             linetype = "dashed", size = 1,
             colour = '#5B544E') +
  geom_vline(xintercept = 760, 
             linetype = "dashed", size = 1,
             colour = '#E32636') +
  scale_x_continuous(limits = c(0, 2486)) +
  labs(x = '11th GART intron', y = ' ') +
  guides(fill=guide_legend(title='Factor'),
         color=guide_legend(title='Factor')) +
  theme_minimal(base_size = 20,
                base_family = 'ubuntu mono') +
  scale_fill_manual(values= wes_palette(n=3, name='Royal1')) +
  scale_color_manual(values=wes_palette(n=3, name='Royal1')) +
  guides(fill = F, color = F)


##### NAP1L ####

napl.z <- subset(z.loc, gene == 'NAP1L',
                 select = c('b_factor',
                            'location'))

napl.z <- subset(napl.z, location<=84+350 & location>=84-350,
                 select = c('b_factor',
                            'location'))

ggplot(napl.z, aes(x = location,
                   fill = b_factor,
                   color = b_factor)) +
  geom_dotplot(alpha = .65) +
  geom_vline(xintercept = 84, 
             linetype = "dashed", size = 1,
             colour = '#5B544E') +
  geom_vline(xintercept = 32, 
             linetype = "dashed", size = 1,
             colour = '#E32636') +
  scale_x_continuous(limits = c(0, 674)) +
  guides(fill = F, color = F) +
  labs(x = '13th NAP1L intron', y = ' ') +
  guides(fill=guide_legend(title='Factor'),
         color=guide_legend(title='Factor')) +
  theme_minimal(base_size = 20,
                base_family = 'ubuntu mono') +
  scale_fill_manual(values=wes_palette(n=3, name='Royal1')) +
  scale_color_manual(values=wes_palette(n=3, name='Royal1')) +
  guides(fill = F, color = F)

#### ZMYM3 ####

zmym.z <- subset(z.loc, gene == 'ZMYM3',
                 select = c('b_factor',
                            'location'))

zmym.z <- subset(zmym.z, location<=251+350 & location>=251-350,
                 select = c('b_factor',
                            'location'))

ggplot(zmym.z, aes(x = location,
                   fill = b_factor,
                   color = b_factor)) +
  geom_dotplot(alpha = .65) +
  geom_vline(xintercept = 251, 
             linetype = "dashed", size = 1,
             colour = '#5B544E') +
  geom_vline(xintercept = 62, 
             linetype = "dashed", size = 1,
             colour = '#E32636') +
  scale_x_continuous(limits = c(0, 674)) +
  guides(fill = F, color = F) +
  labs(x = '7th ZMYM3 intron', y = ' ') +
  guides(fill=guide_legend(title='Factor'),
         color=guide_legend(title='Factor')) +
  theme_minimal(base_size = 20,
                base_family = 'ubuntu mono') +
  scale_fill_manual(values=wes_palette(n=3, name='Royal1')) +
  scale_color_manual(values=wes_palette(n=3, name='Royal1')) +
  guides(fill = F, color = F)

##### CSTF3 #####

cstf.z <- subset(z.loc, gene == 'CSTF3',
                 select = c('b_factor',
                            'location'))

cstf.z <- subset(cstf.z, location<=12607+500 & location>=12607-500,
                 select = c('b_factor',
                            'location'))

ggplot(cstf.z, aes(x = location,
                   fill = b_factor,
                   color = b_factor)) +
  geom_dotplot(alpha = .65) +
  geom_vline(xintercept = 12607, 
             linetype = "dashed", size = 1,
             colour = '#5B544E') +
  geom_vline(xintercept = 2301, 
             linetype = "dashed", size = 1,
             colour = '#E32636') +
  scale_x_continuous(limits = c(0, 33248)) +
  guides(fill = F, color = F) +
  labs(x = '3th CSTF3 intron', y = ' ') +
  guides(fill=guide_legend(title='Factor'),
         color=guide_legend(title='Factor')) +
  theme_minimal(base_size = 20,
                base_family = 'ubuntu mono') +
  scale_fill_manual(values=wes_palette(n=3, name='Royal1')) +
  scale_color_manual(values=wes_palette(n=3, name='Royal1')) +
  guides(fill = F, color = F)

##### PCIF1 #####

pcif.z <- subset(z.loc, gene == 'CSTF3',
                 select = c('b_factor',
                            'location'))

pcif.z <- subset(pcif.z, location<=1018+350 & location>=1018-350,
                 select = c('b_factor',
                            'location'))

ggplot(pcif.z, aes(x = location,
                   fill = b_factor,
                   color = b_factor)) +
  geom_dotplot(alpha = .65) +
  geom_vline(xintercept = 1018, 
             linetype = "dashed", size = 1,
             colour = '#5B544E') +
  geom_vline(xintercept = 267, 
             linetype = "dashed", size = 1,
             colour = '#E32636') +
  scale_x_continuous(limits = c(0, 1396)) +
  guides(fill = F, color = F) +
  labs(x = '2th PCIF1 intron', y = ' ') +
  guides(fill=guide_legend(title='Factor'),
         color=guide_legend(title='Factor')) +
  theme_minimal(base_size = 20,
                base_family = 'ubuntu mono')  +
  scale_fill_manual(values=wes_palette(n=3, name='Royal1')) +
  scale_color_manual(values=wes_palette(n=3, name='Royal1')) +
  guides(fill = F, color = F)




