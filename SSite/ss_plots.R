##### INIT #####

require(ggplot2)
require(gridExtra)

setwd('/home/astria/Bio/Ctools/SiteSet/SSite')

in_df <- read.csv('ss_df.csv', header = T)
res_df <- read.csv('ss_df.csv', header = T)


##### SUBSETTING ######



in_df$pas_loc <- seq(from = 0, to = 0, len = 188)  # gl(1, 1, 188, 'common')

in_df$pas_loc[13] <- 1  # NAP1L
in_df$pas_loc[26] <- 1
in_df$pas_loc[33] <- 1  # ZMYM3
in_df$pas_loc[57] <- 1
in_df$pas_loc[76] <- 1  # PCIF1
in_df$pas_loc[92] <- 1
in_df$pas_loc[109] <- 1  # CSTF3
in_df$pas_loc[129] <- 1
in_df$pas_loc[157] <- 1  # GART
in_df$pas_loc[178] <- 1

in_df$pas_loc <- factor(in_df$pas_loc, levels = c(0, 1))
levels(in_df$pas_loc) <- c('Common', 'APA') 


f_ss <- subset(in_df, 
               end == 5, 
               select = c('score', 'int_num', 'gene', 'pas_loc'))
t_ss <- subset(res_df, 
               end == 3, 
               select = c('score', 'int_num', 'gene', 'pas_loc'))

out_f_ss <- subset(f_ss,                             # 5'ss df exclude APA int val 
                   pas_loc == 0, 
                   select = c('score', 'int_num', 'gene'))


##### STAT CALC #####

n <- c(length(res_df$score[res_df$pas_loc == 0 & res_df$end == 5]),
         length(res_df$score[res_df$pas_loc == 1 & res_df$end == 5]))
int_type <- c('const', 'pas')

f_mean <- c(mean(res_df$score[res_df$pas_loc == 0 & res_df$end == 5]),
            mean(res_df$score[res_df$pas_loc == 1 & res_df$end == 5]))
f_sd <- c(sd(res_df$score[res_df$pas_loc == 0 & res_df$end == 5]),
        sd(res_df$score[res_df$pas_loc == 1 & res_df$end == 5]))
f_se <- c(f_sd[1] / length(res_df$score[res_df$pas_loc == 0 & res_df$end == 5]),
          f_sd[2] / length(res_df$score[res_df$pas_loc == 1 & res_df$end == 5]))
f_stat <- data.frame(n, f_mean, f_sd, f_se, int_type)

t_mean <- c(mean(res_df$score[res_df$pas_loc == 0 & res_df$end == 3]),
            mean(res_df$score[res_df$pas_loc == 1 & res_df$end == 3]))
t_sd <- c(sd(res_df$score[res_df$pas_loc == 0 & res_df$end == 3]),
          sd(res_df$score[res_df$pas_loc == 1 & res_df$end == 3]))
t_se <- c(t_sd[1] / length(res_df$score[res_df$pas_loc == 0 & res_df$end == 3]),
          t_sd[2] / length(res_df$score[res_df$pas_loc == 1 & res_df$end == 3]))
t_stat <- data.frame(n, t_mean, t_sd, t_se, int_type)

wilcox.test(f_ss$score~f_ss$pas_loc)
wilcox.test(t_ss$score~t_ss$pas_loc)


# Shapiro Test

shapiro.test(f_ss$score)  # include PAS int
qqnorm(f_ss$score)


shapiro.test(out_f_ss$score)  # exclude PAS int
qqnorm(out_f_ss$score)


##### PLOTS #####



ggplot(f_ss, aes( x = as.factor(pas_loc),
                  y = score,
                  fill = as.factor(pas_loc))) +
  geom_boxplot(alpha = .35) +
  labs(y = 'Score',
       x = 'Intron type') +
  guides(fill = F, color = F) +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono')

ggplot(f_ss, aes(x = gene, 
                 y = score, 
                 fill = as.factor(gene))) +
  geom_boxplot(alpha = .35) +
  labs(y = 'Score',
       x = 'Gene') +
  guides(fill = F, color = F) +
  theme_minimal(base_size = 12,
                base_family = 'ubuntu mono')+
  scale_fill_brewer(palette="Dark2") +
  scale_color_brewer(palette="Dark2")



