library(vegan)
library(pairwiseAdonis)

##Effect of Sex
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000-sex")

map_sex <- read.table('mapping_dm_sex.txt', header=T)
map_sex_no <- read.table('mapping_dm_sex_no.txt', header=T)

#beta diversity
uw_dm_sex <- as.dist(read.table('uw-distance-matrix-sex.tsv', header=T))
w_dm_sex <- as.dist(read.table('w-distance-matrix-sex.tsv', header=T))
w_dm_sex_no <- as.dist(read.table('w-distance-matrix-sex-no.txt', header=T))

adonis2(uw_dm_sex~Sex, data=map_sex, permutations=5000)
adonis2(w_dm_sex~Sex, data=map_sex, permutations=5000)
adonis2(w_dm_sex_no~Sex, data=map_sex_no, permutations=5000)

adonis2(uw_dm_sex~Age*Sex, data=map_sex, permutations=5000)
adonis2(w_dm_sex~Age*Sex, data=map_sex, permutations=5000)
pairwise.adonis2(w_dm_sex~Age_sex, data=map_sex)

##Effect of Age
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000")

map <- read.table('mapping_dm.txt', header=T)
map_no <- read.table('mapping_dm_no_out.txt', header=T)

uw_dm <- as.dist(read.table('uw-distance-matrix.tsv', header=T))
w_dm <- as.dist(read.table('w-distance-matrix.tsv', header=T))
w_dm_no <- as.dist(read.table('w-distance-matrix-no-out.txt', header=T))

#beta diversity
adonis2(uw_dm~Age, data=map, permutations=5000)
adonis2(w_dm~Age, data=map, permutations=5000)
adonis2(w_dm_no~Age, data=map_no, permutations=5000)

#alpha diversity
library(car)
library(ggplot2)

alpha<-read.table('alpha-summary.txt', header=T)
alpha_sex<-read.table('alpha-summary_sex.txt', header=T)

shan_age<-lm(shannon_entropy~Age, data=alpha)
Anova(shan_age)

faith_age<-lm(faith_pd~Age, data=alpha)
Anova(faith_age)

obs_age<-lm(observed_features~Age, data=alpha)
Anova(obs_age)

shan_sex<-lm(shannon_entropy~Sex, data=alpha_sex)
Anova(shan_sex)

#ANCOMBC
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis")

asv = read.table("feature-table-for-ancom.txt", header=T, check.names=FALSE)
metadata = read.table("mapping_for_ancom.txt", header=T)
taxonomy = read.table("taxonomy-abc.txt", header=T)

asv_matrix = asv %>% column_to_rownames("sampleid") %>% as.matrix()
tax_matrix = taxonomy %>% column_to_rownames("feature") %>% as.matrix()
meta = metadata %>% column_to_rownames("sampleid")

ASV<-otu_table(asv_matrix, taxa_are_rows = TRUE)
TAX<-tax_table(tax_matrix)
samples<-sample_data(meta)

asv_phylo = phyloseq(ASV, TAX, samples)


age_asv = ancombc2(data=asv_phylo, fix_formula="Age",
                      p_adj_method = "fdr",
                      group = "Age", global = T)

res_age_asv<-age_asv$res
res_global_age_asv<-age_asv$res_global

write_csv(res_age_asv, "Diff_abund_age_asv.csv")
write_csv(res_global_age_asv, "Diff_abund_age_global.csv")

##Effect of Diet and Month, Season

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000")

map <- read.table('mapping_dm.txt', header=T)
map_no <- read.table('mapping_dm_no_out.txt', header=T)

uw_dm <- as.dist(read.table('uw-distance-matrix.tsv', header=T))
w_dm <- as.dist(read.table('w-distance-matrix.tsv', header=T))
w_dm_no <- as.dist(read.table('w-distance-matrix-no-out.txt', header=T))

#beta diversity
adonis2(uw_dm~Age+Diet*Season+Group, data=map, permutations=5000)
adonis2(w_dm~Age+Diet*Season+Group, data=map, permutations=5000)
adonis2(w_dm_no~Age+Diet*Season+Group, data=map_no, permutations=5000)


#alpha diversity
alpha<-read.table('alpha-summary.txt', header=T)

shan_diet<-lm(shannon_entropy~Diet*Season, data=alpha)
Anova(shan_diet)

faith_diet<-lm(faith_pd~Diet*Season, data=alpha)
Anova(faith_diet)

obs_diet<-lm(observed_features~Diet*Season, data=alpha)
Anova(obs_diet)

shan_group<-lm(shannon_entropy~Diet+Group, data=alpha)
Anova(shan_group)

faith_group<-lm(faith_pd~Diet+Group, data=alpha)
Anova(faith_group)

obs_group<-lm(observed_features~Diet+Group, data=alpha)
Anova(obs_group)

#ANCOMBC for diet
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis")

asv = read.table("feature-table-for-ancom.txt", header=T, check.names=FALSE)
metadata = read.table("mapping_for_ancom.txt", header=T)
taxonomy = read.table("taxonomy-abc.txt", header=T)

asv_matrix = asv %>% column_to_rownames("sampleid") %>% as.matrix()
tax_matrix = taxonomy %>% column_to_rownames("feature") %>% as.matrix()
meta = metadata %>% column_to_rownames("sampleid")

ASV<-otu_table(asv_matrix, taxa_are_rows = TRUE)
TAX<-tax_table(tax_matrix)
samples<-sample_data(meta)

asv_phylo = phyloseq(ASV, TAX, samples)


diet_asv = ancombc2(data=asv_phylo, fix_formula="Diet",
                   p_adj_method = "fdr",
                   group = "Diet", global = T)

res_diet_asv<-diet_asv$res
res_global_diet_asv<-diet_asv$res_global

write_csv(res_diet_asv, "Diff_abund_diet_asv.csv")
write_csv(res_global_diet_asv, "Diff_abund_diet_global.csv")

#ANCOMBC for season
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis")

asv = read.table("feature-table-for-ancom.txt", header=T, check.names=FALSE)
metadata = read.table("mapping_for_ancom.txt", header=T)
taxonomy = read.table("taxonomy-abc.txt", header=T)

asv_matrix = asv %>% column_to_rownames("sampleid") %>% as.matrix()
tax_matrix = taxonomy %>% column_to_rownames("feature") %>% as.matrix()
meta = metadata %>% column_to_rownames("sampleid")

ASV<-otu_table(asv_matrix, taxa_are_rows = TRUE)
TAX<-tax_table(tax_matrix)
samples<-sample_data(meta)

asv_phylo = phyloseq(ASV, TAX, samples)


season_asv = ancombc2(data=asv_phylo, fix_formula="Season",
                    p_adj_method = "fdr",
                    group = "Season", global = T)

res_season_asv<-season_asv$res
res_global_season_asv<-season_asv$res_global

write_csv(res_season_asv, "Diff_abund_season_asv.csv")
write_csv(res_global_season_asv, "Diff_abund_season_global.csv")

##Effect of Season and Group for Wild

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000-wild")

map_wild <- read.table('dm_mapping_wild.txt', header=T)

uw_dm_wild <- as.dist(read.table('uw-distance-matrix-wild.tsv', header=T))
w_dm_wild <- as.dist(read.table('w-distance-matrix-wild.tsv', header=T))

#beta diversity
adonis2(uw_dm_wild~Age+Season+Group, data=map_wild, permutations=5000)
adonis2(w_dm_wild~Age+Season+Group, data=map_wild, permutations=5000)

#alpha diversity
alpha_wild<-read.table('shannon-diversity.tsv', header=T)

shan_season<-lm(shannon_entropy~Season, data=alpha_wild)
Anova(shan_season)

shan_group<-lm(shannon_entropy~Group, data=alpha_wild)
Anova(shan_group)

##Effect of Season and Group for Supplemented

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000-supp")

map_supp <- read.table('dm_mapping_supp.txt', header=T)
map_supp_no <- read.table('dm_mapping_supp_no.txt', header=T)

uw_dm_supp <- as.dist(read.table('uw-distance-matrix-supp.tsv', header=T))
w_dm_supp <- as.dist(read.table('w-distance-matrix-supp.tsv', header=T))
w_dm_supp_no <- as.dist(read.table('w-distance-matrix-supp-no.txt', header=T))

#beta diversity
adonis2(uw_dm_supp~Age+Season+Group, data=map_supp, permutations=5000)
adonis2(w_dm_supp~Age+Season+Group, data=map_supp, permutations=5000)
adonis2(w_dm_supp_no~Age+Season+Group, data=map_supp_no, permutations=5000)

#alpha diversity
alpha_supp<-read.table('shannon-diversity.tsv', header=T)

shan_season<-lm(shannon_entropy~Season, data=alpha_supp)
Anova(shan_season)

shan_group<-lm(shannon_entropy~Group, data=alpha_supp)
Anova(shan_group)

##PLOTS
#alpha

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000")

#alpha diversity
alpha<-read.table('alpha-summary.txt', header=T)

ggplot(alpha, aes(x=Diet, y=shannon_entropy)) + geom_boxplot() +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
ggplot(alpha, aes(x=Season, y=shannon_entropy)) + geom_boxplot(aes(fill=Diet)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
ggplot(alpha, aes(x=Age, y=shannon_entropy)) + geom_boxplot()
ggplot(alpha, aes(x=Group, y=shannon_entropy)) + geom_boxplot()

ggplot(alpha, aes(x=Season, y=shannon_entropy)) + geom_boxplot(aes(fill=Diet)) +
  +   theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#wild
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000-wild")

alpha<-read.table('shannon-diversity.tsv', header=T)
ggplot(alpha, aes(x=Group, y=shannon_entropy)) + geom_boxplot() +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#supp
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000-supp")

alpha<-read.table('shannon-diversity.tsv', header=T)
ggplot(alpha, aes(x=Group, y=shannon_entropy)) + geom_boxplot() +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#beta
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Bicca_Marques_howler/analysis/core-metrics-results-8000")

map <- read.table('mapping_dm.txt', header=T)
map_no <- read.table('mapping_dm_no_out.txt', header=T)

uw_dm <- as.dist(read.table('uw-distance-matrix.tsv', header=T))
w_dm_no <- as.dist(read.table('w-distance-matrix-no-out.txt', header=T))

u_mds<-metaMDS(uw_dm, trymax=100)
u_mds.points<-merge(x = u_mds$points, y = map, by.x = "row.names", by.y = "sampleid")
unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Diet, shape=Season), size=4)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
unmds

unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(fill=Season), color="black", pch=21, size=5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
unmds

unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(fill=Group), color="black", pch=21, size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black')))
unmds

unmds <- ggplot(u_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(fill=Age), color="black", pch=21, size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
unmds


w_mds<-metaMDS(w_dm_no, trymax=100)
w_mds.points<-merge(x = w_mds$points, y = map_no, by.x = "row.names", by.y = "sampleid")
wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(color=Diet, shape=Season), size=4)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds

wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(fill=Age), color='black', pch=21, size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds

wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(fill=Month), color='black', pch=21, size=5)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds

wnmds <- ggplot(w_mds.points, aes(x = MDS1, y = MDS2)) +  
  geom_point(aes(fill=Group), color='black', pch=21, size=5)+
  scale_fill_brewer(palette="Accent")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
wnmds

##Taxa summary for diet
asv_phylo %>%
  group_by(Diet) %>%
  dplyr::summarize(Mean)



