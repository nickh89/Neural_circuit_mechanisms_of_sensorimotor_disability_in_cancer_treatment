
####### kinematics errors ########
source("~/Dropbox/papers_dropbox/circuits_2021/Code/Final/SummarySE.R")
kinematics_errors <- read_excel("~/Dropbox (GaTech)/papers_dropbox/cancer_Research_2020/behavioral_study/database/kinematics_errors.xlsx")
# 
#   kinematics_errors %>% filter(miss_type_fore_hind == 2) %>%
#   ggboxplot(x = "miss_distance_second_cm", y = "miss_distance_first_cm",
#             palette = c("#707070", "#943bd3"),
#             add = "jitter",boxwex = 10,  fill = "miss_distance_second_cm",add.params = list(size = 0.5, jitter = 0.1))+
#     stat_compare_means(label.y = 4, method = "anova")+
#     labs(title="Hind limb error distance", y = "error (cm)", x = "First vs Second errors")
#   

fig6_f<-kinematics_errors %>% filter(miss_type_fore_hind == 2) %>%
  ggplot( aes(x = miss_distance_second_cm, y = (miss_distance_first_cm), fill=miss_distance_second_cm))+
  geom_boxplot(width = 0.2, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.05))+
  scale_fill_manual(values = c("#808080", "#943bd3" ))+
  theme_classic()+
  labs(title="Hind limb error distance", y = "error (cm)", x = "First vs Second errors")+
  stat_compare_means(label.y = 4, method = "anova")

kinematics_errors %>% filter(miss_type_fore_hind == 2) %>%
  summarySE(measurevar = "miss_distance_first_cm", groupvars = c("miss_distance_second_cm"),na.rm=TRUE)



####### replacement errors ########
replacement_errors <- read_excel("~/Dropbox (GaTech)/papers_dropbox/cancer_Research_2020/behavioral_study/database/replacement_errors.xlsx")

fig6_e<-replacement_errors %>% 
  ggplot( aes(x = treatment, y = (replace_perc_total_error_removed), fill=treatment))+
  geom_boxplot(width = 0.2, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.05))+
  scale_fill_manual(values = c("#808080", "#943bd3" ))+
  theme_classic()+
  labs(title="Hindlimb Replacement %", y = "Replacement (%)", x = "Treatment")+
  stat_compare_means(label.y = 1.01, method = "anova")+
  ylim(0,1)

summarySE(replacement_errors, measurevar = "replace_perc_total_error_removed", groupvars = c("treatment"),na.rm=TRUE)
rm(summarySE,kinematics_errors,replacement_errors)
