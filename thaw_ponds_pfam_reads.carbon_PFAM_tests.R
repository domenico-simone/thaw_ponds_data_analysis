carbon_PFAM_oxic_old <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "old") %>%
  subset(layer != "hypo") %>%
  dplyr::select(total_mapped_reads_percent)

carbon_PFAM_anoxic_old <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "old") %>%
  subset(layer == "hypo") %>%
  dplyr::select(total_mapped_reads_percent)

t_test_old_oxic_vs_anoxic <- t.test(carbon_PFAM_oxic_old, carbon_PFAM_anoxic_old)

t_test_old_oxic_vs_anoxic

# Welch Two Sample t-test
# 
# data:  carbon_PFAM_oxic_old and carbon_PFAM_anoxic_old
# t = -12.824, df = 12.376, p-value = 1.634e-08
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.03344556 -0.02375891
# sample estimates:
#   mean of x  mean of y 
# 0.06608594 0.09468817 

aov_carbon_PFAM_oxic_old <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "old") %>%
  subset(layer != "hypo") %>%
  dplyr::select(age, total_mapped_reads_percent)

aov_carbon_PFAM_oxic_medium <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "medium") %>%
  subset(layer != "hypo") %>%
  dplyr::select(age, total_mapped_reads_percent)

aov_carbon_PFAM_oxic_emerge <- data.frame(sample_name=row.names(all.final.matrix), total_mapped_reads=rowSums(all.final.matrix)) %>%
  left_join(metadata.dataset_size) %>%
  mutate("total_mapped_reads_percent"=total_mapped_reads*100/total_reads) %>%
  left_join(metadata, by = ("sample_name")) %>%
  subset(age == "emerge") %>%
  subset(layer != "hypo") %>%
  dplyr::select(age, total_mapped_reads_percent)

aov.dataframe <- rbind(aov_carbon_PFAM_oxic_old, aov_carbon_PFAM_oxic_medium, aov_carbon_PFAM_oxic_emerge)
aov.dataframe$age <- base::as.factor(aov.dataframe$age)

ggboxplot(aov.dataframe, x = "age", y = "total_mapped_reads_percent", color = "age")

# ANOVA test with assumption of equal variances
summary(aov(total_mapped_reads_percent ~ age, data = aov.dataframe))
# ANOVA test with no assumption of equal variances
oneway.test(total_mapped_reads_percent ~ age, data = aov.dataframe)

plot(aov(total_mapped_reads_percent ~ age, data = aov.dataframe))


#########################################################
