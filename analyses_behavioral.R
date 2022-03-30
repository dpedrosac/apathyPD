# This is code to run analyses on PET data of PD-patients 
# Code developed by David Pedrosa

# Version 1.0 # 2022-03-20, added separate file for meta analyses

# ==================================================================================================
## Specify packages of interest and load them automatically if needed
packages = c("readxl", "tableone", "ggplot2", "tidyverse", "rstatix", "ggpubr",
				"MatchIt") 											# packages needed

## Load or install necessary packages
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# ==================================================================================================
## In case of multiple people working on one project, this helps to create an automatic script
username = Sys.info()["login"]
if (username == "dpedr") {
wdir = "D:/PET-MCI"
} else if (username == "david") {
wdir = "D:/PET-MCI"
}
setwd(wdir)

# ==================================================================================================
## Read and prepare data, so that further analyses are possible
df 				<- read_excel(file.path(wdir, 'data', 'raw_data.xlsx'))
colnames(df) 	<- c(	"record_id", "diagnosis", "age", "gender", "educational_years", "employed", 
						"disease_duration", "mmst", "bdi", paste0("aes", 1:18), "aes_total", 
						"aes_cognitive", "aes_emotional", "aes_behavior", "aes_other", 
						"aes_active_initiation", "aes_passive_engagement", "aes_insight", 
						"aes_social", "aes_curious", "FDG_PET", "aes_cogn_perc", 
						"aes_emot_perc", "aes_behav_perc", "aes_other_perc", "aes_active_perc", 
						"aes_passive_perc", "aes_insight_perc", "aes_social_perc", "aes_curious_perc", 
						"cog_dom_exec_z", "cog_dom_memo_z", "cog_dom_att_z", "cog_dom_lang_z", 
						"cog_dom_visu_z", "objective_cog_z", "panda_total") # rename to unequivocal names
						
df <- df %>% mutate(diagnosis = dplyr::case_when(diagnosis==1 ~"PD without MCI", # change variables' coding
												 diagnosis==2 ~ "PD with MCI", 
												 diagnosis==3 ~"Healthy controls")) %>%
			 mutate(gender = dplyr::case_when(gender==1 ~"male", # change variables' coding
												 gender==2 ~ "female")) %>%
			 mutate(employed = dplyr::case_when(employed==0 ~"no", # change variables' coding
												 employed==1 ~ "yes")) %>%
			 mutate_at(vars(diagnosis, gender, employed), factor) # convert variables to factors
		
# ==================================================================================================
## Create TableOne to visualise results for all three groups of participants

NumVars 				<- c("age", "educational_years", "disease_duration", "mmst", "bdi", 
								"aes_total", "panda_total")
catVars 				<- c("gender", "employed")
myVars 					<- c(catVars, NumVars)
labels_TableOne 		<- c(	"Gender", "Employed", "Age", "Educational years", "Disease duration", 
								"Mini-Mental Status Examination (MMST)", "Apathy evaluation score (AES)", 
								"Parkinson Neuropsychometrie Dementia Assessment (PANDA)")

tableOne <- CreateTableOne(vars = myVars, data = df, factorVars = catVars, strata="diagnosis")
tableOne <- print(tableOne, nonnormal=c("gender", "employed", "diagnosis"))
write.csv(tableOne, file=file.path(wdir, "results", "TableOne_PET-MCI.csv"))

# ==================================================================================================
## Statistiscal analyses for AES-, MMST- and PANDA-total (and sub-) scores

test_name = c("aes_total", "panda_total", "mmst")
plot_list <- list() 
for (i in 1:3){ # loop through screening tests and total AES scores
	df_cognitive <- df %>% select(test_name[i], diagnosis) %>%  
						reorder_levels(diagnosis, order = c("PD without MCI", "PD with MCI", "Healthy controls"))
	colnames(df_cognitive) <- c("measure", "diagnosis")
	res.aov <- df_cognitive %>% anova_test(measure ~ diagnosis)
	pwc <- df_cognitive %>% tukey_hsd(measure ~ diagnosis)

	pwc <- pwc %>% add_xy_position(x = "diagnosis")
	plot_list[[i]] <- ggboxplot(df_cognitive, x = "diagnosis", y = "measure") +
	  stat_pvalue_manual(pwc, hide.ns = TRUE) +
	  labs(
		subtitle = get_test_label(res.aov, detailed = TRUE),
		caption = get_pwc_label(pwc)
		) + 
		theme_minimal()
	#TODO: a) Better labeling is needed for x and y axes! b) subplots should be smaller
	}
grid.arrange(grobs=plot_list,ncol=2)

# Plot group differences in AEs subscores:
df_subscores <- df %>% select(ends_with("perc"), diagnosis) %>%
	  pivot_longer(!diagnosis, names_to = "aes", names_prefix="aes_", values_to = "count")

stat.test <- df_subscores %>% group_by(aes) %>%
  t_test(count ~ diagnosis, p.adjust.method = "bonferroni") %>%
  #adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_y_position(step.increase = 0.1) %>%
  add_x_position(x = "aes", dodge = 1)
  
ggbarplot(df_subscores, x = "aes", y = "count", add = "mean_sd", fill="diagnosis", 
			position = position_dodge(.84)) +
	ylim(-.5, 1) +
	theme_minimal() + 
	stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0, hide.ns = TRUE)
#TODO: a) Change the x position for the sugnificance statement b) change coordinate system
	
# ==================================================================================================
## Analyze subscores in AES scores
	
# Look for differences in subtypes
# 1 Chi-Square test for both groups and stacked bar chart
### UNFINISHED!!!
# Create a summary of available results
df_subgroups1 		<- df %>% select(record_id, diagnosis, aes_active_initiation, aes_passive_engagement, aes_insight, aes_social, aes_curious) %>%  
						mutate_at(vars(diagnosis, record_id), factor) 
summary_subgroup1 	<- df_subgroups1 %>%
								pivot_longer(!record_id, 
										names_to = "aes_subscore", 
										values_to = "count")


summary_care 	<- df_careCovid %>% filter(timing=="before") %>% drop_na() %>% group_by(ratings)  %>% summarize(before=sum(ratings))
temp_data 		<- df_careCovid %>% filter(timing=="during") %>% drop_na()%>% group_by(ratings)  %>% summarize(during=sum(ratings))
summary_care 		<- merge(summary_care, temp_data, by="ratings", all = T)
summary_care[is.na(summary_care)]=0.001

colnames(summary_care) 	<- c("ratings", "before", "during")
summary_care$before 	<- summary_care$before/sum(summary_care$before)
summary_care$during 	<- summary_care$during/sum(summary_care$during)
summary_care 			<- dplyr::mutate(summary_care, ID = row_number()) %>% 
								pivot_longer(cols = c("before", "during"), 
										names_to = "timing", 
										values_to = "percentage")

summary_care		<- summary_care %>% filter(ratings>0)
symmetry_test(ratings~as.factor(timing) | ID, data=df_careCovid, distribution="exact")
stat.test 		<- df_careCovid  %>%
  sign_test(ratings ~ timing) %>%
  add_significance() %>%
  add_xy_position(x = "ratings", dodge = 0.8)
stat.test

# Start plotting results for changes in perceived care during pandemic
bar_width 	<- .7 
y_start 	<- summary_care %>% filter(timing=="before") %>% summarize(y_start=cumsum(percentage))
y_end 		<- summary_care %>% filter(timing=="during") %>% summarize(y_end=cumsum(percentage))
df_lines 	<- data.frame(y_start=y_start, y_end=y_end, x_start=rep(1+bar_width/1.99,4), x_end=rep(2-bar_width/1.99,4))

summary_care$timing <- as.factor(summary_care$timing)
levels(summary_care$timing) <- c("before", "during") # explicit definition of order necessary as otherwise data is sorted alphabetically (cf. https://stackoverflow.com/questions/5208679/order-bars-in-ggplot2-bar-graph)


# 2. Run analyses with matched groups (cf. propensity score matching via {MatchIt})	
df_matching <- df %>% filter(diagnosis!="Healthy controls") %>% arrange(diagnosis) %>% 
				droplevels() # creates matched groups
match.it <- matchit(diagnosis ~ age + gender, data = df_matching, 
						method="optimal", 
						ratio=1)
matched_dataMCI = match.data(match.it, group="all") # matched data for PD-patients w/ MCI
tableOneMatched <- CreateTableOne(vars = myVars, data = matched_dataMCI, factorVars = catVars, strata="diagnosis")
print(tableOneMatched)

#TODO: Create an analyses with subgroups for the matched subjects!