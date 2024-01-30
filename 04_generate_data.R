##################################################################
###
### run the function to generate data
###
##################################################################

source("03_fun_generate_data.R")

dat <- ageperiod_surv_foi_sim_data(
            beta0_survival_sus = beta0_survival_sus,
            beta0_survival_inf = beta0_survival_inf,
            foi_age_effect = foi_age_effect,
            age_effect = age_effect_true,
            period_effect = period_effect_true,
            nT_age = nT_age,
            nT_period = nT_period,
            processnum = processnum
            )

n_fit <- nrow(dat$d_fit)
n_fit_sus_cens_posttest <- nrow(dat$d_fit_sus_cens_posttest)
n_fit_sus_cens_postno <- nrow(dat$d_fit_sus_cens_postno)
n_fit_sus_mort_posttest <- nrow(dat$d_fit_sus_mort_posttest)
n_fit_sus_mort_postno <- nrow(dat$d_fit_sus_mort_postno)
n_fit_icap_cens <- nrow(dat$d_fit_icap_cens)
n_fit_icap_mort <- nrow(dat$d_fit_icap_mort)
n_fit_rec_neg_cens_posttest <- nrow(dat$d_fit_rec_neg_cens_posttest)
n_fit_rec_neg_cens_postno <- nrow(dat$d_fit_rec_neg_cens_postno)
n_fit_rec_neg_mort <- nrow(dat$d_fit_rec_neg_mort)
n_fit_rec_pos_cens <- nrow(dat$d_fit_rec_pos_cens)
n_fit_rec_pos_mort <- nrow(dat$d_fit_rec_pos_mort)
n_fit_idead <- nrow(dat$d_fit_idead)

### test that these data combined are the same
### dimension as the overall generated data
# n_fit_sus_cens_posttest+
# n_fit_sus_cens_postno+
# n_fit_sus_mort_posttest+
# n_fit_sus_mort_postno+
# n_fit_icap_cens+
# n_fit_icap_mort+
# n_fit_rec_neg_cens_posttest+
# n_fit_rec_neg_cens_postno+
# n_fit_rec_neg_mort+
# n_fit_rec_pos_cens+
# n_fit_rec_pos_mort+
# n_fit_idead

# datatypes_out <- data.frame(datatype = c("overall","sus_cens_posttest",
# "sus_cens_postno",
# "sus_mort_posttest",
# "sus_mort_postno",
# "icap_cens",
# "icap_mort",
# "rec_neg_cens_posttest",
# "rec_neg_cens_postno",
# "rec_neg_mort",
# "rec_pos_cens",
# "rec_pos_mort",
# "idead"),
# number_samples = c(n_fit,
# n_fit_sus_cens_posttest,
# n_fit_sus_cens_postno,
# n_fit_sus_mort_posttest,
# n_fit_sus_mort_postno,
# n_fit_icap_cens,
# n_fit_icap_mort,
# n_fit_rec_neg_cens_posttest,
# n_fit_rec_neg_cens_postno,
# n_fit_rec_neg_mort,
# n_fit_rec_pos_cens,
# n_fit_rec_pos_mort,
# n_fit_idead))

# datatypes_out
# write.csv(datatypes_out, file = paste0("datatypes_out_",processnum,".csv"))


##############################################################
### there are too many deer infected at to capture
### because we aren't accounting for 
### disease-associated mortality in the generating function
##############################################################
# icap_eval_df <- data.frame(left_age = dat$left_age[dat$pos1==TRUE])
# icap_eval_plot <- ggplot(data = icap_eval_df) +
#       geom_histogram(aes(x = left_age),bins = 50) +
#       theme_bw() +
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )


# ggsave(paste0("figures/icap_eval_left_age.png"),
#       icap_eval_plot,
#       height = 6,
#       width = 8)

# icap_eval_df$weights <- icap_eval_df$left_age/sum(icap_eval_df$left_age)
# icap_eval_df$weights2 <- (1 - 1/icap_eval_df$left_age)/sum((1 - 1/icap_eval_df$left_age))

# icap_weights_plot <- ggplot(data = icap_eval_df) +
#       geom_point(aes(x = left_age,y = weights), size = 1.3, alpha = .5) +
#       theme_bw() +
#       xlab("Age at Entry") + 
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )

# icap_weights2_plot <- ggplot(data = icap_eval_df) +
#       geom_point(aes(x = left_age,y = weights2), size = 1.3, alpha = .5) +
#       theme_bw() +
#       xlab("Age at Entry") + 
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )

# ggsave(paste0("figures/icap_eval_weights.png"),
#       icap_weights_plot,
#       height = 6,
#       width = 8)



# ggsave(paste0("figures/icap_eval_weights2.png"),
#       icap_weights_plot,
#       height = 6,
#       width = 8)


# n_remove <- sum(dat$pos1) - round(sum(dat$pos1)/5)
# indx_remove <- sample(which(dat$pos1==TRUE),n_remove,replace=FALSE) 
# n_ind <- dat$n_ind - n_remove
# n_ind <- dat$n_ind


############################
# evaluating fawn data
############################
d_fit <- dat$d_fit
names(d_fit)
d_surv_fawn <- d_fit[d_fit$e_age < 52,]
nrow(d_surv_fawn)
#745 total number of captures as fawn
nrow(d_fit[d_fit$e_age < 16,])
#343 captured as neonates
nrow(d_surv_fawn) - nrow(d_fit[d_fit$e_age < 16,])
#310 captured less than 1 year of age, but not as neonates

# how many fawns are cwd+ at capture?
d_surv_fawn$cwd_cap
sum(d_surv_fawn$cwd_cap)
d_surv_fawn$s_age[d_surv_fawn$cwd_cap == 1]
names(d_fit)
d_surv_fawn$data_type[d_surv_fawn$cwd_cap == 1]


# how many fawns are cwd+ at mortality?
d_surv_fawn <- d_surv_fawn[!is.na(d_surv_fawn$s_age),]
nrow(d_surv_fawn)#573 died, 
745-573
#172 right censored

d_surv_fawn <- d_surv_fawn[!is.na(d_surv_fawn$cwd_mort),]
nrow(d_surv_fawn)
#201 fawns at capture tested for cwd at mortality
#number of fawns that were tested at mort that actually died as fawns
length(d_surv_fawn$cwd_mort[d_surv_fawn$s_age < 52])

d_surv_fawn$cwd_mort[d_surv_fawn$s_age < 52]
#0 fawns that died as fawns and were tested for cwd, then tested +

d_fit_fawn <- d_fit[!is.na(d_fit$inf_age),]
sum(d_fit_fawn$inf_age < 52)

d_fit_fawn$data_type[d_fit_fawn$inf_age < 52]



############################
# evaluating yearling data
############################
d_fit <- dat$d_fit
names(d_fit)
d_surv_yrl <- d_fit[d_fit$e_age >= 52 & d_fit$e_age < 104,]

#there were 104 individuals captured as yearlings
nrow(d_surv_yrl)

#there were 104 individuals captured as yearlings
# write.csv(data.frame(table(d_surv_yrl$data_type)),file = "results/data_types_yearlings.csv")


##########################################
# evaluating deer that got infected data
##########################################
d_inf <- d_fit[!is.na(d_fit$inf_age),]
nrow(d_inf)
table(d_inf$data_type)
# write.csv(data.frame(table(d_inf$data_type)),file = "results/data_types_inf.csv")
age_inf <- data.frame(table(d_inf$inf_age))
names(age_inf) <- c("age_infection","num")
head(age_inf)
age_inf$age_infection <- as.numeric(as.character(age_inf$age_infection))
# age_inf_plot <- ggplot(data = age_inf,aes(x = age_infection,y = num)) +
#     geom_bar(stat = "identity", fill = "gray") +
#     ylab("number observed")+
#     xlab("true age of infection")+
#     geom_vline(xintercept = seq(52,962,by = 52), linetype = "dotted")+
#      scale_x_continuous(breaks = seq(0, max(age_inf$age_infection), by = 52),
#                      labels = seq(0, max(age_inf$age_infection), by = 52))+
#     theme_bw()
# age_inf_plot
# ggsave("figures/age_inf_plot.png",age_inf_plot, height = 6, width = 10)    



# Assuming your data frame is named age_inf
age_inf_annual <- age_inf %>%
  mutate(age_year_interval = cut(age_infection, breaks = seq(0, max(age_infection) + 52, by = 52), labels = FALSE)) %>%
  group_by(age_year_interval) %>%
  summarise(sum_num = sum(num))

# age_inf_annual
# age_inf_annual_plot <- ggplot(data = age_inf_annual,aes(x = age_year_interval,y = sum_num)) +
#     xlab("Age(Year) of Infection")+
#     ylab("Number Observed")+
#     geom_bar(stat = "identity", fill = "gray") +
#     theme_bw()
# age_inf_annual_plot
# ggsave("figures/age_inf_annual_plot.png",age_inf_annual_plot, height = 6, width = 10)    
