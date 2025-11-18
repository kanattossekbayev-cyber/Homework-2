############################################################
# 0. Setup: working directory, package, data loading
############################################################

setwd("C:/Users/User/Desktop/HWF")

# install.packages("survival")  # один раз, потом можно закомментировать
library(survival)

load("covid_wtc_2020_pub.Rdata")

df <- wtc_covid_pub

############################################################
# 1. Construct survival time and event variables
############################################################

start_date <- as.Date("2020-01-01")
end_date   <- as.Date("2020-08-31")

df$event_raw <- df$covid_status

df_surv <- subset(df, event_raw %in% c(0, 1))

df_surv$event <- ifelse(df_surv$event_raw == 1, 1, 0)

df_surv$time <- ifelse(
  df_surv$event == 1 & !is.na(df_surv$symptoms_date),
  as.numeric(df_surv$symptoms_date - start_date),
  as.numeric(end_date - start_date)
)

############################################################
# 2. Age categories (for KM) and continuous age (for Cox)
############################################################

df_surv$age_cat_f <- factor(
  df_surv$age_cat,
  levels = c(1, 2, 3, 4, 5),
  labels = c("≤39", "40–49", "50–59", "60–69", "70+")
)

df_surv$age_mid <- ifelse(df_surv$age_cat == 1, 35,
                          ifelse(df_surv$age_cat == 2, 45,
                                 ifelse(df_surv$age_cat == 3, 55,
                                        ifelse(df_surv$age_cat == 4, 65, 75))))

df_surv$sex_f  <- factor(df_surv$sex)
df_surv$race_f <- factor(df_surv$race)
df_surv$empl_f <- factor(df_surv$employment)

############################################################
# 3. Kaplan–Meier plot by age category -> SAVE TO PNG FILE
############################################################

fit_age <- survfit(Surv(time, event) ~ age_cat_f, data = df_surv)

# Сохраняем график в файл, чтобы не было ошибок margins
png("km_age_category.png", width = 1200, height = 900, res = 150)

par(mfrow = c(1, 1))
par(mar = c(5, 5, 4, 2) + 0.1)

plot(
  fit_age,
  col = 1:5,
  lty = 1,
  xlab = "Days since 1 Jan 2020",
  ylab = "COVID-free survival probability",
  main = "Time to COVID-19 by Age Category"
)

legend(
  "bottomleft",
  legend = levels(df_surv$age_cat_f),
  col = 1:5,
  lty = 1,
  bty = "n",
  title = "Age category"
)

dev.off()  # закрываем PNG, график записан в km_age_category.png

############################################################
# 4. Cox proportional hazards model (continuous age) + table
############################################################

cox_model <- coxph(
  Surv(time, event) ~ 
    age_mid + sex_f + race_f + empl_f + morb_obesity,
  data = df_surv,
  na.action = na.omit
)

cox_sum <- summary(cox_model)

cox_table <- data.frame(
  Variable = rownames(cox_sum$coefficients),
  HR       = round(cox_sum$coefficients[, "exp(coef)"], 2),
  Lower95  = round(cox_sum$conf.int[, "lower .95"], 2),
  Upper95  = round(cox_sum$conf.int[, "upper .95"], 2),
  p_value  = signif(cox_sum$coefficients[, "Pr(>|z|)"], 3)
)

print(cox_table)

############################################################
# End of script
############################################################
getwd()
list.files()
