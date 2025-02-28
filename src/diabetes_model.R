# ------------------------------------------------------------
# Set Working Directory & Load Necessary Packages
# ------------------------------------------------------------
library('haven')
library('tidyverse')
library('tableone')
library('dplyr')
library('pROC')
library('stats')    
library('mfp')      
library('lspline')  
library('rms')
library('GGally')
library('tidyverse')
library('broom')
library('ResourceSelection')
library('car')

# Set working directory
#setwd("~/workingDirectory")

#--------------------------------------------
#Data Exploration
#--------------------------------------------

## Load dataframe
data <- read_dta("dataset.dta")

## Calculate size of dataset
n_distinct(data$patient_id)

## Check missing values by column
colSums(is.na(data))

## Check for 0 values by column
colSums(data == 0)

# Proportion of with positive diabetes diagnosis
table(data$diabetes)
prop.table(table(data$diabetes))

# Unique rows in dataset
subset_data <- data[, 2:ncol(data)]
num_unique_rows <- nrow(unique(subset_data))
print(num_unique_rows)

#--------------------------------------------
# Continuous Variables
#--------------------------------------------

vars_cont <- c("glucose", "blood_pressure", "skin_thickness", "insulin", "bmi", "diabetes_genetic_score", "age")

## Check distribution of continuous variable - report mean and SD or med and IQR?
for (col in vars_cont) {
  
  # Histogram
  hist(data[[col]], main = paste("Histogram of", col), 
       xlab = col, col = "blue", border = "black")
  
  # Q-Q Plot
  qqnorm(data[[col]], main = paste("Normal Q-Q Plot of", col))
  qqline(data[[col]], col = "red")
  
  # Boxplot
  boxplot(data[[col]], main = paste("Boxplot of", col), 
          ylab = col, col = "lightgray")
}

## Tabulate mean and SD of continuous variables (to be presented if normally distributed)
tab_cont <- CreateTableOne(vars = vars_cont, data = data, strata = 'diabetes', test = TRUE)
print(tab_cont, nonnormal = vars_cont)

# Get min and max values
data_no_diabetes <- data %>% filter(diabetes == 0)
summary(data_no_diabetes)

data_diabetes <- data %>% filter(diabetes == 1)
summary(data_diabetes)

#--------------------------------------------
# Categorical Variables
#--------------------------------------------

# Create a factorized pregnancies column with a 5+ category
data <- data %>%
  mutate(pregnancies_factor = factor(pregnancies))

vars_cat <- c("pregnancies_factor")

## Tabulate counts and percentages of categorical variables
tab_cat <- CreateTableOne(vars = vars_cat, data = data, strata = 'diabetes', test = TRUE, factorVars = vars_cat)
print(tab_cat, showAllLevels = TRUE)

# Get min and max values
summary(data_no_diabetes)
summary(data_diabetes)

#--------------------------------------------
#Investigating missing data
#--------------------------------------------

#---- Addressing Insulin Missing Values------

# Calculate the percentage of 0 values
zero_percentage <- (sum(data$insulin == 0, na.rm = TRUE) / nrow(data)) * 100

# Display the result
cat("Percentage of data points with a 0 value for insulin:", zero_percentage, "%\n")

# Replace 0 values in the 'insulin' column with NA
data$insulin <- ifelse(data$insulin == 0, NA, data$insulin)

#--------Missingness Mechanism---------

#Create a new column in which TRUE if insulin is missing, and FALSE if not
data$miss_ins <- is.na(data$insulin)

#Check count and % of missing values
table(data$miss_ins)
prop.table(table(data$miss_ins))

miss_model <- glm(miss_ins ~ pregnancies + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age + factor(diabetes),data = data, family=binomial())
summary(miss_model)

tidy(miss_model, exponentiate = TRUE, conf.int = TRUE)

#----- Removal of Insulin---------

#Create model with insulin removed
no_ins_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age, data = data, 
                    family = binomial)
summary(no_ins_model)

#Compare AIC
AIC(no_ins_model)

# Generate ROC curves
roc2 <- roc(data$diabetes, predict(no_ins_model, data = data, type = "response"))

# Calculate AUC values
auc2 <- auc(roc2)
print(auc2)

#------- Complete Records Analysis---------

# Create full model
complete_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + insulin + bmi + diabetes_genetic_score + age, data = data, 
                      family = binomial)
summary(complete_model)

complete_data <- data[complete.cases(data), ]  # Keep rows without NA

#Generate ROC curve
roc1 <- roc(complete_data$diabetes, 
            predict(complete_model, newdata = complete_data, type = "response"))

# Calculate AUC value
auc1 <- auc(roc1)
print(auc1)

#OR and CIs
tidy(complete_model, exponentiate = TRUE, conf.int = TRUE)

#--------Median Imputation-------

# Calculate the median of non-NA insulin values
insulin_median <- median(data$insulin, na.rm = TRUE)
data$insulin_median_imputed <- ifelse(is.na(data$insulin), insulin_median, data$insulin)

#Create model 
med_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + insulin_median_imputed + bmi + diabetes_genetic_score + age, data = data, 
                 family = binomial)

#Compare AIC
AIC(med_model)

# Generate ROC curves
roc2 <- roc(data$diabetes, predict(med_model, data = data, type = "response"))

# Calculate AUC values
auc2 <- auc(roc2)
print(auc2)

#OR and CIs
tidy(med_model, exponentiate = TRUE, conf.int = TRUE)

#--------Mean Imputation-------

# Calculate the median of non-zero insulin values
insulin_mean <- mean(data$insulin, na.rm = TRUE)
data$insulin_mean_imputed <- ifelse(is.na(data$insulin), insulin_mean, data$insulin)

#Create model
mean_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + insulin_mean_imputed + bmi + diabetes_genetic_score + age, data = data, 
                  family = binomial)

#Compare AIC
AIC(mean_model)

# Generate ROC curves
roc2 <- roc(data$diabetes, predict(mean_model, data = data, type = "response"))

# Calculate AUC values
auc2 <- auc(roc2)
print(auc2)

#OR and CIs
tidy(mean_model, exponentiate = TRUE, conf.int = TRUE)

#-------Mean Reression Imputation---------

#Linear regression to predict insulin
ins_predict <- lm(insulin ~ factor(pregnancies) + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age, data=data)
summary(ins_predict)

data$ins_regressimp<-data$insulin #Create new column copying insulin data
data$ins_hat<-predict(ins_predict, newdata = data) # predict insulin
data$ins_regressimp<- ifelse(is.na(data$ins_regressimp),
                             data$ins_hat,
                             data$ins_regressimp)

#Create model 
mean_reg_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + ins_regressimp + bmi + diabetes_genetic_score + age, data = data, 
                      family = binomial)

#Compare AIC
AIC(mean_reg_model)

# Generate ROC curves
roc2 <- roc(data$diabetes, predict(mean_reg_model, data = data, type = "response"))

# Calculate AUC values
auc2 <- auc(roc2)
print(auc2)

#OR and CIs
tidy(mean_reg_model, exponentiate = TRUE, conf.int = TRUE)

#-------Mean Regression Imputation (incl outcome)---------

#Linear regression to predict insulin
ins_predict_dia <- lm(insulin ~ factor(pregnancies) + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age + diabetes, data=data)
summary(ins_predict_dia)

data$dia_ins_regressimp<-data$insulin #Create new column copying insulin data
data$dia_ins_hat<-predict(ins_predict_dia, newdata = data) # predict insulin
data$dia_ins_regressimp<- ifelse(is.na(data$dia_ins_regressimp),
                                 data$dia_ins_hat,
                                 data$dia_ins_regressimp)

#Create model 
dia_mean_reg_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + dia_ins_regressimp + bmi + diabetes_genetic_score + age, data = data, 
                          family = binomial)

#Compare AIC
AIC(dia_mean_reg_model)

# Generate ROC curves
roc2 <- roc(data$diabetes, predict(dia_mean_reg_model, data = data, type = "response"))

# Calculate AUC values
auc2 <- auc(roc2)
print(auc2)

#OR and CIs
tidy(dia_mean_reg_model, exponentiate = TRUE, conf.int = TRUE)

#----- Create Clean Dataset---------

# Specify the columns to keep
columns <- c("patient_id", "pregnancies", "glucose", "blood_pressure", 
             "skin_thickness", "bmi", "diabetes_genetic_score", "age", "diabetes")

clean_data <- data[columns]

#--------------------------------------------
#Assessing Logistic Regression Assumptions
#--------------------------------------------

set.seed(123)  # Set seed for reproducibility

clean_data$S <- 0  # Create a new column 'S' with default value 0
test_indices <- sample(1:nrow(clean_data), size = 0.5 * nrow(clean_data)) # Randomly select 50% of the rows for testing
clean_data$S[test_indices] <- 1 # Set 'S' to 1 for the test data 

# Verify the distribution
cat("Training set size:", sum(clean_data$S == 0), "\n")
cat("Testing set size:", sum(clean_data$S == 1), "\n")

# Filter the data for S == 0
train_data <- clean_data %>% filter(S == 0)

#--------------------------------------------
# Absence of multicollinearity
#--------------------------------------------

# Create correlation matrix
cor(train_data[, c("pregnancies", "glucose", "blood_pressure", "skin_thickness", "diabetes_genetic_score", "age", "bmi")])

#--------------------------------------------
# Assess Linearity in Logit for Continuous Variables
#--------------------------------------------

# Filter the data for S == 0
train_data <- clean_data %>% filter(S == 0)

# Fit the logistic regression model
model_00 <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age, data = train_data, 
                family = binomial)

# Predict the probability (p) of diabetes positivity
probabilities <- predict(model_00, type = "response")

# Create dataframe of continuous variables
lin_check <- train_data %>%
  dplyr::select(pregnancies, glucose, blood_pressure, skin_thickness, bmi, diabetes_genetic_score, age)

# Bind the logit and tidying the data for plot
lin_check <- lin_check %>%
  mutate(logit = log(probabilities/(1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)

#Create plots
ggplot(lin_check, aes(predictor.value, logit))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() + 
  facet_wrap(~predictors, scales = "free_x")

#--------------------------------------------
# Variable Transformations
#--------------------------------------------

#--------------------------------------------
# Create Basic Model
#-------------------------------------------

### Basic Model
basic_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                   family = binomial)
summary(basic_model)

AIC(basic_model)

#-------------------------------------------
# Treating Pregnancy variable
#-------------------------------------------

fac_preg_model <- glm(diabetes ~ as_factor(pregnancies) + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                      family = binomial)
summary(fac_preg_model)

#-------- Compare AIC--------------
AIC(basic_model, fac_preg_model)

#-------- Compare ROC/AUC 
#------- Likelihood ratio test-------
lrtest(fac_preg_model, basic_model)

#-------- Create factors --------------

# Create a factorized pregnancies column with a 5+ category
clean_data <- clean_data %>%
  mutate(pregnancies_5_factor = factor(ifelse(pregnancies >= 5, "5+", pregnancies)))
# Ensure the levels are ordered correctly
clean_data$pregnancies_5_factor <- factor(clean_data$pregnancies_5_factor, 
                                          levels = c("0", "1", "2", "3", "4", "5+"))
# Check the distribution of the factorized column
table(clean_data$pregnancies_5_factor)

# Create a factorized pregnancies column with a 4+ category
clean_data <- clean_data %>%
  mutate(pregnancies_4_factor = factor(ifelse(pregnancies >= 4, "4+", pregnancies)))
# Ensure the levels are ordered correctly
clean_data$pregnancies_4_factor <- factor(clean_data$pregnancies_4_factor, 
                                          levels = c("0", "1", "2", "3", "4+"))
# Check the distribution of the factorized column
table(clean_data$pregnancies_4_factor)

#-------- Create models--------------

#Create model with 5+ category
fac_5_preg_model <- glm(diabetes ~pregnancies_5_factor+ glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                        family = binomial)
summary(fac_5_preg_model)

#Create model with 4+ category
fac_4_preg_model <- glm(diabetes ~pregnancies_4_factor+ glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                        family = binomial)
summary(fac_4_preg_model)

#-------- Compare AIC--------------
AIC(basic_model, fac_5_preg_model, fac_4_preg_model)

#-------- Compare ROC/AUC --------------

# Filter the data for S == 0
train_data <- clean_data %>% filter(S == 0)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(basic_model, newdata = train_data, type = "response"))
roc2 <- roc(train_data$diabetes, predict(fac_5_preg_model, newdata = train_data, type = "response"))
roc3 <- roc(train_data$diabetes, predict(fac_4_preg_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
auc3 <- auc(roc3)
print(auc1)
print(auc2)
print(auc3)

#-------- Clean data accordingly --------------
clean_data <- clean_data %>%
  dplyr::mutate(pregnancies_fac = pregnancies_4_factor) %>%
  dplyr::select(-pregnancies_5_factor, -pregnancies_4_factor)

#-------------------------------------------
# Treating Age variable 
#-------------------------------------------

#---------Log Transformation----------

log_age_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + log(age), data = clean_data %>% filter(S == 0), 
                     family = binomial)
summary(log_age_model)

#-------- Compare AIC--------------
AIC(basic_model, log_age_model)

#-------- Compare ROC/AUC --------------
train_data <- clean_data %>% filter(S == 0)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(basic_model, newdata = train_data, type = "response"))
roc2 <- roc(train_data$diabetes, predict(log_age_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
print(auc1)
print(auc2)

#-------- Assess new logit plot --------------

#Plot Log odds vs Age
train_data$age_logit <- predict(basic_model, newdata = train_data, type = "link")
ggplot(train_data, aes(x = age, y = age_logit)) +
  geom_point(alpha = 0.5, color = "black") + 
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Predicted Log-Odds vs. Age", 
       x = "Transformed Variable", 
       y = "Predicted Log-Odds") +
  theme_minimal()

#Plot Log odds vs Log Age
train_data$log_age_logit <- predict(log_age_model, newdata = train_data, type = "link")
ggplot(train_data, aes(x = log(age), y = log_age_logit)) +
  geom_point(alpha = 0.5, color = "black") + 
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Predicted Log-Odds vs. Log Age", 
       x = "Log Age", 
       y = "Predicted Log-Odds") +
  theme_minimal()

#---------Age group categorisation----------

# Create age categories
clean_data <- clean_data %>%
  mutate(age_group = case_when(
    age >= 18 & age <30 ~ 1,
    age >= 30 & age < 45 ~ 2, 
    age >= 45 & age < 60 ~ 3, 
    age >= 60 ~ 4))

clean_data$age_group <- factor(clean_data$age_group, 
                               levels = c(1, 2, 3, 4), 
                               labels = c("18-29", "30-44", "44-59", "60+"))

table(clean_data$age_group)

group_age_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + bmi + diabetes_genetic_score + as_factor(age_group), data = clean_data %>% filter(S == 0), 
                       family = binomial)
summary(group_age_model)

#-------- Compare AIC--------------
AIC(basic_model, group_age_model)

## Log decreases AIC significantly

#-------- Compare ROC/AUC --------------

train_data <- clean_data %>% filter(S == 0)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(basic_model, newdata = train_data, type = "response"))
roc2 <- roc(train_data$diabetes, predict(group_age_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
print(auc1)
print(auc2)

#-------------------------------------------
# Treating Blood Pressure Variable
#-------------------------------------------

#----- Polynomial Models------
sq_bp_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + I(blood_pressure^2) + skin_thickness + bmi + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                   family = binomial)
summary(sq_bp_model)

cub_bp_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + I(blood_pressure^2) + I(blood_pressure^3) + skin_thickness + bmi + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                    family = binomial)
summary(cub_bp_model)

#-------- Compare AIC (poly models) --------------
AIC(basic_model, sq_bp_model, cub_bp_model)

#-------- Compare ROC/AUC (poly models) --------------

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(basic_model, newdata = train_data, type = "response"))
roc2 <- roc(train_data$diabetes, predict(sq_bp_model, newdata = train_data, type = "response"))
roc3 <- roc(train_data$diabetes, predict(cub_bp_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
auc3 <- auc(roc3)
print(auc1)
print(auc2)
print(auc3)

#------- Likelihood ratio test-------
#Quadratic model
lrtest(sq_bp_model, basic_model)
#Cubic model
lrtest(cub_bp_model, basic_model)

#------ Categorise blood pressure-----

# Define the function to categorise blood pressure
categorize_blood_pressure <- function(diastolic) {
  if (diastolic < 60) {
    return("Low")
  } else if (diastolic >= 60 & diastolic <= 80) {
    return("Normal")
  } else if (diastolic > 80 & diastolic <= 89) {
    return("High (Stage 1)")
  } else if (diastolic >= 90 & diastolic < 120) {
    return("High (Stage 2)")
  } else if (diastolic >= 120) {
    return("Hypertensive Crisis")
  } else {
    return(NA)
  }
}

# Apply the function to the blood_pressure column and create a new column cat_bp
clean_data$cat_1_bp <- sapply(clean_data$blood_pressure, categorize_blood_pressure)

clean_data$cat_1_bp <- factor(clean_data$cat_1_bp, 
                              levels = c("Normal", "Low","High (Stage 1)", 
                                         "High (Stage 2)", "Hypertensive Crisis"))

# Check the distribution of the categorised column
table(clean_data$cat_1_bp)

# Define the function to categorise blood pressure
categorize_blood_pressure_2 <- function(diastolic) {
  if (diastolic < 60) {
    return("Low")
  } else if (diastolic >= 60 & diastolic <= 80) {
    return("Normal")
  } else if (diastolic > 80 & diastolic <= 89) {
    return("High (Stage 1)")
  } else if (diastolic >= 90) {
    return("High (Stage 2)")
  } else {
    return(NA)
  }
}

# Apply the function to the blood_pressure column and create a new column cat_bp
clean_data$cat_2_bp <- sapply(clean_data$blood_pressure, categorize_blood_pressure_2)

clean_data$cat_2_bp <- factor(clean_data$cat_2_bp, 
                              levels = c("Normal", "Low","High (Stage 1)", 
                                         "High (Stage 2)"))

# Check the distribution of the categorised column
table(clean_data$cat_2_bp)


cat_1_bp_model <- glm(diabetes ~ pregnancies + glucose + cat_1_bp + skin_thickness + bmi + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                      family = binomial)
summary(cat_1_bp_model)

cat_2_bp_model <- glm(diabetes ~ pregnancies + glucose + cat_2_bp + skin_thickness + bmi + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                      family = binomial)
summary(cat_2_bp_model)

#-------- Compare AIC (categorised models) --------------
AIC(basic_model, cat_1_bp_model, cat_2_bp_model)

#-------- Compare ROC/AUC (categorised models) --------------
# Filter the data for S == 0
train_data <- clean_data %>% filter(S == 0)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(basic_model, newdata = train_data, type = "response"))
roc2 <- roc(train_data$diabetes, predict(cat_1_bp_model, newdata = train_data, type = "response"))
roc3 <- roc(train_data$diabetes, predict(cat_2_bp_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
auc3 <- auc(roc3)
print(auc1)
print(auc2)
print(auc3)

#------- Likelihood ratio test-------
lrtest(cat_1_bp_model, basic_model)
lrtest(cat_2_bp_model, basic_model)

#-------- Clean data accordingly --------------
clean_data <- clean_data %>%
  dplyr::select(-cat_1_bp, -cat_2_bp)

#------ Fractional Polynomials -----

fp_bp_model <- mfp(
  formula = diabetes ~ pregnancies + glucose + fp(blood_pressure) + skin_thickness + bmi + diabetes_genetic_score + age,
  family = binomial,
  data = clean_data %>% filter(S == 0)
)
summary(fp_bp_model)

#-------- Compare AIC (fractional polynomial model) --------------
AIC(basic_model, fp_bp_model)

#-------- Compare ROC/AUC (fractional polynomial model) --------------

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(basic_model, newdata = train_data, type = "response"))
roc2 <- roc(train_data$diabetes, predict(fp_bp_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
print(auc1)
print(auc2)

#------- Likelihood ratio test-------
lrtest(fp_bp_model, basic_model)

#------ Linear Splines -----

#------- Quartile Knots---------

# Calculate quartiles for blood pressure
quartiles <- quantile(train_data$blood_pressure, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# Print quartiles for reference
print(quartiles)

# Fit logistic regression model with linear splines at quartile knots
spl_quart_p_model <- glm(
  diabetes ~ lspline(blood_pressure, knots = c(quartiles), marginal = FALSE) +
    pregnancies + glucose + skin_thickness + bmi + diabetes_genetic_score + age,
  family = binomial,
  data = clean_data %>% filter(S == 0)
)
summary(spl_quart_p_model)

#-------- Compare AIC (linear spline model) --------------
AIC(basic_model, spl_quart_p_model)

#------- Logit Visualised Knots ---------

#Plot Log odds vs Age
train_data$bp_logit <- predict(basic_model, newdata = train_data, type = "link")
ggplot(train_data, aes(x = blood_pressure, y = bp_logit)) +
  geom_point(alpha = 0.5, color = "black") + 
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Predicted Log-Odds vs. Blood Pressure", 
       x = "Transformed Variable", 
       y = "Predicted Log-Odds") +
  theme_minimal()

# Fit logistic regression model with custom knots
spline_model_custom <- glm(
  diabetes ~ lspline(blood_pressure, knots = c(50, 60, 75, 100), marginal = FALSE) +
    pregnancies + glucose + skin_thickness + bmi + diabetes_genetic_score + age,
  family = binomial,
  data = clean_data %>% filter(S == 0)
)

# View summary of the model
summary(spline_model_custom)

#-------- Compare AIC (linear spline model) --------------
AIC(basic_model, spline_model_custom)

#------- Cubic Splines --------

# Create a datadist object
ddist <- datadist(clean_data %>% filter(S == 0))
options(datadist = "ddist")  # Specify the datadist object

# Fit the cubic spline model with 4 knots
spline_model_cubic <- ols(
  diabetes ~ rcs(blood_pressure, 4) +
    pregnancies + glucose + skin_thickness + bmi + diabetes_genetic_score + age,
  data = clean_data %>% filter(S == 0)
)

# Print the model summary
print(spline_model_cubic)

#-------- Compare AIC (cubic spline model) --------------
AIC(basic_model, spline_model_cubic)

#--------------------------------------------
# Create BMI categories
#-------------------------------------------

# Create BMI categories
clean_data <- clean_data %>%
  mutate(bmi_group = case_when(
    bmi < 25 ~ 1, # Normal weight
    bmi >= 25 & bmi < 30 ~ 2, # Overweight
    bmi >= 30 & bmi < 40 ~ 3, # Obese
    bmi >= 40 & bmi < 50 ~ 4, # Morbidly obese
    bmi >= 50 ~ 5)) # Super Morbidly obese


table(clean_data$bmi_group)

clean_data$bmi_group <- factor(clean_data$bmi_group, 
                               levels = c(1, 2, 3, 4, 5), 
                               labels = c("Normal Weight", "Overweight", "Obese", "Morbidly Obese", "Super Morbidly Obese"))

fac_bmi_model <- glm(diabetes ~ pregnancies + glucose + blood_pressure + skin_thickness + as_factor(bmi_group) + diabetes_genetic_score + age, data = clean_data %>% filter(S == 0), 
                     family = binomial)
summary(fac_bmi_model)

#-------- Compare AIC--------------
AIC(basic_model, fac_bmi_model)

#-------- Compare ROC/AUC --------------

# Filter the data for S == 0
train_data <- clean_data %>% filter(S == 0)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(basic_model, newdata = train_data, type = "response"))
roc2 <- roc(train_data$diabetes, predict(fac_bmi_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
print(auc1)
print(auc2)

#------- Likelihood ratio test-------
lrtest(fac_bmi_model, basic_model)

#--------------------------------------------
# Interaction Terms
#--------------------------------------------

#---------Transformed Basic model--------------

trans_model <- glm(diabetes ~ pregnancies_fac + glucose + blood_pressure + skin_thickness + as_factor(bmi_group) + diabetes_genetic_score + as_factor(age_group), data = clean_data %>% filter(S == 0), 
                   family = binomial)
summary(trans_model)

#---------Create models with Interaction Terms----------

gl_dia_model <- glm(diabetes ~ glucose*diabetes_genetic_score + pregnancies_fac + glucose + blood_pressure + skin_thickness + as_factor(bmi_group) + diabetes_genetic_score + as_factor(age_group), data = clean_data %>% filter(S == 0), 
                    family = binomial)
summary(gl_dia_model)
## Interaction term is statistically significant p<0.001

gl_bmi_model <- glm(diabetes ~ glucose*as_factor(bmi_group) + pregnancies_fac + glucose + blood_pressure + skin_thickness + as_factor(bmi_group) + diabetes_genetic_score + as_factor(age_group), data = clean_data %>% filter(S == 0), 
                    family = binomial)
summary(gl_bmi_model)
## No terms statistically significant

gl_age_model <- glm(diabetes ~ glucose*as_factor(age_group) + pregnancies_fac + glucose + blood_pressure + skin_thickness + as_factor(bmi_group) + diabetes_genetic_score + as_factor(age_group), data = clean_data %>% filter(S == 0), 
                    family = binomial)
summary(gl_age_model)
## Interaction term is statistically significant p<0.01

gen_age_model <- glm(diabetes ~ diabetes_genetic_score*as_factor(age_group) + pregnancies_fac + glucose + blood_pressure + skin_thickness + as_factor(bmi_group) + diabetes_genetic_score + as_factor(age_group), data = clean_data %>% filter(S == 0), 
                     family = binomial)
summary(gen_age_model)
## Not Statistically significant

gen_bmi_model <- glm(diabetes ~ diabetes_genetic_score*as_factor(bmi_group) + pregnancies_fac + glucose + blood_pressure + skin_thickness + as_factor(bmi_group) + diabetes_genetic_score + as_factor(age_group), data = clean_data %>% filter(S == 0), 
                     family = binomial)
summary(gen_bmi_model)
## Not Statistically significant

#-------- Compare AIC --------------
#Compare AIC of statistically significant glucose interaction terms
AIC(trans_model, gl_dia_model, gl_bmi_model, gl_age_model, gen_age_model, gen_bmi_model)

#-------- Compare ROC/AUC (categorised models) --------------
# Filter the data for S == 0
train_data <- clean_data %>% filter(S == 0)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(trans_model, newdata = train_data, type = "response"))
roc2 <- roc(train_data$diabetes, predict(gl_dia_model, newdata = train_data, type = "response"))
roc3 <- roc(train_data$diabetes, predict(gl_bmi_model, newdata = train_data, type = "response"))
roc4 <- roc(train_data$diabetes, predict(gl_age_model, newdata = train_data, type = "response"))
roc5 <- roc(train_data$diabetes, predict(gen_age_model, newdata = train_data, type = "response"))
roc6 <- roc(train_data$diabetes, predict(gen_bmi_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
auc2 <- auc(roc2)
auc3 <- auc(roc3)
auc4 <- auc(roc4)
auc5 <- auc(roc5)
auc6 <- auc(roc6)
print(auc1)
print(auc2)
print(auc3)
print(auc4)
print(auc5)
print(auc6)

#------- Likelihood ratio test-------
lrtest(gl_bmi_model, trans_model)

#--------------------------------------------
# Final Model Design
#--------------------------------------------

#--------Null Model---------

# Fit a null model
null_model <- glm(diabetes ~ 1, data = clean_data %>% filter(S == 0), family = binomial)
summary(null_model)

AIC(null_model)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(null_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
print(auc1)

# Log-likelihood
log_likelihood <- logLik(null_model)
print(log_likelihood)

#--------Full Model---------

#Fit full model
full_model <- glm(diabetes ~ pregnancies_fac + glucose + blood_pressure + skin_thickness + diabetes_genetic_score + factor(age_group) + factor(bmi_group), data = clean_data %>% filter(S == 0), 
                  family = binomial)
summary(full_model)

AIC(full_model)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(full_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
print(auc1)

# Log-likelihood
log_likelihood <- logLik(full_model)
print(log_likelihood)

#--------Less insignificant variables Model---------

#Fit model
model_3 <- glm(diabetes ~ glucose + skin_thickness + diabetes_genetic_score + factor(age_group) + factor(bmi_group), data = clean_data %>% filter(S == 0), 
               family = binomial)
summary(model_3)

AIC(model_3)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(model_3, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
print(auc1)

# Log-likelihood
log_likelihood <- logLik(model_3)
print(log_likelihood)

#--------Add Gluc*BMI ---------

#Fit model
model_4 <- glm(diabetes ~ glucose + skin_thickness + diabetes_genetic_score + factor(age_group) + factor(bmi_group) + glucose*factor(bmi_group), data = clean_data %>% filter(S == 0), 
               family = binomial)
summary(model_4)

AIC(model_4)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(model_4, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
print(auc1)

# Log-likelihood
log_likelihood <- logLik(model_4)
print(log_likelihood)

#--------------------------------------------
# Influential Data Points
#--------------------------------------------

# Filter the data for S == 0
train_data <- clean_data %>% filter(S == 0)

#-------Leverage vs Residuals-----------

train_data <- train_data %>%
  mutate(
    leverage = hatvalues(model_4),                            # Leverage (hat values)
    residuals_raw = residuals(model_4, type = "response"),    # Raw residuals
    rss = sum(residuals_raw^2),                               # Residual sum of squares
    norm_resid = residuals_raw / sqrt(rss),                   # Normalized residuals
    norm_resid_sq = norm_resid^2                              # Squared normalized residuals
  )

# Plot leverage vs squared normalized residuals
ggplot(train_data, aes(x = norm_resid_sq, y = leverage)) +
  geom_point() +
  labs(title = "Leverage vs. Squared Standardized Residuals for model_4",
       x = "Normalized Residuals Squared",
       y = "Leverage") +
  theme_minimal()

#--------Influence Plot-----------

#Influence plot showing leverage and Cook's distance
influencePlot(model_4, data = train_data, main = "Influence Plot", id.method = "identify", id.n = 20)

#--------- dfBetas - Skin thickness ---------

# Calculate DFBETAS for model_4
dfbetas_values <- dfbetas(model_4)

# Access DFBETAS for the variable "skin_thickness"
train_data$dfbeta_skin_thickness <- dfbetas_values[, "skin_thickness"]

# Generate absolute DFBETAS
train_data$abs_dfbeta_skin_thickness <- abs(train_data$dfbeta_skin_thickness)

# Histogram of absolute DFBETAS
ggplot(train_data, aes(x = abs_dfbeta_skin_thickness)) +
  geom_histogram(binwidth = 0.005, fill = "blue", color = "black") +
  labs(title = "Frequency Plot of Absolute DFBETA for Skin Thickness",
       x = "Absolute DFBETA (Skin Thickness)",
       y = "Frequency") +
  theme_minimal()

# Determine the cut-off for DFBETAS
cutoff <- 2 / sqrt(nrow(train_data))
cat("Cut-off for DFBETAS:", cutoff, "\n")

#---------Re-fit final model with clean data-----------

# Remove patient 2360 from the clean_data dataset
clean_data <- clean_data %>% filter(!patient_id %in% c(2360, 2337))
train_data <- clean_data %>% filter(S == 0)

#Fit model
final_model<- glm(diabetes ~ glucose + skin_thickness + diabetes_genetic_score + factor(age_group) + factor(bmi_group) + glucose*factor(bmi_group), data = clean_data %>% filter(S == 0), 
                  family = binomial)
summary(final_model)

AIC(final_model)

# Generate ROC curves
roc1 <- roc(train_data$diabetes, predict(final_model, newdata = train_data, type = "response"))

# Calculate AUC
auc1 <- auc(roc1)
print(auc1)

# Log-likelihood
log_likelihood <- logLik(final_model)
print(log_likelihood)

#--------------------------------------------
# Model Interpretation
#--------------------------------------------

#Print model summary
summary(final_model)

#Analyse odds ratios
odds_ratios <- exp(coef(final_model))
print(odds_ratios)

#Analyse odds ratios
tidy_model <- tidy(final_model, exponentiate = TRUE, conf.int = TRUE)
print(tidy_model)

#--------Understand Effect of Interaction Term---------

#------Plot Glucose vs BMI Category--------

# Define the range of glucose values
glucose_vals <- seq(50, 200, by = 5)  # Glucose range


# Define BMI categories as a factor with labels
bmi_categories <- factor(
  c(1, 2, 3, 4, 5), 
  levels = c(1, 2, 3, 4, 5),
  labels = c("Normal Weight", "Overweight", "Obese", "Morbidly Obese", "Super Morbidly Obese")
)

# Create a prediction grid with labeled BMI categories
prediction_grid <- expand.grid(
  glucose = glucose_vals,
  bmi_group = bmi_categories  # Use labeled BMI categories
)


# Create a data frame for predictions
prediction_grid <- expand.grid(
  glucose = glucose_vals,
  bmi_group = bmi_categories  # Use numeric BMI categories
)

# Add other variables with constant/default values
prediction_grid$skin_thickness <- mean(train_data$skin_thickness, na.rm = TRUE)
prediction_grid$diabetes_genetic_score <- mean(train_data$diabetes_genetic_score, na.rm = TRUE)

# Add the age variable with the reference category (20–30 years)
prediction_grid$age_group <- "18-29"

# Predict probabilities from the model (replace final_model with your fitted model)
prediction_grid$probability <- predict(final_model, newdata = prediction_grid, type = "response")

# Plot Glucose vs Probability for numeric BMI categories
library(ggplot2)
ggplot(prediction_grid, aes(x = glucose, y = probability, color = as.factor(bmi_group))) +
  geom_line(size = 1) +
  labs(
    title = "Glucose vs Probability of Diabetes at Fixed BMI categories",
    x = "Glucose (mg/dl)",
    y = "Probability of Diabetes",
    color = "BMI Category"
  ) +
  theme_minimal()


#--------------------------------------------
# Model Validation
#--------------------------------------------

# Predict risk of death for entire dataset (S=0 (train) and S=1 (test))

## Predict the linear predictor (log-odds) for all individuals
clean_data$log_odds <- predict(final_model, newdata = clean_data, type = "link")

## Predict the probabilities (i.e., risk) for all individuals
clean_data$probs <- predict(final_model, newdata = clean_data, type = "response")

## Summarise predicted probabilities by diabested status
clean_data %>% 
  group_by(diabetes) %>%
  summarise(mean_prob = mean(probs, na.rm = TRUE),
            sd_prob   = sd(probs, na.rm = TRUE), #standard deviation 
            min_prob  = min(probs, na.rm = TRUE), #min
            max_prob  = max(probs, na.rm = TRUE)) #max

# Create a boxplot of predicted probabilities by the 'dead' variable
boxplot(probs ~ as_factor(diabetes), data = clean_data,
        xlab = "Diabetes Status", ylab = "Predicted Probabilities",
        main = "Boxplot of Predicted Probabilities by Diabetes Status")

#---------ROC Curves----------

# ROC curve for S == 0
roc_0 <- roc(clean_data$diabetes[clean_data$S == 0], clean_data$probs[clean_data$S == 0])

# Plot ROC curve with truncated x-axis for specificity
plot(roc_0, main = "ROC Curve for S=0", print.thres = T, print.auc
     = TRUE, xlim = c(0, 1), ylim = c(0, 1))

# ROC curve for S == 1
roc_1 <- roc(clean_data$diabetes[clean_data$S == 1], clean_data$probs[clean_data$S == 1])

# Plot ROC curve with truncated x-axis for specificity
plot(roc_1, main = "ROC Curve for S=1", print.thres = F, print.auc
     = T, xlim = c(0, 1), ylim = c(0, 1), )

#Combine into single plot
# Plot ROC curve for S == 0
plot(roc_0, main = "ROC Curves for S=0 and S=1", 
     xlab = "Specificity", ylab = "Sensitivity", 
     col = "blue", print.auc = TRUE, 
     xlim = c(0, 1), ylim = c(0, 1), print.auc.y = 0.4)

# Add the ROC curve for S == 1 to the same plot
plot(roc_1, col = "red", add = TRUE, print.auc = TRUE,
     print.auc.y = 0.3)

# Add a legend
legend("bottomright", legend = c("S = 0", "S = 1"), 
       col = c("blue", "red"), lwd = 2)

#---------Sensitivity and Specificity----------

#-------------Contingency Table S=0-------------------

# Filter data for S = 0
filtered_data <- clean_data[clean_data$S == 0, ]

# Define a threshold for classification
threshold <- 0.5  

# Generate predicted classes based on the threshold
filtered_data$predicted_class <- ifelse(filtered_data$probs > threshold, 1, 0)

# Create a 2x2 contingency table
contingency_table <- table(Predicted = filtered_data$predicted_class, Actual = filtered_data$diabetes)

# Print the contingency table
print("2x2 Contingency Table for S = 0:")
print(contingency_table)

# Extract values from the contingency table
true_positive <- ifelse("1" %in% rownames(contingency_table) & "1" %in% colnames(contingency_table), contingency_table["1", "1"], 0)
false_negative <- ifelse("0" %in% rownames(contingency_table) & "1" %in% colnames(contingency_table), contingency_table["0", "1"], 0)
true_negative <- ifelse("0" %in% rownames(contingency_table) & "0" %in% colnames(contingency_table), contingency_table["0", "0"], 0)
false_positive <- ifelse("1" %in% rownames(contingency_table) & "0" %in% colnames(contingency_table), contingency_table["1", "0"], 0)

# Calculate metrics
sensitivity <- true_positive / (true_positive + false_negative)  # TP / (TP + FN)
specificity <- true_negative / (true_negative + false_positive)  # TN / (TN + FP)
PPV <- true_positive / (true_positive + false_positive)    # TP / (TP + FP)
NPV <- true_negative / (true_negative + false_negative)  # TN / (TN + FN)

# Print metrics
print(paste("Sensitivity:", round(sensitivity, 3)))
print(paste("Specificity:", round(specificity, 3)))
print(paste("PPV:", round(PPV, 3)))
print(paste("NPV:", round(NPV, 3)))

#-------------Contingency Table S=1-------------------

# Filter data for S = 0
filtered_data <- clean_data[clean_data$S == 1, ]

# Define a threshold for classification
threshold <- 0.5  

# Generate predicted classes based on the threshold
filtered_data$predicted_class <- ifelse(filtered_data$probs > threshold, 1, 0)

# Create a 2x2 contingency table
contingency_table <- table(Predicted = filtered_data$predicted_class, Actual = filtered_data$diabetes)

# Print the contingency table
print("2x2 Contingency Table for S = 0:")
print(contingency_table)

# Extract values from the contingency table
true_positive <- ifelse("1" %in% rownames(contingency_table) & "1" %in% colnames(contingency_table), contingency_table["1", "1"], 0)
false_negative <- ifelse("0" %in% rownames(contingency_table) & "1" %in% colnames(contingency_table), contingency_table["0", "1"], 0)
true_negative <- ifelse("0" %in% rownames(contingency_table) & "0" %in% colnames(contingency_table), contingency_table["0", "0"], 0)
false_positive <- ifelse("1" %in% rownames(contingency_table) & "0" %in% colnames(contingency_table), contingency_table["1", "0"], 0)

# Calculate metrics
sensitivity <- true_positive / (true_positive + false_negative)  # TP / (TP + FN)
specificity <- true_negative / (true_negative + false_positive)  # TN / (TN + FP)
PPV <- true_positive / (true_positive + false_positive)    # TP / (TP + FP)
NPV <- true_negative / (true_negative + false_negative)  # TN / (TN + FN)

# Print metrics
print(paste("Sensitivity:", round(sensitivity, 3)))
print(paste("Specificity:", round(specificity, 3)))
print(paste("PPV:", round(PPV, 3)))
print(paste("NPV:", round(NPV, 3)))

#--------Analyse model performance across dif threshholds (S=0) ----------

# Filter the data for S = 0
filtered_data <- clean_data[clean_data$S == 0, ]

thresholds <- seq(0.1, 0.9, by = 0.1)
results <- data.frame()

for (t in thresholds) {
  # Filter data for S = 0
  filtered_data <- clean_data[clean_data$S == 0, ]
  
  # Generate predicted classes based on the threshold
  filtered_data$predicted_class <- ifelse(filtered_data$probs > t, 1, 0)
  
  # Create contingency table
  contingency_table <- table(Predicted = filtered_data$predicted_class, Actual = filtered_data$diabetes)
  
  # Extract values from contingency table
  true_positive <- contingency_table[2, 2]
  false_negative <- contingency_table[1, 2]
  true_negative <- contingency_table[1, 1]
  false_positive <- contingency_table[2, 1]
  
  # Calculate metrics
  sensitivity <- true_positive / (true_positive + false_negative)  # TP / (TP + FN)
  specificity <- true_negative / (true_negative + false_positive)  # TN / (TN + FP)
  PPV <- true_positive / (true_positive + false_positive)    # TP / (TP + FP)
  NPV <- true_negative / (true_negative + false_negative)  # TN / (TN + FN)
  
  # Append results
  results <- rbind(results, data.frame(
    Threshold = t,
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = PPV,
    NPV = NPV
  ))
}

# Print the results
print(results)

#--------Analyse model performance across dif threshholds (S=1) ----------

thresholds <- seq(0.1, 0.9, by = 0.1)
results <- data.frame()

for (t in thresholds) {
  # Filter data for S = 0
  filtered_data <- clean_data[clean_data$S == 1, ]
  
  # Generate predicted classes based on the threshold
  filtered_data$predicted_class <- ifelse(filtered_data$probs > t, 1, 0)
  
  # Create contingency table
  contingency_table <- table(Predicted = filtered_data$predicted_class, Actual = filtered_data$diabetes)
  
  # Extract values from contingency table
  true_positive <- contingency_table[2, 2]
  false_negative <- contingency_table[1, 2]
  true_negative <- contingency_table[1, 1]
  false_positive <- contingency_table[2, 1]
  
  # Calculate metrics
  sensitivity <- true_positive / (true_positive + false_negative)  # TP / (TP + FN)
  specificity <- true_negative / (true_negative + false_positive)  # TN / (TN + FP)
  PPV <- true_positive / (true_positive + false_positive)    # TP / (TP + FP)
  NPV <- true_negative / (true_negative + false_negative)  # TN / (TN + FN)
  
  # Append results
  results <- rbind(results, data.frame(
    Threshold = t,
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = PPV,
    NPV = NPV
  ))
}

# Print the results
print(results)

#---------Hosmer-Lemeshow (S=0)  ---------

# Create a Hosmer-Lemeshow goodness of fit table for S=0
hl_data_0 <- na.omit(data.frame(diabetes = clean_data$diabetes[clean_data$S == 0], fitted_values = clean_data$probs[clean_data$S == 0]))

# Perform Hosmer-Lemeshow test for S == 0
hl_test_0 <- hoslem.test(hl_data_0$diabetes, hl_data_0$fitted_values, g = 10)
print(hl_test_0)

# Tabulate observed vs expected counts by deciles for S == 0
hl_table_S0 <- clean_data %>%
  filter(S == 0) %>%
  mutate(decile = ntile(probs, 10)) %>%     # Create decile groups
  group_by(decile) %>%
  summarise(
    observed_0 = sum(diabetes == 0),
    observed_1 = sum(diabetes == 1),
    expected_0 = sum(1 - probs),
    expected_1 = sum(probs),
    total_observed = observed_0 + observed_1,
    total_expected = expected_0 + expected_1
  )
print("Hosmer-Lemeshow Table for S = 0")
print(hl_table_S0)

#---------Graph Hosmer Lemeshow (S=0)----------

# Create decile groups for predicted probabilities (S == 0 and S == 1)
clean_data <- clean_data %>%
  mutate(
    # Create deciles for S == 0
    decile_s0 = if_else(S == 0, 
                        as.integer(cut(probs, 
                                       breaks = quantile(probs[S == 0], probs = seq(0, 1, by = 0.1)), 
                                       include.lowest = TRUE, labels = FALSE)) - 1, 
                        NA_integer_),
    # Create deciles for S == 1
    decile_s1 = if_else(S == 1, 
                        as.integer(cut(probs, 
                                       breaks = quantile(probs[S == 1], probs = seq(0, 1, by = 0.1)), 
                                       include.lowest = TRUE, labels = FALSE)) - 1, 
                        NA_integer_)
  )

# Summarize data for S == 0
data_s0 <- clean_data %>%
  filter(S == 0) %>%
  group_by(decile_s0) %>%
  summarise(
    sum_probs = sum(probs),  # Sum of predicted probabilities
    sum_diabetes = sum(diabetes)  # Count of observed diabetes cases
  ) %>%
  pivot_longer(cols = c(sum_probs, sum_diabetes), names_to = "variable", values_to = "value")

# Bar plot for S == 0
library(ggplot2)
ggplot(data_s0, aes(x = factor(decile_s0), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(
    values = c("sum_probs" = "blue", "sum_diabetes" = "red")) +
  labs(x = "Predicted Probability Deciles (S=0)", y = "Sum",
       title = "Sum of Predicted Probabilities and Observed Diabetes Cases (S=0)", fill = "") +
  theme_minimal()

#---------Hosmer-Lemeshow (S=1) ---------

# Create a Hosmer-Lemeshow goodness of fit table for S=1
hl_data_1 <- na.omit(data.frame(diabetes = clean_data$diabetes[clean_data$S == 1], fitted_values = clean_data$probs[clean_data$S == 1]))

# Perform Hosmer-Lemeshow test for S == 1
hl_test_1 <- hoslem.test(hl_data_1$diabetes, hl_data_1$fitted_values, g = 10)
print(hl_test_1)

# Tabulate observed vs expected counts by deciles for S == 1
hl_table_S1 <- clean_data %>%
  filter(S == 1) %>%
  mutate(decile = ntile(probs, 10)) %>%     # Create decile groups
  group_by(decile) %>%
  summarise(
    observed_0 = sum(diabetes == 0),
    observed_1 = sum(diabetes == 1),
    expected_0 = sum(1 - probs),
    expected_1 = sum(probs),
    total_observed = observed_0 + observed_1,
    total_expected = expected_0 + expected_1
  )
print("Hosmer-Lemeshow Table for S = 1")
print(hl_table_S1)

#---------Graph Hosmer Lemeshow (S=1)--------

# Summarize and reshape data for S == 1 (for plotting)
data_s1 <- clean_data %>%
  filter(S == 1) %>%
  group_by(decile_s1) %>%  # Use the correct decile column
  summarise(
    sum_probs = sum(probs, na.rm = TRUE),        # Sum of predicted probabilities
    sum_diabetes = sum(diabetes, na.rm = TRUE)  # Sum of observed diabetes cases
  ) %>%
  # Re-shape the data from wide to long
  pivot_longer(cols = c(sum_probs, sum_diabetes), names_to = "variable", values_to = "value")

# Bar plot for S == 1
library(ggplot2)
ggplot(data_s1, aes(x = factor(decile_s1), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(
    values = c("sum_probs" = "blue", "sum_diabetes" = "red")) +
  labs(x = "Predicted Probability Groups (S=1)", y = "Sum",
       title = "Sum of Predicted Probabilities and Diabetes (S=1)", fill = "") +
  theme_minimal()

#-------- Calibration Plot (S=0) ----------

# Step 1: Prepare the data
calibration_data <- clean_data %>%
  filter(S == 0) %>%                # Filter for S == 0
  mutate(decile = ntile(probs, 10)) %>% # Create deciles based on predicted probabilities
  group_by(decile) %>%
  summarise(
    mean_predicted = mean(probs),         # Mean predicted probability
    mean_observed = mean(diabetes),      # Proportion of observed diabetes cases
    count = n()                          # Count of individuals in each decile
  )

# Step 2: Calculate Calibration in the Large
mean_predicted <- mean(clean_data$probs[clean_data$S == 0])  # Mean predicted probability
observed_event_rate <- mean(clean_data$diabetes[clean_data$S == 0])  # Observed event rate
calibration_in_the_large <- mean_predicted - observed_event_rate  # Calibration in the large

# Step 3: Create the calibration plot
ggplot(calibration_data, aes(x = mean_predicted, y = mean_observed)) +
  geom_point(size = 3) +                  # Add points for each decile
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + # Perfect calibration line
  labs(
    title = "Calibration Plot (S=0)",
    x = "Mean Predicted",
    y = "Mean Observed"
  ) +
  theme_minimal() +
  annotate(
    "text", x = 0.5, y = 0.1, 
    label = paste0(
      "Calibration in the Large: ", round(calibration_in_the_large, 3), 
      "\nMean Predicted: ", round(mean_predicted, 3), 
      "\nObserved Event Rate: ", round(observed_event_rate, 3)
    ), 
    size = 4, hjust = 0, color = "blue"
  )

#-------- Calibration Plot (S=1) ----------

# Step 1: Prepare the data
calibration_data_1 <- clean_data %>%
  filter(S == 1) %>%                # Filter for S == 0
  mutate(decile = ntile(probs, 10)) %>% # Create deciles based on predicted probabilities
  group_by(decile) %>%
  summarise(
    mean_predicted = mean(probs),         # Mean predicted probability
    mean_observed = mean(diabetes),      # Proportion of observed diabetes cases
    count = n()                          # Count of individuals in each decile
  )

# Step 2: Calculate Calibration in the Large
mean_predicted <- mean(clean_data$probs[clean_data$S == 1])  # Mean predicted probability
observed_event_rate <- mean(clean_data$diabetes[clean_data$S == 1])  # Observed event rate
calibration_in_the_large <- mean_predicted - observed_event_rate  # Calibration in the large

# Step 3: Create the calibration plot
ggplot(calibration_data_1, aes(x = mean_predicted, y = mean_observed)) +
  geom_point(size = 3) +                  # Add points for each decile
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + # Perfect calibration line
  labs(
    title = "Calibration Plot (S=1)",
    x = "Mean Predicted",
    y = "Mean Observed"
  ) +
  theme_minimal() +
  annotate(
    "text", x = 0.5, y = 0.1, 
    label = paste0(
      "Calibration in the Large: ", round(calibration_in_the_large, 3), 
      "\nMean Predicted: ", round(mean_predicted, 3), 
      "\nObserved Event Rate: ", round(observed_event_rate, 3)
    ), 
    size = 4, hjust = 0, color = "blue"
  )

#-------- Calibration Plot (full dataset) ----------

# Step 1: Prepare the data
calibration_data_2 <- clean_data %>%
  mutate(decile = ntile(probs, 10)) %>% # Create deciles based on predicted probabilities
  group_by(decile) %>%
  summarise(
    mean_predicted = mean(probs),         # Mean predicted probability
    mean_observed = mean(diabetes),      # Proportion of observed diabetes cases
    count = n()                          # Count of individuals in each decile
  )

# Step 2: Create the calibration plot
ggplot(calibration_data_2, aes(x = mean_predicted, y = mean_observed)) +
  geom_point(size = 3) +                  # Add points for each decile
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + # Perfect calibration line
  labs(
    title = "Calibration Plot (full dataset)",
    x = "Mean Predicted",
    y = "Mean Observed"
  ) +
  theme_minimal()

#------------Sub-Group Analysis---------

#-----------BMI (S=0)----------

# Filter the data for S == 0
filtered_data <- clean_data %>% filter(S == 0)

# Initialize an empty data frame to store results
results <- data.frame(
  bmi_group = numeric(),
  sensitivity = numeric(),
  specificity = numeric()
)

# Loop through each BMI group
for (group in unique(filtered_data$bmi_group)) {
  
  # Subset data for the current BMI group
  subgroup_data <- filtered_data %>% filter(bmi_group == group)
  
  # Generate predictions
  predictions <- predict(final_model, newdata = subgroup_data, type = "response")
  
  # Classify predictions into binary outcomes (threshold = 0.5, modify if needed)
  predicted_class <- ifelse(predictions > 0.5, 1, 0)
  
  # Actual outcomes
  actual <- subgroup_data$diabetes  # Assuming 'diabetes' is the binary outcome (0 = no diabetes, 1 = diabetes)
  
  # Create confusion matrix
  confusion_matrix <- table(Predicted = predicted_class, Actual = actual)
  
  # Calculate sensitivity and specificity
  sensitivity <- ifelse(
    "1" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix),
    confusion_matrix["1", "1"] / sum(confusion_matrix[, "1"], na.rm = TRUE),
    NA
  )
  
  specificity <- ifelse(
    "0" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix),
    confusion_matrix["0", "0"] / sum(confusion_matrix[, "0"], na.rm = TRUE),
    NA
  )
  
  # Append results to the data frame
  results <- rbind(
    results,
    data.frame(
      bmi_group = group,
      sensitivity = sensitivity,
      specificity = specificity
    )
  )
}

# Display results
print(results)

#-----------BMI (S=1)----------

# Filter the data for S == 1
filtered_data <- clean_data %>% filter(S == 1)

# Initialize an empty data frame to store results
results <- data.frame(
  bmi_group = numeric(),
  sensitivity = numeric(),
  specificity = numeric()
)

# Loop through each BMI group
for (group in unique(filtered_data$bmi_group)) {
  
  # Subset data for the current BMI group
  subgroup_data <- filtered_data %>% filter(bmi_group == group)
  
  # Generate predictions
  predictions <- predict(final_model, newdata = subgroup_data, type = "response")
  
  # Classify predictions into binary outcomes (threshold = 0.5, modify if needed)
  predicted_class <- ifelse(predictions > 0.5, 1, 0)
  
  # Actual outcomes
  actual <- subgroup_data$diabetes  # Assuming 'diabetes' is the binary outcome (0 = no diabetes, 1 = diabetes)
  
  # Create confusion matrix
  confusion_matrix <- table(Predicted = predicted_class, Actual = actual)
  
  # Calculate sensitivity and specificity
  sensitivity <- ifelse(
    "1" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix),
    confusion_matrix["1", "1"] / sum(confusion_matrix[, "1"], na.rm = TRUE),
    NA
  )
  
  specificity <- ifelse(
    "0" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix),
    confusion_matrix["0", "0"] / sum(confusion_matrix[, "0"], na.rm = TRUE),
    NA
  )
  
  # Append results to the data frame
  results <- rbind(
    results,
    data.frame(
      bmi_group = group,
      sensitivity = sensitivity,
      specificity = specificity
    )
  )
}

# Display results
print(results)

#-----------Age S=0--------------

# Filter the data for S == 0
filtered_data <- clean_data %>% filter(S == 0)

# Initialize an empty data frame to store results
results <- data.frame(
  age_group = numeric(),
  sensitivity = numeric(),
  specificity = numeric()
)

# Loop through each age group
for (group in unique(filtered_data$age_group)) {
  
  # Subset data for the current age group
  subgroup_data <- filtered_data %>% filter(age_group == group)
  
  # Generate predictions
  predictions <- predict(final_model, newdata = subgroup_data, type = "response")
  
  # Classify predictions into binary outcomes (threshold = 0.5, modify if needed)
  predicted_class <- ifelse(predictions > 0.5, 1, 0)
  
  # Actual outcomes
  actual <- subgroup_data$diabetes  # Assuming 'diabetes' is the binary outcome (0 = no diabetes, 1 = diabetes)
  
  # Create confusion matrix
  confusion_matrix <- table(Predicted = predicted_class, Actual = actual)
  
  # Calculate sensitivity and specificity
  sensitivity <- ifelse(
    "1" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix),
    confusion_matrix["1", "1"] / sum(confusion_matrix[, "1"], na.rm = TRUE),
    NA
  )
  
  specificity <- ifelse(
    "0" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix),
    confusion_matrix["0", "0"] / sum(confusion_matrix[, "0"], na.rm = TRUE),
    NA
  )
  
  # Append results to the data frame
  results <- rbind(
    results,
    data.frame(
      age_group = group,
      sensitivity = sensitivity,
      specificity = specificity
    )
  )
}

# Display results
print(results)  

#-----------Age S=1--------------

# Filter the data for S == 1
filtered_data <- clean_data %>% filter(S == 1)

# Initialize an empty data frame to store results
results <- data.frame(
  age_group = numeric(),
  sensitivity = numeric(),
  specificity = numeric()
)

# Loop through each age group
for (group in unique(filtered_data$age_group)) {
  
  # Subset data for the current age group
  subgroup_data <- filtered_data %>% filter(age_group == group)
  
  # Generate predictions
  predictions <- predict(final_model, newdata = subgroup_data, type = "response")
  
  # Classify predictions into binary outcomes (threshold = 0.5, modify if needed)
  predicted_class <- ifelse(predictions > 0.5, 1, 0)
  
  # Actual outcomes
  actual <- subgroup_data$diabetes  # Assuming 'diabetes' is the binary outcome (0 = no diabetes, 1 = diabetes)
  
  # Create confusion matrix
  confusion_matrix <- table(Predicted = predicted_class, Actual = actual)
  
  # Calculate sensitivity and specificity
  sensitivity <- ifelse(
    "1" %in% rownames(confusion_matrix) && "1" %in% colnames(confusion_matrix),
    confusion_matrix["1", "1"] / sum(confusion_matrix[, "1"], na.rm = TRUE),
    NA
  )
  
  specificity <- ifelse(
    "0" %in% rownames(confusion_matrix) && "0" %in% colnames(confusion_matrix),
    confusion_matrix["0", "0"] / sum(confusion_matrix[, "0"], na.rm = TRUE),
    NA
  )
  
  # Append results to the data frame
  results <- rbind(
    results,
    data.frame(
      age_group = group,
      sensitivity = sensitivity,
      specificity = specificity
    )
  )
}

# Display results
print(results)  