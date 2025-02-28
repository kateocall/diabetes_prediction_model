# Diabetes Prediction Model Repository

## Model Development, Interpretation and Evaluation
The model is a logistic regression model which predicts if a patient has diabetes.

The R script (src/diabetes_model.R) contains the code for exploratory data analysis, model development, model evaluation and model interpretation. 

The model was internally validated as follows:
- Discrimination and Calibration: ROC curve, Hosmer-Lemeshow test
- Sensitivity Analysis: Confusion Matrix, Subgroup Analysis

## Dataset
The dataset (data/dataset.dta) for this project is a modified version of the Pima Indian Diabetes dataset from Kaggle. 

Original dataset: https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database
