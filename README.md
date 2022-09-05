# Project Bank-Data
Built-in decision tree in scikit-learn cannot handle categorical data well. However, a dataset can have both categorical and numeric features.
For categorical features, a feature can have more than 2 possible values. So a node which separates data based on a categorical feature can have more than 2 children.
In case of numeric feature, we won’t allow more than 2 children. To split based on a numeric feature, we sort the data according to the values of the numeric feature, and
then try splitting on every possible midpoint.
1. Implement a function named entropy(S,targetAttribute) to compute Entropy of a set S of examples with respect to a target attribute. 
2. Implement a function named informationGain(S,attribute) to compute Information gain for the set of examples S, if we split a S based on the given attribute.
3. Break ties (if any) arbitrarily. Implement decision tree algorithm ID3 (the one we discussed in class, uses entropy and information gain using the above functions.

# Housing Data
Boston Dataset
Contains information about different houses in Boston. There are 506 samples and 13 feature variables in this dataset.
https://archive.ics.uci.edu/ml/machine-learning-databases/housing/
Loading and Analysis of Dataset (EDA)
        Your dataset is split into two files, combine them into a CSV. Then load them using pandas and do the following data analysis (atleast)

Describe
Unique and Null values
HeatMap
KDE-plot
Box plot

Train-Test Split (70/30, 80/20, 75/25)
Training the models
Evaluation of Models
Confusion Matrix (Is it possible?)
R^2 value (Use Adjusted R^2 value)
MAE value
MSE value
RMSE value

Report
Report comparing all models 

R^2 Value

MAE value

MSE Value

RMSE value

KNN

SVM

LR

Report for different test-train split (same table - different format)
K-Fold Cross Validation (Bonus)
