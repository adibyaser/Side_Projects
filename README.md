# Project Bank-Data
Built-in decision tree in scikit-learn cannot handle categorical data well. However, a dataset can have both categorical and numeric features.
For categorical features, a feature can have more than 2 possible values. So a node which separates data based on a categorical feature can have more than 2 children.
In case of numeric feature, we won’t allow more than 2 children. To split based on a numeric feature, we sort the data according to the values of the numeric feature, and
then try splitting on every possible midpoint.
1. Implement a function named entropy(S,targetAttribute) to compute Entropy of a set S of examples with respect to a target attribute. 
2. Implement a function named informationGain(S,attribute) to compute Information gain for the set of examples S, if we split a S based on the given attribute.
3. Break ties (if any) arbitrarily. Implement decision tree algorithm ID3 (the one we discussed in class, uses entropy and information gain using the above functions.
