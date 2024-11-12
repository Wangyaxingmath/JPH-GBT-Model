# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 17:00:53 2024

@author: yxwang
"""

import numpy as np
import scipy as sp
import pandas as pd
import ast

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error





df_train = pd.read_csv('/home/yxwang/GEPH-GBT/train.csv')
df_train['Feature'] = df_train['Feature'].apply(ast.literal_eval)
X_train = np.stack(df_train['Feature'].values)
y_train = df_train['Label'].values


df_test = pd.read_csv('/home/yxwang/GEPH-GBT/test.csv')
df_test['Feature'] = df_test['Feature'].apply(ast.literal_eval)
X_test = np.stack(df_test['Feature'].values)
y_test = df_test['Label'].values




params = {'n_estimators': 60000, 'max_depth': 7, 'min_samples_split':2,
                'learning_rate': 0.001, 'loss': 'squared_error',
                'max_features':'sqrt','subsample':0.7}

result = []
for i in range(10):
    clf = GradientBoostingRegressor(**params)
    clf = clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    rmse = pow(mse, 0.5)
    pearcorr = sp.stats.pearsonr(y_test, y_pred)
    temp = [i, rmse, pearcorr[0]]
    result.append(temp)
    
print(result)


    





    

    
    
