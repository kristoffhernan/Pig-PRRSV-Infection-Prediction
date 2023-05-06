import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import splev, splrep
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import GridSearchCV


# Splits data into those with and without na given a gene
def na_split(gene, gd_df):
    eg = gd_df[[gene, 'dpi']]
    eg = eg.reset_index()

    # df with x as days and y and individual
    df = eg.pivot(index='index', columns='dpi')[gene].reset_index(drop=True).drop(11, axis=1)

    # df with day and ge columns
    melt = pd.melt(df, value_vars=df.columns, var_name='day', value_name='ge')

    # splitting into two df, with and without NA values
    w_na = melt[melt.isna().any(axis=1)].sort_values('day')
    wo_na = melt[melt.notna().all(axis=1)].sort_values('day')
    
    return(w_na, wo_na)


# includes two different imputing methods to estimate the dpi or gene expression one of the two
def impute_methods(w_na, wo_na):

    x = np.array(wo_na['day'])
    y = np.array(wo_na['ge'])
    x_test = np.array(w_na['day'])

    # Linear interpolation regression
    interp_func = interp1d(x, y, kind='linear', bounds_error=False, assume_sorted=False)

    # evaluate the interpolation function at new x-values
    interp_y = interp_func(x_test)
    
    # impute values to df
    interp_vals = {'day': x_test, 'ge':interp_y}
    
    # combine imputed df with df wo na's
    interp_df = pd.concat([pd.DataFrame(interp_vals), wo_na], axis=0).sort_values('day')
    
    
    
    # polynomial regression
    coeffs = np.polyfit(x, y, 3)

    # define polynomial function
    poly_func = np.poly1d(coeffs)

    # interpolate values for new x values
    poly_y = poly_func(x_test)
    
    # impute values to df
    poly_vals = {'day': x_test, 'ge':poly_y}
    
    # combine imputed df with df wo na's
    poly_df = pd.concat([pd.DataFrame(poly_vals), wo_na], axis=0).sort_values('day')
    
    
    # KNN regression
    # kneighbors only accepts 2d arrays, so may into array
    # then reshape it, -1 is 1 shape dimension, then 1 is 1 column
    x = np.array(wo_na['day']).reshape(-1, 1)
    y = np.array(wo_na['ge']).reshape(-1, 1)

    x_test = np.array(w_na['day']).reshape(-1, 1)

    
    # Nearest Neighbors regression
    # hyperparameters to search over
    param_grid = {'n_neighbors': [1, 3, 5, 7, 9, 10],
                 'p': [1,2]}

    # KNN regressor 
    knn = KNeighborsRegressor(n_jobs=2)

    # create a grid search object to search over the hyperparameters
    # by default, gridsearchcv tries to maximize whatever is inside of the scoring parameter
    # so by putting in the -mse, it will maximize it, essentially minimising the mse which is what we wants
    grid_search = GridSearchCV(knn, param_grid, cv=5, scoring='neg_mean_squared_error', n_jobs=2)

    # fit the grid search
    grid_search.fit(x, y)

    # append the best score
    knn_mse = -grid_search.best_score_

    # get the best hyperparameters from the grid search
    best_params = grid_search.best_params_

    # create a new KNN regressor with the best hyperparameters
    knn_func = KNeighborsRegressor(n_neighbors=best_params['n_neighbors'], 
                                   p=best_params['p'], n_jobs=2)
    knn_func.fit(x, y)

    # impute based on x test
    knn_y = knn_func.predict(x_test)
   
    # impute values to df, use ravel to unlist to be 1d array
    knn_vals = {'day': x_test.ravel(), 'ge':knn_y.ravel()}

    # combine imputed df with df wo na's
    knn_df = pd.concat([pd.DataFrame(knn_vals), wo_na], axis=0).sort_values('day') 
    
    
    return(poly_df, poly_func, interp_df, interp_func, knn_df, knn_func, knn_mse)