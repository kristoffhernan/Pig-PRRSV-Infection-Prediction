import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import splev, splrep


def na_split(gene, base_df):
    eg = base_df[[gene, 'dpi']]
    eg = eg.reset_index()

    # df with x as days and y and individual
    df = eg.pivot(index='index', columns='dpi')[gene].reset_index(drop=True).drop(11, axis=1)

    # df with day and ge columns
    melt = pd.melt(df, value_vars=df.columns, var_name='day', value_name='ge')

    # splitting into two df, with and without NA values
    w_na = melt[melt.isna().any(axis=1)].sort_values('day')
    wo_na = melt[melt.notna().all(axis=1)].sort_values('day')
    
    return(w_na, wo_na)



def impute_methods(w_na, wo_na):

    x = np.array(wo_na['day'])
    y = np.array(wo_na['ge'])
    x_test = np.array(w_na['day'])

    # interpolation regression
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
    
    
    return(poly_df, poly_y, interp_df, interp_y)