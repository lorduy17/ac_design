import numpy as np
def coefficients(W_empty,W_gross):
    """
    Determine equation of empty weight allowed since a data base of similar
    A/C's related Gross Weight and empty weight. For a power law and regression
    lineal. The equation is determined by the least square method.
    Input:
    W_empty: asarray(dimensions=(n,),dtype=float),
        empty weight of similar A/C's in lb.
    W_gross: asarray(dimensions=(n,),dtype=float),
        gross weight of similar A/C's in lb.
    Output:
    A: float,
        coefficient of the power law.
    B: float,
        exponent of the power law.
    """
    # y = Ax^b
    y = W_empty
    x = W_gross
    if len(x) != len(y):
        raise ValueError('x and y must have the same length.')
    if isinstance(x, list) and isinstance(y, list):
        x = np.asarray(x)
        y = np.asarray(y)
    x_log = np.log(x)
    y_log = np.log(y)
    # square method Y_pred = b0 + b1*X
    x_mean = np.mean(x_log)
    y_mean = np.mean(y_log)
    B = np.sum((x_log-x_mean)*(y_log-y_mean))/np.sum((x_log-x_mean)**2)
    A = y_mean - B*x_mean
    return A,B
    


