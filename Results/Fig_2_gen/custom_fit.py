import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t

def jackknife(xa, med):
    # Now, should actually jackknife across anneal times.
    polyjack = []
    for i in range(len(xa)):
        inds = list(range(len(xa)))
        inds.remove(i)
        polyjack.append(np.polyfit(np.array(xa)[inds], np.array(med)[inds], 1))

    # Length of half error bar from jackknifing
    slopes = np.array(polyjack)[:, 0]
    jk_ci = np.std(slopes) * np.sqrt(len(xa) - 1) * 2
    ret = (np.mean(slopes), jk_ci)
    print(f'Jackknife slope is ({ret[0]:.3f},{ret[1]:.3f})')
    return np.mean(slopes), jk_ci


def custom_fit(x, y, weights=None):

    weights = np.ones(len(x))
    # Define the linear model
    def linear_model(x, a, b):
        return a * x + b

    # Function to compute the Jacobian matrix for a linear model
    def compute_jacobian(f, xdata, popt):
        """ Compute the Jacobian matrix for the fitted function f at xdata and popt. """
        eps = np.sqrt(np.finfo(float).eps)
        jacobian = np.zeros((len(xdata), len(popt)))
        for i, p in enumerate(popt):
            p_eps = np.copy(popt)
            p_eps[i] += eps
            jacobian[:, i] = (f(xdata, *p_eps) - f(xdata, *popt)) / eps
        return jacobian

    # Initial guess for the parameters
    initial_guess = [1, 0]

    # Perform the fit
    popt, _ = curve_fit(linear_model, x, y, sigma=1/np.sqrt(weights), p0=initial_guess, absolute_sigma=True)

    # Compute the Jacobian matrix at the fitted parameters
    jacobian = compute_jacobian(linear_model, x, popt)

    # Calculate the covariance matrix from the Jacobian
    residuals = y - linear_model(x, *popt)
    residual_variance = np.var(residuals, ddof=len(popt))
    covariance_matrix = residual_variance * np.linalg.inv(np.dot(jacobian.T, jacobian))

    # Calculate the degrees of freedom
    dof = max(0, len(x) - len(popt))

    # Calculate the standard errors from the covariance matrix
    standard_errors = np.sqrt(np.diag(covariance_matrix))

    # Calculate the t-value for 95% confidence interval
    confidence_level = 0.95
    alpha = 1.0 - confidence_level
    t_val = t.ppf(1.0 - alpha / 2.0, dof)

    # Compute the confidence intervals for each parameter
    ci_from_jacobian = np.zeros((len(popt), 2))
    for i, p in enumerate(popt):
        ci_from_jacobian[i, 0] = p - t_val * standard_errors[i]
        ci_from_jacobian[i, 1] = p + t_val * standard_errors[i]


    # alternatively, do jackknifing
    ci_from_jk = np.zeros((2, 2))
    popt_jk, ci_from_jkn = jackknife(x, y)
    ci_from_jk[0, 0] = popt_jk - ci_from_jkn
    ci_from_jk[0, 1] = popt_jk + ci_from_jkn
    popt_jk = np.array([popt_jk, 0])

    return popt, ci_from_jacobian, popt_jk, ci_from_jk

# Print the results
#print('Fitted Parameters:', popt)
#print('Standard Errors (from Jacobian):', standard_errors)
#print(f'95% Confidence Intervals (from Jacobian):\n', ci_from_jacobian)
