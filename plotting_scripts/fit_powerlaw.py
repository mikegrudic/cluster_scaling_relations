import emcee
from multiprocessing import Pool
from scipy.optimize import minimize
import numpy as np
from scipy.stats import theilslopes
from numba import jit

def log_likelihood(params,r,errors,*indep_vars,priors=True):
    if errors is None: errors = np.zeros_like(r)
    coeffs = params[1:-1] # slopes/power-law indices
    scatter = params[-1] 
    if priors:
        if scatter <= 0 or scatter>10: return -np.inf # prior on size of scatter in dex
        if np.any(np.abs(coeffs) > 10): return -np.inf # prior on slopes (nothing steeper than 10)
        if np.any(np.abs(params[0]) > 100): return -np.inf
    norm = params[0] # normalization of relation
    r_med = norm + np.inner(coeffs, np.c_[indep_vars]) # median value predicted by relation
    dlogr = r - r_med # differences from median - should be log-normal
    lnprob = np.log( 1/np.sqrt(2*np.pi*(errors**2 + scatter**2)) * np.exp(-dlogr**2 / (2 * (errors**2 +scatter**2)))) # likelihoods of deviations, given the scatter
    #lnprob += np.log(weights)
    if np.any(np.invert(np.isfinite(lnprob))): return -np.inf
    if np.any(np.isnan(lnprob)): return -np.inf
    return np.sum((lnprob)[np.isfinite(lnprob)])

def log_likelihood_variable_calibration(params,r,errors,calibration_sigma,*indep_vars,priors=True):
    if errors is None: errors = np.zeros_like(r)
    calibration = params[-2]
    coeffs = params[1:-2] # slopes/power-law indices
    scatter = params[-1] 
    if priors:
        if scatter <= 0 or scatter>10: return -np.inf # prior on size of scatter in dex
        if np.any(np.abs(coeffs) > 10): return -np.inf # prior on slopes (nothing steeper than 10)
        if np.any(np.abs(params[0]) > 100): return -np.inf
    norm = params[0] # normalization of relation
   # print(norm.shape,coeffs.shape,np.c_[indep_vars].shape)
    r_med = calibration * (norm + np.inner(coeffs, np.c_[indep_vars])) # median value predicted by relation
    dlogr = r - r_med # differences from median - should be log-normal
    lnprob = np.log( 1/np.sqrt(2*np.pi*(errors**2 + scatter**2)) * np.exp(-dlogr**2 / (2 * (errors**2 + scatter**2)))) # likelihoods of deviations, given the scatter
    lnprob += np.log(1/np.sqrt(2*np.pi*calibration_sigma**2) * np.exp(-calibration**2/(2*calibration_sigma**2))) # lognormal calibration to allow for variable slope
    #lnprob += np.log(weights)
    if np.any(np.invert(np.isfinite(lnprob))): return -np.inf
    if np.any(np.isnan(lnprob)): return -np.inf
    return np.sum((lnprob)[np.isfinite(lnprob)])


def multilinear_fit(dep_var, *indep_vars):
    # return the generalization of m and b in y = m x + b for the optimal multilinear fit of dep_var to indep_vars
    X = np.c_[indep_vars]
    y = np.copy(dep_var)
    finite = np.all(np.isfinite(X),axis=1)*np.isfinite(y)
    X, y = X[finite], y[finite]
   # print(X.shape, y.shape)
    X_mean = np.mean(X,axis=0)
    y_mean = np.mean(dep_var)
    y -= y_mean
    X -= X_mean
    coeffs, rms_error = np.linalg.lstsq(X,y,rcond=None)[:2]
    #y_pred = np.inner(coeffs, X)
    return coeffs, y_mean - np.inner(coeffs, X_mean), np.sqrt(rms_error/len(y))

def fit_powerlaw_model_with_scatter(dep_var,*indep_vars,errors=None,  chainlength=10000, return_likelihood=False, return_DIC=False,calibration_sigma=0):
    if errors is None: errors = np.zeros_like(dep_var)
    ndim = 2 + len(indep_vars) # normalization + slope for each variable + scatter
    if calibration_sigma != 0: ndim += 1
    nwalkers = 20
    
    coeffs, b, rms_error = multilinear_fit(dep_var, *indep_vars)    
    ypred = np.inner(coeffs,np.c_[indep_vars]) + b
    scatter_guess = np.std(dep_var-ypred)
    if(np.abs(b)>10): 
        b = 0; coeffs[:]=0    
    
    if calibration_sigma==0:
        p0 = np.array([b, *coeffs, scatter_guess]) + 1e-2 * np.random.normal(size=(nwalkers,ndim))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood,args=(dep_var,errors,*indep_vars))
    else:
        p0 = np.array([b, calibration_sigma, *coeffs, scatter_guess]) + 1e-2 * np.random.normal(size=(nwalkers,ndim))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood_variable_calibration,args=(dep_var,errors,calibration_sigma,*indep_vars))
    burnin = max(200,chainlength//10)
    sampler.run_mcmc(p0, chainlength,progress=False)
    samples = sampler.get_chain(discard=burnin,flat=True,thin=100)

    if not return_likelihood or return_DIC: return samples
    stuff_toreturn = [samples,]
    if return_likelihood:
        med = samples[sampler.get_log_prob(discard=burnin,flat=True,thin=100).argmax()]
        if calibration_sigma==0:
            func =  lambda x: -log_likelihood(x,dep_var,errors,*indep_vars,priors=False) 
        else:
            func = lambda x: -log_likelihood_variable_calibration(x,dep_var,errors,calibration_sigma,*indep_vars,priors=False)
        sol = minimize(func,med,tol=1e-3)
        stuff_toreturn.append(-sol.fun)
        #return samples, -sol.fun #sampler.get_log_prob(discard=burnin,flat=True,thin=100).max()
    if return_DIC:
        logprobs = sampler.get_log_prob(discard=burnin,flat=True,thin=100)
        if calibration_sigma==0:
            DIC = 2*np.average(-2*logprobs) + 2*log_likelihood(np.mean(samples,axis=0),dep_var,errors,*indep_vars)
        else:
            DIC = 2*np.average(-2*logprobs) + 2*log_likelihood_variable_calibration(np.mean(samples,axis=0),dep_var,errors,calibration_sigma,*indep_vars)
        stuff_toreturn.append(DIC)
    return stuff_toreturn



def log_likelihood_variable_scatter(params,r,time,*indep_vars,errors=None):
    if errors is None: errors = np.ones_like(r)
    coeffs = params[1:-2] # slopes/power-law indices
    scatter_intercept, scatter_slope = params[-2], params[-1]
    scatter = scatter_intercept + scatter_slope * time
    if np.any(scatter <= 0) or np.any(scatter>10): return -np.inf # prior on size of scatter in dex
    if np.any(np.abs(coeffs) > 10): return -np.inf # prior on slopes (nothing steeper than 10)
    if np.any(np.abs(params[0]) > 10): return -np.inf
    norm = params[0] # normalization of relation
    #print(coeffs.shape, np.c_[indep_vars].shape)
    r_med = norm + np.inner(coeffs, np.c_[indep_vars]) # median value predicted by relation
    dlogr = r - r_med # differences from median - should be log-normal
    lnprob = np.log( 1/np.sqrt(2*np.pi*(errors**2 + scatter**2)) * np.exp(-dlogr**2 / (2 * (errors**2 + scatter**2)))) # likelihoods of deviations, given the scatter
    if np.any(np.invert(np.isfinite(lnprob))): return -np.inf
    if np.any(np.isnan(lnprob)): return -np.inf
    return np.sum(lnprob)

def fit_powerlaw_model_with_variable_scatter(dep_var,time,*indep_vars,chainlength=10000,errors=None):
    if errors is None: errors = np.zeros_like(dep_var)
    ndim = 3 + len(indep_vars) # normalization + slope for each variable + scatter slope + scatter intercept
    nwalkers = 40
    
    coeffs, b, rms_error = multilinear_fit(dep_var, *indep_vars)    
    ypred = np.inner(coeffs,np.c_[indep_vars]) + b
    scatter_intercept_guess = np.std(dep_var-ypred)
    scatter_slope_guess = 0.
    if(np.abs(b)>10): 
        b = 0; coeffs[:]=0

    p0 = np.array([b, *coeffs, scatter_intercept_guess, scatter_slope_guess]) + 1e-2 * np.random.normal(size=(nwalkers,ndim))
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood_variable_scatter,args=(dep_var,time, errors, *indep_vars))
    burnin = max(200,chainlength//10)
    sampler.run_mcmc(p0, chainlength,progress=False)
    samples = sampler.get_chain(discard=burnin,flat=True,thin=100)
    return samples

@jit
def log_likelihood_orth(params,x,y,dx,dy):
    theta,bperp,scatter = params
    if scatter < 0: return -np.inf
    if theta < -np.pi/2: return -np.inf
    if theta > np.pi/2: return -np.inf

    slope = np.tan(theta)
    intercept = bperp / np.cos(theta)
    scatter = 10**logscatter
    sigma_i = np.sqrt((slope**2 * dx**2 + dy**2)/(1+slope**2))
    sigma_eff2 = sigma_i**2 + scatter**2
    delta_i = (y - (slope*x + intercept))/np.sqrt(1+slope**2)
    logL = np.log((2*np.pi*sigma_eff2)**-0.5 * np.exp(-delta_i**2/(2*sigma_eff2))) #-0.5*np.log(2*np.pi*(sigma_i**2 + scatter_orth**2)) - delta_i**2/(2*(sigma_i**2 + scatter_orth**2))
    if not np.all(np.isfinite(logL)): return  -np.inf
    return logL.sum()


def orthogonal_linear_regression(x,y,dx=None,dy=None,chainlength=1000):
    """Hogg 2010 orthogonal linear regression assuming uncorrelated errors in x and y variables; returns MCMC samples"""
    if dx is None: dx = np.zeros_like(x)
    if dy is None: dy = np.zeros_like(y)

    ndim, nwalkers = 3, 100
    slope, intercept = theilslopes(y,x)[:2]
    logscatter = np.log10(np.std(x*slope + intercept - y))

    theta = np.arctan(slope)
    bperp = np.cos(theta) * intercept

    p0 = np.array([theta,bperp,logscatter]) + 1e-2*np.random.normal(size=(nwalkers,ndim))
    p0[:,0] = (p0[:,0]+np.pi/2)%np.pi - np.pi/2
    burnin = max(200,chainlength//10)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood_orth,args=(x,y,dx,dy))
    sampler.run_mcmc(p0, chainlength,progress=False)
    samples = sampler.get_chain(discard=burnin,flat=True,thin=100)
    slope = np.tan(samples[:,0])
    intercept = samples[:,1] / np.cos(samples[:,0])
    scatter = 10**samples[:,2] * np.sqrt(1+slope**2)
#    samples[:,2] = 10**samples[:,2]
    return np.c_[slope,intercept,scatter]