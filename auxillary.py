import itertools, typing
from copy import deepcopy
import scipy
import time
import numpy as np
                      
    
class ConvergenceError(Exception):
    "Raised when the input value is less than 18"
    pass
    
    
def Newton_solver(C: np.array, moments: np.array, conv_tolerance: float, max_iter: int=10**4, p0: np.ndarray = None, info: bool=False) -> typing.Tuple[np.ndarray, dict]:
    """
    
    """
    assert C.shape[0] == moments.shape[0], "number of moment/marginal constraints must equal number of rows of C"
    details = dict()  
    #
    tolerance = conv_tolerance * np.amax(np.abs(C), axis=1)
    if info:
        print("using tolerance of", np.amax(tolerance))
    
    dimS = C.shape[1]

    if p0 is None:
        p_old = np.ones(dimS) / dimS
    else:
        assert len(p0.shape) == 1 and p0.shape[0] == dimS, "base distribution ill-defined"
        p_old = p0.copy()                                                                   
    if info:
        print("starting ITEM from p =", p_old)
        details['error p'] = []
    
    t0 = time.time()
    for counter in range(1, max_iter+1):
        try:
            inv_grad = np.linalg.inv(p_old[np.newaxis] * C @ C.T)  
        except np.linalg.LinAlgError:
            if info:
                print("Jacobian not invertible or poorly nomerically stable")
                print("is C of full row-rank?", np.linalg.matrix_rank(C)==C.shape[0])
                print(" - marginals/moments estimate vs. empirical:", np.dot(C, p_old), "vs.", moments)
            raise ConvergenceError("Cannot invert Jacobian matrix at " + str(counter) + ". iteration of ITEM.")
        
        update = C.T @ inv_grad @ (C @ p_old - moments)
        p = p_old * np.exp(-update) 
        
        error_p = np.abs(p - p_old)                      
        if info:                          
            details['error p'].append(error_p.max())   

        if np.all(error_p < conv_tolerance):            
            moment_error = np.abs(C@p - moments)         
            if not np.all(moment_error < tolerance):
                if info: print("p has converged within specified tolerance, max|error p| = " + str(np.amax(error_p)) + ", but error on moments is still " + str(np.amax(moment_error)) + ", above tolerance of ", list(tolerance))
                details['status'] = 'converged, but not updating' 
            else:
                details['status'] = 'converged'
            #
            if info: print("--> convergence of Newton: max |C p-m| =", np.amax(moment_error), "after", counter, "iterations with max p-convergence error", np.amax(error_p))
                
            details['duration'] = round(time.time() - t0, 6)
            details['number of iterations'] = counter

            return p, details                            
        
        if np.isnan(p).any():
            if info:
                print("p =", p)
            raise ConvergenceError("ITEM de-railed at " + str(counter) + ". iteration. Potentially due to the presence zero-probability microstates")
        
        vanishing_p = p < conv_tolerance               
        if np.any(vanishing_p):
            p[vanishing_p] = 0.
    
        p_old = p
    
    raise RuntimeError("ITEM could not converge within set tolerance " + str(conv_tolerance) + " reaching max number of iterations = " + str(max_iter))