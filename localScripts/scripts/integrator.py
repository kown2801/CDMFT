import numpy as np
discard = 1e-8

#func is the 2D function we wish to integrate
#x and Y are the upper limits of the integral in the x and y directions (the lower bounds are -X and -Y)
#nMin, nMax and errorLevel are parameters of the EulerMaclaurin2D integration
def integrate(func,X,Y,nMin,nMax,errorLevel):
    steps = 1 << (nMin - 1)
    error = 0
    result = 0
    result,error,steps = do_while_func(func,X,Y,steps,result)
    nbsteps = 0
    while(np.sum(abs(error)) > errorLevel*np.sum(abs(result)) and steps < (1 << nMax) and np.sum(abs(result)) > discard):
        result,error,steps = do_while_func(func,X,Y,steps,result)
    return result,error

def do_while_func(func,X,Y,steps,result):
    steps*=2
    error = result
    #First we compute the function on the 2D grid
    kx = -X*(2*np.arange(steps) - steps + 1)/(steps-1)
    ky = -Y*(2*np.arange(steps) - steps + 1)/(steps-1)
    KY,KX = np.meshgrid(kx,ky)
    this_call = func(KX,KY)
    #We create the integration mask for the EulerMaclaurin2D integration
    fact_matrix = np.zeros(this_call.shape) + 1
    fact_matrix[:,0]/=2
    fact_matrix[:,-1]/=2
    fact_matrix[0,:]/=2
    fact_matrix[-1,:]/=2
    #We compute the result
    result = this_call*fact_matrix/((steps - 2.)*(steps - 2.) + 2.*(steps - 2.) + 1.)
    result = np.sum(result,axis=(0,1))
    #Compute the error to the last iteration
    error -= result
    return result,error,steps

#Old integration function that worked but used too much for loops. The new one takes advantage of speedups coming from numpy array operations.
"""       
def do_while_func(func,X,Y,steps,result):
    steps*=2
    error = result
    result = 0
    for i in range(steps):
        for j in range(steps):
            kx = -X*(2.*i - steps + 1.)/(steps - 1)
            ky = -Y*(2.*j - steps + 1.)/(steps - 1)
            fact = 1
            if(i%(steps - 1) == 0):
                fact/=2
            if(j%(steps - 1) == 0):
                fact/=2
            result += fact*func(kx,ky)
    result *= 1./((steps - 2.)*(steps - 2.) + 2.*(steps - 2.) + 1.)
    error-= result
    return result,error,steps
"""