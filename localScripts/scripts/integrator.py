discard = 1e-8

def integrate(func,X,Y,nMin,nMax,errorLevel):
    steps = 1 << (nMin - 1)
    error = 0
    result = 0
    result,error,steps = do_while_func(func,X,Y,steps,result)
    nbsteps = 0
    while(abs(error) > errorLevel*abs(result) and steps < (1 << nMax) and abs(result) > discard):
        result,error,steps = do_while_func(func,X,Y,steps,result)
    return result,error
          
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


from scipy.integrate import  dblquad

def easy_integrate(func,X,Y):
    return dblquad(func, -X, X, lambda x: -Y, lambda x: Y)