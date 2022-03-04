def rosen(x_in,args):
    """The Rosenbrock function"""
    x=x_in[0]   
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)
