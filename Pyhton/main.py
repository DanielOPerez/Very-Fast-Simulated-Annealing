import numpy as np
from vfsa_mod import vfsa
from cost_fun_mod import cost_fun
from rosenbrock_mod import rosen

var1_op={'name': 'coordenada-x',
         'npar': 5,
         'xa': -10,
         'xb': 10,
         'overlap': 'y'}

var_op_list=[var1_op]



result=vfsa(rosen,
            var_op=var_op_list,
            itmax=30000,
            iop_prt=-100,
            quench=1.0)

#print(result.cost)
#print(result.x[0])
