#!/usr/bin/env python  
# -*- coding: utf-8 -*- 
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.
#
#----------------------------------------------------------------------------
# Created By  : Daniel O. Perez
# Created Date: November of 2019
# email: perez.daniel.omar@gmail.com
# ---------------------------------------------------------------------------

import numpy as np
from vfsa_mod import vfsa
from cost_fun_mod import cost_fun

# Parameter example
# var_op={'name': 'var1' <-- string (optional): name of the variable
#         'npar': 14,    <-- integer: number of parameters
#         'xa': 0        <-- float: lower limit
#         'xb': 100      <-- float: upper limit
#         'overlap': '   <-- string: if the variables can overlap each other. (y,n,yes,no,Yes,No,YES,NO)
#         'dx': 0.1      <-- float: minimun distance between variables, mandatory if overlap is "no"
#         'x0': x01}     <-- vector (optional): initial values of the variables. Must be of dimension npar.

def rosen(x_in,args):
    """The Rosenbrock function"""
    x=x_in[0]   
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)


var1_op={'name': 'x-coor',
         'npar': 5,
         'xa': -10,
         'xb': 10,
         'overlap': 'y'}

var_op_list=[var1_op]

result=vfsa(rosen,
            var_op=var_op_list,
            itmax=30000, #maximun number of iterations
            iop_prt=-100, 
            quench=1.0)

print(result.cost)
print(result.x[0])
