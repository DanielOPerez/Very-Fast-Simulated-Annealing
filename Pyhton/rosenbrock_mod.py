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

def rosen(x_in,args):
    """The Rosenbrock function"""
    x=x_in[0]   
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)
