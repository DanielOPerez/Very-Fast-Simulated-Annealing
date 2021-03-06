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

module cost_fun_mod
export cost_fun

function cost_fun(var_in::Array,
                  args_in)

    x=var_in[1]
    return sum(100.0*(x[2:end]-x[1:end-1].^2.0).^2.0+(1.0.-x[1:end-1]).^2.0)

end

end
