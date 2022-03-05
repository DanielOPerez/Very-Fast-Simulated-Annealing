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

module vfsa_mod
import Random.rand
export vfsa
using Printf


struct Variable
    overlap::String
    npar::Int64
    xa::Float64
    xb::Float64
    dx::Float64
    tmp::Array
    opt::Array
    old::Array
end

function vfsa(func_in::Function,
              var_op_in::Tuple;
              arg_op_in::Tuple=(),
              itmax::Int64=10000,
              itmax_con::Int64=100,
              thres::Float64=0.1,
              xmisfit::Float64=1.0e-3,
              quench::Float64=1.0,
              tscale::Float64=1.0,
              trs::Float64=1.0e-5,
              tas::Int64=100,
              ifreq_re::Int64=0,
              ifreq_tscale::Int64=1,
              iop_prt::Int64=-100,
              fileout::String="out_vfsa90")

   
    
    mpar=length(var_op_in)
    #creo un arreglo var, de tipo Variable, de dimension mpar
    var=Array{Variable,1}(undef,mpar)
    #lleno el arreglo var objetos de la clase Variable
    icont=0
    mnpar=0
    for i in 1:mpar
        #si overlap no existe lo pongo "yes" por defecto
        if(~haskey(var_op_in[i],"overlap")); var_op_in[i]["overlap"]="yes" end
        #defino dx
        dx=0.0
        if((var_op_in[i]["overlap"]=="no") & haskey(var_op_in[i],"dx"))
           dx=var_op_in[i]["dx"]
        end        
        #defino el modelo inicial, si no existe
        if(haskey(var_op_in[i],"x0"))
            x0=var_op_in[i]["x0"]
        else
            x0=collect(range(var_op_in[i]["xa"],
                             stop=var_op_in[i]["xb"],
                             length=var_op_in[i]["npar"]+2))[2:end-1]

        end
       
        
        var[i]=Variable(var_op_in[i]["overlap"],
                        var_op_in[i]["npar"],
                        var_op_in[i]["xa"],
                        var_op_in[i]["xb"],
                        dx,
                        x0[1:end],  #hay que ponerlo asi para que inicialice
                        x0[1:end],  #y no que todos apunten a x0 
                        x0[1:end])

        mnpar=mnpar+var[i].npar
        if (var[i].xa==var[i].xb); icont+=1 end

    end



    if (quench>=1.0)
        qd=quench/(mnpar-icont)
        c=-log(trs)*exp(-log(tas)*qd)
        qd_cost=qd
        c_cost=c*tscale
    elseif ((quench>0.0) & (quench<1.0))
        qd=1.0/(mnpar-icont)        
        c=-log(trs)*exp(-log(tas)*qd)
        qd_cost=1.0
        c_cost=quench
    else 
        qd=1.0
        c=-quench
        qd_cost=qd
        c_cost=c*tscale
    end

    
    k=0
    acc_dn=0
    acc_up=0
    rej=0
    delta_cost=1.e150

    #cost of initial parameters is initalized
    cost_0=func_in([var[i].tmp for i in 1:mpar],arg_op_in)
   
    # ----------------------------------------------------
    # Variable init.
    # ----------------------------------------------------
    t_min=1.0e-150  #minimum temperature
    t0_cost=cost_0  #initial cost temperature
    cost_opt=cost_0 #initial optimum cost value
    t_cost=t0_cost  #cost temperature
    var_t0=1.0      #initial variables temperature
    var_t=var_t0    #variables temperature  
    # ----------------------------------------------------

    @printf("%i ( %i %i ) %f %e \n", k,acc_dn,acc_up,cost_opt,t_cost)
    #MAIN LOOP STARTS
    flag_cut=true
    while flag_cut

        #--------------------------------------------------------------------------------
        for i in 1:mpar
            #se guarda temporariamente el valor de 
            #la variable sin perturbar para poder
            #ser utilizada por el criterio de Metropolis.
            var[i].old[1:end]=var[i].tmp[1:end]
            #----------------------------------------------------------------------------
            # New config
            #----------------------------------------------------------------------------
            ##without overap
            if (var[i].overlap == "no")
                #if j=0
                var[i].tmp[1]=perturb(var[i].tmp[1],
                                      var_t,
                                      var[i].xa,
                                      var[i].tmp[2]-var[i].dx)[1]
                #if j from 1 to npar-1
                for j in 2:var[i].npar-1
                    var[i].tmp[j]=perturb(var[i].tmp[j],
                                          var_t,
                                          var[i].tmp[j-1]+var[i].dx,
                                          var[i].tmp[j+1]-var[i].dx)[1]
                end
                # if j is the last
                var[i].tmp[end]=perturb(var[i].tmp[end],
                                        var_t,
                                        var[i].tmp[end-1]+var[i].dx,
                                        var[i].xb)[1]
            else
                #With overlap, whole array is perturbed
                var[i].tmp[1:end]=perturb(var[i].tmp[1:end],
                                          var_t,
                                          var[i].xa,
                                          var[i].xb)[1:end]
            end
                        
            #----------------------------------------------------------------------------
            #----------------------------------------------------------------------------
        end #<- for i in 1:mpar
        #--------------------------------------------------------------------------------
       
        
        k+=1 #iteration
                
        ##new cost temperature.
        if (t_cost <= t_min)
            t_cost=t_min
        else
            t_cost=t0_cost*exp(-c_cost*((k-1)^qd_cost))
        end

        
        ##new parameter temperature. 
        if (var_t <= t_min)
            var_t=t_min
        else
            var_t=var_t0*exp(-c*((k-1)^qd))
        end

        cost_tmp=func_in([var[i].tmp for i in 1:mpar],arg_op_in)
        delta_cost=cost_tmp-cost_0
      
        
        #Metropolis
        if(delta_cost<0.0)
            acc_dn+=1
            cost_0=cost_tmp
            if(cost_tmp<cost_opt)
                cost_opt=cost_tmp
                for i in 1:mpar
                    var[i].opt[1:end]=var[i].tmp[1:end] #copiar por valor y no por referencia
                end
            end
        else
            prob=exp(-delta_cost/t_cost)
            b=rand()
            if(b < prob)
                acc_up+=1
                cost_0=cost_tmp
                if(cost_tmp<cost_opt)
                    cost_opt=cost_tmp
                    for i in 1:mpar
                        var[i].opt[1:end]=var[i].tmp[1:end] #copiar por valor y no por referencia
                    end
                end
            else
                for i in 1:mpar
                    var[i].tmp[1:end]=var[i].old[1:end]  #copiar por valor y no por referencia
                end
                rej+=1
            end
        end
        
        #if(mod(k,iop_prt)==0); println(k," ",acc_dn," ",acc_up," ",cost_opt," ",t_cost," ",var_t) end
        if((mod(k,iop_prt)==0) & (iop_prt<0))
            @printf("%i ( %i %i ) %f %e \n", k,acc_dn,acc_up,cost_opt,t_cost)
        end
        
        if(k>itmax); flag_cut=false end
        
    end #<-do while


        
end

function perturb(x_in,t_in,xa,xb)
    flag_cut_p=true
    while flag_cut_p
        u=rand(size(x_in,1))
        global perturbed=x_in.+t_in*(xb-xa)*sign.(u.-0.5).*((1.0+1.0/t_in).^abs.(2.0*u.-1.0).-1.0)
        flag_cut_p=(any(x -> (x<xa), perturbed) | any(x -> (x>xb), perturbed))
    end
    return perturbed
end


end


