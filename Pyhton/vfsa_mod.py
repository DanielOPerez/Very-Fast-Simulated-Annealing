import numpy as np
import sys
import copy

def vfsa(cost_fun,
         var_op,
         args=(),
         itmax=1000,
         itmax_con=100,
         thres=0.1,
         xmisfit=1.0e-3,
         quench=1.0,
         tscale=1.0,
         trs=1.0e-5,
         tas=100,
         ifreq_re=0,
         ifreq_tscale=1,
         iop_prt=100,
         fileout='out_vfsa90'):

    class variable:
        ""
        def __init__(self,npar):
            self.name=''
            self.overlap=''
            self.npar=0
            self.xa=0.0
            self.xb=0.0
            self.dx=0.0
            self.tmp=np.zeros(npar)
            self.opt=np.zeros(npar)
            self.old=np.zeros(npar)

    class optimize_result:
        ""
        def __init__(self ):
            self.cost_opt=0
            self.x=[]

    def perturb(x_in,t_in,xa,xb,npar):
        "Perturb function"
        flag_cut=True
        while (flag_cut):
            u=np.random.uniform(size=npar)
            perturbed=x_in+np.sign(u-0.5)*t_in*((1.0+1.0/t_in)**np.abs(2.0*u-1.0)-1.0)*(xb-xa)
            flag_cut=(any(perturbed<xa) or any(perturbed>xb))
            
        return perturbed
    
    #==================================================================================
    # Check all variable options dictionaries and see if mandatory keys exists
    # Optional keys are completed if not exists
    #==================================================================================

    
    #Mandatory keys
    var_op_m_keys={'npar', 'xa', 'xb', 'overlap'} #esto es un set, no una lista

    #Overla values
    overlap_yes=('y','Y','yes','YES','Yes')
    overlap_no=('n','N','no','NO','No')
    
    #number of different types of variables
    mpar=len(var_op)

    #check all variables 
    for i in range(mpar): 

        #Optional 'name' key
        if 'name' not in var_op[i]:
            var_op[i].update(name='number '+str(i+1))
       
        #Check mandatory keys
        for key in var_op_m_keys:
            if key not in var_op[i]: 
                print("The {0} is not in defined in variable {1}. STOP."
                      .format(key,var_op[i].get('name')))
                sys.exit()

        #If npar is <=0 or not integer abort
        if var_op[i].get('npar') <= 0 or not isinstance(var_op[i].get('npar'),int):
            print("In variable {0} npar must be an integer greater than 0. STOP."
                  .format(var_op[i].get('name')))
            sys.exit()
                              
        #Check if 'overlap' value is correct            
        if var_op[i].get('overlap') not in overlap_yes+overlap_no:
            print("Invalid overlap value in variable {0}. STOP."
                  .format(var_op[i].get('name')))
            sys.exit()            
        #If overlap equal to 'n' check 'dx' 
        elif var_op[i].get('overlap') in overlap_no:
            if 'dx' not in var_op[i] or var_op[i].get('dx') < 0.0:
                print("The dx is not defined or is negative in variable {0}. STOP."
                      .format(var_op[i].get('name')))
                sys.exit()
                #var_op[i].update(dx=0.0)
        else: 
            var_op[i].update(dx=0.0)
                
        #If optional xma and xmb keys do not exists then a value is asigned
        # if 'xma' not in var_op[i]:
        #     var_op[i].update(xma=-1000*var_op[i].get('xa'))
        # if 'xmb' not in var_op[i]:
        #     var_op[i].update(xmb=1000*var_op[i].get('xb'))   
            

        #Check dimension of initial model, if  not exists is created
        if 'x0' in var_op[i]:
            if np.size(var_op[i].get('x0')) !=  var_op[i].get('npar'):
                print("Invalid dimension of initial model in variable {0}."
                      .format(var_op[i].get('name')))
                sys.exit()
        if 'x0' not in var_op[i]:
            var_op[i].update(x0=np.linspace(
                var_op[i].get('xa'),var_op[i].get('xb'),var_op[i].get('npar')+2)[1:-1])
  
    #==================================================================================
    #==================================================================================

    #The object var is created and initialized
    var = []
    for i in range(mpar):
        var.append(variable(var_op[i].get('npar')))
        var[i].name=var_op[i].get('name')
        var[i].overlap=var_op[i].get('overlap')
        var[i].npar=var_op[i].get('npar')
        var[i].xa=var_op[i].get('xa')
        var[i].xb=var_op[i].get('xb')
        var[i].dx=var_op[i].get('dx')
        var[i].tmp=var_op[i].get('x0')
        var[i].opt=copy.copy(var[i].tmp)
      
    
    mnpar=np.sum([var[i].npar for i in range(mpar)])
    icont=0
    for i in range(mpar):
        if(var[i].xa==var[i].xb):
            icont+=1
                        
    if(quench>=1.0):
        qd=quench/(mnpar-icont)
        c=-np.log(trs)*np.exp(-np.log(tas)*qd)
        qd_cost=qd
        c_cost=c*tscale
    elif(quench>0.0 and quench<1.0):
        qd=1.0/(mnpar-icont)
        c=-np.log(trs)*np.exp(-np.log(tas)*qd)
        qd_cost=1.0
        c_cost=quench
    else:
        qd=1.0
        c=-quench
        qd_cost=qd
        c_cost=c*tscale


    
    k=0
    acc_dn=0
    acc_up=0
    rej=0
    delta_cost=1.e150
            
    #cost of initial parameters is initalized
    cost_0=cost_fun([var[i].tmp for i in range(mpar)],args)
 
    #===================================================
    # Variable init.
    #===================================================
    t_min=1.0e-150  #minimum temperature
    t0_cost=cost_0  #initial cost temperature
    cost_opt=cost_0 #initial optimum cost value
    t_cost=t0_cost  #cost temperature
    var_t0=1.0      #initial variables temperature
    var_t=var_t0    #variables temperature  
    #===================================================
    print("{:7d}  ( {:4d}  {:4d} )  {:7f}   {:6e}"
          .format(k,acc_dn,acc_up,cost_opt,t_cost))


    
    #MAIN LOOP STARTS
    flag_cut=True
    while(flag_cut):

        
        #================================================================================
        for i in range(mpar):
            #se guarda temporariamente el valor de 
            #la variable sin perturbar para poder
            #ser utilizada por el criterio de Metropolis.
            var[i].old=list(var[i].tmp)
           
            #============================================================================
            # New config
            #============================================================================
            ##without overap
            if var[i].overlap in overlap_no:
                #if j=0
                var[i].tmp[0]=float(perturb(var[i].tmp[0],var_t,
                                            var[i].xa,var[i].tmp[1]-var[i].dx,1))
                #if j from 1 to npar-1
                for j in range(1,var[i].npar-1):
                    var[i].tmp[j]=float(perturb(var[i].tmp[j],var_t,
                                                var[i].tmp[j-1]+var[i].dx, 
                                                var[i].tmp[j+1]-var[i].dx,
                                                1))
                # if j is the last
                var[i].tmp[-1]=float(perturb(var[i].tmp[-1],var_t,
                                             var[i].tmp[-2]+var[i].dx,var[i].xb,1))  
            else:
                #With overlap, whole array is perturbed
                var[i].tmp=perturb(var[i].tmp,
                                   var_t,
                                   var[i].xa,
                                   var[i].xb,
                                   var[i].npar)
                
            #============================================================================
            #============================================================================
            
        

        #================================================================================

        k+=1 #iteration
        ##new cost temperature.
        
        if (t_cost <= t_min):
            t_cost=t_min
        else:
            t_cost=t0_cost*np.exp(-c_cost*((k-1)**qd_cost))
            
            
        ##new parameter temperature. 
        if (var_t <= t_min):
            var_t=t_min
        else:
            var_t=var_t0*np.exp(-c*((k-1)**qd))
            
       
        cost_tmp=cost_fun([var[i].tmp for i in range(mpar)],args)
        delta_cost=cost_tmp-cost_0
        
            
        
        #Metropolis
        if(delta_cost<0.0):
            acc_dn+=1
            cost_0=cost_tmp
            if(cost_tmp<cost_opt):
                cost_opt=cost_tmp
                for i in range(mpar):
                    var[i].opt=copy.copy(var[i].tmp) #copiar por valor y no por referencia
             
        else:
            prob=np.exp(-delta_cost/t_cost)
            b=np.random.uniform()
            if(b < prob):
                acc_up+=1
                cost_0=cost_tmp
                if(cost_tmp<cost_opt):
                    for i in range(mpar):
                        var[i].opt=copy.copy(var[i].tmp) #copiar por valor y no por referencia
                    cost_opt=cost_tmp
            else:
                for i in range(mpar):
                    var[i].tmp=copy.copy(var[i].old)
                rej+=1
       
        if(np.mod(k,iop_prt)==0 and iop_prt<0):
            print("{:7d}  ( {:4d}  {:4d} )  {:7f}   {:6e} {:6e}"
                  .format(k,acc_dn,acc_up,cost_opt,t_cost,var_t))
        if(k>itmax):flag_cut=False
    
        
    result=optimize_result()
    result.cost=cost_opt
    result.x=[var[i].opt for i in range(mpar)]
    
    return result

