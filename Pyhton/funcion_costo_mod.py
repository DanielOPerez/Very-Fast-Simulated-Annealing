import numpy as np

def funcion_costo(var_in,args_in):

    pos_spikes=var_in[0]
    amp_spikes=var_in[1]
    
    traza_obs=args_in[0]
    wav=args_in[1]
    ns=args_in[2]

    #ojo hay que pasar los pos_spikes a enteros
    #para eso uso np.int 
    ref=np.zeros(ns)
    for i in range(np.size(pos_spikes)):
        ref[np.int(pos_spikes[i])]=amp_spikes[i]

    traza_calc=np.convolve(wav,ref,mode='same')
    
    cost=np.sqrt(np.sum(np.power(traza_calc-traza_obs,2)))

    return cost
