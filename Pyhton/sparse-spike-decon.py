import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from vfsa_mod import vfsa
from funcion_costo_mod import funcion_costo
import time


#Genero el dato observado
ns=250
dt=0.004
#reflectividad
ref=np.zeros(ns)
spikes_pos=[20,50,80,120,200]
spikes_amp=[0.1,-0.1,-0.2,0.3,-0.3]
for i in range(np.size(spikes_pos)):
    ref[spikes_pos[i]]=spikes_amp[i]
#ondicula
wav=signal.ricker(101, 5)
#traza
traza=np.convolve(wav,ref,mode='same')
#traza con ruido
snr=10
traza_n=traza+np.random.normal(0,
                               np.max(np.abs(traza))/snr,
                               size=ns)

#===========================================================
# se usaran dos grupos de parametros ambos se pueden cruzar
# voy a buscar 5 spikes entre la muestra 0 y la 250
# con amplitudes entre -0.5 y 0.5

nspikes=5

#ejemplo de parametros
#var_op={'name': 'var1' <-- string, opcional, nombre de la variable (ni se para que lo puse)
#        'npar': 14,    <-- entero, obligatorio, numero de parametros
#        'xa': 0        <-- real, obligatorio, limite inferior
#        'xb': 100      <-- real, obligatorio, limite superior
#        'overlap': '   <-- stringr, obligatorio (y,n,yes,no,Yes,No,YES,NO)
#        'dx': 0.1      <-- real, obligatorio si overlap es no
#        'x0': x01}     <-- arreglo, opcional: modelo inicial de esa variable


#opciones de los parametros de la posicion (diccionario)
var1_op={'name': 'posicion spikes',
         'npar': nspikes,
         'xa': 0,
         'xb': ns,
         'overlap': 'y'}

#opciones de los parametros de la amplitud (diccionario)
var2_op={'npar': nspikes,
         'xa': -0.5,
         'xb': 0.5,
         'overlap': 'y'}

#datos que necesita la funcion de costo (lista)
args_in=[traza_n,wav,ns]

start_time = time.time()

#se llama al vfsa
resultado=vfsa(funcion_costo,       #obligatorio: funcion a minimizar
               [var1_op,var2_op],   #obligatorio: lista con las op de los param.
               args=args_in,        #opcional, depende de la funcion a minimizar
               itmax=10000,         #opcional (por defecto 1000)
               iop_prt=100)        #opcional (por defecto 100)

#Ver vfsa_mod.py para otras opciones del vfsa, todavia no funcionan todas
#La funcion a minimizar tiene que ser de la forma f(var_in,args_in), donde
#var_in y args_in van a ser listas. 

# el resultado vuelve en una estructura:
#resulado.cost: valor optimo de la fucion costo
#resultado.x arreglo de las variables optimas

#===========================================================

print("--- %s seconds ---" % (time.time() - start_time))


# print(resultado.cost)
# print('Posiciones',resultado.x[0])
# print('Amplitudes',resultado.x[1])

# #grafico

# ref_out=np.zeros(ns)
# for i in range(np.size(resultado.x[0])):
#     ref_out[np.int(resultado.x[0][i])]=resultado.x[1][i]

# traza_calc=np.convolve(wav,ref_out,mode='same')


# plt.figure(1)
# plt.plot(range(0,ns),traza_n,range(0,ns),traza_calc)
# plt.show()

# plt.figure(2)
# plt.plot(range(0,ns),ref,range(0,ns),ref_out)
# plt.show()
