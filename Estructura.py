import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

# Configuración global de la fuente
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Latin Modern Roman']
plt.rcParams['mathtext.fontset'] = 'cm'  # Fuente cm, es similar y serif
#plt.rcParams['font.size'] = 70  # Tamaño general para etiquetas y leyenda

#EoS-----------------------------------------------------------------
# Cargar los datos
data = np.loadtxt('EoS.dat')

# Crear funciones de interpolación cúbica
f_eos = interp1d(data[:, 1], data[:, 2], kind='cubic')

# Puntos nuevos para interpolación 
x_nuevos_eos = np.linspace(data[:, 1].min(), data[:, 1].max(), 500)

# Interpolaciones poblaciones
y_eos = f_eos(x_nuevos_eos)


#Masa-Radio-----------------------------------------------------------------
data = np.loadtxt('Radio-Masa.dat')

# Extraer columnas
radio = data[:, 0]
masa = data[:, 1]
#Figura
fig, ax = plt.subplots(1,2, layout = 'constrained', figsize = (15, 6))

ax[0].plot(x_nuevos_eos, y_eos, color='#FF6666', label=r'$\epsilon(p)$', linewidth=2)

ax[0].set_xlabel(r"Densidad de energía [MeV fm$^{-3}$]", fontsize=25)
ax[0].text(0,1200, r"(a)", fontsize=20, ha="center")
ax[0].set_ylabel(r"Presión [MeV fm$^{-3}$]", fontsize = 25)

ax[0].tick_params(axis='both', which = 'major', labelsize=20) 

ax[0].legend(loc = 'best', frameon = False, fontsize = 20)


# Graficar línea que conecta los puntos directamente
ax[1].plot(radio, masa, color='#FF6666', label=r'$M(r)$', linewidth=2)

# Etiquetas de ejes
ax[1].set_xlabel(r"Radio [km]", fontsize=25)
ax[1].text(8.3,1.4, r"(b)", fontsize=20, ha="center")
ax[1].set_ylabel(r"$M/M_\odot$", fontsize=25)
ax[1].set_xlim(7.7,20)
# Estilo de ticks
ax[1].tick_params(axis='both', which='major', labelsize=20)

# Leyenda
ax[1].legend(loc='best', frameon=False, fontsize=20)
plt.savefig('Estructura_SORP.pdf')
plt.show()

plt.close()
