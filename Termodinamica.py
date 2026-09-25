import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Configuración global de la fuente
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Latin Modern Roman']
plt.rcParams['mathtext.fontset'] = 'cm'  # Fuente cm, es similar y serif
#plt.rcParams['font.size'] = 70  # Tamaño general para etiquetas y leyenda

# Cargar datos
data = np.loadtxt('masas.dat')
radio_masa = np.loadtxt('bulk_mass.dat')

# Extraer columnas
x_eos = data[:, 0]
y_eos = data[:, 1]

radio = radio_masa[:, 0]
masa = radio_masa[:, 1]

# Crear figura y subplots
fig, ax = plt.subplots(1,2, layout = 'constrained', figsize = (15, 6))

# Gráfico 1: Masa vs Radio (bulk_mass.dat)
ax[0].plot(radio, masa, color='#FF6666', label=r'$K(\rho)$', linewidth=2)
ax[0].axvline(0.16, ls='-.', lw=1, c='gray')
ax[0].axhline(0, ls='-.', lw=1, c='gray')

# Etiquetas y estilos
ax[0].set_xlabel(r"Densidad bariónica [fm$^{-3}$]", fontsize=25)
ax[0].set_ylabel(r"K [MeV]", fontsize=25)
ax[0].set_xlim(0,1.5)
#ax[0].set_ylim(-20,1800)
ax[0].text(0.07, 6000, r"(a)", fontsize=20, ha="center")
ax[0].tick_params(axis='both', which='major', labelsize=20)
ax[0].legend(loc='best', frameon=False, fontsize=18)

# Gráfico 2: Datos directos sin interpolar (masas.dat)
ax[1].plot(x_eos, y_eos, color='#FF6666', label=r'm$_\text{eff}/m$', linewidth=2, linestyle='-')
ax[1].axvline(0.16, ls='-.', lw=1, c='gray')

# Etiquetas y estilos
ax[1].set_xlabel(r'Densidad bariónica [fm$^{-3}$]', fontsize=25)  # Cambia el label según convenga
ax[1].set_ylabel(r'm$_\text{eff}/m$', fontsize=25)
ax[1].text(0.07, 0.83, r"(b)", fontsize=20, ha="center")
ax[1].set_xlim(0,1.5)
ax[1].tick_params(axis='both', which='major', labelsize=20)
ax[1].legend(loc='best', frameon=False, fontsize=18)

plt.savefig('Termodinamica_SORP.pdf')
plt.show()

plt.close()
