import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

# Configuración global de la fuente
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Latin Modern Roman']
plt.rcParams['mathtext.fontset'] = 'cm'  # Fuente cm, es similar y serif

# POBLACIONES -----------------------------------------------------------------
# Cargar los datos
data_pob = np.loadtxt('poblaciones.dat')

# Igualar p y e⁻ para x < 0.13
mascara_pob = data_pob[:, 0] < 0.13
data_pob[mascara_pob, 1] = data_pob[mascara_pob, 3]

# Crear funciones de interpolación cúbica
f_poblaciones1 = interp1d(data_pob[:, 0], data_pob[:, 1], kind='linear')
f_poblaciones2 = interp1d(data_pob[:, 0], data_pob[:, 2], kind='linear')
f_poblaciones3 = interp1d(data_pob[:, 0], data_pob[:, 3], kind='linear')
f_poblaciones4 = interp1d(data_pob[:, 0], data_pob[:, 4], kind='linear')
f_poblaciones5 = interp1d(data_pob[:, 0], data_pob[:, 5], kind='linear')

# Puntos nuevos para interpolación 
x_nuevos_poblaciones = np.linspace(data_pob[:, 0].min(), data_pob[:, 0].max(), 500)

# Interpolaciones poblaciones
y_poblaciones_1 = f_poblaciones1(x_nuevos_poblaciones)
y_poblaciones_2 = f_poblaciones2(x_nuevos_poblaciones)
y_poblaciones_3 = f_poblaciones3(x_nuevos_poblaciones)
y_poblaciones_4 = f_poblaciones4(x_nuevos_poblaciones)
y_poblaciones_5 = f_poblaciones5(x_nuevos_poblaciones)

# CAMPOS MEDIOS ---------------------------------------------------------------
# Cargar los datos
data_mom = np.loadtxt('momentos.dat')

# Igualar rho y pi⁻ para x < 0.13
mascara_mom = data_mom[:, 0] < 0.13
data_mom[mascara_mom, 8] = data_mom[mascara_mom, 6]

# Crear funciones de interpolación cúbica
f_cmedios_1 = interp1d(data_mom[:, 0], data_mom[:, 4], kind='cubic')
f_cmedios_2 = interp1d(data_mom[:, 0], data_mom[:, 5], kind='cubic')
f_cmedios_3 = interp1d(data_mom[:, 0], data_mom[:, 6], kind='cubic')
f_cmedios_4 = interp1d(data_mom[:, 0], data_mom[:, 8], kind='cubic')

# Puntos nuevos para interpolación 
x_nuevos = np.linspace(data_mom[:, 0].min(), data_mom[:, 0].max(), 500)

# Interpolaciones campos medios
y_cmedios_1 = f_cmedios_1(x_nuevos)
y_cmedios_2 = f_cmedios_2(x_nuevos)
y_cmedios_3 = f_cmedios_3(x_nuevos)
y_cmedios_4 = f_cmedios_4(x_nuevos)

# GRAFICADO -------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, layout='constrained', figsize=(15, 6))

# Poblaciones
ax[0].plot(x_nuevos_poblaciones, y_poblaciones_1, color='#FF6666', label=r'$p$', linewidth=1.5)
ax[0].plot(x_nuevos_poblaciones, y_poblaciones_2, color="#009235", label=r'$n$', linewidth=1.5)
ax[0].plot(x_nuevos_poblaciones, y_poblaciones_3, color='#6666FF', label=r'$e^{-}$', linewidth=1.5)
ax[0].plot(x_nuevos_poblaciones, y_poblaciones_4, color='#66CCFF', label=r'${\mu}^{-}$', linewidth=1.5)
ax[0].axhline(1, ls='-.', lw=1, c='gray')

ax[0].set_xlim(0, 1.5)
ax[0].set_ylim(0.0001, 3)
ax[0].set_yscale('log')

ax[0].set_xlabel(r"Densidad bariónica [fm$^{-3}$]", fontsize=25)
ax[0].text(0.1, 1.8, r"(a)", fontsize=20, ha="center")
ax[0].set_ylabel(r"Poblaciones relativas", fontsize=25)

ax[0].tick_params(axis='both', which='major', labelsize=20)
ax[0].legend(loc='best', frameon=False, fontsize=20)

# Campos medios
ax[1].plot(x_nuevos, y_cmedios_1, color='#FF6666', label=r'$\langle \sigma \rangle$', linewidth=1.5)
ax[1].plot(x_nuevos, y_cmedios_2, color="#009235", label=r'$\langle \omega \rangle$', linewidth=1.5)
ax[1].plot(x_nuevos, y_cmedios_3, color='#6666FF', label=r'$\langle \rho \rangle$', linewidth=1.5)
ax[1].plot(x_nuevos, y_cmedios_4, color='#66CCFF', label=r'$\langle{\pi}^{-}\rangle$', linewidth=1.5)
ax[1].axhline(0, ls='-.', lw=1, c='gray')

ax[1].set_xlim(0, 1.5)
#ax[1].set_ylim(-2, 1)

ax[1].set_xlabel(r"Densidad bariónica [fm$^{-3}$]", fontsize=25)
ax[1].text(0.1, 0.75, r"(b)", fontsize=20, ha="center")
ax[1].set_ylabel(r"Campos Medios [MeV]", fontsize=25)

ax[1].tick_params(axis='both', which='major', labelsize=20)
ax[1].legend(loc='best', frameon=False, fontsize=20)

plt.savefig('Sistema_SORP.pdf')
plt.show()
plt.close(fig)