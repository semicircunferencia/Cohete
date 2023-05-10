import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "DejaVu Sans"
#fig = plt.figure(figsize=(11.,8.))      # figura (ancho, largo)
ax = plt.subplot()      # subfigura

# Datos
data = np.loadtxt('H_modificado.dat')
tiempo = np.linspace[0,30,6e-6] # cojo el tiempo de la primera columna
H = data

# configurar ejes
ax.set_ylabel('H modificado (meses$^-2$)', fontname='DejaVu Sans', fontsize='12')
ax.set_xlabel('Tiempo (días)', fontname='DejaVu Sans', fontsize='12')

#Cambiar ticks
#for label in ax.get_xticklabels():
    #label.set_fontproperties('Times New Roman')
#plt.xticks(fontsize='15')
#for label in ax.get_yticklabels():
    #label.set_fontproperties('Times New Roman')
#plt.yticks(fontsize='15')

# Creación de la gráfica
ax.plot(tiempo, H, linestyle='-', marker='', markersize=4, color='#B4045F', label="H modificado", linewidth=1.0)  #marker=puntos


# Guardar la gráfica
plt.savefig('grafica_energia.pdf',dpi=300, bbox_inches = "tight")

# Mostrarla en pantalla
#plt.show()