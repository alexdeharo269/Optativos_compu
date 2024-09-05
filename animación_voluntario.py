# Importa los módulos necesarios
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
import numpy as np

# Parámetros
# ========================================
file_in = r"C:\Users\Ale\Desktop\UGR FISICA\FISICA_COMPUTACIONAL\Alejandro de Haro García\Planetas\Animacion\planets_data_seed_1567234587_h_0.5_vr_0.050_radii_10000.0.dat" # Nombre del fichero de datos
file_out = "planetas" # Nombre del fichero de salida (sin extensión)

# Límites de los ejes X e Y
x_min = -50
x_max = 50
y_min = -50 
y_max = 50

interval = 30 # Tiempo entre fotogramas en milisegundos
show_trail = False # Muestra la "estela" del planeta
trail_width = 0.2 # Ancho de la estela
save_to_file = True # False: muestra la animación por pantalla,
                     # True: la guarda en un fichero 
dpi = 150 # Calidad del vídeo de salida (dots per inch)

# Radio del planeta, en las mismas unidades que la posición
base_radius =1.812/5  # Base radius as specified

# Lectura del fichero de datos
# ========================================
# Lee el fichero a una cadena de texto
with open(file_in, "r") as f:
    data_str = f.read()

# Inicializa la lista con los datos de cada fotograma.
# frames_data[j] contiene los datos del fotograma j-ésimo
frames_data = list()

# Itera sobre los bloques de texto separados por líneas vacías
# (cada bloque corresponde a un instante de tiempo)
for frame_data_str in data_str.strip().split("\n\n"):
    # Inicializa la lista con la posición de cada planeta
    frame_data = list()

    # Itera sobre las líneas del bloque
    # (cada línea da la posición de un planeta)
    for planet_pos_str in frame_data_str.strip().split("\n"):
        # Lee los datos de la línea
        planet_data = np.fromstring(planet_pos_str, sep=",")
        # Si la línea no está vacía y el planeta es real, añade planet_data a la lista de posiciones del fotograma
        if planet_data.size > 0 and planet_data[0] == 1:
            frame_data.append(planet_data[1:])

    # Añade los datos de este fotograma a la lista
    frames_data.append(frame_data)

# El número de planetas es el número de líneas en cada bloque
# Lo calculamos del primer bloque
nplanets = len(frames_data[0])


# Lectura del fichero de datos
# ========================================
# Lee el fichero a una cadena de texto
with open(file_in, "r") as f:
    data_str = f.read()

# Inicializa la lista con los datos de cada fotograma.
# frames_data[j] contiene los datos del fotograma j-ésimo
frames_data = list()

# Itera sobre los bloques de texto separados por líneas vacías
# (cada bloque corresponde a un instante de tiempo)
for frame_data_str in data_str.strip().split("\n\n"):
    # Inicializa la lista con la posición de cada planeta
    frame_data = list()

    # Itera sobre las líneas del bloque
    # (cada línea da la posición de un planeta)
    for planet_pos_str in frame_data_str.strip().split("\n"):
        # Lee los datos de la línea
        planet_data = np.fromstring(planet_pos_str, sep=",")
        # Si la línea no está vacía, añade planet_data a la lista de posiciones del fotograma
        if planet_data.size > 0:
            frame_data.append(planet_data)

    # Añade los datos de este fotograma a la lista si hay planetas reales
    if frame_data:
        frames_data.append(frame_data)

# El número de planetas es el número de líneas en cada bloque
# Lo calculamos del primer bloque
nplanets = len(frames_data[0])

# Creación de la animación/gráfico
# ========================================
# Crea los objetos figure y axis
fig, ax = plt.subplots()

# Define el rango de los ejes
ax.axis("equal")  # Misma escala para ejes X e Y
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# Representa el primer fotograma
# Pinta un punto en la posición de cada planeta y guarda el objeto asociado
# al punto en una lista
planet_points = list()
planet_trails = list()
for planet_data in frames_data[0]:
    _, x, y, radius, rock = planet_data
    color_planet = 'blue' if rock == 0 else 'red'
    planet_point = Circle((x, y), radius, color=color_planet)
    ax.add_artist(planet_point)
    planet_points.append(planet_point)

    # Inicializa las estelas (si especificado en los parámetros)
    if show_trail:
        planet_trail, = ax.plot(
            x, y, "-", linewidth=trail_width,
            color=planet_point.get_facecolor())
        planet_trails.append(planet_trail)


# Función que actualiza la posición de los planetas en la animación
def update(j_frame, frames_data, planet_points, planet_trails, show_trail):
    # Actualiza la posición del correspondiente a cada planeta
    for j_planet, planet_data in enumerate(frames_data[j_frame]):
        _, x, y, radius, rock = planet_data
        planet_points[j_planet].center = (x, y)
        planet_points[j_planet].radius = radius
        if show_trail:
            old_x, old_y = planet_trails[j_planet].get_data()
            new_x = np.append(old_x, x)
            new_y = np.append(old_y, y)
            planet_trails[j_planet].set_data(new_x, new_y)

    return planet_points + planet_trails


def init_anim():
    # Clear trails
    if show_trail:
        for j_planet in range(nplanets):
            planet_trails[j_planet].set_data(list(), list())

    return planet_points + planet_trails


# Calcula el nº de frames
nframes = len(frames_data)

# Si hay más de un instante de tiempo, genera la animación
if nframes > 1:
    # Info sobre FuncAnimation: https://matplotlib.org/stable/api/animation_api.html
    animation = FuncAnimation(
        fig, update, init_func=init_anim,
        fargs=(frames_data, planet_points, planet_trails, show_trail),
        frames=nframes, blit=True, interval=interval)

    # Muestra por pantalla o guarda según parámetros
    if save_to_file:
        animation.save("{}.gif".format(file_out), dpi=dpi)
    else:
        plt.show()
# En caso contrario, muestra o guarda una imagen
else:
    # Muestra por pantalla o guarda según parámetros
    if save_to_file:
        fig.savefig("{}.pdf".format(file_out))
    else:
        plt.show()