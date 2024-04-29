import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import pdist, squareform

# FUNCKA
#######################################################################################
def atoms_plot_3D(data):
    points = data

    # Obliczenie odległości między wszystkimi parami punktów
    distances = squareform(pdist(points))

    # Wartość maksymalna, która będzie uznana za wystarczająco bliską
    max_distance = 1.8  # Tutaj możesz ustawić odpowiednią wartość

    # Tworzenie listy przechowującej połączenia między punktami
    connections = []

    # Szukanie połączeń
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            if distances[i, j] <= max_distance:
                connections.append((points[i], points[j]))

    # Tworzenie wykresu 3D
    fig = plt.figure(figsize=(4,4), dpi=200)
    ax = fig.add_subplot(111, projection='3d')

    r = 3.0

    # Make data
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))

    # Plot the surface
    ax.plot_surface(x, y, z, zorder=-1, alpha=0.3)

    # Dodawanie linii łączących punkty
    for connection in connections:
        ax.plot3D(*zip(*connection), color='black', linewidth=1.5, alpha=1, zorder=-1)

    # Dodawanie punktów
    ax.scatter(points[:,0], points[:,1], points[:,2], color="tab:orange", zorder=1)

    # Ustawienie etykiet osi
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_aspect('equal')

##################################################################################################
data = np.loadtxt("data.dat", dtype=float, delimiter=" ")
print(data)

X = data[:,0]
Y = data[:,1]
Z = data[:,2]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Plot the surface
ax.plot(X, Y, Z, "ko")

# Set an equal aspect ratio
ax.set_aspect('equal')

plt.show()

zad3 = np.loadtxt("atoms_avogadro_zad3.xyz", dtype=float, delimiter=" ", skiprows=1, usecols=range(1,4))
zad4 = np.loadtxt("atoms_avogadro_zad4.xyz", dtype=float, delimiter=" ", skiprows=1, usecols=range(1,4))
zad5 = np.loadtxt("atoms_avogadro_zad5.xyz", dtype=float, delimiter=" ", skiprows=1, usecols=range(1,4))

atoms_plot_3D(zad5)
atoms_plot_3D(zad4)
atoms_plot_3D(zad3)
atoms_plot_3D(data)

plt.show()