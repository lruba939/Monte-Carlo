import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import pdist, squareform

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

# Dodawanie punktów
ax.scatter(points[:,0], points[:,1], points[:,2], color="red")

# Dodawanie linii łączących punkty
for connection in connections:
    ax.plot3D(*zip(*connection), color='black', linewidth=1.5, alpha=0.8)

# Ustawienie etykiet osi
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_aspect('equal')

# Wyświetlenie wykresu
plt.show()