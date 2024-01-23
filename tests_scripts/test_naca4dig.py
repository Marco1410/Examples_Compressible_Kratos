import numpy as np
import matplotlib.pyplot as plt

def naca4_digit_series(number, c, n=100):
    # Desglosa el número en sus dígitos
    m = number // 1000 / 100.0
    p = (number // 100 % 10) / 10.0
    t = (number % 100) / 100.0

    # Crea el rango de puntos en el eje x
    x = np.linspace(0, c, n)

    # Calcula el grosor
    yt = 5 * t * (0.2969 * np.sqrt(x / c) - 0.1260 * (x / c) - 0.3516 * (x / c)**2 + 0.2843 * (x / c)**3 - 0.1015 * (x / c)**4)

    # Calcula la línea media del perfil
    yc = np.where((x / c) <= p, m * (2 * p * (x / c) - (x / c)**2) / p**2, m * (1 - 2 * p + 2 * p * (x / c) - (x / c)**2) / (1 - p)**2)

    # Calcula los puntos superiores e inferiores
    xu = x - yt * np.sin(np.arctan2(np.diff(yc, prepend=0), np.diff(x, prepend=0)))
    yu = yc + yt * np.cos(np.arctan2(np.diff(yc, prepend=0), np.diff(x, prepend=0)))
    xl = x + yt * np.sin(np.arctan2(np.diff(yc, prepend=0), np.diff(x, prepend=0)))
    yl = yc - yt * np.cos(np.arctan2(np.diff(yc, prepend=0), np.diff(x, prepend=0)))

    return xu, yu, xl, yl, x, yc

# Parámetros del perfil NACA
naca_number = 2412  # Número NACA de cuatro dígitos
chord_length = 1    # Longitud de la cuerda

# Genera los puntos del perfil NACA
x_upper, y_upper, x_lower, y_lower, x_camber, y_camber = naca4_digit_series(naca_number, chord_length)

# Gráfica del perfil NACA
plt.figure(figsize=(10, 5))
plt.plot(x_upper, y_upper, label="Superficie Superior")
plt.plot(x_lower, y_lower, label="Superficie Inferior")
plt.plot(x_camber, y_camber, label="Línea Media", linestyle='--')
plt.title(f"Perfil NACA {naca_number}")
plt.xlabel("Cuerda (x/c)")
plt.ylabel("Altura (y/c)")
plt.axis('equal')
plt.legend()
plt.grid(True)
plt.show()

