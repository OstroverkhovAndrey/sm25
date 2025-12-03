
# Визуализация работы последовательной и OpenMP программы

import argparse
import numpy as np
import matplotlib.pyplot as plt
import ast
from mpl_toolkits.mplot3d import Axes3D

# Границы прямоугольника
xmin, xmax = -4.0, 4.0
ymin, ymax = -1.0, 3.0

def visvisualization_2d(legend, matrix, cmap):

    Nx = len(matrix[0])
    Ny = len(matrix)

    x = np.linspace(xmin, xmax, Nx)
    y = np.linspace(ymin, ymax, Ny)
    X, Y = np.meshgrid(x, y)

    Z = matrix

    plt.figure(figsize=(7, 4))
    im = plt.imshow(Z, origin='lower', extent=[xmin, xmax, ymin, ymax], cmap=cmap, interpolation='nearest')
    plt.colorbar(im, label='значение')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(legend)
    plt.tight_layout()
    plt.show()


def visvisualization_3d(legend, matrix, cmap):

    Nx = len(matrix[0])
    Ny = len(matrix)

    x = np.linspace(xmin, xmax, len(matrix[0]))
    y = np.linspace(ymin, ymax, len(matrix))
    X, Y = np.meshgrid(x, y)

    Z = np.array(matrix)

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    surface = ax.plot_surface(X, Y, Z, cmap=cmap, rcount=100, ccount=100)
    ax.set_box_aspect([np.ptp(x), np.ptp(y), 2])
    ax.set_zlim(Z.min(), Z.min() + 0.5)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(legend)
    fig.colorbar(surface, ax=ax, shrink=0.5, aspect=5, pad=0.15)
    plt.show()


def visvisualization(legend, matrix, mod, cmap):

    if mod == '2d':
        visvisualization_2d(legend, matrix, cmap)
    elif mod == '3d':
        visvisualization_3d(legend, matrix, cmap)

def parse_file(filename):

    with open(filename) as f:
        lines = [line.rstrip() for line in f if line.strip()]

    legend = lines[0]
    i = 1
    matrix = " ".join(lines[1:])
    matrix = ast.literal_eval(matrix)
    return legend, matrix


def main():
    parser = argparse.ArgumentParser(description="Визуализация работы последовательной и OpenMP программы")
    parser.add_argument('filename', type=str, help='Путь до файл с матрицей значений')
    parser.add_argument('mod', type=str, nargs='?', default='2d', help='Тип визулизации: "2d" или "3d"')
    parser.add_argument('cmap', type=str, nargs='?', default='inferno', help='Цвет "inferno" или "viridis"')
    args = parser.parse_args()

    legend, matrix = parse_file(args.filename)
    visvisualization(legend, matrix, args.mod, args.cmap)


if __name__ == '__main__':
    main()
