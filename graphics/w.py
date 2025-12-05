
# Визуализация работы последовательной и OpenMP программы
# Пример запуска python3 w.py ./txt/M_80__N_120__MPI_16.txt 3d inferno --mpi

import argparse
import numpy as np
import matplotlib.pyplot as plt
import ast
from mpl_toolkits.mplot3d import Axes3D
import re

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


def parse_file_sequential(filename):

    with open(filename) as f:
        lines = [line.rstrip() for line in f if line.strip()]

    legend = lines[0]
    i = 1
    matrix = " ".join(lines[1:])
    matrix = ast.literal_eval(matrix)
    return legend, matrix


# парсим строку с координатами
def parse_coords(header):
    match = re.search(r'coords=\((\d+),\s*(\d+)\)', header)
    if match:
        return int(match.group(1)), int(match.group(2))
    else:
        raise ValueError(f'Ошибка распознавания строки: {header}')


def parse_file_mpi(filename):
    blocks = {}

    with open(filename) as f:
        lines = [line.rstrip() for line in f if line.strip()]

    legend = lines[0]
    i = 1
    while i < len(lines):
        if lines[i].startswith('coords='):
            coords = parse_coords(lines[i])

            block_lines = []
            i += 1
            while i < len(lines) and not lines[i].startswith('coords='):
                block_lines.append(lines[i])
                i += 1
            block_str = " ".join(block_lines)
            block_str = block_str.replace(",]", "]")
            block = ast.literal_eval(block_str)
            block = [ row[1:-1] for row in block[1:-1] ]

            blocks[coords] = block
        else:
            i += 1
    block_size_x = max(i[0] for i in blocks) + 1
    block_size_y = max(i[1] for i in blocks) + 1
    
    matrix = []
    for y in range(block_size_y):
        matrix.append(np.hstack(([blocks[i] for i in blocks if i[1] == y])))
    matrix = np.vstack((matrix))

    return legend, matrix


def parse_file(filename, mpi):
    
    if mpi:
        return parse_file_mpi(filename)
    else:
        return parse_file_sequential(filename)


def main():
    parser = argparse.ArgumentParser(description="Визуализация работы последовательной и OpenMP программы")
    parser.add_argument('filename', type=str, help='Путь до файл с матрицей значений')
    parser.add_argument('--mpi', action='store_true', help=' Если флаг установлен, матрица значение получена с помощью MPI и ее надо собрать по кусочкам')
    parser.add_argument('mod', type=str, nargs='?', default='2d', help='Тип визулизации: "2d" или "3d"')
    parser.add_argument('cmap', type=str, nargs='?', default='inferno', help='Цвет "inferno" или "viridis"')
    args = parser.parse_args()

    legend, matrix = parse_file(args.filename, args.mpi)
    visvisualization(legend, matrix, args.mod, args.cmap)


if __name__ == '__main__':
    main()
