
# сравниваем два файла с результирующими матрицами
# опциональные флаги --mpi1 --mpi2 указыают что соответствующая матрица получена в резальтае работы mpi программы и ее надо собрать по кусочкам
# пример запуска: python3 compare.py ./../sequential/txt/M_80__N_120.txt ./M_80__N_120__MPI_2.txt --mpi2

import argparse
import numpy as np
import matplotlib.pyplot as plt
import ast
from mpl_toolkits.mplot3d import Axes3D
import re

def parse_file_sequential(filename):

    with open(filename) as f:
        lines = [line.rstrip() for line in f if line.strip()]

    legend = lines[0]
    i = 1
    matrix = " ".join(lines[1:])
    matrix = ast.literal_eval(matrix)
    matrix = [ row[1:-1] for row in matrix[1:-1] ]
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


def compare_matrix(legend1, matrix1, legend2, matrix2):
    print("legend1: ", legend1)
    print("legend2: ", legend2)

    if len(matrix1) != len(matrix2):
        print("not equal, size_y")
        return
    for j in range(len(matrix1)):
        if len(matrix1) != len(matrix2):
            print(f"not equal, y={j}, size_x")
            return
        for i in range(len(matrix1[j])):
            if matrix1[j][i] != matrix2[j][i]:
                print(f"not equal, y={j}, x={i}, value")
                return
    print("equal")


def main():
    parser = argparse.ArgumentParser(description="Визуализация работы последовательной и OpenMP программы")
    parser.add_argument('filename_1', type=str, help='Путь до первого файла с матрицей значений')
    parser.add_argument('--mpi1', action='store_true', help=' Если флаг установлен, первая матрица значение получена с помощью MPI и ее надо собрать по кусочкам')
    parser.add_argument('filename_2', type=str, help='Путь до второго файла с матрицей значений')
    parser.add_argument('--mpi2', action='store_true', help=' Если флаг установлен, вторая матрица значение получена с помощью MPI и ее надо собрать по кусочкам')
    args = parser.parse_args()

    legend1, matrix1 = parse_file(args.filename_1, args.mpi1)
    legend2, matrix2 = parse_file(args.filename_2, args.mpi2)
    compare_matrix(legend1, matrix1, legend2, matrix2)



if __name__ == '__main__':
    main()
