import re
import ast
import numpy as np

# собираем матрицу из блоков

def parse_coords(header):
    match = re.search(r'coords=\((\d+),\s*(\d+)\)', header)
    if match:
        return int(match.group(2)), int(match.group(1))
    else:
        raise ValueError(f'Ошибка распознавания строки: {header}')    

def assemble_block_matrix(filename):
    blocks = []
    block_row_col = set()
    block_shape = None

    with open(filename) as f:
        lines = [line.rstrip() for line in f if line.strip()]

    i = 0
    while i < len(lines):
        if lines[i].startswith('coords='):
            coords = parse_coords(lines[i])
            # Собираем многострочный список
            block_lines = []
            i += 1
            while i < len(lines) and not lines[i].startswith('coords='):
                block_lines.append(lines[i])
                i += 1
            block_str = " ".join(block_lines)
            # Убираем лишние запятые внутри списка
            block_str = block_str.replace(",]", "]")
            block = ast.literal_eval(block_str)
            # Удаляем граничные элементы
            block = [ row[1:-1] for row in block[1:-1] ]
            if block_shape is None:
                block_shape = (len(block), len(block[0]))
            blocks.append({'coords': coords, 'data': block})
            block_row_col.add(coords)
        else:
            i += 1

    n_blocks_rows = max(c[0] for c in block_row_col) + 1
    n_blocks_cols = max(c[1] for c in block_row_col) + 1
    block_rows, block_cols = block_shape
    total_rows = n_blocks_rows * block_rows
    total_cols = n_blocks_cols * block_cols
    full_matrix = np.zeros((total_rows, total_cols), dtype=float)

    for b in blocks:
        i_block, j_block = b['coords']
        data = np.array(b['data'])
        row_start = i_block * block_rows
        row_end = row_start + block_rows
        col_start = j_block * block_cols
        col_end = col_start + block_cols
        full_matrix[row_start:row_end, col_start:col_end] = data

    return full_matrix

filename = "all.txt"
matrix = assemble_block_matrix(filename)


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

m, n = 800, 1200

x = np.linspace(-4, 4, m)
y = np.linspace(-1, 3, n)
X, Y = np.meshgrid(x, y)

Z = matrix

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, Y, Z, cmap='inferno', edgecolor='k')

ax.set_xlabel('x')
ax.set_ylabel('y')
#ax.set_zlabel('f(x, y)')
ax.set_zlim(Z.min(), Z.min() + 1)
ax.set_title('Значение w')

fig.colorbar(surf, ax=ax, shrink=0.6, aspect=10, label='значение')
plt.tight_layout()
plt.show()