
# Строит графики ускорений разных версий программ
# Пример запуска: python3 boost.py [ openmp | mpi | mpi_openmp ]
import matplotlib.pyplot as plt
import argparse

def openmp():
    d1 = {2: 1.740, 4: 3.471, 8: 6.587, 16: 7.164}
    d2 = {4: 3.500, 8: 6.840, 16: 7.533, 32: 8.011}

    xticks = [2, 4, 8, 16, 32]

    def plot_runtime(ax, data: dict, title: str, color):
        xs = sorted(data.keys())
        ys = [data[t] for t in xs]
        ax.plot(xs, ys, marker='o', linewidth=2, color=color)
        ax.set_title(title)
        ax.set_ylabel('Ускорение')
        ax.grid(True, alpha=0.3)
        ax.set_xticks(xticks)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 9), sharex=True)

    plot_runtime(axes[0], d1, '400 X 600', 'green')
    plot_runtime(axes[1], d2, '800 X 1200', 'blue')

    axes[-1].set_xlabel('Число нитей')

    fig.suptitle('Зависимость ускорения от числа нитей', y=0.98)
    fig.tight_layout()
    plt.show()


def mpi():
    d1 = {2: 2.061, 4: 3.325, 8: 8.242, 16: 16.000}
    d2 = {4: 4.130, 8: 8.280, 16: 13.105, 32: 32.126}

    xticks = [2, 4, 8, 16, 32]

    def plot_runtime(ax, data: dict, title: str, color):
        xs = sorted(data.keys())
        ys = [data[t] for t in xs]
        ax.plot(xs, ys, marker='o', linewidth=2, color=color)
        ax.set_title(title)
        ax.set_ylabel('Ускорение')
        ax.grid(True, alpha=0.3)
        ax.set_xticks(xticks)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 9), sharex=True)

    plot_runtime(axes[0], d1, '400 X 600', 'green')
    plot_runtime(axes[1], d2, '800 X 1200', 'blue')

    axes[-1].set_xlabel('Число процессов')

    fig.suptitle('Зависимость ускорения от числа процессов', y=0.98)
    fig.tight_layout()
    plt.show()


def mpi_openmp():
    d2 = {1: 1.274, 2: 2.513, 4: 4.946, 8: 9.503}
    d3 = {1: 2.570, 2: 5.000, 4: 9.655, 8: 19.138}

    xticks = [1, 1, 2, 4, 8]

    def plot_runtime(ax, data: dict, title: str, color):
        xs = sorted(data.keys())
        ys = [data[t] for t in xs]
        ax.plot(xs, ys, marker='o', linewidth=2, color=color)
        ax.set_title(title)
        ax.set_ylabel('Ускорение')
        ax.grid(True, alpha=0.3)
        ax.set_xticks(xticks)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 9), sharex=True)

    plot_runtime(axes[0], d2, '400 X 600 (2 mpi процесса)', 'green')
    plot_runtime(axes[1], d3, '800 X 1200 (4 mpi процесса)', 'blue')

    axes[-1].set_xlabel('Число нитей')

    fig.suptitle('Зависимость ускорения от числа нитей', y=0.98)
    fig.tight_layout()
    plt.show()

def impl(mod):
    if mod == 'openmp':
        openmp()
    elif mod == 'mpi':
        mpi()
    elif mod == 'mpi_openmp':
        mpi_openmp()
    else:
        print('Unsupported mod')

def main():
    parser = argparse.ArgumentParser(description="Визуализация ускорений прграмм")
    parser.add_argument('mod', type=str, help='Какой именно график хотим построить openmp, mpi или mpi_openmp')
    args = parser.parse_args()

    impl(args.mod)


if __name__ == '__main__':
    main()