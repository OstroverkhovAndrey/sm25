import matplotlib.pyplot as plt

d1 = {1: 1.000, 4: 1.680, 16: 1.059}
d2 = {1: 1.000, 2: 2.008, 4: 3.974, 8: 7.280, 16: 8.583}
d3 = {1: 1.000, 4: 3.999, 8: 7.522, 16: 8.623, 32: 9.185}

xticks = [1, 2, 4, 8, 16, 32]

def plot_runtime(ax, data: dict, title: str, color):
    xs = sorted(data.keys())
    ys = [data[t] for t in xs]
    ax.plot(xs, ys, marker='o', linewidth=2, color=color)
    ax.set_title(title)
    ax.set_ylabel('Ускорение')
    ax.grid(True, alpha=0.3)
    ax.set_xticks(xticks)

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 9), sharex=True)

plot_runtime(axes[0], d1, '40 X 40', 'red')
plot_runtime(axes[1], d2, '400 X 600', 'green')
plot_runtime(axes[2], d3, '800 X 1200', 'blue')

axes[-1].set_xlabel('Число нитей')

fig.suptitle('Зависимость ускорения от числа нитей', y=0.98)
fig.tight_layout()
plt.show()