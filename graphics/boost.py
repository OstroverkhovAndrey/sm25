import numpy as np
import matplotlib.pyplot as plt

# Фиксированный набор потоков на оси X
threads = np.array([1, 2, 4, 8, 16, 32], dtype=int)

cases = {
    "40 X 40": {
        "value": {1: 1.0, 2: 2.3, 4: 3.2,8: 4.1, 16: 4.3},
        "color": "red",
    },
    "400 X 600": {
        "value": {1: 1.0, 2: 2.3, 4: 3.2,8: 4.1, 16: 4.3, 32: 7.1},
        "color": "green",
    },
    "800 X 1200": {
        "value": {1: 1.0, 2: 2.3, 4: 3.2,8: 4.1, 16: 4.3, 32: 7.1},
        "color": "blue",
    },
}

def speedup_array(time_by_threads, threads_axis):
    # Вернёт массив ускорений согласно threads_axis, для отсутствующих нитей — NaN
    y = []
    for t in threads_axis:
        if t in time_by_threads:
            y.append(time_by_threads[t])
        else:
            y.append(np.nan)
    return np.array(y, dtype=float)

titles = list(cases.keys())

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 10), sharex=True, sharey=True)

for ax, title in zip(axes, titles):

    s = speedup_array(cases[title]['value'], threads)
    ax.plot(threads, s, marker="o", color=cases[title]["color"])  # NaN создаёт разрывы линии

    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_ylabel("Ускорение")
    ax.set_xticks(threads)
    ax.legend(loc="upper left")

axes[-1].set_xlabel("Число нитей")

fig.suptitle("Зависимость ускорения от числа нитей", y=0.98)
fig.tight_layout()
plt.show()