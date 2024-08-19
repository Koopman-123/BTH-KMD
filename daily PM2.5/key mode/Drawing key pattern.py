import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("mode134_dataset.csv", header=None)

cmap = plt.get_cmap('Spectral')
reversed_cmap = cmap.reversed()

fig, ax = plt.subplots(figsize=(12, 7))

X, Y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))

cax = ax.pcolormesh(X, Y, data, cmap=reversed_cmap, shading='auto', vmin=-10, vmax=10)

cbar = fig.colorbar(cax, pad=0.08)

cbar.ax.tick_params(labelsize=32)

dates = pd.date_range(start='2019-01-01', end='2021-12-31')

ax.set_xticks(np.arange(0, data.shape[1], 100))
ax.set_xticklabels(dates[::100].strftime('%Y-%m-%d'), rotation=45, fontsize=32)
ax.set_yticks(np.arange(0, data.shape[0], 10))
ax.set_yticklabels(np.arange(0, data.shape[0], 10), fontsize=32)

ax2 = ax.twinx()

yticks = [23, 49, 54, 59, 64, 76]
ylabels = ['1', '2', '3', '4', '5', '6']
ax2.set_yticks(yticks)
ax2.set_yticklabels(ylabels, fontsize=32)
# ax2.set_ylabel('Region', fontsize=28)

for tick in yticks:
    ax2.axhline(y=tick, color='black', linewidth=1, linestyle='--', xmin=0, xmax=1)

ax.tick_params(axis='both', which='major', labelsize=32)
ax2.tick_params(axis='both', which='major', labelsize=32)

# ax.set_title('Mode3 (353day)', fontsize=36)
# ax.set_title('Mode6 (185day)', fontsize=36)
# ax.set_title('Mode9 (124day)', fontsize=36)
# ax.set_title('Mode48 (22day)', fontsize=36)
ax.set_title('Mode134 (8day)', fontsize=36)

plt.subplots_adjust(bottom=0.25)

plt.show()
