import matplotlib
from matplotlib.text import TextPath
import numpy as np
import sys
import subprocess

matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print('')
    print('Manhattan plot for GWAS -log(p)value.')
    print(
        'Usage: python manhattanplot.py [p-value file] [pvalue column],[chromosome column],[position column],[snp column] [threshold] [number of lines in header]')
    print('')
    sys.exit()

### in our case: python manhattan_plot.py Netblotch.txt 3,1,2,0 3 1
colors = ['steelblue', 'lightblue']
header = int(sys.argv[4])
pval = int(sys.argv[2].split(',')[0])
chrom = int(sys.argv[2].split(',')[1])
position = int(sys.argv[2].split(',')[2])
snp = int(sys.argv[2].split(',')[3])
threshold = float(sys.argv[3])
hsnps = []

with open(sys.argv[1], 'r') as f:
    data = f.readlines()

chroms = list(set(i.split()[chrom] for i in data[header:]))
chroms.sort(key=int)
plotdata = []
chromindices = []
for i in chroms:
    plotdata.append([])
    chromindices.append([])

for j, i in enumerate(data):
    k = i.split()
    try:
        if float(k[pval]) < threshold:
            plotdata[chroms.index(k[chrom])].append([k[pval], k[chrom], k[position], k[snp], j])
    except ValueError:
        continue

fig = plt.figure()
manhattanplot = fig.add_subplot(111, facecolor='white')
fig.set_facecolor('none')
fig.set_figheight(6)
fig.set_figwidth(12)

indexmax = 0
for j, i in enumerate(plotdata):
    if len(i) == 0:
        chromindices[j].append(0)
        chromindices[j].append(0)
        continue
    pvalues = np.absolute(np.log10([float(k[0]) for k in i]))
    indices = [k[-1] for k in i]
    color = colors[j % len(colors)]
    chromindices[j].append(np.min(indices))
    chromindices[j].append(np.max(indices))
    indexmax = np.max(indices + [indexmax])
    manhattanplot.scatter(indices, pvalues, marker='o', s=1, color=color)

plt.axhline(y=3., color='blue', alpha=0.5)

plt.xlim(-10, indexmax + 10)
plt.ylim(abs(np.log10(threshold)))

plt.tick_params(axis='x', which='both', top='off', direction='out')
plt.tick_params(axis='y', which='both', direction='out')
plt.xticks([i[0] + (i[1] - i[0]) / 2 for i in chromindices], chroms, fontsize=8, rotation=0)

plt.title('Manhattan plot of GWAS -log(p)value')
plt.xlabel('Chromosomes')
plt.ylabel('$-log_{10}(P)$')

plt.savefig('manhattan-plot-' + sys.argv[1].split('/')[-1].split('.')[0] + '.png', bbox_inches='tight', facecolor=fig.get_facecolor(),
            dpi=600, edgecolor='none')
