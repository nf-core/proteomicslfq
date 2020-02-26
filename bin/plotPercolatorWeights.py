#!/usr/bin/python

import csv
import os
import numpy as np
import matplotlib.pyplot as plt

mean_weight = None

files = [f for f in os.listdir('.') if os.path.isfile(f)]

count = 0

filename = args[0]
f = open(filename, 'r')
reader = csv.reader(f, delimiter='\t')
headers = next(reader, None)
converters = [str.strip] + [float] * (len(headers) - 1)
w1 = next(reader, None)
next(reader, None)
next(reader, None) # 2nd header
w2 = next(reader, None)
next(reader, None)
next(reader, None) # 3rd header
w3 = next(reader, None)

f1 = np.asarray([float(x) for x in w1])
f2 = np.asarray([float(x) for x in w2])
f3 = np.asarray([float(x) for x in w3])

mean_weight = (f1+f2+f3)/3.0

# set width of bar
barWidth = 0.25

# Set position of bar on X axis
r1 = np.asarray(np.arange(len(f1)))
r2 = np.asarray([x + barWidth for x in r1])
r3 = np.asarray([x + barWidth for x in r2])

# Make the plot
plt.bar(r1, mean_weight, color='#7f6d5f', width=barWidth, edgecolor='white', label='var1')
#plt.bar(r2, f2, color='#557f2d', width=barWidth, edgecolor='white', label='var2')
#plt.bar(r3, f3, color='#2d7f5e', width=barWidth, edgecolor='white', label='var3')

# Add xticks on the middle of the group bars
plt.xlabel('group', fontweight='bold')
plt.xticks([r + barWidth for r in range(len(f1))], headers, rotation=90)

# Create legend & Show graphic
plt.legend()
plt.show()