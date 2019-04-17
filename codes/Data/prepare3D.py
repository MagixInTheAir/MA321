import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

index = []
precisions = []
times = []
iterations = []
errors = []

with open('sor_traite.txt') as csv_file:
    for i, line in enumerate(csv_file):
        if (i % 2) == 0:
            line = line.split(',')
            if int(line[1]) == 100:
                if float(line[2]) == -5:
                    errors.append(float(line[3]))
                    times.append(float(line[4]))
                    iterations.append(int(line[5]))


with open('sor_for3D.txt', mode='w') as data_traitees:
    data_writer = csv.writer(data_traitees, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    newrows = []
    for i in range(len(errors)):
        newrow = [errors[i], times[i], iterations[i]]
        newrows.append(newrow)
        data_writer.writerow(newrows[i])
    print(newrows)