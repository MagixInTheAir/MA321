import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

# 5, 10, 25, 50, 100
index = []
matrix_sizes = []
errors_size5 = []
errors_size10 = []
errors_size25 = []
errors_size50 = []
errors_size100 = []
times = []
times_size5 = []
times_size10 = []
times_size25 = []
times_size50 = []
times_size100 = []

with open('gauss-seidel_traite.txt') as csv_file:
    for i, line in enumerate(csv_file):
        if (i % 2) == 0:
            line = line.split(',')
            print(line[2])
            if float(line[2]) == -5:
                if int(line[1]) == 5:
                    errors_size5.append(float(line[3]))
                    times_size5.append(float(line[4]))
                if int(line[1]) == 10:
                    errors_size10.append(float(line[3]))
                    times_size10.append(float(line[4]))
                if int(line[1]) == 25:
                    errors_size25.append(float(line[3]))
                    times_size25.append(float(line[4]))
                if int(line[1]) == 50:
                    errors_size50.append(float(line[3]))
                    times_size50.append(float(line[4]))
                if int(line[1]) == 100:
                    errors_size100.append(float(line[3]))
                    times_size100.append(float(line[4]))

    meanResolError5 = np.mean(errors_size5)
    meanResolError10 = np.mean(errors_size10)
    meanResolError25 = np.mean(errors_size25)
    meanResolError50 = np.mean(errors_size50)
    meanResolError100 = np.mean(errors_size100)
    meanResolError = [meanResolError5, meanResolError10, meanResolError25, meanResolError50, meanResolError100]
    print(meanResolError)
