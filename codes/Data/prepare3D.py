import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def prepare3D(fileToOpen, fileToGenerate):
    index = []
    precisions = []
    times = []
    iterations = []
    errors = []
    with open(fileToOpen) as csv_file:
        for i, line in enumerate(csv_file):
            if (i % 2) == 0:
                line = line.split(',')
                if int(line[1]) == 10000:
                    if float(line[2]) == -15:
                        errors.append(float(line[3]))
                        times.append(float(line[4]))
                        iterations.append(int(line[5]))


    with open(fileToGenerate, mode='w') as data_traitees:
        data_writer = csv.writer(data_traitees, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        newrows = []
        for i in range(len(errors)):
            newrow = [errors[i], times[i], iterations[i]]
            newrows.append(newrow)
            data_writer.writerow(newrows[i])
        print(newrows)
