import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

# 5, 10, 25, 50, 10

def prepareBars(fileToOpen, fileToGenerate):
    errors_size5 = []
    errors_size10 = []
    errors_size25 = []
    errors_size50 = []
    errors_size100 = []
    errors_size250 = []
    errors_size500 = []
    errors_size1000 = []
    errors_size5000 = []
    errors_size10000 = []
    times_size5 = []
    times_size10 = []
    times_size25 = []
    times_size50 = []
    times_size100 = []
    times_size250 = []
    times_size500 = []
    times_size1000 = []
    times_size5000 = []
    times_size10000 = []
    iter_size5 = []
    iter_size10 = []
    iter_size25 = []
    iter_size50 = []
    iter_size100 = []
    iter_size250 = []
    iter_size500 = []
    iter_size1000 = []
    iter_size5000 = []
    iter_size10000 = []
    with open(fileToOpen) as csv_file:
        for i, line in enumerate(csv_file):
            if (i % 2) == 0:
                line = line.split(',')
                print(line[2])
                if float(line[2]) == -15:
                    if int(line[1]) == 5:
                        errors_size5.append(float(line[3]))
                        times_size5.append(float(line[4]))
                        iter_size5.append(int(line[5]))
                    if int(line[1]) == 10:
                        errors_size10.append(float(line[3]))
                        times_size10.append(float(line[4]))
                        iter_size10.append(int(line[5]))
                    if int(line[1]) == 25:
                        errors_size25.append(float(line[3]))
                        times_size25.append(float(line[4]))
                        iter_size25.append(int(line[5]))
                    if int(line[1]) == 50:
                        errors_size50.append(float(line[3]))
                        times_size50.append(float(line[4]))
                        iter_size50.append(int(line[5]))
                    if int(line[1]) == 100:
                        errors_size100.append(float(line[3]))
                        times_size100.append(float(line[4]))
                        iter_size100.append(int(line[5]))
                    if int(line[1]) == 250:
                        errors_size250.append(float(line[3]))
                        times_size250.append(float(line[4]))
                        iter_size250.append(int(line[5]))
                    if int(line[1]) == 500:
                        errors_size500.append(float(line[3]))
                        times_size500.append(float(line[4]))
                        iter_size500.append(int(line[5]))
                    if int(line[1]) == 1000:
                        errors_size1000.append(float(line[3]))
                        times_size1000.append(float(line[4]))
                        iter_size1000.append(int(line[5]))
                    if int(line[1]) == 5000:
                        errors_size5000.append(float(line[3]))
                        times_size5000.append(float(line[4]))
                        iter_size5000.append(int(line[5]))
                    if int(line[1]) == 1000:
                        errors_size10000.append(float(line[3]))
                        times_size10000.append(float(line[4]))
                        iter_size10000.append(int(line[5]))

        meanResolError5 = np.mean(errors_size5)
        meanResolError10 = np.mean(errors_size10)
        meanResolError25 = np.mean(errors_size25)
        meanResolError50 = np.mean(errors_size50)
        meanResolError100 = np.mean(errors_size100)
        meanResolError250 = np.mean(errors_size250)
        meanResolError500 = np.mean(errors_size500)
        meanResolError1000 = np.mean(errors_size1000)
        meanResolError5000 = np.mean(errors_size5000)
        meanResolError10000 = np.mean(errors_size10000)
        meanResolError = [meanResolError5, meanResolError10, meanResolError25, meanResolError50, meanResolError100,
                          meanResolError250, meanResolError500, meanResolError1000, meanResolError5000, meanResolError10000]
        meanResolTime5 = np.mean(times_size5)
        meanResolTime10 = np.mean(times_size10)
        meanResolTime25 = np.mean(times_size25)
        meanResolTime50 = np.mean(times_size50)
        meanResolTime100 = np.mean(times_size100)
        meanResolTime250 = np.mean(times_size250)
        meanResolTime500 = np.mean(times_size500)
        meanResolTime1000 = np.mean(times_size1000)
        meanResolTime5000 = np.mean(times_size5000)
        meanResolTime10000 = np.mean(times_size10000)
        meanIter5 = np.mean(iter_size5)
        meanIter10 = np.mean(iter_size10)
        meanIter25 = np.mean(iter_size25)
        meanIter50 = np.mean(iter_size50)
        meanIter100 = np.mean(iter_size100)
        meanIter250 = np.mean(iter_size250)
        meanIter500 = np.mean(iter_size500)
        meanIter1000 = np.mean(iter_size1000)
        meanIter5000 = np.mean(iter_size5000)
        meanIter10000 = np.mean(iter_size10000)
        meanResolTime = [meanResolTime5, meanResolTime10, meanResolTime25, meanResolTime50, meanResolTime100,
                         meanResolTime250, meanResolTime500, meanResolTime1000, meanResolTime5000, meanResolTime10000]
        meanIter = [meanIter5, meanIter10, meanIter25, meanIter50, meanIter100,
                    meanIter250, meanIter500, meanIter1000, meanIter5000, meanIter10000]
        print(meanResolError)
        print(meanResolTime)
        print(meanIter)

    with open(fileToGenerate, mode='w') as data_traitees:
        data_writer = csv.writer(data_traitees, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        data_writer.writerow(meanResolError)
        data_writer.writerow(meanResolTime)
        data_writer.writerow(meanIter)
