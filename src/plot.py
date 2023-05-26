# Copyright (C) BGI-Reasearch - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
import os, sys
import pandas as pd
import numpy as np
import csv
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def main():
    if len(sys.argv) != 4:
        print("enter <saturationFile> <ratio> <outdir>")
        sys.exit(-1)

    saturationFile = sys.argv[1]
    ratio = float(sys.argv[2])
    outdir = sys.argv[3]
    getSaturationFig(saturationFile, ratio, outdir)

def getSaturationFig(saturationFile, ratio, outdir, binSize=200, readsScale=1):
    os.makedirs(outdir, exist_ok=True)
    sadf = pd.read_csv(saturationFile, sep="\t", quoting=csv.QUOTE_NONE)
    sadf['bar_x_total'] = sadf['bar_x']/ratio
    sadf['bin_x_total'] = sadf['bin_x']/ratio

    labelsize = 14
    fontsize = 18
    fig = plt.figure(figsize=(12,5),dpi=100)
    ax = plt.subplot(1, 2, 1)
    plt.tick_params(labelsize=labelsize)
    requiredSa = 0.8
    ax.plot(sadf['bar_x_total']/readsScale, sadf['bar_y1'], marker='o')
    ax.set_xlabel("Total reads number of sampling", fontsize=fontsize)
    ax.set_ylabel("Sequencing saturation", fontsize=fontsize)
    plt.axhline(requiredSa, color='green', lw=1, alpha=0.7)
    ax.grid()
    ax = plt.subplot(1, 2, 2)
    plt.tick_params(labelsize=labelsize)
    ax.plot(sadf['bar_x_total']/readsScale, sadf['bar_y2'], marker='o')
    ax.set_xlabel("Total reads number of sampling", fontsize=fontsize)
    ax.set_ylabel("Median genes per barcode", fontsize=fontsize)
    ax.grid()
    figFilePath = os.path.join(outdir, "plot_1x1_saturation.png")
    plt.tight_layout()
    plt.savefig(figFilePath, format="png", bbox_inches='tight')
    plt.clf()
    
    #plot saturation figure of bin200
    fig=plt.figure(figsize=(16,5),dpi=100)
    ax = plt.subplot(1, 3, 1)
    plt.tick_params(labelsize=labelsize)
    ax.plot(sadf['bin_x_total']/readsScale, sadf['bar_y1'], marker='o')
    ax.set_xlabel("Total reads number of sampling", fontsize=fontsize)
    ax.set_ylabel("Sequencing saturation", fontsize=fontsize)
    plt.axhline(requiredSa, color='green', lw=1, alpha=0.7)
    ax.grid()

    ax = plt.subplot(1, 3, 2)
    plt.tick_params(labelsize=labelsize)
    ax.plot(sadf['bin_x_total']/readsScale, sadf['bin_y2'], marker='o')
    ax.set_xlabel("Total reads number of sampling", fontsize=fontsize)
    ax.set_ylabel("Median genes per bin", fontsize=fontsize)
    ax.grid()

    xData = sadf['bin_x_total']
    yData = sadf['bin_umi']
    threshold = 5000
    maxY = max(yData)
    Rsquared, fittedParameters = _cur_fit(xData, yData)
    bstr=str(round(fittedParameters[1], 2)) if fittedParameters[1]<0 else "+{0}".format(round(fittedParameters[1], 2))
    cstr=str(round(fittedParameters[2], 2)) if fittedParameters[2]<0 else "+{0}".format(round(fittedParameters[2], 2))
    dstr="%.2e"%fittedParameters[3]
    rstr="R\u00b2={:0.3f}".format(Rsquared)
    labelstr='y={0}*log({1}x{2}){3}\n{4}'.format(round(fittedParameters[0],2), dstr, bstr, cstr, rstr)
    maxX = _reFunc(threshold, *fittedParameters)

    ax = plt.subplot(1, 3, 3)
    plt.tick_params(labelsize=labelsize)
    ax.plot(xData/readsScale, yData, marker='o')
    innerfontsize = 10
    if (maxY < threshold and Rsquared >= 0.9 and maxX < 100000000000):
        #maxX = _reFunc(threshold, *fittedParameters)
        xModel = np.linspace(min(xData), max(xData))
        yModel = _func(xModel, *fittedParameters)
        xModel1 = np.linspace(max(xData), maxX)
        yModel1 = _func(xModel1, *fittedParameters)
        ax.plot(xModel/readsScale, yModel, color="red", lw=1)
        ax.plot(xModel1/readsScale, yModel1, color="red", linestyle='dashed', lw=0.5)
        plt.axvline(maxX/readsScale, color='green', lw=1, alpha=0.7)
        plt.text(int(maxX/readsScale),threshold-100,(int(maxX/readsScale),threshold),color='b', fontsize=innerfontsize)
    ax.set_xlabel("Total reads number of sampling", fontsize=fontsize)
    ax.set_ylabel("unique reads number per bin", fontsize=fontsize)
    ax.grid()
    #ax.legend()
    plt.axhline(threshold, color='green', lw=1, alpha=0.7)
    plt.text(int(max(xData)/readsScale/1.8),max(yData),(int(max(xData)/readsScale),max(yData)),color='b', fontsize=innerfontsize)
    plt.text(min(xData)*2, min(yData)+100, labelstr, fontsize=innerfontsize)
    figFilePath = os.path.join(outdir, "plot_{0}x{0}_saturation.png".format(binSize))
    plt.tight_layout()
    plt.savefig(figFilePath, format="png", bbox_inches='tight')

def _cur_fit(xData, yData):
    initialParameters = np.array([1.0, 1.0, 1.0, 1.0])
    Rsquared = 0
    try:
        fittedParameters, pcov = curve_fit(_func, xData, yData, initialParameters)
        modelPredictions = _func(xData, *fittedParameters)
        absError = modelPredictions - yData
        SE = np.square(absError) # squared errors
        MSE = np.mean(SE) # mean squared errors
        RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
        Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    except:
        fittedParameters = initialParameters
        Rsquared = 0
    return Rsquared, fittedParameters

def _func(x, a, b, c, d):
    result = a*np.log(d * x + b) + c
    return result

def _reFunc(y, a, b, c, d):
    result = (np.power(np.e, (y-c)/a) - b) / d
    return result

if __name__ == "__main__":
    main()
