from importlib.metadata import files
from idna import check_initial_combiner
from ROOT import TF1, TH1F, gPad, TCanvas, TFile, TGraph, kRed, kBlack, TMultiGraph, TLegend, kBlue, kGreen, kMagenta, kOrange
import numpy as np
import array as arr
import codecs
import math
import scipy 
import json
import os 
import natsort
import random as rand
from multiprocessing import Pool
import time

def split(string): #Scoate toate literele in ordine dintr-un string
        return [char for char in string]

def mergeLists(l1, l2): #Uneste doualiste, primul element cu primul element, al doile cu al doilea etc.
    return list(map(lambda x, y:(x,y), l1, l2)) 



def getInfo(files): #Scoate informatiile de care ai nevoie din fisierel e tip info pe care le scoate Janus
    infofiles = files
    holdDelay = np.array([])
    with codecs.open(infofiles, 'r', 'utf-8', errors='ignore') as file:
            for i, line in enumerate(file):
                line = line.strip()
                info = line.split()
                if(len(info) == 0):
                 i = i+1
                else:
                    if info[0] == 'LG_ShapingTime':
                        holdDelay = np.append(holdDelay, info[1])
    return holdDelay


def getPeak(hist, start, end): #Scoate niste parametrii preliminari pentru a fita gausiene
    mean = 0
    height = 0
    fhwm = 0
    sigma = 0
    baseline = np.array([])
    par = np.array([]), 
    print ('start = ', start, "end = ", end)
    for i in range(start, start + 10):
        baseline = np.append(baseline, hist.GetBinContent(i))
    bslValue = np.sum(baseline) / len(baseline)
    offset = 50
    hist.GetXaxis().SetRange(start, end)
    treshold = float(hist.GetMaximum())*0.8
    print('treshold= ', treshold)
    for i in range(start, end):
        dif1 = offset + hist.GetBinContent(i) - bslValue
        dif2 = offset + hist.GetBinContent(i+1) - bslValue +  0.03*(offset + hist.GetBinContent(i+1) - bslValue)
        # print('dif1 = ', dif1, "dif2 = ", dif2)
        if dif2 < dif1 and dif1 > treshold and dif2 > treshold:
            height = hist.GetBinContent(i)
            mean = hist.GetBinCenter(i)
            par = np.append(par, height)
            par = np.append(par, mean)
            break
    print("mean = ", mean, 'height = ', height)
    for i in range(hist.FindBin(mean), end):
        h = hist.GetBinContent(i)
        c = hist.GetBinCenter(i)
        uplim = height/2 + 0.15 * (height/2)
        downlim = height/2 - 0.15 * (height/2)
        if downlim < h and h < uplim:
            fhwm = (c - mean) * 2
            sigma = fhwm/2.355 
            par = np.append(par, sigma)
            par = np.append(par, fhwm) 
            break
    print('fhwm = ', fhwm, 'sigma = ', sigma)
    return par

def fit(hist, par, start, end): #fiteaza gausiana cu un vector de parametrii dati
    fit = TF1("fitfunc", "gaus(0) + pol1(3)", start, end)
    print(start, end)
    fit.SetParameters(par[0], par[1], par[2])
    print(fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2))
    fit.SetParLimits(0, par[0] - par[0]*0.75, par[0] + par[0]*0.75)
    fit.SetParLimits(1, par[1] - par[1]*0.75, par[1] + par[1]*0.75)
    fit.SetParLimits(2, par[2] - par[2]*0.75, par[2] + par[2]*0.75)
    hist.Fit(fit, 'RNQ')
    return fit

def integrate(hist, start, end): #O mica functie de integrare pentru binii histogramei
    integration = 0
    for i in range(hist.FindBin(start), hist.FindBin(end)):
        integration += hist.GetBinContent(i)
    return integration

def gaussianArea(func): #Aria gausienei
    mu = func.GetParameter(1)
    sigma = func.GetParameter(2)
    return mu * sigma / 0.3989

def plot(datafiles, infofiles): #Ploeaza histogramele generate din fisierele de tip list pe care le scoate Janus
    holdDelay = getInfo(infofiles)
    histsCh0 = []
    histsCh32 =[]
    Ch0 = TCanvas('c1', "Channel 0")
    Ch32 = TCanvas('c2', 'Channel 32')
    for  i in range(len(datafiles)):
        print(i)
        with codecs.open(datafiles[i], 'r', 'utf-8', errors='ignore') as file:
            # if i == 18:
            #     continue
            histCh0 = TH1F(f"h0_{i}", f"Histogram Channel 0_{i}_ShapingTime_" + holdDelay[i] + "ns", 2**13, 0, 2**13-1)
            histCh32 = TH1F(f"h32_{i}", f"Histogram Channel 32_{i}_ShapingTime" + holdDelay[i] + "ns", 2**13, 0, 2**13-1)
            for i, line in enumerate(file):
                if (i>8):
                    line = line.strip()
                    event = line.split()
                    if event[-3] == '00':
                        histCh0.Fill(int(event[-2]))
                    if event[-3] == '32':
                        histCh32.Fill(int(event[-2]))
            histsCh0.append(histCh0)
            histsCh32.append(histCh32)
            # print(f'Max hist {i} = ', histCh0.GetMaximum(), f'Max hist {i} = ', histCh32.GetMaximum())
    print(len(histCh32), " ", len(histCh0))
    

    file = TFile('hists2.root', 'RECREATE')

    for hist in histsCh0:
        if hist.Integral():
            hist.Write()
    for hist in histsCh32:
        if hist.Integral():
            hist.Write()

    for i in range(len(histCh0)):   
        Ch0.cd()
        histsCh0[i].Draw('hist')
        Ch32.cd()
        histsCh32[i].SetLineColor(2)
        histsCh32[i].Draw('hist')
        gPad.WaitPrimitive('ggg')


def loadHists(files): #Incarca histogamele din fisierele de date
    datafiles = files[0][0]
    infofiles = files[0][1]
    hists = []
    info = getInfo(infofiles)
    with codecs.open(datafiles, 'r', 'utf-8', errors='ignore') as file:
            histCh0 = TH1F("h0_" + info[0], "Histogram Channel 0__ShapingTime_" + info[0] + "ns", 2**13, 0, 2**13-1)
            histCh32 = TH1F("h32_" + info[0],"Histogram Channel 32__ShapingTime" + info[0] + "ns", 2**13, 0, 2**13-1)
            for i, line in enumerate(file):
                if (i>8):
                    line = line.strip()
                    event = line.split()
                    if event[-3] == '00':
                        histCh0.Fill(int(event[-2]))
                    if event[-3] == '32':
                        histCh32.Fill(int(event[-2]))
    hists.append(histCh0)
    hists.append(histCh32)
    return hists


def peakPositions(par): #
    mean = []
    for i in range(len(par)):
        mean.append(par[i])
    return mean 

def calibration(hist,par): #returneaza parametri de calibrare pentru Histograma
    fitPar = np.array([])
    energies = np.array([661.657, 1173.237, 1332.501])
    graph = TGraph()
    print(par)
    for i in range(len(energies)):
        print(par[i], energies[i])
        graph.AddPoint(float(par[i]), float(energies[i]))
    cal = TF1('calibration', 'pol1(0)', par[0], par[2])
    graph.Fit(cal, 'R')
    # graph.Draw('AP')
    # gPad.WaitPrimitive('ggg')
    fitPar = np.append(fitPar, (float(cal.GetParameter(0)), float(cal.GetParameter(1))))
    return fitPar

def uncalAnalysis(files): #Scoate parametrii de fit din histogramele necalibrate 

    datafiles = files[0][0]
    infofiles = files[0][1]

    fileNumber = split(infofiles)
    shapingTime = getInfo(infofiles)

    parCh0 = np.array([])
    parmsCh0Co1 = np.array([])
    parmsCh0Co2 = np.array([])
    parCh32 = np.array([])
    parmsCh32Co1 = np.array([])
    parmsCh32Co2 = np.array([])

    with codecs.open(datafiles, 'r', 'utf-8', errors='ignore') as file:
            pCh0 = np.array([])
            pCh0Co1 = np.array([])
            pCh0Co2 = np.array([])
            pCh32 = np.array([])
            pCh32Co1 = np.array([])
            pCh32Co2 = np.array([])

            histCh0 = TH1F("h0_" + shapingTime[0], "Histogram Channel 0__ShapingTime"  + shapingTime[0] + "ns", 2**13, 0, 2**13-1)
            histCh32 = TH1F("h32_" + shapingTime[0], "Histogram Channel 32__ShapingTime" + shapingTime[0] + "ns", 2**13, 0, 2**13-1)
            for j, line in enumerate(file):
                if (j>8):
                    line = line.strip()
                    event = line.split()
                    if event[-3] == '00':
                        histCh0.Fill(int(event[-2]))
                    if event[-3] == '32':
                        histCh32.Fill(int(event[-2]))
            histCh0Co1 = histCh0.Clone()
            histCh0Co2 = histCh0.Clone()
            histCh32Co1 = histCh32.Clone()
            histCh32Co2 = histCh32.Clone()

            par0 = getPeak(histCh0, histCh0.FindBin(config["ch0"][str(fileNumber[13])]['start']), histCh0.FindBin(config["ch0"][str(fileNumber[13])]['end']))
            parCh0Co1 = getPeak(histCh0Co1, histCh0.FindBin(config["ch0"][str(fileNumber[13])]['startCo1']), histCh0.FindBin(config["ch0"][str(fileNumber[13])]['endCo1']))
            parCh0Co2 = getPeak(histCh0Co2, histCh0.FindBin(config["ch0"][str(fileNumber[13])]['startCo2']), histCh0.FindBin(config["ch0"][str(fileNumber[13])]['endCo2']))

            par32 = getPeak(histCh32, histCh32.FindBin(config["ch32"][str(fileNumber[13])]['start']), histCh32.FindBin(config["ch32"][str(fileNumber[13])]['end']))
            parCh32Co1 = getPeak(histCh32Co1, histCh0.FindBin(config["ch32"][str(fileNumber[13])]['startCo1']), histCh0.FindBin(config["ch32"][str(fileNumber[13])]['endCo1']))
            parCh32Co2 = getPeak(histCh32Co2, histCh0.FindBin(config["ch32"][str(fileNumber[13])]['startCo2']), histCh0.FindBin(config["ch32"][str(fileNumber[13])]['endCo2']))
 

            fitFunc0 = fit(histCh0, par0, config["ch0"][str(fileNumber[13])]['start'], config["ch0"][str(fileNumber[13])]['end'])
            fitfuncCh0Co1 = fit(histCh0Co1, parCh0Co1, config["ch0"][str(fileNumber[13])]['startCo1'], config["ch0"][str(fileNumber[13])]['endCo1'])
            fitfuncCh0Co2 = fit(histCh0Co2, parCh0Co2, config["ch0"][str(fileNumber[13])]['startCo2'], config["ch0"][str(fileNumber[13])]['endCo2'])

            fitFunc32 = fit(histCh32, par32, config["ch32"][str(fileNumber[13])]['start'], config["ch32"][str(fileNumber[13])]['end'])
            fitfuncCh32Co1 = fit(histCh32Co1, parCh32Co1, config["ch32"][str(fileNumber[13])]['startCo1'], config["ch32"][str(fileNumber[13])]['endCo1'])
            fitfuncCh32Co2 = fit(histCh32Co2, parCh32Co2, config["ch32"][str(fileNumber[13])]['startCo2'], config["ch32"][str(fileNumber[13])]['endCo2'])

            # Ch0 = TCanvas('c1', "Channel 0")
            # Ch32 = TCanvas('c2', 'Channel 32')

            # Ch0.cd()
            # histCh0.Draw('hist')
            # histCh0.Fit(fitFunc0, "R")
            # fitFunc0.Draw("same")

            # histCh0Co1.Draw('same')
            # histCh0Co1.Fit(fitfuncCh0Co1, "R")
            # fitfuncCh0Co1.Draw('same')

            # histCh0Co2.Draw('same')
            # histCh0Co2.Fit(fitfuncCh0Co2, "R")
            # fitfuncCh0Co2.Draw('same')

            # Ch32.cd()
            # histCh32.Draw('hist')
            # histCh32.Fit(fitFunc32, "R")
            # fitFunc32.Draw('same')

            # histCh32Co1.Draw('same')
            # histCh32Co1.Fit(fitfuncCh32Co1, "R")
            # fitfuncCh32Co1.Draw('same')

            # histCh32Co2.Draw('same')
            # histCh32Co2.Fit(fitfuncCh32Co2, "R")
            # fitfuncCh32Co2.Draw('same')    
           
            # gPad.WaitPrimitive('ggg')
            
            # rez0 = fitFunc0.GetParameter(2) 
            # rezCh0Co1 =  fitfuncCh0Co1.GetParameter(2)
            # rezCh0Co2 =   fitfuncCh0Co2.GetParameter(2)
            # rez32 = fitFunc32.GetParameter(2) 
            # rezCh32Co1 =  fitfuncCh32Co1.GetParameter(2)
            # rezCh32Co2 =  fitfuncCh32Co2.GetParameter(2)

            peakCsCh0 = fitFunc0.GetParameter(1)
            peakCo1Ch0 = fitfuncCh0Co1.GetParameter(1)
            peakCo2Ch0 = fitfuncCh0Co2.GetParameter(1)
            peakCsCh32 = fitFunc32.GetParameter(1)
            peakCo1Ch32 = fitfuncCh32Co1.GetParameter(1)
            peakCo2Ch32= fitfuncCh32Co2.GetParameter(1)

            pCh0 = np.append(pCh0,( peakCsCh0))
            pCh0Co1 = np.append(pCh0Co1, ( peakCo1Ch0))
            pCh0Co2 = np.append(pCh0Co2, ( peakCo2Ch0))

            print(pCh0)

            pCh32 = np.append(pCh32, ( peakCsCh32))
            pCh32Co1 = np.append(pCh32Co1, ( peakCo1Ch32))
            pCh32Co2 = np.append(pCh32Co2, ( peakCo2Ch32))

            print(pCh32)
            
            parCh0 = np.append(parCh0, pCh0)
            parmsCh0Co1 = np.append(parmsCh0Co1, pCh0Co1)
            parmsCh0Co2 = np.append(parmsCh0Co2, pCh0Co2)

            parCh32 = np.append(parCh32, pCh32)
            parmsCh32Co1 = np.append(parmsCh32Co1, pCh32Co1)
            parmsCh32Co2 = np.append(parmsCh32Co2, pCh32Co2)
        
    print(parCh32)
    print(parCh0)
    print(parmsCh0Co1)
    print(parmsCh0Co2)
    print(parmsCh32Co1)
    print(parmsCh32Co2)

    return parCh0, parmsCh0Co1, parmsCh0Co2, parCh32, parmsCh32Co1, parmsCh32Co2 

def histCalibration(files): #Calibreaza histogramele si returneaza parametrii de calibrare doriti
   

    datafiles = files[0][0]
    infofiles = files[0][1]
    
    print(datafiles)
    print(infofiles)

    hists = loadHists(files)
    histsCh0 = hists[0]
    histsCh32 = hists[1]

    parameters = uncalAnalysis(files)
    info = getInfo(infofiles)
    
    file = TFile("calibratedHists.root", "RECREATE")

    finalParameters = np.array([])

    pCh0Cs = np.array([])
    pCh0Co1 = np.array([])
    pCh0Co2 = np.array([])
    pCh32Cs = np.array([])
    pCh32Co1 = np.array([])
    pCh32Co2 = np.array([])

    histCh0Cal = TH1F("h0_" + info[0], "Histogram Channel 0__ShapingTime_" + info[0] + "ns_Calibrated", 2**13, 0, 2**13-1)
    histCh32Cal = TH1F("h32_" +info[0], "Histogram Channel 32__ShapingTime" + info[0] + "ns_Calibrated", 2**13, 0, 2**13-1)

    ch0Data = [parameters[0], parameters[1], parameters[2]]
    ch32Data =[parameters[3], parameters[4], parameters[5]]

    parCh0 = calibration(histsCh0, ch0Data)
    parCh32 = calibration(histsCh32, ch32Data)

    for j in range(0, histsCh0.GetNbinsX()):
         for k in range(int(histsCh0.GetBinContent(j))):
                histCh0Cal.Fill((j + rand.uniform(-0.5, 0.5))*parCh0[1] + parCh0[0])
    # histCh0Cal.Draw('hist')
    # gPad.WaitPrimitive('ggg')
    histCh0Cal.Write()
  
    for j in range(0, histsCh32.GetNbinsX()):
        for k in range(int(histsCh32.GetBinContent(j))):
                histCh32Cal.Fill((j + rand.uniform(-0.5, 0.5))*parCh32[1] + parCh32[0])
    # histCh32Cal.Draw('hist')
    # gPad.WaitPrimitive('ggg')

    histCh32Cal.Write()

    parCh0Cs = getPeak(histCh0Cal, histCh0Cal.FindBin(config['startCsCal']), histCh0Cal.FindBin(config['endCsCal']))
    parCh0Co1 = getPeak(histCh0Cal, histCh0Cal.FindBin(config['startCo1Cal']), histCh0Cal.FindBin(config['endCo1Cal']))
    parCh0Co2 = getPeak(histCh0Cal, histCh0Cal.FindBin(config['startCo2Cal']), histCh0Cal.FindBin(config['endCo2Cal']))

    parCh32Cs = getPeak(histCh32Cal, histCh32Cal.FindBin(config['startCsCal']), histCh32Cal.FindBin(config['endCsCal']))
    parCh32Co1 = getPeak(histCh32Cal, histCh32Cal.FindBin(config['startCo1Cal']), histCh32Cal.FindBin(config['endCo1Cal']))
    parCh32Co2 = getPeak(histCh32Cal, histCh32Cal.FindBin(config['startCo2Cal']), histCh32Cal.FindBin(config['endCo2Cal']))

    histCh0Cal.GetXaxis().SetRange(0,0)
    fitFunc0 = fit(histCh0Cal, parCh0Cs, config['startCsCal'], config['endCsCal'])
    fitfuncCh0Co1 = fit(histCh0Cal, parCh0Co1, config['startCo1Cal'], config['endCo1Cal'])
    fitfuncCh0Co2 = fit(histCh0Cal, parCh0Co2, config['startCo2Cal'], config['endCo2Cal'])

    histCh32Cal.GetXaxis().SetRange(0,0)
    fitFunc32 = fit(histCh32Cal, parCh32Cs, config['startCsCal'], config['endCsCal'])
    fitfuncCh32Co1 = fit(histCh32Cal, parCh32Co1, config['startCo1Cal'], config['endCo1Cal'])
    fitfuncCh32Co2 = fit(histCh32Cal, parCh32Co2, config['startCo2Cal'], config['endCo2Cal'])      

    # Ch0 = TCanvas('c1', "Channel 0")
    # Ch32 = TCanvas('c2', 'Channel 32')

    # Ch0.cd()
    # histCh0Cal.Draw('hist')
    # histCh0Cal.Fit(fitFunc0, "R+")
    # # gPad.WaitPrimitive("ggg")
    # fitFunc0.Update()
    # fitFunc0.Draw("same")
    # histCh0Cal.Fit(fitfuncCh0Co1, "R+")
    # # gPad.WaitPrimitive("ggg")
    # fitfuncCh0Co1.Update()
    # fitfuncCh0Co1.Draw("same")
    # histCh0Cal.Fit(fitfuncCh0Co2, "R+")
    # # gPad.WaitPrimitive("ggg")

    # fitfuncCh0Co2.Update()
    # fitfuncCh0Co2.Draw("same")
    # # gPad.WaitPrimitive("ggg")

    # Ch32.cd()
    # histCh32Cal.Draw('hist')
    # histCh32Cal.Fit(fitFunc32, "R+")
    # # gPad.WaitPrimitive("ggg")
    # fitFunc32.Update()
    # fitFunc32.Draw('same')
    # histCh32Cal.Fit(fitfuncCh32Co1,"R+")
    # # gPad.WaitPrimitive("ggg")
    # fitfuncCh32Co1.Update()
    # fitfuncCh32Co1.Draw("same")
    # histCh32Cal.Fit(fitfuncCh0Co2,"R=")
    # # gPad.WaitPrimitive("ggg")
    # fitfuncCh0Co2.Update()
    # fitfuncCh32Co2.Draw("same")
    # # gPad.WaitPrimitive("ggg")
             
    rez0 = fitFunc0.GetParameter(2)
    rezCh0Co1 =  fitfuncCh0Co1.GetParameter(2)
    rezCh0Co2 =   fitfuncCh0Co2.GetParameter(2)
    rez32 = fitFunc32.GetParameter(2)
    rezCh32Co1 =  fitfuncCh32Co1.GetParameter(2)
    rezCh32Co2 =  fitfuncCh32Co2.GetParameter(2)
    print(rez0, rez32)

    pCh0Cs = np.append(pCh0Cs,(rez0))
    pCh0Co1 = np.append(pCh0Co1, (rezCh0Co1))
    pCh0Co2 = np.append(pCh0Co2, (rezCh0Co2))

    pCh32Cs = np.append(pCh32Cs, (rez32))
    pCh32Co1 = np.append(pCh32Co1, (rezCh32Co1))
    pCh32Co2 = np.append(pCh32Co2, (rezCh32Co2))

    finalParameters = np.append(finalParameters, (pCh0Cs, pCh0Co1, pCh0Co2, pCh32Cs, pCh32Co1, pCh32Co2))
    return finalParameters 


def plotGraphs(infoFilesLocation, parameters): #Ploteaza mai multe grafice 
  

    infoFiles = [os.path.relpath(x) for x in os.scandir(infoFilesLocation)]
    infoFiles = natsort.natsorted(infoFiles)
    info = np.array([])
    for i in range(len(infoFiles)):
        with codecs.open(infoFiles[i], 'r', 'utf-8', errors='ignore') as file:
            for i, line in enumerate(file):
                line = line.strip()
                infoL = line.split()
                if(len(infoL) == 0):
                 i = i+1
                else:
                    if infoL[0] == 'LG_ShapingTime':
                        info = np.append(info, infoL[1])
    
    c1 = TCanvas('graphs', 'Energy Rezsolution vs Shaping Time Ch0')
    mg1 = TMultiGraph("graphs", "Energy Resolution vs Shaping Time")
    c2 = TCanvas('graphs32', 'Energy Resolution vs Shaping Time Ch32')
    graph1 = TGraph()
    graph2 = TGraph()
    graph3 = TGraph()
    graph4 = TGraph()
    graph5 = TGraph()
    graph6 = TGraph()

    graph1.SetName('Energy Resolution vs Shaping Time Ch0Cs')
    graph1.GetXaxis().SetTitle("Shaping Time")
    graph1.GetYaxis().SetTitle("Energy Resolution")
    graph1.SetMarkerSize(1.2)
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerColor(kBlack)

    graph2.SetName('Energy Resolution vs Shaping Time Ch32Cs')
    graph2.SetMarkerSize(1.2)
    graph2.SetMarkerStyle(20)
    graph2.SetMarkerColor(kRed)

    graph3.SetName('Energy Resolution vs Shaping Time Ch0Co1')
    graph3.SetMarkerSize(1.2)
    graph3.SetMarkerStyle(20)
    graph3.SetMarkerColor(kBlue)

    graph4.SetName('Energy Resolution vs Shaping Time Ch0Co2')
    graph4.SetMarkerSize(1.2)
    graph4.SetMarkerStyle(20)
    graph4.SetMarkerColor(kGreen)

    graph5.SetName('Energy Resolution vs Shaping Time Ch32Co1')
    graph5.SetMarkerSize(1.2)
    graph5.SetMarkerStyle(20)
    graph5.SetMarkerColor(kMagenta)

    graph6.SetName('Energy Resolution vs Shaping Time Ch32Co2')
    graph6.SetMarkerSize(1.2)
    graph6.SetMarkerStyle(20)
    graph6.SetMarkerColor(kOrange)

    for i in range(len(info)):
        graph1.AddPoint(float(info[i]), float(parameters[i][0]))
        graph2.AddPoint(float(info[i]), float(parameters[i][1]))
        graph3.AddPoint(float(info[i]), float(parameters[i][2]))
        graph4.AddPoint(float(info[i]), float(parameters[i][3]))
        graph5.AddPoint(float(info[i]), float(parameters[i][4]))
        graph6.AddPoint(float(info[i]), float(parameters[i][5]))

    file = TFile('graphs.root', 'RECREATE')
    legend = TLegend(0.1,0.7,0.48,0.9)
    legend.SetHeader("Legend", "C")
    legend.AddEntry(graph1, "Energy Resolution vs Shaping Time Ch0Cs", 'p')
    legend.AddEntry(graph2, 'Energy Resolution vs Shaping Time Ch32Cs', 'p')
    legend.AddEntry(graph3, "Energy Resolution vs Shaping Time Ch0Co1", 'p')
    legend.AddEntry(graph4, 'Energy Resolution vs Shaping Time Ch0Co2', 'p')
    legend.AddEntry(graph5, "Energy Resolution vs Shaping Time Ch32Co1", 'p')
    legend.AddEntry(graph6, 'Energy Resolution vs Shaping Time Ch32Co2', 'p')


    mg1.Add(graph2)
    mg1.Add(graph1)
    mg1.Add(graph3)
    mg1.Add(graph4)
    mg1.Add(graph5)
    mg1.Add(graph6)

    mg1.Draw('ap')
    legend.Draw()
    c1.Write()
    gPad.WaitPrimitive('ggg')


# plot(config['dataFiles'],config['infoFiles'])

configfile = 'config.json'
with open(configfile) as datafile:
    config = json.load(datafile)

infoFiles = [os.path.relpath(x) for x in os.scandir(config['infoFiles'])]
infoFiles = natsort.natsorted(infoFiles)

dataFiles = [os.path.relpath(x) for x in os.scandir(config['dataFiles'])]
dataFiles = natsort.natsorted(dataFiles)
parallelFiles = mergeLists(dataFiles, infoFiles)


parallelFiles = [[parallelFile] for parallelFile in parallelFiles]
print(parallelFiles)
time_start = time.perf_counter()

with Pool() as pool:
    parameters = pool.map(histCalibration, parallelFiles)
end_time = time.perf_counter()
plotGraphs(config["infoFiles"], parameters)
print("Program finished in {} seconds".format(end_time - time_start))



    

