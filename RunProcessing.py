from ROOT import TF1, TH1F, gPad, TCanvas, TFile, TGraph, kRed, kBlack, TMultiGraph, TLegend, kBlue, kGreen, kMagenta, kOrange
import numpy as np
import array as arr
import codecs
import math
import scipy 
import json
import os 
import natsort



configfile = 'config.json'
with open(configfile) as datafile:
    config = json.load(datafile)

def getInfo(files): #Scoate informatiile de care ai nevoie din fisierel e tip info pe care le scoate Janus
    infoFiles = [os.path.relpath(x) for x in os.scandir(files)]
    infoFiles = natsort.natsorted(infoFiles)
    holdDelay = np.array([])
    for i in range(len(infoFiles)):
        with codecs.open(infoFiles[i], 'r', 'utf-8', errors='ignore') as file:
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
        dif2 = offset + hist.GetBinContent(i+1) - bslValue +  0.05*(offset + hist.GetBinContent(i+1) - bslValue)
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
    fit.SetParameters(par[0], par[1], par[2])
    fit.SetParLimits(0, par[0] - par[0]*0.75, par[0] + par[0]*0.75)
    fit.SetParLimits(1, par[1] - par[1]*0.75, par[1] + par[1]*0.75)
    fit.SetParLimits(2, par[2] - par[2]*0.75, par[2] + par[2]*0.75)
    hist.Fit(fit, 'RN')
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
    dataFiles = [os.path.relpath(x) for x in os.scandir(datafiles)]
    dataFiles = natsort.natsorted(dataFiles)
    Ch0 = TCanvas('c1', "Channel 0")
    Ch32 = TCanvas('c2', 'Channel 32')
    for  i in range(len(dataFiles)):
        print(i)
        with codecs.open(dataFiles[i], 'r', 'utf-8', errors='ignore') as file:
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

def analyseFiles(datafiles, infofiles): #Scoate parametrii doriti din fisierele de tip list de la Janus
    holdDelay = getInfo(infofiles)
    parCh0 = np.array([])
    # parmsCh0Co1 = np.array([])
    # parmsCh0Co2 = np.array([])
    parCh32 = np.array([])
    # parmsCh32Co1 = np.array([])
    # parmsCh32Co2 = np.array([])

    histsCh0 = []
    histsCh32 =[]
    dataFiles = [os.path.relpath(x) for x in os.scandir(datafiles)]
    dataFiles = natsort.natsorted(dataFiles)
    for  i in range(len(dataFiles)):
        print(i)
        if i ==18:
            continue
        with codecs.open(dataFiles[i], 'r', 'utf-8', errors='ignore') as file:
            pCh0 = np.array([])
            # pCh0Co1 = np.array([])
            # pCh0Co2 = np.array([])
            pCh32 = np.array([])
            # pCh32Co1 = np.array([])
            # pCh32Co2 = np.array([])

            histCh0 = TH1F(f"h0_{i}", f"Histogram Channel 0_{i}_HoldDelay"  + holdDelay[i], 2**13, 0, 2**13-1)
            histCh32 = TH1F(f"h32_{i}", f"Histogram Channel 32_{i}_HoldDelay" + holdDelay[i], 2**13, 0, 2**13-1)
            for i, line in enumerate(file):
                if (i>8):
                    line = line.strip()
                    event = line.split()
                    if event[-3] == '00':
                        histCh0.Fill(int(event[-2]))
                    if event[-3] == '32':
                        histCh32.Fill(int(event[-2]))
            # histCh0Co1 = histCh0.Clone()
            # histCh0Co2 = histCh0.Clone()
            # histCh32Co1 = histCh32.Clone()
            # histCh32Co2 = histCh32.Clone()

            par0 = getPeak(histCh0, histCh0.FindBin(config['start']), histCh0.FindBin(config['end']))
            par32 = getPeak(histCh32, histCh32.FindBin(config['start']), histCh32.FindBin(config['end']))
            # parCh0Co1 = getPeak(histCh0Co1, histCh0.FindBin(config['startCo1']), histCh0.FindBin(config['endCo1']))
            # parCh0Co2 = getPeak(histCh0Co2, histCh0.FindBin(config['startCo2']), histCh0.FindBin(config['endCo2']))
            # parCh32Co1 = getPeak(histCh32Co1, histCh0.FindBin(config['startCo1']), histCh0.FindBin(config['endCo1']))
            # parCh32Co2 = getPeak(histCh32Co2, histCh0.FindBin(config['startCo2']), histCh0.FindBin(config['endCo2']))
 
            print(par0)
            print(par32)
            # print(parCh0Co1)
            # print(parCh0Co2)
            # print(parCh32Co1)
            # print(parCh32Co2)

            fitFunc0 = fit(histCh0, par0, config['start'], config['end'])
            fitFunc32 = fit(histCh32, par32, config['start'], config['end'])
            # fitfuncCh0Co1 = fit(histCh0Co1, parCh0Co1, config['startCo1'], config['endCo1'])
            # fitfuncCh0Co2 = fit(histCh0Co2, parCh0Co2, config['startCo2'], config['endCo2'])
            # fitfuncCh32Co1 = fit(histCh32Co1, parCh32Co1, config['startCo1'], config['endCo1'])
            # fitfuncCh32Co2 = fit(histCh32Co2, parCh32Co2, config['startCo2'], config['endCo2'])

            totalCounts0 = integrate(histCh0, config['specStart'], config['specEnd'])
            totalCounts32 = integrate(histCh0, config['specStart'], config['specEnd'])

            Ch0 = TCanvas('c1', "Channel 0")
            Ch32 = TCanvas('c2', 'Channel 32')

            Ch0.cd()
            histCh0.Draw('hist')
            histCh0.Fit(fitFunc0, "R")
            fitFunc0.Draw("same")

            # histCh0Co1.Draw('same')
            # histCh0Co1.Fit(fitfuncCh0Co1, "R")
            # fitfuncCh0Co1.Draw('same')

            # histCh0Co2.Draw('same')
            # histCh0Co2.Fit(fitfuncCh0Co2, "R")
            # fitfuncCh0Co2.Draw('same')

            Ch32.cd()
            histCh32.Draw('hist')
            histCh32.Fit(fitFunc32, "R")
            fitFunc32.Draw('same')

            # histCh32Co1.Draw('same')
            # histCh32Co1.Fit(fitfuncCh32Co1, "R")
            # fitfuncCh32Co1.Draw('same')

            # histCh32Co2.Draw('same')
            # histCh32Co2.Fit(fitfuncCh32Co2, "R")
            # fitfuncCh32Co2.Draw('same')    
           
            gPad.WaitPrimitive('ggg')
            
            rez0 = fitFunc0.GetParameter(2) 
            # rezCh0Co1 =  fitfuncCh0Co1.GetParameter(2)
            # rezCh0Co2 =   fitfuncCh0Co2.GetParameter(2)
            rez32 = fitFunc32.GetParameter(2) 
            # rezCh32Co1 =  fitfuncCh32Co1.GetParameter(2)
            # rezCh32Co2 =  fitfuncCh32Co2.GetParameter(2)

           

           

            pCh0 = np.append(pCh0, rez0)
            # pCh0Co1 = np.append(pCh0Co1, rezCh0Co1)
            # pCh0Co2 = np.append(pCh0Co2, rezCh0Co2)

            print(pCh0)

            pCh32 = np.append(pCh32, rez32)
            # pCh32Co1 = np.append(pCh32Co1, rezCh32Co1)
            # pCh32Co2 = np.append(pCh32Co2, rezCh32Co2)

            print(pCh32)
            
            parCh0 = np.append(parCh0, pCh0)
            # parmsCh0Co1 = np.append(parmsCh0Co1, pCh0Co1)
            # parmsCh0Co2 = np.append(parmsCh0Co2, pCh0Co2)

            parCh32 = np.append(parCh32, pCh32)
            # parmsCh32Co1 = np.append(parmsCh32Co1, pCh32Co1)
            # parmsCh32Co2 = np.append(parmsCh32Co2, pCh32Co2)

            histsCh0.append(histCh0)
            histsCh32.append(histCh32)
        
    print(parCh32)
    print(parCh0)
    # print(parmsCh0Co1)
    # print(parmsCh0Co2)
    # print(parmsCh32Co1)
    # print(parmsCh32Co2)

    return parCh0, parCh32
# , parmsCh0Co1, parmsCh0Co2, parmsCh32Co1, parmsCh32Co2 

def plotGraphs(parms1, parms2, parms3): #Ploteaza mai multe grafice 
  
    c1 = TCanvas('graphs', 'Hold Delay vs Sigma graphs')
    mg1 = TMultiGraph("graphs", "Hold Delay vs Sigma")
    # c2 = TCanvas('graphs32', 'Hold Delay vs Sigma graphs Ch32')
    graph1 = TGraph()
    graph2 = TGraph()
    # graph3 = TGraph()
    # graph4 = TGraph()
    # graph5 = TGraph()
    # graph6 = TGraph()

    graph1.SetName('Hold Delay vs Sigma Ch0Cs')
    graph1.GetXaxis().SetTitle("Hold Delay")
    graph1.GetYaxis().SetTitle("Sigma")
    graph1.SetMarkerSize(1.2)
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerColor(kBlack)

    graph2.SetName('Hold Delay vs Sigma Ch32Cs')
    graph2.SetMarkerSize(1.2)
    graph2.SetMarkerStyle(20)
    graph2.SetMarkerColor(kRed)

    # graph3.SetName('Hold Delay vs Sigma Ch0Co1')
    # graph3.SetMarkerSize(1.2)
    # graph3.SetMarkerStyle(20)
    # graph3.SetMarkerColor(kBlue)

    # graph4.SetName('Hold Delay vs Sigma Ch0Co2')
    # graph4.SetMarkerSize(1.2)
    # graph4.SetMarkerStyle(20)
    # graph4.SetMarkerColor(kGreen)

    # graph5.SetName('Hold Delay vs Sigma Ch32Co1')
    # graph5.SetMarkerSize(1.2)
    # graph5.SetMarkerStyle(20)
    # graph5.SetMarkerColor(kMagenta)

    # graph6.SetName('Hold Delay vs Sigma Ch32Co2')
    # graph6.SetMarkerSize(1.2)
    # graph6.SetMarkerStyle(20)
    # graph6.SetMarkerColor(kOrange)

    for i in range(len(parms1)):
        graph1.AddPoint(float(parms1[i]), float(parms2[i]))
        graph2.AddPoint(float(parms1[i]), float(parms3[i]))
        # graph3.AddPoint(float(parms1[i]), float(parms4[i]))
        # graph4.AddPoint(float(parms1[i]), float(parms5[i]))
        # graph5.AddPoint(float(parms1[i]), float(parms6[i]))
        # graph6.AddPoint(float(parms1[i]), float(parms7[i]))
    file = TFile('graphs.root', 'RECREATE')
    legend = TLegend(0.1,0.7,0.48,0.9)
    legend.SetHeader("Legend", "C")
    legend.AddEntry(graph1, "Hold Delay vs Sigma Ch0Cs", 'p')
    legend.AddEntry(graph2, 'Hold Delay vs Sigma Ch32Cs', 'p')
    # legend.AddEntry(graph3, "Hold Delay vs Sigma Ch0Co1", 'p')
    # legend.AddEntry(graph4, 'Hold Delay vs Sigma Ch0Co2', 'p')
    # legend.AddEntry(graph5, "Hold Delay vs Sigma Ch32Co1", 'p')
    # legend.AddEntry(graph6, 'Hold Delay vs Sigma Ch32Co2', 'p')


    mg1.Add(graph2)
    mg1.Add(graph1)
    # mg1.Add(graph3)
    # mg1.Add(graph4)
    # mg1.Add(graph5)
    # mg1.Add(graph6)

    mg1.Draw('ap')
    legend.Draw()
    c1.Write()
    gPad.WaitPrimitive('ggg')
# print(getInfo(config['infoFiles']))
plot(config['dataFiles'],config['infoFiles'])
# print(getInfo(config['infoFiles']))
# info = getInfo2(config['infoFiles'])
# parameters = analyseFiles(config['dataFiles'], config['infoFiles'])
# plotGraphs(info, parameters[0], parameters[1])




 



    

