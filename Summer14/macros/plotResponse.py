#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

import ROOT
#from ROOT import TMath

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);

############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-i','--input',action="store",type="string",dest="input",default="outtre.root")
parser.add_option('-o','--outdir',action="store",type="string",dest="outdir",default="plots")
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=40)
parser.add_option('-r',action="store",type="float",dest="radius",default=0.8)
parser.add_option('--minPt',action="store",type="float",dest="minPt",default=25.)
parser.add_option('--maxPt',action="store",type="float",dest="maxPt",default=200.)

(options, args) = parser.parse_args()

############################################################
def makeTexts():
    # text                                                                                                                                                                                                                 
    latex1 = ROOT.TLatex(0.20,0.89,("Anti-kT (R=%.1f)"%(options.radius)))
    latex1.SetNDC()
    latex1.SetTextSize(0.03)
    latex2 = ROOT.TLatex(0.20,0.84,("n_{PU} = "+str(options.nPU)))
    latex2.SetNDC()
    latex2.SetTextSize(0.03)
    latex3 = ROOT.TLatex(0.20,0.79,("%.0f GeV < p_{T} < %.0f GeV "%(options.minPt, options.maxPt)))
    latex3.SetNDC()
    latex3.SetTextSize(0.03)
    return latex1, latex2, latex3


def makeTrendResponse(f, types, xvar, yvar, styles, rebin, outdir):

    legend = ROOT.TLegend(0.7,0.7,0.93,0.9);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);

    xtitle = 'N_{PU}'
    if xvar == 'pt':
        xtitle = 'p^{T} (GeV)'
    if xvar == 'eta':
        xtitle = '#eta'

    for ivar,var in enumerate(yvar):
        c1 = ROOT.TCanvas(var+'_mean_vs_'+xvar,var+'_mean_vs_'+xvar,700,700);
        c2 = ROOT.TCanvas(var+'_rms_vs_'+xvar,var+'_rms_vs_'+xvar,700,700);
        h = []
        graphmean = []
        graphrms = []

        maxmean = []
        maxrms  = []
        minmean = []
        minrms = []
        
        if 'pt' not in var:
            meanytitle = '<%s - gen m> (GeV)'%var
            rmsytitle = 'RMS(%s - gen m) (GeV)'%var
        else:
            meanytitle = '<(%s - ptgen)/ptgen>'%var
            rmsytitle = 'RMS((%s - ptgen)/ptgen)'%var

        n = 0
        for typ, suff in types.iteritems():

            h.append( f.Get(suff+'/h'+var+'_response_vs_'+xvar+'_'+suff))
            
            graphmean.append(ROOT.TGraphErrors())
            graphrms.append(ROOT.TGraphErrors())
            graphmean[n].SetLineColor(styles[typ][0])
            graphrms[n].SetLineColor(styles[typ][0])
            graphmean[n].SetLineWidth(2)
            graphrms[n].SetLineWidth(2)
            graphmean[n].SetMarkerColor(styles[typ][0])
            graphrms[n].SetMarkerColor(styles[typ][0])
            graphmean[n].SetMarkerStyle(21)
            graphrms[n].SetMarkerStyle(21)
            graphmean[n].SetName(suff)
            graphrms[n].SetName(suff)

            if (ivar == 0):
                legend.AddEntry(graphmean[n],typ,'L')

            j = 0
            px = (h[n].ProjectionX('px')).Clone('px') 
            for bin in range(1,h[n].GetNbinsX()+1,rebin): # rebin
                firstbin = bin
                lastbin  = bin+(rebin-1)
                py       = (h[n].ProjectionY('py',firstbin,lastbin)).Clone('py')
                mean     = py.GetMean()
                meanerr  = py.GetMeanError()
                rms      = py.GetRMS()
                rmserr   = py.GetRMSError()
                if (py.GetEntries()>0):
                    j = j + 1 
                    #px.GetXaxis().SetRange(firstbin,lastbin)    
                    #x = px.GetMean()  
                    #xx = (px.GetXaxis().GetBinUpEdge(lastbin) - px.GetXaxis().GetBinLowEdge(firstbin))/2
                    x  = (px.GetXaxis().GetBinUpEdge(lastbin) + px.GetXaxis().GetBinLowEdge(firstbin))/2
                    xx = (px.GetXaxis().GetBinUpEdge(lastbin) - px.GetXaxis().GetBinLowEdge(firstbin))/2
                    graphmean[n].SetPoint(j,x,mean)
                    graphmean[n].SetPointError(j,xx,meanerr)
                    #graphrms[n] .SetPoint(j,x, rms)
                    #graphrms[n] .SetPointError(j,xx,rmserr)
                    if ('pt' in var):
                        graphrms[n] .SetPoint(j,x, rms/(1+mean)) # take response corrected RMS (N.B. pt histograms are filled with pt/genpt-1 
                        graphrms[n] .SetPointError(j,xx,rmserr/(1+mean))   
                    else:
                        graphrms[n] .SetPoint(j,x, rms) # NB: mass plots are filled with m -genm
                        graphrms[n] .SetPointError(j,xx,rmserr)

            #imaxmean = ROOT.TMath.LocMax(graphmean[n].GetN(),graphmean[n].GetY())
            #iminmean = ROOT.TMath.LocMin(graphmean[n].GetN(),graphmean[n].GetY())
            #imaxrms  = ROOT.TMath.LocMax(graphrms[n].GetN(),graphrms[n].GetY())
            #iminrms  = ROOT.TMath.LocMin(graphrms[n].GetN(),graphrms[n].GetY())            
            
            #maxmean.append(graphmean[n].GetY()[imaxmean])
            #minmean.append(graphmean[n].GetY()[iminmean])
            #maxrms.append(graphrms[n].GetY()[imaxrms])
            #minrms.append(graphrms[n].GetY()[iminrms])            

            n = n + 1

        # set axis range and plot
        for i in range(0,n):
            
            if (i == 0):
                #graphmean[i].SetMinimum(min(minmean)-10)
                #graphmean[i].SetMaximum(max(maxmean)+10)
                #graphrms[i].SetMinimum(min(minrms)-10)
                #graphrms[i].SetMaximum(max(maxrms)+10)

                if 'pt' in xvar:
                    graphmean[i].GetHistogram().GetXaxis().SetRangeUser(options.minPt,options.maxPt)
                    graphrms[i].GetHistogram().GetXaxis().SetRangeUser(options.minPt,options.maxPt)

                if ('pt' in var):
                    graphmean[i].SetMinimum(-0.5)
                    graphmean[i].SetMaximum(0.5)
                    graphrms[i].SetMinimum(0)
                    graphrms[i].SetMaximum(0.5)
                else:
                    graphmean[i].SetMinimum(-20)
                    graphmean[i].SetMaximum(30)
                    graphrms[i].SetMinimum(0)
                    graphrms[i].SetMaximum(30)

                graphmean[i].GetHistogram().SetXTitle(xtitle)
                graphrms[i].GetHistogram().SetXTitle(xtitle)
                graphmean[i].GetHistogram().SetYTitle(meanytitle)
                graphrms[i].GetHistogram().SetYTitle(rmsytitle)

                c1.cd()
                graphmean[i].Draw("ap")
                c2.cd()
                graphrms[i].Draw("ap")
            else:
                c1.cd()
                graphmean[i].Draw("p*same")
                c2.cd()
                graphrms[i].Draw("p*same")
            
                
        # save plots
        text1, text2, text3 = makeTexts()
        for c in c1,c2:
            c.SetGridx()
            c.SetGridy()
            c.cd()
            text1.Draw()
            text2.Draw()
            text3.Draw()
            legend.Draw()
            c.Modified()
            c.Update()

        for p in '.pdf', '.png','.root':
            c1.SaveAs(outdir+'/'+c1.GetName()+p)
            c2.SaveAs(outdir+'/'+c2.GetName()+p)


if __name__ == '__main__':

    filename = options.input
    outdir   = options.outdir
    try:
        os.mkdir(outdir)
    except:
        print 'Cannot create output directory: directory already exists'
        sys.exit()

    docmssw = False
    
    types = {'PUPPI':'puppi','PFlow':'pf','PFlowCHS':'pfchs'}
    if (docmssw):
        types = {'PUPPI':'puppi','PFlow':'pf','PFlowCHS':'pfchs','PF-CMSSW':'pfcmssw'}

    styles = {} # color, linestyle, line width
    styles['GEN'] = [ROOT.kBlack, 1, 2]
    styles['PUPPI'] = [ROOT.kGreen+1, 1, 2]
    styles['PFlow'] = [ROOT.kBlue, 1, 2]
    styles['PFlowCHS'] = [ROOT.kMagenta, 1, 2]
    styles['PF-CMSSW'] = [ROOT.kOrange, 1, 2]

    
    # -- open file
    f = ROOT.TFile.Open(filename);

    # -- make plots 
    masses = ['mraw','m','mtrim','mtrimsafe','mconst']

    makeTrendResponse(f, types, 'npu', masses, styles, 5 , options.outdir)
    makeTrendResponse(f, types, 'eta', masses, styles, 5 , options.outdir)
    makeTrendResponse(f, types, 'pt' , masses, styles, 10, options.outdir)

    pts = ['ptraw','pt','ptcorr']
    makeTrendResponse(f, types, 'npu', pts, styles, 5 , options.outdir)
    makeTrendResponse(f, types, 'eta', pts, styles, 5 , options.outdir)
    makeTrendResponse(f, types, 'pt' , pts, styles, 10, options.outdir)


    raw_input('ok?')

        
        
