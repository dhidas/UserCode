# A version of plotMttbar.py that has code for plotting the top-like histograms.  This file probably requires 
# modification before use.

from __future__ import division

from tdrStyle import *
from ROOT import *

canvases = []
scanvases = []
counter = 0

def plotMttbar():
    tdrstyle = setTDRStyle();
    gStyle.SetHatchesSpacing(1.0);
    # lumi = 36.135;#3.06 + 1.37 + 2.53 + 0.26;
    lumi = 203.8;
    oldLumi = 203.8;
    # oldLumi = 36.135;
    scale = lumi / oldLumi;
    qcdScale = {'default':1.72, 'withMETAndAsymJets': 3.03};
    #    data = TFile.Open("data2.9pb_fullSetOfVars.root");
    #data =	TFile.Open("Run2010_35pb_PFElectron_PF2PATJets_PFMET.root");
    data =	TFile.Open("Run_203.8pb_PFElectron_PF2PATJets_PFMET.root");
    ttbar =	TFile.Open("TT_Tune_203.8pb_PFElectron_PF2PATJets_PFMET.root");
    WprimeTToTTD_M1000 = TFile.Open("WprimeTToTTD_M1000_203.8pb_PFElectron_PF2PATJets_PFMET.root");
    WprimeTToTTD_M600 = TFile.Open("WprimeTToTTD_M600_203.8pb_PFElectron_PF2PATJets_PFMET.root");
    # ttbar =	TFile.Open("TTJets_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # wjets =	TFile.Open("WJets_35pb_PFElectron_PF2PATJets_PFMET.root");
    # zjets =	TFile.Open("DYJets_35pb_PFElectron_PF2PATJets_PFMET.root");
    # bce1 =	TFile.Open("QCD_Pt-20to30_BCtoE_35pb_PFElectron_PF2PATJets_PFMET.root");
    # bce2 =	TFile.Open("QCD_Pt-30to80_BCtoE_35pb_PFElectron_PF2PATJets_PFMET.root");
    # bce3 =	TFile.Open("QCD_Pt-80to170_BCtoE_35pb_PFElectron_PF2PATJets_PFMET.root");
    # enri1 =	TFile.Open("QCD_Pt-20to30_EMEnriched_35pb_PFElectron_PF2PATJets_PFMET.root");
    # enri2 =	TFile.Open("QCD_Pt-30to80_EMEnriched_35pb_PFElectron_PF2PATJets_PFMET.root");
    # enri3 =	TFile.Open("QCD_Pt-80to170_EMEnriched_35pb_PFElectron_PF2PATJets_PFMET.root");
    # zjets = TFile.Open("DYJets_35pb_PFElectron_PF2PATJets_PFMET.root");
    # pj1 = TFile.Open("pj1_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # pj2 = TFile.Open("pj2_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # pj3 = TFile.Open("pj3_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # tW = TFile.Open("tW_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # tchan = TFile.Open("tchan_36.145pb_PFElectron_PF2PATJets_PFMET.root");
#    vqq = TFile.Open("vqq_7.22pb_V4PFCalo.root__fullSetOfVars.root");
#    Zprime500 = TFile.Open("/storage/workspace/BristolAnalysisTools/outputfiles/Zprime_M500GeV_W5GeV_36.135pb.root");
#    Zprime750 = TFile.Open("/storage/workspace/BristolAnalysisTools/outputfiles/Zprime_M750GeV_W7500MeV_36.135pb.root");
#    Zprime1000 = TFile.Open("/storage/workspace/BristolAnalysisTools/outputfiles/Zprime_M1TeV_W10GeV_36.135pb.root");
#    Zprime1250 = TFile.Open("/storage/workspace/BristolAnalysisTools/outputfiles/Zprime_M1250GeV_W12500MeV_36.135pb.root");
#    Zprime1500 = TFile.Open("/storage/workspace/BristolAnalysisTools/outputfiles/Zprime_M1500GeV_W15GeV_36.135pb.root");

    hists = [];
    # hists.append('mttbar_QCDEnriched')
    # hists.append('mttbar_conversions')
    # hists.append('mttbar_conversions_withMETCut')
    # hists.append('mttbar_conversions_withMETAndAsymJets')
    # hists.append('mttbar_conversions_withAsymJetsCut')
    
    hists.append("mttbar");
    # hists.append("tPrimeMass");
    # hists.append("tPrimepT");
    # hists.append("tPrimeHT");
    # hists.append("tPrime_pt");
    # hists.append("tPrime_px");
    # hists.append("tPrime_py");
    # hists.append("tPrime_pz");
    # hists.append("mLeptonicTopLone");
    # hists.append("mHadronicTopLone");
    # hists.append("mLepTopLoneNM");
    # hists.append("mHadTopLoneNM");
    # # hists.append("chiHadronicTop");
    # hists.append("chiLeptonicTop");
    # hists.append("chiGlobal");
    # hists.append("chiTotal");
    # hists.append("lepTopDeltaR");
    # hists.append("hadTopDeltaR");
    # hists.append("lepTopDeltaRNM");
    # hists.append("hadTopDeltaRNM");
    # hists.append("lepTopLep");
    # hists.append("hadTopHad");
    # hists.append("numHadTopMCMatches");
    # hists.append("numLepTopMCMatches");
    # hists.append("numHadTopCorrectID");
    # hists.append("numLepTopCorrectID");
    # hists.append("numHadTopMCMatchesNM");
    # hists.append("numLepTopMCMatchesNM");
    # hists.append("numHadTopCorrectIDNM");
    # hists.append("numLepTopCorrectIDNM");
    # hists.append("mttbar_2ndSolution");
    # hists.append("mttbar_3rdSolution");
    # hists.append("mttbar_allSolutions");
##    
    # hists.append("mttbar_withMETCut");
    # hists.append("mttbar_2ndSolution_withMETCut");
    # hists.append("mttbar_3rdSolution_withMETCut");
    # hists.append("mttbar_allSolutions_withMETCut");
#    
    # hists.append("mttbar_withMETAndAsymJets");
    # hists.append("mttbar_2ndSolution_withMETAndAsymJets");
    # hists.append("mttbar_3rdSolution_withMETAndAsymJets");
    # hists.append("mttbar_allSolutions_withMETAndAsymJets");
#    
    # hists.append("mttbar_withAsymJetsCut");
    # hists.append("mttbar_2ndSolution_withAsymJetsCut");
    # hists.append("mttbar_3rdSolution_withAsymJetsCut");
    # hists.append("mttbar_allSolutions_withAsymJetsCut");
#    
    # hists.append("ttbar_pt");
    # hists.append("ttbar_pt_2ndSolution");
    # hists.append("ttbar_pt_3rdSolution");
    # hists.append("ttbar_pt_allSolutions");
#    
    # hists.append("ttbar_pt_withMETCut");
    # hists.append("ttbar_pt_2ndSolution_withMETCut");
    # hists.append("ttbar_pt_3rdSolution_withMETCut");
    # hists.append("ttbar_pt_allSolutions_withMETCut");
#    
    # hists.append("ttbar_pt_withMETAndAsymJets");
    # hists.append("ttbar_pt_2ndSolution_withMETAndAsymJets");
    # hists.append("ttbar_pt_3rdSolution_withMETAndAsymJets");
    # hists.append("ttbar_pt_allSolutions_withMETAndAsymJets");
#    
    # hists.append("ttbar_pt_withAsymJetsCut");
    # hists.append("ttbar_pt_2ndSolution_withAsymJetsCut");
    # hists.append("ttbar_pt_3rdSolution_withAsymJetsCut");
    # # hists.append("ttbar_pt_allSolutions_withAsymJetsCut");
    # hists.append("angleTops");
    # hists.append("angleTops_withMETCut");
    # hists.append("angleTops_withMETAndAsymJets");
    # hists.append("angleTops_withAsymJetsCut");
#    
    hists.append("pt_leadingTop");
    # hists.append("pt_leadingTop_withMETCut");
    # hists.append("pt_leadingTop_withMETAndAsymJets");
    # hists.append("pt_leadingTop_withAsymJetsCut");
#    
    hists.append("pt_NextToLeadingTop");
    # hists.append("pt_NextToLeadingTop_withMETCut");
    # hists.append("pt_NextToLeadingTop_withMETAndAsymJets");
    # hists.append("pt_NextToLeadingTop_withAsymJetsCut");
#    
    hists.append("mLeptonicTop");
    hists.append("mHadronicTop");
    hists.append("mAllTop");
#    
    # hists.append("ttbar_px");
    # hists.append("ttbar_py");
    # hists.append("ttbar_pz");
#    
    # hists.append("m3");
    # hists.append("HT");
    # hists.append("MET");
    # hists.append("leadingJetMass");
    # hists.append("mtW");
    # hists.append("neutrino_pz");
    # suffixes = ["allBtags",
    suffixes = ["0orMoreBtag",
        "1orMoreBtag",
        "2orMoreBtags",
        "3orMoreBtags",'4orMoreBtags' ]
    hists = [hist + '_' + suffix for hist in hists for suffix in suffixes]
    
    # jetBinned = ["allJets",
        # "1orMoreJets",
        # "2orMoreJets",
        # "3orMoreJets" , "4orMoreJets"]
    # jetBinnedhists = []
    # jetBinnedhists.append("QCDest_CombRelIso")
    # jetBinnedhists.append("QCDest_CombRelIso_1btag")
    # jetBinnedhists.append("QCDest_CombRelIso_2btag")
    # jetBinnedhists.append("QCDest_CombRelIso_controlRegion")
    # jetBinnedhists.append("QCDest_CombRelIso_controlRegion_1btag")
    # jetBinnedhists.append("QCDest_CombRelIso_controlRegion_2btag")
    # 
    # jetBinnedhists.append("QCDest_PFIsolation")
    # jetBinnedhists.append("QCDest_PFIsolation_WithMETCut")
    # jetBinnedhists.append("QCDest_PFIsolation_WithMETCutAndAsymJetCuts")
    # jetBinnedhists.append("QCDest_PFIsolation_WithAsymJetCuts")
    
    # jetBinnedhists.append("QCDest_PFIsolation_1btag")
    # jetBinnedhists.append("QCDest_PFIsolation_2btag")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion_1btag")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion_2btag")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion2")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion2_WithMETCut")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion2_WithMETCutAndAsymJetCuts")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion2_WithAsymJetCuts")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion2_1btag")
    # jetBinnedhists.append("QCDest_PFIsolation_controlRegion2_2btag")
    
    # jetBinnedhists = [hist + '_' + suffix for hist in jetBinnedhists for suffix in jetBinned]
    # hists.extend(jetBinnedhists)
    gcd = gROOT.cd

    for histname in hists:
        currentSelection = 'default'
        if 'withMETAndAsymJets' in histname:
            currentSelection = 'withMETAndAsymJets'
        gcd()
        print histname
        hist_data =  data.Get(histname);

#        hist_data2;
#        if (histname == "mttbar_rebinned")
#            hist_data2 =  data2.Get(histname);
#        else
#            hist_data2 = TH1F(*hist_data);
#        hist_data.Sumw2();
#        hist_data2.Sumw2();
        hist_ttbar =  ttbar.Get(histname);
        hist_wpm1000 =  WprimeTToTTD_M1000.Get(histname);
        hist_wpm600 =  WprimeTToTTD_M600.Get(histname);
        # hist_wjets =  wjets.Get(histname);
        # hist_zjets =  zjets.Get(histname);

        # hist_bce1 =  bce1.Get(histname);
        # hist_bce2 =  bce2.Get(histname);
        # # hist_bce3 =  bce3.Get(histname);
        # hist_enri1 =  enri1.Get(histname);
        # hist_enri2 =  enri2.Get(histname);
        # hist_enri3 =  enri3.Get(histname);
        # hist_pj1 =  pj1.Get(histname);
        # hist_pj2 =  pj2.Get(histname);
        # hist_pj3 =  pj3.Get(histname);
        # hist_singleTop = tW.Get(histname)
        # hist_singleTop.Add(tchan.Get(histname))

#        hist_Zprime500 =  Zprime500.Get(histname);
#        hist_Zprime750 =  Zprime750.Get(histname);
#        hist_Zprime1000 =  Zprime1000.Get(histname);
#        hist_Zprime1250 =  Zprime1250.Get(histname);
#        hist_Zprime1500 =  Zprime1500.Get(histname);

        hist_ttbar.Scale(scale);
        hist_wpm1000.Scale(scale);
        hist_wpm600.Scale(scale);
        # hist_wjets.Scale(scale);
        # hist_zjets.Scale(scale);

        # hist_bce1.Scale(scale);
        # hist_bce2.Scale(scale);
        # hist_bce3.Scale(scale);

        # hist_enri1.Scale(scale);
        # hist_enri2.Scale(scale);
        # hist_enri3.Scale(scale);

        # hist_pj1.Scale(scale);
        # hist_pj2.Scale(scale);
        # hist_pj3.Scale(scale);

#        hist_Zprime500.Scale(scale);
#        hist_Zprime750.Scale(scale);
#        hist_Zprime1000.Scale(scale);
#        hist_Zprime1250.Scale(scale);
#        hist_Zprime1500.Scale(scale);

        # hist_qcd = hist_bce1.Clone("qcd")#TH1F(*hist_bce1);
        # hist_qcd.Add(hist_bce2);
        # hist_qcd.Add(hist_bce3);
        # hist_qcd.Add(hist_enri1);
        # hist_qcd.Add(hist_enri2);
        # hist_qcd.Add(hist_enri3);
        # hist_qcd.Add(hist_pj1);
        # hist_qcd.Add(hist_pj2);
        # hist_qcd.Add(hist_pj3);
        # hist_qcd.Scale(qcdScale[currentSelection]);
        #        ndata = hist_data.Integral();
        #        ntop = hist_ttbar.Integral();
        #        nwj = hist_wjets.Integral();
        #        nzj = hist_zjets.Integral();
        # nqcd = hist_qcd.Integral();
        #        sumMC = ntop + nwj + nzj + nqcd;
        #        cout << ndata << " " << sumMC << endl;
        #                        hist_wjets.Scale(ndata / sumMC);
        #                        hist_ttbar.Scale(ndata / sumMC);
        #                        hist_zjets.Scale(ndata / sumMC);
        #                        hist_qcd.Scale(ndata / sumMC);
        # mttbars = ['mttbar_' + suffix for suffix in suffixes]
        # mttbars2 = ['mttbar_withMETAndAsymJets_' + suffix for suffix in suffixes]
#        if histname in mttbars or histname in mttbars2:
#            print "taking QCD shape from DATA"
#            name = histname.replace('mttbar', 'mttbar_conversions')
#            hist_qcd = data.Get(name)
#            if( hist_qcd.Integral() > 0):
#                hist_qcd.Scale(nqcd/hist_qcd.Integral())
        # hist_mc = hist_qcd.Clone("all_mc")
        # hist_mc = hist_ttbar.Clone("all_mc")
        # hist_mc.Add(hist_ttbar);
        # hist_mc.Add(hist_zjets);
        # hist_mc.Add(hist_wjets);
        # hist_mc.Add(hist_singleTop);

        rebin = 1;
        Urange = (0,5000)
        Yrange = (0,0)
        if ("mttbar" in histname):
            hist_data.SetXTitle("M_{t#bar{t}}/GeV");
            hist_data.SetYTitle("Events/(100 GeV)");
            rebin = 100;
            Urange = (0, 2200)
        elif ("tPrimeMass" in histname):
            hist_data.SetXTitle("T' mass/GeV");
            hist_data.SetYTitle("Events/(100 GeV)");
            rebin = 10;
            Urange = (0, 2000)
        elif ("tPrimepT" in histname):
            hist_data.SetXTitle("T' pT/GeV");
            hist_data.SetYTitle("Events/(8 GeV)");
            rebin = 5;
            Urange = (0, 200)
        elif ("tPrimeHT" in histname):
            hist_data.SetXTitle("T' HT");
            hist_data.SetYTitle("Events");
            rebin = 5;
            Urange = (0, 1.5)
        elif ("m3" in histname):
            hist_data.SetXTitle("M3/GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 1500)
        elif (histname == "electron_et"):
            hist_data.SetXTitle("electron p_{T}/GeV");
            hist_data.SetYTitle("Events/(5 GeV)");
            rebin = 5;
        elif ("ttbar_pt" in histname):
            hist_data.SetXTitle("p_{T} of t#bar{t} system/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 500)
        elif ("ttbar_px" in histname):
            hist_data.SetXTitle("p_{x} of t#bar{t} system/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 500)
        elif ("ttbar_py" in histname):
            hist_data.SetXTitle("p_{y} of t#bar{t} system/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 500)
        elif ("ttbar_pz" in histname):
            hist_data.SetXTitle("p_{z} of t#bar{t} system/GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 2000)
        elif ("tPrime_pt" in histname):
            hist_data.SetXTitle("p_{T} of T' system/GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 20;
            Urange = (0, 200)
        elif ("tPrime_px" in histname):
            hist_data.SetXTitle("p_{x} of T' system/GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 20;
            Urange = (0, 200)
        elif ("tPrime_py" in histname):
            hist_data.SetXTitle("p_{y} of T' system/GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 20;
            Urange = (0, 200)
        elif ("tPrime_pz" in histname):
            hist_data.SetXTitle("p_{z} of T' system/GeV");
            hist_data.SetYTitle("Events/(100 GeV)");
            rebin = 20;
            Urange = (0, 1800)
        elif ("HT" in histname):
            hist_data.SetXTitle("#Sigma p_{T}/GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 2000)
        elif (histname == "numberOfJets"):
            hist_data.SetXTitle("number of jets");
            hist_data.SetYTitle("Events");
        elif (histname == "numberOfBJets"):
            hist_data.SetXTitle("number of b-tagged jets (SSVHE medium)");
            hist_data.SetYTitle("Events");
        elif ('MET_' in histname):
            hist_data.SetXTitle("MET/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 2;
            Urange = (0, 500)
        elif ("leadingJetMass" in histname):
            hist_data.SetXTitle("leading jet mass/GeV");
            hist_data.SetYTitle("Events/(5 GeV)");
            rebin = 5;
            Urange = (0, 150)
        elif ("mtW" in histname):
            hist_data.SetXTitle("transverse W-boson mass/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
        elif ("electronD0" in histname):
            hist_data.SetXTitle("electron d_{0} / cm");
            hist_data.SetYTitle("Events/(0.001 cm)");
            rebin = 10;
        elif ("angleTops" in histname):
            hist_data.SetXTitle("angle between top quarks");
            hist_data.SetYTitle("Events/(0.1 rad)");
            rebin = 10;
        elif ("neutrino_pz" in histname):
            hist_data.SetXTitle("neutrino p_{Z} /GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (-500, 500)
        elif ('mHadronicTopLone' in histname or 'mHadTopLoneNM' in histname):
            hist_data.SetXTitle("top mass /GeV");
            hist_data.SetYTitle("Events/(25 GeV)");
            rebin = 2;
            Urange = (0, 400)
        elif ('mLeptonicTopLone' in histname or 'mLepTopLoneNM' in histname):
            hist_data.SetXTitle("top mass /GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 2;
            Urange = (0, 450)
        elif ('chiHadronicTop' in histname or 'chiLeptonicTop' in histname):
            hist_data.SetXTitle("top #chi^{2} /GeV");
            hist_data.SetYTitle("Events/(1)");
            rebin = 1;
            Urange = (0, 30)
        elif ('chiGlobal' in histname or 'chiTotal' in histname):
            hist_data.SetXTitle("Event #chi^{2} /GeV");
            hist_data.SetYTitle("Events/(100)");
            rebin = 20;
        elif ('numHadTopMCMatches' in histname or 'numLepTopMCMatches' in histname or 'numHadTopMCMatchesNM' in histname or 'numLepTopMCMatchesNM' in histname):
            hist_data.SetXTitle("Number of MC-Truth Matches");
            hist_data.SetYTitle("Events");
            Yrange = (0, 150)
        elif ('numHadTopCorrectID' in histname or 'numLepTopCorrectID' in histname or 'numLepTopCorrectIDNM' in histname or 'numHadTopCorrectIDNM' in histname):
            hist_data.SetXTitle("Number of Reco Candidates Matched to MC-Truth Particles");
            hist_data.SetYTitle("Events");
            Yrange = (0, 150)
        elif ('lepTopDeltaR' in histname or 'hadTopDeltaR' in histname or 'lepTopDeltaRNM' in histname or 'hadTopDeltaRNM' in histname):
            hist_data.SetXTitle("#Delta R of MC truth top");
            hist_data.SetYTitle("Events");
            rebin = 5;
        elif ('lepTopLep' in histname or 'hadTopHad' in histname):
            hist_data.SetXTitle("top MC truth identity");
            hist_data.SetYTitle("Events");
            Urange = (0, 2)
        elif ('mHadronicTop' in histname or 'mLeptonicTop' in histname or 'mAllTop' in histname):
            hist_data.SetXTitle("top mass /GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
        elif ('pt_leadingTop' in histname or 'pt_NextToLeadingTop' in histname):
            hist_data.SetXTitle("top p_{T} /GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 10;
            Urange = (0, 700)
        elif('QCDest_CombRelIso' in histname):
            hist_data.SetXTitle("relative isolation");
            hist_data.SetYTitle("Events/(0.1)");
            rebin = 10;
            Urange = (0, 2)
        elif('QCDest_PFIsolation' in histname):
            hist_data.SetXTitle("ParticleFlow isolation");
            hist_data.SetYTitle("Events/(0.1)");
            rebin = 10
            Urange = (0, 2)
        

        hist_data.Rebin(rebin);
        hist_ttbar.Rebin(rebin);
        hist_wpm1000.Rebin(rebin);
        hist_wpm600.Rebin(rebin);
        # hist_wjets.Rebin(rebin);
        # hist_zjets.Rebin(rebin);
        # hist_qcd.Rebin(rebin);
        # hist_singleTop.Rebin(rebin)
        
        hist_data.SetAxisRange(Urange[0], Urange[1]);
        hist_ttbar.SetAxisRange(Urange[0], Urange[1]);
        hist_wpm1000.SetAxisRange(Urange[0], Urange[1]);
        hist_wpm600.SetAxisRange(Urange[0], Urange[1]);
        # hist_wjets.SetAxisRange(Urange[0], Urange[1]);
        # hist_zjets.SetAxisRange(Urange[0], Urange[1]);
        # hist_qcd.SetAxisRange(Urange[0], Urange[1]);
        # hist_singleTop.SetAxisRange(Urange[0], Urange[1]);

        # hist_data.SetMarkerStyle(8);
        # hist_data.SetMarkerSize(2.0);
        hist_data.SetLineWidth(2);
        
        hist_ttbar.SetFillStyle(1001);
        hist_ttbar.SetFillColor(kRed + 1);
        # hist_wjets.SetFillStyle(1001);
        # hist_wjets.SetFillColor(kGreen - 3);
        # hist_zjets.SetFillStyle(1001);
        # hist_zjets.SetFillColor(kAzure - 2);
        # hist_qcd.SetFillStyle(1001);
        # hist_qcd.SetFillColor(kYellow);
        # hist_singleTop.SetFillStyle(1001);
        # hist_singleTop.SetFillColor(kMagenta)
#        nbins = hist_qcd.GetXaxis().GetNbins();
#        binwidth = (hist_qcd.GetXaxis().GetXmax() - hist_qcd.GetXaxis().GetXmin()) / nbins;
#        for i in range(1, nbins + 1):
#            yvalue = hist_qcd.GetBinContent(i);
#            #            float xvalue = hist_qcd.GetBinCenter(i);
#            #            float ymin = yvalue - yvalue*0.5;
#            #            float ymax = yvalue + yvalue*0.5;
#            #            float xmin = xvalue - 0.5*binwidth;
#            #            float xmax = xvalue + 0.5*binwidth;
#            error = yvalue * 0.5;
#            hist_mc.SetBinError(i, error);
#            #            qcdUncert.SetPointError(i, xmin, xmax, ymin, ymax);
#        
#        qcdUncert = TGraphAsymmErrors(hist_mc);

        linewidth = 6;
        hist_wpm1000.SetLineColor(kOrange + 1);
        hist_wpm1000.SetLineWidth(linewidth);
        hist_wpm600.SetLineColor(kGray + 3);
        hist_wpm600.SetLineWidth(linewidth);
#        hist_Zprime500.SetLineColor(kTeal - 9);
#        hist_Zprime500.SetLineWidth(3);
#        hist_Zprime500.SetFillStyle(0);
#        hist_Zprime500.SetFillColor(kWhite);
#
#        hist_Zprime750.SetLineColor(kBlue - 6);
#        hist_Zprime750.SetLineWidth(3);
#        hist_Zprime750.SetFillStyle(0);
#        hist_Zprime750.SetFillColor(kWhite);
#
#        hist_Zprime1000.SetLineColor(28);
#        hist_Zprime1000.SetLineWidth(3);
#        hist_Zprime1000.SetFillStyle(0);
#        hist_Zprime1000.SetFillColor(kWhite);
#
#        hist_Zprime1250.SetLineColor(kCyan - 5);
#        hist_Zprime1250.SetLineWidth(3);
#        hist_Zprime1250.SetFillStyle(0);
#        hist_Zprime1250.SetFillColor(kWhite);
#
#        hist_Zprime1500.SetLineColor(kOrange + 1);
#        hist_Zprime1500.SetLineWidth(3);
#        hist_Zprime1500.SetFillStyle(0);
#        hist_Zprime1500.SetFillColor(kWhite);

#        qcdUncert.SetFillColor(kGray + 3);
#        qcdUncert.SetFillStyle(3003);

        leg = TLegend(0.75, 0.75, 0.98, 0.98);
        # leg = TLegend(0.696, 0.35, 0.94, 0.92);
        leg.SetBorderSize(0);
        leg.SetLineStyle(0);
        leg.SetTextFont(42);
        leg.SetFillStyle(0);

        # leg.AddEntry(hist_data, "data", "P");
        leg.AddEntry(hist_data, "data", "L");
        #        leg.AddEntry(hist_data2, "data(no HLT)", "P");
        leg.AddEntry(hist_ttbar, "t#bar{t}", "f");
        leg.AddEntry(hist_wpm1000, "W' 1TeV (50 pb)", "L");
        leg.AddEntry(hist_wpm600, "W' 600GeV (50 pb)", "L");
        # leg.AddEntry(hist_wjets, "W#rightarrowl#nu", "f");
        # leg.AddEntry(hist_zjets, "Z/#gamma*#rightarrowl^{+}l^{-}", "f");
        # leg.AddEntry(hist_qcd, "QCD");
        # leg.AddEntry(hist_qcd, "QCD/#gamma + jets");
        # leg.AddEntry(hist_singleTop, "Single-Top")
#        leg.AddEntry(hist_Zprime500, "Z' 0.5TeV (50pb)");
#        leg.AddEntry(hist_Zprime750, "Z' 0.75TeV (50pb)");
#        leg.AddEntry(hist_Zprime1000, "Z' 1TeV (50pb)");
#        leg.AddEntry(hist_Zprime1250, "Z' 1.25TeV (50pb)");
#        leg.AddEntry(hist_Zprime1500, "Z' 1.5TeV (50pb)");

        
        canvases.append(TCanvas("cname" + histname, histname, 1200, 900))
        canvases[-1].cd().SetRightMargin(0.04);
        hs = THStack("MC", "MC");
        # hs.Add(hist_qcd);
        # hs.Add(hist_zjets);
        # hs.Add(hist_wjets);
        # hs.Add(hist_singleTop);
        hs.Add(hist_ttbar);
        max = 0
        if hs.GetMaximum() > hist_data.GetMaximum():
            max = hs.GetMaximum()*1.1
        else:
            max = hist_data.GetMaximum()*1.1
        
        hist_data.GetYaxis().SetRangeUser(0, max);
        hist_data.SetStats(1);
        # hist_data.Draw('axis');
        hist_data.Draw('hist error');
        hs.Draw("hist same");
        hist_wpm1000.Draw("same");
        hist_wpm600.Draw("same");
#        hist_Zprime500.Draw("same");
#        hist_Zprime750.Draw("same");
#        hist_Zprime1000.Draw("same");
#        hist_Zprime1250.Draw("same");
#        hist_Zprime1500.Draw("same");
        #        qcdUncert.Draw("1 same");
        #        hist_data2.Draw("error same");
        hist_data.Draw("hist error same");
        leg.Draw();

        text1 = TLatex(3.570061, 23.08044, "CMS Preliminary");
        text1.SetNDC();
        text1.SetTextAlign(13);
        text1.SetX(0.38);
        text1.SetY(0.928);
                #text1.SetLineWidth(2);
        text1.SetTextFont(42);
        text1.SetTextSizePixels(24);# dflt=28
        text1.Draw();

        text2 = TLatex(3.570061, 23.08044, "~%.1f pb^{-1} at #sqrt{s} = 7 TeV" % lumi);
        text2.SetNDC();
        text2.SetTextAlign(13);
        text2.SetX(0.38);
        text2.SetY(0.88);
        #text2.SetLineWidth(2);
        text2.SetTextFont(42);
        text2.SetTextSizePixels(24);# dflt=28
        text2.Draw();
        canvases[-1].SaveAs('plots/' + histname + '.png')
        # canvases[-1].SaveAs('plots/' + histname + '.root')
        # canvases[-1].SaveAs('/storage/results/' + histname + '.png')

        cu_hist_data = getCumulativePlot(hist_data, "data");
        cu_hist_ttbar = getCumulativePlot(hist_ttbar, "ttbar");
        # cu_hist_wjets = getCumulativePlot(hist_wjets, "wjets");
        # cu_hist_zjets = getCumulativePlot(hist_zjets, "zjets");
        # cu_hist_qcd = getCumulativePlot(hist_qcd, "qcd");
        # cu_hist_singleTop = getCumulativePlot(hist_singleTop, "singleTop");
##        cu_hist_Zprime500 = getCumulativePlot(hist_Zprime500, "Zprime500");
##        cu_hist_Zprime750 = getCumulativePlot(hist_Zprime750, "Zprime750");
##        cu_hist_Zprime1000 = getCumulativePlot(hist_Zprime1000, "Zprime1000");
##        cu_hist_Zprime1250 = getCumulativePlot(hist_Zprime1250, "Zprime1250");
##        cu_hist_Zprime1500 = getCumulativePlot(hist_Zprime1500, "Zprime1500");
        cu_hist_data.SetTitleSize(0.045, "X");
        cu_hist_data.SetTitleSize(0.045, "Y");
        cu_hist_data.SetYTitle("Integrated Events");
        # cu_hist_data.SetYTitle("Integrated Events/(50 GeV)");
        cu_hist_ttbar.SetYTitle("Integrated Events/(50 GeV)");
##        
#
        cu_hist_data.SetAxisRange(Urange[0], Urange[1]);
        cu_hist_ttbar.SetAxisRange(Urange[0], Urange[1]);
        if (Yrange[1] > 0):
            cu_hist_data.SetAxisRange(Yrange[0], Yrange[1], "Y");
            cu_hist_ttbar.SetAxisRange(Yrange[0], Yrange[1], "Y");
        # cu_hist_wjets.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_zjets.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_qcd.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_singleTop.SetAxisRange(Urange[0], Urange[1]);
        
        # cu_hs = THStack("cu_MC", "cu_MC");
        # cu_hs.Add(cu_hist_qcd);
        # cu_hs.Add(cu_hist_zjets);
        # cu_hs.Add(cu_hist_wjets);
        # cu_hs.Add(cu_hist_ttbar);
        
        scanvases.append(TCanvas("cu_cname" + histname, histname + "(cu)", 1200, 900))
        scanvases[-1].cd().SetRightMargin(0.04);
        cu_hist_ttbar.Draw("error");
        cu_hist_data.Draw("hist error");
        # cu_hs.Draw("hist same");
###        cu_hist_Zprime500.Draw("same");
###        cu_hist_Zprime750.Draw("same");
###        cu_hist_Zprime1000.Draw("same");
###        cu_hist_Zprime1250.Draw("same");
###        cu_hist_Zprime1500.Draw("same");
##        #        cu_hist_data2.Draw("error same");
        cu_hist_data.Draw("hist error same");
        cu_hist_ttbar.Draw("same");
        leg.Draw();
##
        text1.Draw();
##
        text2.Draw();
        # scanvases[-1].SaveAs('/storage/results/' + histname + '_integrated.png')
        scanvases[-1].SaveAs('plots/' + histname + '_integrated.png')
    


# tdrGrid: Turns the grid lines on (true) or off (false)

#def tdrGrid(bool gridOn):
#    tdrStyle.SetPadGridX(gridOn);
#    tdrStyle.SetPadGridY(gridOn);
#

# fixOverlay: Redraws the axis

def fixOverlay():
    gPad.RedrawAxis();


def getCumulativePlot(initial, type):
    global counter
    counter = counter + 1;
    name = initial.GetName()
    name = "cu_" + name + "_" + type + str(counter);
    title = initial.GetTitle()
    title = "cu_" + title + "_" + type;
    xaxis = initial.GetXaxis().GetTitle();
    yaxis = initial.GetYaxis().GetTitle();
    nBins = initial.GetNbinsX();
    cu = TH1F(name, title, nBins, initial.GetXaxis().GetXmin(), initial.GetXaxis().GetXmax());
    for bin in range(1,nBins+1):
        cu.SetBinContent(bin, initial.Integral(bin, nBins));
    
    cu.SetFillStyle(initial.GetFillStyle());
    cu.SetFillColor(initial.GetFillColor());
    cu.SetLineColor(initial.GetLineColor());
    # cu.SetMarkerSize(initial.GetMarkerSize());
    # cu.SetMarkerStyle(initial.GetMarkerStyle());
    # cu.SetMarkerColor(initial.GetMarkerColor());
    cu.SetLineWidth(initial.GetLineWidth());
    cu.GetXaxis().SetTitle(xaxis);
    cu.GetYaxis().SetTitle(yaxis);
    return cu;

if __name__ == "__main__":
    gROOT.SetBatch(True)
    gROOT.ProcessLine('gErrorIgnoreLevel = 1001;')
    plotMttbar()
#    print "press enter to quit"
#    a = raw_input()
