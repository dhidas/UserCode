from __future__ import division

from tdrStyle import *
from math import sqrt
from ROOT import *

canvases = []
scanvases = []
counter = 0
asymrange = 0.5

def plotMttbar():
    tdrstyle = setTDRStyle();
    gStyle.SetHatchesSpacing(1.0);
    gStyle.SetOptStat("emri"); # Print integral & usual stuff
    lumi = 5000;
    oldLumi = 5000;
    scale = lumi / oldLumi;
    qcdScale = {'default':1.72, 'withMETAndAsymJets': 3.03};
    #    data = TFile.Open("data2.9pb_fullSetOfVars.root");

    # data =  TFile.Open("Run_4000pb_PFElectron_PF2PATJets_PFMET.root");
    ttbar =  TFile.Open("TTJets_5000pb_PFElectron_PF2PATJets_PFMET.root");
    # ttbar =  TFile.Open("TT_Tune_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # wjets =  TFile.Open("WJetsToLNu_5000pb_PFElectron_PF2PATJets_PFMET.root");
    # wjets =  TFile.Open("WToENu_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # zjets =  TFile.Open("DYJetsToLL_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # bce1 =  TFile.Open("QCD_Pt-20to30_BCtoE_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # bce2 =  TFile.Open("QCD_Pt-30to80_BCtoE_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # bce3 =  TFile.Open("QCD_Pt-80to170_BCtoE_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # enri1 =  TFile.Open("QCD_Pt-20to30_EMEnriched_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # enri2 =  TFile.Open("QCD_Pt-30to80_EMEnriched_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # enri3 =  TFile.Open("QCD_Pt-80to170_EMEnriched_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # zjets = TFile.Open("DYJets_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # Zprime500 = TFile.Open("Zprime_M500GeV_W50GeV_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # Zprime1500 = TFile.Open("Zprime_M1500GeV_W15GeV_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # Zprime750 = TFile.Open("Zprime_M750GeV_W7500MeV_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # Zprime1000 = TFile.Open("Zprime_M1000GeV_W10GeV_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # Zprime2000 = TFile.Open("Zprime_M2000GeV_W20GeV_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # Zprime4000 = TFile.Open("Zprime_M4000GeV_W400GeV_4000pb_PFElectron_PF2PATJets_PFMET.root");
    # Tprime300 = TFile.Open("TprimeToWb_M300_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # Tprime350 = TFile.Open("TprimeToWb_M350_35.9pb_PFElectron_PF2PATJets_PFMET.root");
    # WprimeTToTTD_M400 = TFile.Open("WprimeTToTTD_M400_5000pb_PFElectron_PF2PATJets_PFMET.root");
    WprimeTToTTD_M600 = TFile.Open("WprimeTToTTD_M600_5000pb_PFElectron_PF2PATJets_PFMET.root");
    # WprimeTToTTD_M800 = TFile.Open("WprimeTToTTD_M800_5000pb_PFElectron_PF2PATJets_PFMET.root");
    # WprimeTToTTD_M1000 = TFile.Open("WprimeTToTTD_M1000_5000pb_PFElectron_PF2PATJets_PFMET.root");
    # WprimeToTTBbar_M1000 = TFile.Open("WprimeToTBbar_M-1000_4000pb_PFElectron_PF2PATJets_PFMET.root");

    # pj1 = TFile.Open("pj1_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # pj2 = TFile.Open("pj2_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # pj3 = TFile.Open("pj3_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # tW = TFile.Open("tW_36.145pb_PFElectron_PF2PATJets_PFMET.root");
    # tchan = TFile.Open("tchan_36.145pb_PFElectron_PF2PATJets_PFMET.root");
#    vqq = TFile.Open("vqq_7.22pb_V4PFCalo.root__fullSetOfVars.root");
#    Zprime1250 = TFile.Open("/storage/workspace/BristolAnalysisTools/outputfiles/Zprime_M1250GeV_W12500MeV_36.135pb.root");
#    Zprime1500 = TFile.Open("/storage/workspace/BristolAnalysisTools/outputfiles/Zprime_M1500GeV_W15GeV_36.135pb.root");

    hists = [];
    # hists.append('mttbar_QCDEnriched')
    # hists.append('mttbar_conversions')
    # hists.append('mttbar_conversions_withMETCut')
    # hists.append('mttbar_conversions_withMETAndAsymJets')
    # hists.append('mttbar_conversions_withAsymJetsCut')
    # hists.append('mttbar_controlRegion')
    # hists.append("jet12vsjet13");
    # hists.append("jet12vsjet23");
    # hists.append("jet13vsjet23");
    # hists.append("mttbar");
    # hists.append("tPrimeMass");
    # hists.append("tPrimepT");
    # hists.append("tPrimeHT");
    # hists.append("tPrime_pt");
    # hists.append("tPrime_px");
    # hists.append("tPrime_py");
    # hists.append("tPrime_pz");
    # hists.append("mLeptonicTopLone");
    # hists.append("mLepTopLoneNM");
    # hists.append("mHadTopLoneNM");
    # hists.append("chiHadronicTop");
    # hists.append("chiLeptonicTop");
    # hists.append("chiGlobal");
    # hists.append("chiTotal");
    # hists.append("lepTopDeltaR");
    # hists.append("hadTopDeltaR");
    # hists.append("lepTopLep");
    # hists.append("hadTopHad");
    # hists.append("numHadTopMCMatches");
    # hists.append("numLepTopMCMatches");
    # hists.append("numHadTopCorrectID");
    # hists.append("numLepTopCorrectID");
    # hists.append("numHadTopMCMatches2top");
    # hists.append("numLepTopMCMatches2top");
    # hists.append("numHadTopCorrectID2top");
    # hists.append("numLepTopCorrectID2top");
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
    hists.append("leadingJetPt");
    hists.append("nextLeadingJetPt");
    hists.append("trailingJetPt");
    # hists.append("pt_leadingTop");
    # hists.append("pt_loneLepTop");
    # hists.append("pt_lepTop");
    # hists.append("pt_hadTop");
    # hists.append("pt_loneHadTop");
    # hists.append("pt_leadingTop_withMETCut");
    # hists.append("pt_leadingTop_withMETAndAsymJets");
    # hists.append("pt_leadingTop_withAsymJetsCut");
#    
    # hists.append("pt_NextToLeadingTop");
    # hists.append("pt_NextToLeadingTop_withMETCut");
    # hists.append("pt_NextToLeadingTop_withMETAndAsymJets");
    # hists.append("pt_NextToLeadingTop_withAsymJetsCut");
#    
    # hists.append("mAllTop");
#    
    # hists.append("ttbar_px");
    # hists.append("ttbar_py");
    # hists.append("ttbar_pz");
#    
    # hists.append("m3");
    # hists.append("HT");
    # hists.append("qkjet_eta");
    # hists.append("qkjet_rapid");
    # hists.append("leadjet_rap");
    # hists.append("lepb_angle");
    # hists.append("angleHadTopD");
    hists.append("wpmassplus");
    hists.append("wpmassminus");
    # hists.append("MET");
    # hists.append("METtopangle_ltop");
    # hists.append("HT3jet_ltop");
    # hists.append("HT4jet_ltop");
    # hists.append("mt2jbothchglep4j");
    # hists.append("mt2jbothchgelec4j");
    # hists.append("mt2jbothchgmu4j");
    # hists.append("mt2jbothchglep5j");
    # hists.append("mt2jbothchgmu5j");
    # hists.append("mt2jbothchgelec5j");
    # hists.append("mt2jpluselec4j");
    # hists.append("mt2jplusmu4j");
    # hists.append("mt2jpluslep4j");
    # hists.append("mt2jpluslep5j");
    # hists.append("mt2jpluselec5j");
    # hists.append("mt2jplusmu5j");
    # hists.append("mt2jminuselec4j");
    # hists.append("mt2jminusmu4j");
    # hists.append("mt2jminuslep4j");
    # hists.append("mt2jminuselec5j");
    # hists.append("mt2jminusmu5j");
    # hists.append("mt2jminuslep5j");
    # hists.append("mt4jbothchglep4j");
    # hists.append("mt4jbothchgelec4j");
    # hists.append("mt4jbothchgmu4j");
    # hists.append("mt4jbothchglep5j");
    # hists.append("mt4jbothchgmu5j");
    # hists.append("mt4jbothchgelec5j");
    # hists.append("mt4jpluselec4j");
    # hists.append("mt4jplusmu4j");
    # hists.append("mt4jpluslep4j");
    # hists.append("mt4jpluslep5j");
    # hists.append("mt4jpluselec5j");
    # hists.append("mt4jplusmu5j");
    # hists.append("mt4jminuselec4j");
    # hists.append("mt4jminusmu4j");
    # hists.append("mt4jminuslep4j");
    # hists.append("mt4jminuselec5j");
    # hists.append("mt4jminusmu5j");
    # hists.append("mt4jminuslep5j");
    # hists.append("leadingJetMass");
    # hists.append("mtW");
    # hists.append("neutrino_pz");
    # hists.append("lone_neutrino_pz");
    suffixes = ["0orMoreBtag",
        "1orMoreBtag",
        "2orMoreBtags",
        "3orMoreBtags",'4orMoreBtags' ]
    hists = [hist + '_' + suffix for hist in hists for suffix in suffixes]

    histslep = [];
    # histslep.append("HT_ltop");
    histslep.append("MET_ltop");
    histslep.append("METpt_ltop");
    # histslep.append("mass4jets");
    # histslep.append("ST_ltop");
    histslep.append("mLeptonicTop");
    histslep.append("mHadronicTop");
    histslep.append("mHadTopLone");
    histslep.append("mLepTopLone");
    histslep.append("mHadTopLonepos");
    histslep.append("mLepTopLonepos");
    histslep.append("mHadTopLoneneg");
    histslep.append("mLepTopLoneneg");
    histslep.append("wpmass");
    histslep.append("wpnegmass");
    histslep.append("wpnegmassbad");
    histslep.append("wpnegmasspos");
    histslep.append("wpnegmassneg");
    histslep.append("wp1topmass");
    histslep.append("wp1topmassbad");
    histslep.append("wp1topmasspos");
    histslep.append("wp1topmassbadpos");
    histslep.append("wp1topmassneg");
    histslep.append("wp1topmassbadneg");
    # histslep.append("wpanglmass");
    histslep.append("tausInJets");
    histslep.append("numTaus");
    histslep.append("deltaPhiTops");
    histslep.append("deltaphilepjetpos");
    histslep.append("deltaphilepdjetpos");
    histslep.append("deltaphilepjetneg");
    histslep.append("deltaphilepdjetneg");
    # histslep.append("numHadTopDauMCMatches");
    # histslep.append("numLepTopDauMCMatches");
    histslep.append("numHadTopDauMCMatches2top");
    histslep.append("numLepTopDauMCMatches2top");
    histslep.append("numGenJets");
    histslep.append("numDiffJets");
    histslep.append("missGenJets");
    histslep.append("missRecoJets");
   # histslep.append("transMasslm1b");
    suffixes = ["elec", "mu", "lepton"]
    histslep = [hist + '_' + suffix for hist in histslep for suffix in suffixes]

    hists += histslep;
    print hists
    hists.append("eleciso");
    hists.append("looseeleciso");
    hists.append("muiso");
    # hists.append("electron_et");
    hists.append("elvertdist");
    hists.append("muvertdist");
    hists.append("jetmultip");
    # hists.append("electminus_pt");
    # hists.append("positron_pt");
    hists.append("electron_pt");
    # hists.append("muon_et");
    hists.append("muon_pt");
    
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
        # hist_data =  data.Get(histname);

#        hist_data2;
#        if (histname == "mttbar_rebinned")
#            hist_data2 =  data2.Get(histname);
#        else
#            hist_data2 = TH1F(*hist_data);
#        hist_data.Sumw2();
#        hist_data2.Sumw2();
        hist_ttbar =  ttbar.Get(histname);
        # hist_wjets =  wjets.Get(histname);
        # hist_zjets =  zjets.Get(histname);

        # hist_bce1 =  bce1.Get(histname);
        # hist_bce2 =  bce2.Get(histname);
        # hist_bce3 =  bce3.Get(histname);
        # hist_enri1 =  enri1.Get(histname);
        # hist_enri2 =  enri2.Get(histname);
        # hist_enri3 =  enri3.Get(histname);
        # # hist_pj1 =  pj1.Get(histname);
        # # hist_pj2 =  pj2.Get(histname);
        # hist_pj3 =  pj3.Get(histname);
        # hist_singleTop = tW.Get(histname)
        # hist_singleTop.Add(tchan.Get(histname))

        # hist_Zprime500 =  Zprime500.Get(histname);
        # hist_Zprime750 =  Zprime750.Get(histname);
        # hist_Zprime1000 =  Zprime1000.Get(histname);
        # hist_Zprime1500 =  Zprime1500.Get(histname);
        # hist_Zprime2000 =  Zprime2000.Get(histname);
        # hist_Zprime4000 =  Zprime4000.Get(histname);
        # hist_Tprime300  =  Tprime300.Get(histname);
        # hist_Tprime350  =  Tprime350.Get(histname);
# #        hist_Zprime1250 =  Zprime1250.Get(histname);
        # hist_wpm400 =  WprimeTToTTD_M400.Get(histname);
        hist_wpm600 =  WprimeTToTTD_M600.Get(histname);
        # hist_wpm800 =  WprimeTToTTD_M800.Get(histname);
        # hist_wpm1000 =  WprimeTToTTD_M1000.Get(histname);
        # hist_wptb1000 =  WprimeToTTBbar_M1000.Get(histname);
        
        # hist_ttbar.Sumw2(); # No good, must be done before filling

        # if ('iso' in histname or 'jetmultip' in histname or 'DeltaR' in histname or ('num' in histname and 'Top' in histname)):
        if ('looseeleciso' in histname or 'tau' in histname or 'DeltaR' in histname or ('num' in histname and 'Top' in histname)):
          numevts = hist_ttbar.Integral();
          if (numevts > 0):
            hist_ttbar.Scale(1.0/numevts);
            # numevts = hist_wpm400.Integral();
            # hist_wpm400.Scale(1.0/numevts);
            numevts = hist_wpm600.Integral();
            hist_wpm600.Scale(1.0/numevts);
            # numevts = hist_wpm800.Integral();
            # hist_wpm800.Scale(1.0/numevts);
            # numevts = hist_wpm1000.Integral();
            # hist_wpm1000.Scale(1.0/numevts);
        else:
          hist_ttbar.Scale(scale);
          # hist_wpm400.Scale(scale);
          hist_wpm600.Scale(scale);
          # hist_wpm1000.Scale(scale);
          # hist_wpm800.Scale(scale);
          # hist_wjets.Scale(scale);
        
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

        # scale = scale * 10.0;
        # hist_Zprime500.Scale(scale);
        # hist_Zprime750.Scale(scale);
        # hist_Zprime1000.Scale(scale);
        # hist_Zprime2000.Scale(scale);
        # hist_Zprime1500.Scale(scale);
        # hist_Zprime4000.Scale(scale);
        # hist_Tprime300.Scale(scale);
        # hist_Tprime350.Scale(scale);
#        hist_Zprime1250.Scale(scale);
        # hist_wptb1000.Scale(scale);


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
        # hist_mc.Add(hist_ttbar);
        # hist_mc.Add(hist_zjets);
        # hist_mc.Add(hist_wjets);
        # hist_mc.Add(hist_singleTop);

        hist_data = hist_ttbar.Clone("test");
        rebin = 1;
        Urange = (0,5000)
        Yrange = (0,0)
        
        if ("mttbar" in histname):
            hist_data.SetXTitle("M_{t#bar{t}}/GeV");
            hist_data.SetYTitle("Events/(100 GeV)");
            rebin = 100;
            Urange = (0, 2000)
        elif ("tPrimeMass" in histname):
            hist_data.SetXTitle("T' mass/GeV");
            hist_data.SetYTitle("Events/(100 GeV)");
            rebin = 200;
            Urange = (0, 2000)
        elif ("wpmass" in histname or "wpnegmass" in histname or "wpanglmass" in histname or "wp1top" in histname):
            hist_data.SetXTitle("W' mass/GeV");
            hist_data.SetYTitle("Events/(60 GeV)");
            rebin = 60;
            Urange = (0, 1500)
        elif ("jet12" in histname or "jet13" in histname):
            hist_data.SetXTitle("Jet mass/GeV");
            hist_data.SetYTitle("Jet mass/GeV)");
        elif ("tPrimepT" in histname):
            hist_data.SetXTitle("T' pT/GeV");
            hist_data.SetYTitle("Events/(8 GeV)");
            rebin = 8;
            Urange = (0, 200)
        elif ("tPrimeHT" in histname):
            hist_data.SetXTitle("T' HT");
            hist_data.SetYTitle("Events");
            rebin = 4;
            Urange = (0, 1.5)
        elif ("m3" in histname):
            hist_data.SetXTitle("M3/GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 1500)
        elif (histname == "electron_et"):
            hist_data.SetXTitle("electron E_{T}/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 300)
        elif (histname == "electron_pt" or histname == "electminus_pt" or histname == "positron_pt"):
            hist_data.SetXTitle("electron p_{T}/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 300)
        elif (histname == "muon_et"):
            hist_data.SetXTitle("muon E_{T}/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 300)
        elif (histname == "muon_pt"):
            hist_data.SetXTitle("muon p_{T}/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 300)
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
        elif ("METpt_ltop" in histname):
            hist_data.SetXTitle("MET p_{T}/GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 20;
            Urange = (0, 400)
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
            rebin = 200;
            Urange = (0, 1800)
        elif ("HT_" in histname):
            hist_data.SetXTitle("#Sigma p_{T}/GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 20;
            Urange = (0, 2000)
        elif ("ST" in histname):
            hist_data.SetXTitle("ST/GeV");
            hist_data.SetYTitle("Events/(30 GeV)");
            rebin = 30;
            Urange = (0, 1000)
        elif ("HT3jet_ltop" in histname):
            hist_data.SetXTitle("#Sigma p_{T} of leading 3 jets/GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 20;
            Urange = (0, 1000)
        elif ("HT4jet_ltop" in histname):
            hist_data.SetXTitle("#Sigma p_{T} of leading 4 jets/GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 1200)
        elif (histname == "numberOfJets"):
            hist_data.SetXTitle("number of jets");
            hist_data.SetYTitle("Events");
        elif (histname == "numberOfBJets"):
            hist_data.SetXTitle("number of b-tagged jets (SSVHE medium)");
            hist_data.SetYTitle("Events");
        elif ('MET_' in histname and 'MET_ltop' not in histname):
            hist_data.SetXTitle("MET/GeV");
            hist_data.SetYTitle("Events/(5 GeV)");
            rebin = 5;
            Urange = (0, 400)
        elif ('MET_ltop' in histname):
            hist_data.SetXTitle("MET/GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 20;
            Urange = (0, 400)
        elif ("leadingJetMass" in histname):
            hist_data.SetXTitle("leading jet mass/GeV");
            hist_data.SetYTitle("Events/(5 GeV)");
            rebin = 5;
            Urange = (0, 150)
        elif ("mtW" in histname):
            hist_data.SetXTitle("transverse W-boson mass/GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
        elif ("both" in histname or "plus" in histname or "minus" in histname):
            hist_data.SetXTitle("transverse mass/GeV");
            hist_data.SetYTitle("Events/(100 GeV)");
            rebin = 100;
            Urange = (0, 2000)
        elif ("transMasslm1b" in histname):
            hist_data.SetXTitle("M_{T} (lep+MET+1b)/GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 700)
        elif ("electronD0" in histname):
            hist_data.SetXTitle("electron d_{0} / cm");
            hist_data.SetYTitle("Events/(0.001 cm)");
            rebin = 10;
        elif ("angleTops" in histname or "deltaPhiTops" in histname):
            hist_data.SetXTitle("angle between top quarks");
            hist_data.SetYTitle("Events/(0.1 rad)");
            rebin = 10;
        elif ( "deltaphilep" in histname):
            hist_data.SetXTitle("angle between lepton & leading jet");
            hist_data.SetYTitle("Events/(0.2 rad)");
            rebin = 20;
        elif ("lepb_angle" in histname):
            hist_data.SetXTitle("angle between electron and b jet");
            hist_data.SetYTitle("Events/(0.1 rad)");
            rebin = 10;
        elif ("angleHadTopD" in histname):
            hist_data.SetXTitle("angle between hadronic top and light jet");
            hist_data.SetYTitle("Events/(0.2 rad)");
            rebin = 20;
        elif ("METtopangle_ltop" in histname):
            hist_data.SetXTitle("angle between MET and top p_{T}");
            hist_data.SetYTitle("Events/(0.1 rad)");
            rebin = 10;
            Urange = (-2, 2)
        elif ("qkjet_eta" in histname):
            hist_data.SetXTitle("#eta of b jet");
            hist_data.SetYTitle("Events/(0.1 rad)");
            Urange = (-3, 3)
            rebin = 10;
        elif ("_rap" in histname):
            hist_data.SetXTitle("Jet Rapidity");
            hist_data.SetYTitle("Events/(0.1 rad)");
            Urange = (-3, 3)
            rebin = 20;
        elif ("neutrino_pz" in histname):
            hist_data.SetXTitle("neutrino p_{Z} /GeV");
            hist_data.SetYTitle("Events/(20 GeV)");
            rebin = 20;
            Urange = (-500, 500)
        elif ('mHadTopLoneNM' in histname):
            hist_data.SetXTitle("top mass /GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 400)
        elif ('mass4jets' in histname):
            hist_data.SetXTitle("4-jet mass /GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 1500)
        elif ('mLepTopLoneNM' in histname):
            hist_data.SetXTitle("top mass /GeV");
            hist_data.SetYTitle("Events/(10 GeV)");
            rebin = 10;
            Urange = (0, 400)
        elif ('chiHadronicTop' in histname or 'chiLeptonicTop' in histname):
            hist_data.SetXTitle("top #chi^{2} /GeV");
            hist_data.SetYTitle("Events/(1)");
            rebin = 1;
            Urange = (0, 30)
        elif ('chiGlobal' in histname or 'chiTotal' in histname):
            hist_data.SetXTitle("Event #chi^{2} /GeV");
            hist_data.SetYTitle("Events/(100)");
            rebin = 100;
        elif ('numHadTopMCMatches' in histname or 'numLepTopMCMatches' in histname):
            hist_data.SetXTitle("Number of Exact Decay-Chain Matches");
            hist_data.SetYTitle("Events");
            Yrange = (0, 1.1)
        elif ('numHadTopDauMCMatches' in histname or 'numLepTopDauMCMatches' in histname):
            hist_data.SetXTitle("Number of Top Daughter Matches");
            hist_data.SetYTitle("Events");
            Yrange = (0, 1.1)
        elif ('numHadTopCorrectID' in histname or 'numLepTopCorrectID' in histname):
            hist_data.SetXTitle("Number of Reco Candidates Matched to MC-Truth Particles");
            hist_data.SetYTitle("Events");
            Yrange = (0, 1.1)
        elif ('lepTopDeltaR' in histname or 'hadTopDeltaR' in histname):
            hist_data.SetXTitle("#Delta R of MC truth top");
            hist_data.SetYTitle("Events");
            rebin = 20;
            # Urange = (0, 4)
        elif ('lepTopLep' in histname or 'hadTopHad' in histname):
            hist_data.SetXTitle("top MC truth identity");
            hist_data.SetYTitle("Events");
            Urange = (0, 2)
        elif ('mHadTopLone' in histname or 'mHadronicTop' in histname or 'mAllTop' in histname):
            hist_data.SetXTitle("top mass /GeV");
            hist_data.SetYTitle("Events/(15 GeV)");
            rebin = 15;
            Urange = (0, 700)
        elif ('mLepTopLone' in histname or 'mLeptonicTop' in histname):
            hist_data.SetXTitle("top mass /GeV");
            hist_data.SetYTitle("Events/(15 GeV)");
            rebin = 15;
            Urange = (0, 700)
        elif ('pt_' in histname and 'Top' in histname):
            hist_data.SetXTitle("top p_{T} /GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 500)
        elif ('trailingJetPt' in histname):
            hist_data.SetXTitle("Jet p_{T} /GeV");
            hist_data.SetYTitle("Events/(25 GeV)");
            rebin = 30;
            Urange = (0, 270)
        elif ('JetPt' in histname):
            hist_data.SetXTitle("Jet p_{T} /GeV");
            hist_data.SetYTitle("Events/(50 GeV)");
            rebin = 50;
            Urange = (0, 600)
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
        elif('iso' in histname):
            hist_data.SetXTitle("ParticleFlow Isolation");
            hist_data.SetYTitle("Events/(0.02)");
            rebin = 4
            Urange = (0, 0.5)
        elif('vertdist' in histname):
            hist_data.SetXTitle("Steps from PV");
            hist_data.SetYTitle("Events");
            # rebin = 200
            Urange = (0, 2)
        elif('jetmultip' in histname):
            hist_data.SetXTitle("Number of jets");
            hist_data.SetYTitle("Events");
            # rebin = 60
            Urange = (0, 11)
        elif('tausInJets' in histname):
            hist_data.SetXTitle("Number of tau jets");
            hist_data.SetYTitle("Events");
            # rebin = 60
            Urange = (0, 3)
        elif('numTaus' in histname):
            hist_data.SetXTitle("Number of taus");
            hist_data.SetYTitle("Events");
            # rebin = 60
            Urange = (0, 15)
        elif('numGenJets' in histname):
            hist_data.SetXTitle("Number of GenJets");
            hist_data.SetYTitle("Events");
            # rebin = 60
            Urange = (3, 15)
        elif('numDiffJets' in histname):
            hist_data.SetXTitle("Difference between Reco Jets & GenJets");
            hist_data.SetYTitle("Events");
            # rebin = 60
            Urange = (-7, 4)
        elif('miss' in histname and 'Jets' in histname):
            if('miss' in histname and 'missGenJets' in histname):
              hist_data.SetXTitle("Number of Unmatched Reco Jets");
            else:
              hist_data.SetXTitle("Number of Unmatched GenJets");
            hist_data.SetYTitle("Events");
            # rebin = 60
            Urange = (0, 8)
       

        hist_data.SetTitleOffset(1.5, "Y");
        if ('jet1' in histname):
                  hist_data.Rebin2D();
                  hist_ttbar.Rebin2D();
                  hist_Tprime300.Rebin2D();
                  hist_Tprime350.Rebin2D();
                  hist_data.SetAxisRange(0, 600, "X");
                  hist_data.SetAxisRange(0, 600, "Y");
                  hist_ttbar.SetAxisRange(0, 600, "X");
                  hist_ttbar.SetAxisRange(0, 600, "Y");
                  # hist_Tprime300.SetAxisRange(0, 600, "X");
                  # hist_Tprime300.SetAxisRange(0, 600, "Y");
                  # hist_Tprime350.SetAxisRange(0, 600, "X");
                  # hist_Tprime350.SetAxisRange(0, 600, "Y");
        else:
                  hist_data.Rebin(rebin);
                  hist_ttbar.Rebin(rebin);
                  # hist_wjets.Rebin(rebin);
                  # hist_zjets.Rebin(rebin);
                  # hist_qcd.Rebin(rebin);
                  # hist_Zprime500.Rebin(rebin);
                  # hist_Zprime1500.Rebin(rebin);
                  # hist_Zprime750.Rebin(rebin);
                  # hist_Zprime1000.Rebin(rebin);
                  # hist_Zprime2000.Rebin(rebin);
                  # hist_Zprime4000.Rebin(rebin);
                  # hist_Tprime300.Rebin(rebin);
                  # hist_Tprime350.Rebin(rebin);
                  # hist_singleTop.Rebin(rebin)
                  # hist_wpm400.Rebin(rebin);
                  hist_wpm600.Rebin(rebin);
                  # hist_wpm800.Rebin(rebin);
                  # hist_wpm1000.Rebin(rebin);
                  # hist_wptb1000.Rebin(rebin);

                  hist_data.SetAxisRange(Urange[0], Urange[1]);
                  hist_ttbar.SetAxisRange(Urange[0], Urange[1]);
                  # hist_wjets.SetAxisRange(Urange[0], Urange[1]);
                  # hist_zjets.SetAxisRange(Urange[0], Urange[1]);
                  # hist_qcd.SetAxisRange(Urange[0], Urange[1]);
                  # hist_Zprime500.SetAxisRange(Urange[0], Urange[1]);
                  # hist_Zprime750.SetAxisRange(Urange[0], Urange[1]);
                  # hist_Zprime1000.SetAxisRange(Urange[0], Urange[1]);
                  # hist_Zprime1500.SetAxisRange(Urange[0], Urange[1]);
                  # hist_Zprime2000.SetAxisRange(Urange[0], Urange[1]);
                  # hist_Zprime4000.SetAxisRange(Urange[0], Urange[1]);
                  # hist_Tprime300.SetAxisRange(Urange[0], Urange[1]);
                  # hist_Tprime350.SetAxisRange(Urange[0], Urange[1]);
                  # hist_singleTop.SetAxisRange(Urange[0], Urange[1]);
                  # hist_wpm400.SetAxisRange(Urange[0], Urange[1]);
                  hist_wpm600.SetAxisRange(Urange[0], Urange[1]);
                  # hist_wpm800.SetAxisRange(Urange[0], Urange[1]);
                  # hist_wpm1000.SetAxisRange(Urange[0], Urange[1]);
                  # hist_wptb1000.SetAxisRange(Urange[0], Urange[1]);


        hist_data.SetMarkerStyle(8);
        # hist_data.SetMarkerSize(1.5);
        hist_data.SetMarkerSize(0.1);
        
        # hist_data.SetFillStyle(1001);
        # hist_data.SetFillColor(kRed);
        hist_ttbar.SetFillStyle(1001);
        hist_ttbar.SetFillColor(kRed);
        # hist_wjets.SetFillStyle(1001);
        # hist_wjets.SetFillColor(4);
        # hist_zjets.SetFillStyle(1001);
        # hist_zjets.SetFillColor(8);
        # hist_qcd.SetFillStyle(1001);
        # hist_qcd.SetFillColor(kYellow);
        # hist_qcd.SetFillColor(12);
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
        # hist_Zprime500.SetLineColor(6);
        # hist_Zprime500.SetLineWidth(linewidth);
        # hist_Zprime500.SetFillStyle(0);
        # hist_Zprime500.SetFillColor(kWhite);

        # hist_Zprime750.SetLineColor(5);
        # hist_Zprime750.SetLineWidth(linewidth);
        # hist_Zprime750.SetFillStyle(0);
        # hist_Zprime750.SetFillColor(kWhite);

        # hist_Zprime1000.SetLineColor(7);
        # hist_Zprime1000.SetLineWidth(linewidth);
        # hist_Zprime1000.SetFillStyle(0);
        # hist_Zprime1000.SetFillColor(kWhite);
        # hist_Zprime2000.SetLineColor(42);
        # hist_Zprime2000.SetLineWidth(linewidth);
        # hist_Zprime2000.SetFillStyle(0);
        # hist_Zprime2000.SetFillColor(kWhite);
        # # hist_Zprime4000.SetLineColor(50);
        # hist_Zprime4000.SetLineWidth(linewidth);
        # hist_Zprime4000.SetFillStyle(0);
        # hist_Zprime4000.SetFillColor(kWhite);

        # hist_Tprime300.SetLineColor(kOrange + 1);
        # hist_Tprime300.SetLineWidth(linewidth);
        # hist_Tprime350.SetLineColor(kGray + 3);
        # hist_Tprime350.SetLineWidth(linewidth);
#
#        hist_Zprime1250.SetLineColor(kCyan - 5);
#        hist_Zprime1250.SetLineWidth(linewidth);
#        hist_Zprime1250.SetFillStyle(0);
#        hist_Zprime1250.SetFillColor(kWhite);
#
        # hist_Zprime1500.SetLineColor(42);
        # hist_Zprime1500.SetLineWidth(linewidth);
        # hist_Zprime1500.SetFillStyle(0);
        # hist_Zprime1500.SetFillColor(kWhite);

        # hist_wpm400.SetLineColor(kGreen + 1);
        # hist_wpm400.SetLineWidth(linewidth);
        hist_wpm600.SetLineColor(kGray + 3);
        hist_wpm600.SetLineWidth(linewidth);
        # hist_wpm800.SetLineColor(kCyan);
        # hist_wpm800.SetLineWidth(linewidth);
        # hist_wpm1000.SetLineColor(kOrange + 1);
        # hist_wpm1000.SetLineWidth(linewidth);
        # hist_wptb1000.SetLineColor(kCyan);
        # hist_wptb1000.SetLineWidth(linewidth);

#        qcdUncert.SetFillColor(kGray + 3);
#        qcdUncert.SetFillStyle(3003);

        leg = TLegend(0.696, 0.35, 0.95, 0.92);
        leg.SetBorderSize(0);
        leg.SetLineStyle(0);
        leg.SetTextFont(42);
        leg.SetFillStyle(0);

        # leg.AddEntry(hist_data, "data", "P");
        #        leg.AddEntry(hist_data2, "data(no HLT)", "P");
        leg.AddEntry(hist_ttbar, "t#bar{t}", "f");
        # leg.AddEntry(hist_wjets, "W#rightarrowl#nu", "f");
        # leg.AddEntry(hist_zjets, "Z/#gamma*#rightarrowl^{+}l^{-}", "f");
        # leg.AddEntry(hist_qcd, "QCD", "F");
        # leg.AddEntry(hist_qcd, "QCD/#gamma + jets");
        # leg.AddEntry(hist_singleTop, "Single-Top")
        # leg.AddEntry(hist_Zprime500, "Z' 0.5TeV (50 pb)", "L");
        # leg.AddEntry(hist_Zprime750, "Z' 0.75TeV", "L");
        # leg.AddEntry(hist_Zprime1000, "Z' 1TeV", "L");
        # leg.AddEntry(hist_Zprime1500, "Z' 1.5TeV (50 pb)", "L");
        # leg.AddEntry(hist_Zprime2000, "Z' 2TeV", "L");
        # leg.AddEntry(hist_Zprime4000, "Z' 4TeV (50 pb)", "L");
        # leg.AddEntry(hist_Tprime300, "T' 300GeV", "L");
        # leg.AddEntry(hist_Tprime350, "T' 350GeV", "L");
#        leg.AddEntry(hist_Zprime1250, "Z' 1.25TeV");
        # if (not 'num' in histname or not 'Top' in histname):
          # leg.AddEntry(hist_wpm400, "W't .4TeV (32 pb)", "L");
        leg.AddEntry(hist_wpm600, "W't .6TeV (8 pb)", "L");
        # if (not 'num' in histname or not 'Top' in histname):
          # leg.AddEntry(hist_wpm800, "W't .8TeV (2.2 pb)", "L");
          # leg.AddEntry(hist_wpm1000, "W't 1TeV (0.72 pb)", "L");
        # leg.AddEntry(hist_wptb1000, "W'-> tb 1TeV (8 pb)", "L");

        
        canvases.append(TCanvas("cname" + histname, histname, 1200, 900))
        
        # if ('iso' in histname or 'jetmultip' in histname or 'num' in histname and 'Top' in histname):
        canvases[-1].SetGrid();
        canvases[-1].cd().SetRightMargin(0.04);
        
        if ('jet1' in histname):
          hist_data.SetMarkerStyle(8);
          hist_ttbar.SetMarkerStyle(21);
          hist_ttbar.SetMarkerColor(kRed);
          hist_Tprime300.SetMarkerColor(kOrange + 1);
          hist_Tprime300.SetMarkerStyle(22);
          hist_Tprime350.SetMarkerColor(kGray + 3);
          hist_Tprime350.SetMarkerStyle(3);
        else:
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
          if hist_wpm600.GetMaximum() > max :
              max = hist_wpm600.GetMaximum() * 1.1
#           if hist_wpm400.GetMaximum() > max :
              # max = hist_wpm400.GetMaximum() * 1.1
          
        if (false and 'wp1topmass' in histname):
          hist_data.GetYaxis().SetRangeUser(0, 310);
        elif (false and 'mLepTopLone' in histname):
          hist_data.GetYaxis().SetRangeUser(0, 510);
        elif (false and 'mHadTopLone' in histname):
          hist_data.GetYaxis().SetRangeUser(0, 400);
        elif ('jetmultip' in histname):
          hist_data.GetYaxis().SetRangeUser(0, 1600);
        elif ('looseeleciso' in histname):
          hist_data.GetYaxis().SetRangeUser(0, 0.2);
        else:
          hist_data.GetYaxis().SetRangeUser(0, max);
        # hist_data.SetStats(1);  # SetOptStat at start turns on stats
        # hist_data.Draw('error');
        # if ('iso' in histname):
          # canvases[-1].SetLogy();
        hist_data.Draw("hist");
        hist_data.SetMarkerStyle(8);
        if ('jet1' in histname):
          hist_ttbar.Draw("same");
        else:
          hs.Draw("hist same");
        # hist_Zprime500.Draw("same");
        # hist_Zprime750.Draw("same");
        # hist_Zprime1000.Draw("same");
        # hist_Zprime2000.Draw("same");
        # hist_Zprime4000.Draw("same");
        # hist_Tprime300.Draw("same");
        # hist_Tprime350.Draw("same");
#        hist_Zprime1250.Draw("same");
        # hist_Zprime1500.Draw("same");
        hist_wpm600.Draw("hist same");
        # if (not 'num' in histname or not 'Top' in histname):
          #hist_wpm400.Draw("hist same");
          # hist_wpm800.Draw("hist same");
          # hist_wpm1000.Draw("hist same");
        # hist_wptb1000.Draw("same");
        #        qcdUncert.Draw("1 same");
        #        hist_data2.Draw("error same");
        # hist_data.Draw("error same");
        # hist_data.Draw("same");
        leg.Draw();

        text1 = TLatex(3.570061, 23.08044, "CMS Preliminary");
        text1.SetNDC();
        text1.SetTextAlign(13);
        text1.SetX(0.38);
        text1.SetY(0.928);
                #text1.SetLineWidth(2);
        text1.SetTextFont(42);
        text1.SetTextSizePixels(24);# dflt=28
        text1.SetTextSize(0.03);
        text1.Draw();

        text2 = TLatex(3.570061, 23.08044, "%d pb^{-1} at #sqrt{s} = 7 TeV" % lumi);
        text2.SetNDC();
        text2.SetTextAlign(13);
        text2.SetX(0.38);
        text2.SetY(0.88);
        #text2.SetLineWidth(2);
        text2.SetTextFont(42);
        text2.SetTextSizePixels(24);# dflt=28
        text2.SetTextSize(0.03);
        text2.Draw();
        # canvases[-1].SaveAs('plots/' + histname + '.png')
        if (false and 'iso' in histname):
          canvases[-1].SaveAs('plots/' + histname + '.root')
        else:
            canvases[-1].SaveAs('plots/' + histname + '.png')
        
        if (histname == "wpnegmass_lepton" or histname == "wp1topmass_lepton" or histname == "wp1topmasspos_lepton" or histname == "wp1topmassneg_lepton"):
            hist_ttbargood = hist_ttbar.Clone("test");
            # hist_ttbargood.Add(hist_wjets);
            hist_wprimegood = hist_wpm600.Clone("test");
        elif (histname == "wpnegmassbad_lepton" or histname == "wp1topmassbad_lepton" or histname == "wp1topmassbadpos_lepton" or histname == "wp1topmassbadneg_lepton"):
            if (histname == "wp1topmassbad_lepton"):
              suffix = "1top";
            elif (histname == "wp1topmassbadpos_lepton"):
              suffix = "1toppos";
            elif (histname == "wp1topmassbadneg_lepton"):
              suffix = "1topneg";
            else:
                suffix = "";
            ttbname = 'ttbarasym' + suffix;
            wpname = 'wprimeasym' + suffix;
            comasymname = 'wpttcombined' + suffix;
            wpgoodbad = "wpnegmassgoodbad" + suffix;
            ttbgoodbad = "ttbargoodbad" + suffix;
            totgoodbad = "wpnegmtotgoodbad" + suffix;
            # hist_ttbar.Add(hist_wjets);
            asymmetryplot(hist_ttbargood, hist_ttbar, ttbname);
            asymmetryplot(hist_wprimegood, hist_wpm600, wpname);
            asymmetry2plot(hist_ttbargood, hist_ttbar, hist_wprimegood, hist_wpm600, comasymname);
            plusminuscmp(hist_wpm600, hist_wprimegood, wpgoodbad, 0);
            plusminuscmp(hist_ttbar, hist_ttbargood, ttbgoodbad, 0);
            hist_wpm600.Add(hist_ttbar);
            hist_wprimegood.Add(hist_ttbargood);
            plusminuscmp(hist_wpm600, hist_wprimegood, totgoodbad, 0);
        elif (histname == "wpnegmasspos"):
            hist_wprimepos = hist_wpm600.Clone("test");
        elif (histname == "wpnegmassneg"):
            plusminuscmp(hist_wprimepos, hist_wpm600, 'wpmassplusminus', 1);


        # canvases[-1].SaveAs('plots/' + histname + '.png')
        # canvases[-1].SaveAs('/storage/results/' + histname + '.png')

        cu_hist_data = getCumulativePlot(hist_data, "data");
        cu_hist_ttbar = getCumulativePlot(hist_ttbar, "ttbar");
        # cu_hist_wjets = getCumulativePlot(hist_wjets, "wjets");
        # cu_hist_zjets = getCumulativePlot(hist_zjets, "zjets");
        # cu_hist_qcd = getCumulativePlot(hist_qcd, "qcd");
        # cu_hist_singleTop = getCumulativePlot(hist_singleTop, "singleTop");
        # cu_hist_wpm400 =  getCumulativePlot(hist_wpm400, "Wprimet400");
        cu_hist_wpm600 =  getCumulativePlot(hist_wpm600, "Wprimet600");
        # cu_hist_wpm800 =  getCumulativePlot(hist_wpm800, "Wprimet800");
        # cu_hist_wpm1000 =  getCumulativePlot(hist_wpm1000, "Wprimet1000");
        # cu_hist_wptb1000 =  getCumulativePlot(hist_wptb1000, "Wprimetb1000");
        # cu_hist_Zprime500 = getCumulativePlot(hist_Zprime500, "Zprime500");
        # cu_hist_Zprime750 = getCumulativePlot(hist_Zprime750, "Zprime750");
        # cu_hist_Zprime1000 = getCumulativePlot(hist_Zprime1000, "Zprime1000");
        # cu_hist_Zprime2000 = getCumulativePlot(hist_Zprime2000, "Zprime2000");
        # cu_hist_Zprime4000 = getCumulativePlot(hist_Zprime4000, "Zprime4000");
        # cu_hist_Tprime300 = getCumulativePlot(hist_Tprime300, "Tprime300");
        # cu_hist_Tprime350 = getCumulativePlot(hist_Tprime350, "Tprime350");
##        cu_hist_Zprime1250 = getCumulativePlot(hist_Zprime1250, "Zprime1250");
##        cu_hist_Zprime1500 = getCumulativePlot(hist_Zprime1500, "Zprime1500");

        cu_hist_data.SetTitleSize(0.04, "X");
        cu_hist_data.SetTitleSize(0.04, "Y");
        cu_hist_data.SetYTitle("Integrated Events");
        cu_hist_data.SetTitleOffset(1.3, "X");
        cu_hist_data.SetTitleOffset(1.7, "Y");
##        
#
        cu_hist_data.SetAxisRange(Urange[0], Urange[1]);
        cu_hist_ttbar.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_wjets.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_zjets.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_qcd.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_singleTop.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_wpm400.SetAxisRange(Urange[0], Urange[1]);
        cu_hist_wpm600.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_wpm800.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_wpm1000.SetAxisRange(Urange[0], Urange[1]);
        # cu_hist_wptb1000.SetAxisRange(Urange[0], Urange[1]);
        if (Yrange[1] > 0):
            cu_hist_data.SetAxisRange(Yrange[0], Yrange[1], "Y");
            cu_hist_ttbar.SetAxisRange(Yrange[0], Yrange[1], "Y");
            cu_hist_wpm600.SetAxisRange(Yrange[0], Yrange[1], "Y");
        
        cu_hs = THStack("cu_MC", "cu_MC");
        # cu_hs.Add(cu_hist_qcd);
        # cu_hs.Add(cu_hist_zjets);
        # cu_hs.Add(cu_hist_wjets);
        cu_hs.Add(cu_hist_ttbar);
        
        scanvases.append(TCanvas("cu_cname" + histname, histname + "(cu)", 1200, 900))

        # if ('iso' in histname or 'jetmultip' in histname or 'num' in histname and 'Top' in histname):
        scanvases[-1].SetGrid();
        if ('muiso' in histname):
          scanvases[-1].SetLogy();
          
        scanvases[-1].cd().SetRightMargin(0.04);
        # cu_hist_data.Draw("error");
        # cu_hist_data.Draw("axis");
        cu_hist_data.Draw();
        cu_hs.Draw("hist same");
        cu_hist_wpm600.Draw("same");
        # if (not 'num' in histname or not 'Top' in histname):
          # cu_hist_wpm400.Draw("same");
          # cu_hist_wpm800.Draw("same");
          # cu_hist_wpm1000.Draw("same");
        # cu_hist_wptb1000.Draw("same");
        # cu_hist_Zprime500.Draw("same");
        # cu_hist_Zprime750.Draw("same");
        # cu_hist_Zprime1000.Draw("same");
        # cu_hist_Zprime2000.Draw("same");
        # cu_hist_Zprime4000.Draw("same");
        # cu_hist_Tprime300.Draw("same");
        # cu_hist_Tprime350.Draw("same");
###        cu_hist_Zprime1250.Draw("same");
###        cu_hist_Zprime1500.Draw("same");
##        #        cu_hist_data2.Draw("error same");
        # cu_hist_data.Draw("error same");
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

def setCorrErrs(hist_good, hist_bad, hist_diff):
    nbins = hist_diff.GetXaxis().GetNbins();
    for i in range(1, nbins):
      goodval = hist_good.GetBinContent(i);
      badval = hist_bad.GetBinContent(i);
      gooderr = hist_good.GetBinError(i);
      baderr = hist_bad.GetBinError(i);
      errsum = goodval + badval;
      if errsum > 0:
        # print 'evts err ', goodval, gooderr;
        goodprod = goodval * baderr;
        badprod = gooderr * badval;
        tot =  (goodprod * goodprod) + (badprod * badprod);
        # Assume 90% correlation
        tot = tot - (1.8 * goodprod * badprod);
        if tot > 0:
          tot = 4.0 * sqrt(tot) / (errsum * errsum);
          hist_diff.SetBinError(i, tot);



def setUncorrErrs(hist_good, hist_bad, hist_diff):
    nbins = hist_diff.GetXaxis().GetNbins();
    for i in range(1, nbins):
      goodval = hist_good.GetBinContent(i);
      badval = hist_bad.GetBinContent(i);
      err = goodval + badval;
      err = err * err * err;
      if err > 0:
        err = 4.0 * sqrt((goodval * badval) / err);
        hist_diff.SetBinError(i, err);


def asymmetryplot(hist_good, hist_bad, histname):
    hist_total = hist_good.Clone("test");
    # hist_total.Sumw2();
    hist_total.Add(hist_bad);
    hist_diff = hist_good.Clone("test");
    # hist_diff.Sumw2();
    hist_diff.Add(hist_bad, -1);
    hist_total.Scale(0.5);
    hist_diff.Divide(hist_total);
    # setCorrErrs(hist_good, hist_bad, hist_diff);
    hist_diff.SetXTitle("W' mass/GeV");
    hist_diff.SetYTitle("Asymmetry");
    hist_diff.SetAxisRange(200, 1100);
    hist_diff.SetLineWidth(4);
    # hist_diff.SetFillStyle(1001);
    hist_diff.SetFillStyle(3001);
    hist_diff.SetFillColor(kRed);
    canvases[-1].SetGridy();
    hist_diff.GetYaxis().SetRangeUser(-asymrange, asymrange);
    hist_diff.Draw("hist E1");
    canvases[-1].SaveAs('plots/' + histname + '.png')


def asymmetry2plot(ttbar_good, ttbar_bad, wp_good, wp_bad, histname):
    hist_total = ttbar_good.Clone("test");
    # hist_total.Sumw2();
    hist_total.Add(ttbar_bad);
    hist_total.Add(wp_good);
    hist_total.Add(wp_bad);
    hist_total.Scale(0.5);
    hist_diff = ttbar_good.Clone("test");
    # hist_diff.Sumw2();
    hist_diff.Add(wp_good);
    hist_diff.Add(ttbar_bad, -1);
    hist_diff.Add(wp_bad, -1);
    hist_diff.Divide(hist_total);
    # nbins = hist_diff.GetXaxis().GetNbins();
    hist_good = ttbar_good.Clone("test");
    hist_good.Add(wp_good);
    hist_bad = ttbar_bad.Clone("test");
    hist_bad.Add(wp_bad);
    # setCorrErrs(hist_good, hist_bad, hist_diff);
    hist_diff.SetXTitle("W' mass/GeV");
    hist_diff.SetYTitle("Asymmetry");
    hist_diff.SetAxisRange(200, 1100);
    hist_diff.SetLineWidth(4);
    # hist_diff.SetFillStyle(1001);
    hist_diff.SetFillStyle(3001);
    hist_diff.SetFillColor(kRed);
    canvases[-1].SetGridy();
    hist_diff.GetYaxis().SetRangeUser(-asymrange, asymrange);
    hist_diff.Draw("hist E1");
    canvases[-1].SaveAs('plots/' + histname + '.png')



def plusminuscmp(hist_pos, hist_neg, hist_name, sign):
    hist_neg.SetXTitle("Reco W' mass/GeV");
    hist_neg.SetYTitle("Events/(60 GeV)");
    hist_pos.SetLineColor(kGray + 3);
    hist_neg.SetLineColor(kRed);
    hist_neg.SetFillStyle(0);
    hist_pos.SetFillStyle(0);
    linewidth = 6;
    hist_pos.SetLineWidth(linewidth);
    hist_neg.SetLineWidth(linewidth);
    leg = TLegend(0.696, 0.35, 0.95, 0.92);
    leg.SetBorderSize(0);
    leg.SetLineStyle(0);
    leg.SetTextFont(42);
    leg.SetFillStyle(0);
    if sign == 1:
      leg.AddEntry(hist_pos, "W't+ .6TeV (8 pb)", "L");
      leg.AddEntry(hist_neg, "W't- .6TeV (8 pb)", "L");
    else:
      leg.AddEntry(hist_pos, "W't bad .6TeV (8 pb)", "L");
      leg.AddEntry(hist_neg, "W't good .6TeV (8 pb)", "L");
    hist_pos.SetStats(1);
    hist_neg.Draw("hist E1");
    hist_pos.Draw("hist same");
    leg.Draw();
    canvases[-1].SaveAs('plots/' + hist_name + '.png')


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
    cu.SetMarkerSize(initial.GetMarkerSize());
    cu.SetMarkerStyle(initial.GetMarkerStyle());
    cu.SetMarkerColor(initial.GetMarkerColor());
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
