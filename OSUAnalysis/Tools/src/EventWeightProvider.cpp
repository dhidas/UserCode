/*
 * CrossSections.cpp
 *
 *  Created on: 13 Aug 2010
 *      Author: kreczko
 */

#include "../interface/EventWeightProvider.h"
#include "TFile.h"
#include <boost/scoped_ptr.hpp>
#include <iostream>

namespace BAT {

boost::array<float, DataType::NUMBER_OF_DATA_TYPES> sevenTeV::getXSections() {
    boost::array<float, DataType::NUMBER_OF_DATA_TYPES> xsection;
    xsection[DataType::DATA] = 0;
    xsection[DataType::ttbar] = 157.5;
    xsection[DataType::ttjets] = 157.5;

    xsection[DataType::Zjets] = 3048.;
		// From twiki.cern.ch/twiki/bin/view/CMS/CrossSectionDetails
		// For DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola

    xsection[DataType::Wjets] = 31314.;
		// From twiki.cern.ch/twiki/bin/view/CMS/CrossSectionDetails
		// for WJetsToLNu_TuneZ2_7TeV-madgraph-tauola
	
    // xsection[DataType::Wjets] = 24640.;
    xsection[DataType::WToENu] = 7899.;

    xsection[DataType::QCD_Flat_15to3000] = 2.213e5;	// Fake number 
    // xsection[DataType::QCD_Flat_15to3000] = 2.213e10;  // Too big
		// Need per-event weighting
		// From twiki.cern.ch/twiki/bin/view/CMS/ProductionSummer2011
		// for QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6

    xsection[DataType::QCD_EMEnriched_Pt20to30] = 0.2355e9 * 0.0073;//xs 0.2355 mb (filter efficiency=0.0073)
    xsection[DataType::QCD_EMEnriched_Pt30to80] = 0.0593e9 * 0.059; //xs 0.0593 mb
    // xsection[DataType::QCD_EMEnriched_Pt80to170] = 0.906e6 * 0.148; //xs 0.906e-3 mb, total 134088 pb
    xsection[DataType::QCD_EMEnriched_Pt80to170] = 139500; // From twiki.cern.ch/twiki/bin/view/CMS/CrossSectionDetails

    xsection[DataType::QCD_MuEnriched_Pt20Pt15] = 84680; // From twiki.cern.ch/twiki/bin/view/CMS/CrossSectionDetails

    xsection[DataType::QCD_BCtoE_Pt20to30] = 0.2355e9 * 0.00046; //xs 0.2355 mb (filter efficiency=0.00046)
    // xsection[DataType::QCD_BCtoE_Pt30to80] = 0.0593e9 * 0.00234; //xs 0.0593 mb, total 138762 pb
    xsection[DataType::QCD_BCtoE_Pt30to80] = 136800;	// From twiki.cern.ch/twiki/bin/view/CMS/CrossSectionDetails
    // xsection[DataType::QCD_BCtoE_Pt80to170] = 0.906e6 * 0.0104; //xs 0.906e-3 mb, total 9422 pb
    xsection[DataType::QCD_BCtoE_Pt80to170] = 9360;	// From twiki.cern.ch/twiki/bin/view/CMS/CrossSectionDetails

    xsection[DataType::PhotonJets_Pt40to100] = 23620.; //pb
    xsection[DataType::PhotonJets_Pt100to200] = 3476.; //pb
    xsection[DataType::PhotonJets_Pt200toInf] = 485.; //pb

    xsection[DataType::WWtoAnything] = 43.; //pb
    xsection[DataType::WZtoAnything] = 18.; //pb
    xsection[DataType::ZZtoAnything] = 5.9; //pb

    xsection[DataType::singleTop_And_W] = 10.6; //xs  11 pb (NLO MCFM) inclusive t,W decay
    xsection[DataType::singleTopTChannel] = 21.53;
    xsection[DataType::singleTopSChannel] = 1.40; //=4.21/3 15Jul

    xsection[DataType::VQQ] = 36.;
    xsection[DataType::Zprime_M500GeV_W5GeV] = 50;
    xsection[DataType::Zprime_M500GeV_W50GeV] = 50;
    xsection[DataType::Zprime_M750GeV_W7500MeV] = 50;
    xsection[DataType::Zprime_M1TeV_W10GeV] = 50;
    xsection[DataType::Zprime_M1TeV_W100GeV] = 50;
    xsection[DataType::Zprime_M1250GeV_W12500MeV] = 50;
    xsection[DataType::Zprime_M1500GeV_W15GeV] = 50;
    xsection[DataType::Zprime_M1500GeV_W150GeV] = 50;
    xsection[DataType::Zprime_M2TeV_W20GeV] = 50;
    xsection[DataType::Zprime_M2TeV_W200GeV] = 50;
    xsection[DataType::Zprime_M3TeV_W30GeV] = 50;
    xsection[DataType::Zprime_M3TeV_W300GeV] = 50;
    xsection[DataType::Zprime_M4TeV_W40GeV] = 50;
    xsection[DataType::Zprime_M4TeV_W400GeV] = 50;

    //cross section S-channel only
    
    //xsection[DataType::WprimeTToTTD_M600] = 8;
    // xsection[DataType::WprimeTToTTD_M800] = 2.2;
    //xsection[DataType::WprimeTToTTD_M1000] = 0.72;
    xsection[DataType::WprimeToTBbar_M1000] = 8;
    xsection[DataType::WprimeTToTTD_M600] = 18.2;
    xsection[DataType::WprimeTToTTD_M800] = 6.5;


    xsection[DataType::WprimeTToTTD_M1000] = 0.72;

    return xsection;
}

EventWeightProvider::EventWeightProvider(float lumiInInversePb, unsigned short tev, std::string pileUpEstimationFile) :
    lumiInInversePb(lumiInInversePb),
    tev(tev),
    useSkimEff(true),
    xsection(),
    estimatedPileUp(getPileUpHistogram(pileUpEstimationFile)),
    pileUpWeights(),
    numberOfEventsWithTooHighPileUp(0),
    numberOfProcessedEvents(),
    numberOfPattplSkimEvents(),
    numberOfNtplSkimEvents(),
    numberOfElectronSkimEvents(),
    numberOfMuonSkimEvents()
{
    generate_flat10_weights();
    if (tev == 7)
        xsection = sevenTeV::getXSections();
    defineNumberOfProducedEvents();
//    defineNumberOfSkimmedEvents();
}

void EventWeightProvider::defineNumberOfProducedEvents() {
    numberOfProcessedEvents[DataType::DATA] = 0;

    
		/*
    numberOfProcessedEvents[DataType::ttjets] = 3688248;  // V8 PATtuple
    numberOfPattplSkimEvents[DataType::ttjets] = 3683795;  // V8 PATtuple
    numberOfNtplSkimEvents[DataType::ttjets] = 1173710;  // V8 PATtuple
    numberOfElectronSkimEvents[DataType::ttjets] = 727781;	  // V8 PATtuple
    numberOfMuonSkimEvents[DataType::ttjets] = 510117;	// V8 PATtuple
		*/
		// PATtuple /TTJets_TuneZ2_7TeV-madgraph-tauola/srappocc-ttbsm_v8_Summer11-PU_S4_START42_V11-v1-5c91b0700768331a44f51c8a9892d391/USER 

		numberOfProcessedEvents[DataType::ttjets] = 3701947 ;  // V9 PATtuple
    numberOfPattplSkimEvents[DataType::ttjets] = 3697476;  // V9 PATtuple
    numberOfNtplSkimEvents[DataType::ttjets] = 1348512;  // V9 PATtuple
    numberOfElectronSkimEvents[DataType::ttjets] = 730452;	  // V9 PATtuple
    numberOfMuonSkimEvents[DataType::ttjets] = 716210;	// V9 PATtuple
		// PATtuple /TTJets_TuneZ2_7TeV-madgraph-tauola/srappocc-ttbsm_v9_Summer11-PU_S4_START42_V11-v1-bf57a985b107a689982b667a3f2f23c7/USER 

		/*
    // numberOfProcessedEvents[DataType::ttjets] = 3701947;	// v2 AOD ntuple original dataset number
    numberOfProcessedEvents[DataType::ttjets] = 3664929;	// v2 AOD ntuple, two jobs failed
    numberOfPattplSkimEvents[DataType::ttjets] = 3664929;	// v2 AOD ntuple
    numberOfNtplSkimEvents[DataType::ttjets] = 1687141;	// v2 AOD ntuple
    numberOfElectronSkimEvents[DataType::ttjets] = 1123778;	// v2 AOD ntuple
    numberOfMuonSkimEvents[DataType::ttjets] = 738210;	// v2 AOD ntuple
*/

    // numberOfProcessedEvents[DataType::ttjets] = 3652589;	// TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v2_AODSIM_CMSSW_4_2_4_Ov3.5 04/08/2011
    // numberOfNtplSkimEvents[DataType::ttjets] = 1696252;	// Previous ntuple collection
    // numberOfElectronSkimEvents[DataType::ttjets] = 1129134;
    // numberOfMuonSkimEvents[DataType::ttjets] = 742786;

    // numberOfProcessedEvents[DataType::ttbar] = 3570035;
    numberOfProcessedEvents[DataType::ttbar] = 1089625; // TT_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2_AODSIM. Checked 21/07/11
    // numberOfProcessedEvents[DataType::Zjets] = 34016401;
    numberOfProcessedEvents[DataType::Zjets] = 32512091;
    numberOfPattplSkimEvents[DataType::Zjets] = 32475188;
    numberOfNtplSkimEvents[DataType::Zjets] = 15161964;
    numberOfElectronSkimEvents[DataType::Zjets] = 7503990;
    numberOfMuonSkimEvents[DataType::Zjets] = 7680017;
		// Previous values before LoosePFlow
    // numberOfNtplSkimEvents[DataType::Zjets] = 16354092;
    // numberOfElectronSkimEvents[DataType::Zjets] = 8390919;
    // numberOfMuonSkimEvents[DataType::Zjets] = 8056270;
		// PATtuple V8 /DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/skhalil-ttbsm_v8_Summer11-PU_S4_START42_V11-v1-87037ef7c828ea57e128f1ace23a632e/USER

    // numberOfProcessedEvents[DataType::Wjets] = 54895290;	// Old AOD dataset
    // numberOfNtplSkimEvents[DataType::Wjets] = 17166386; // previous
    // numberOfElectronSkimEvents[DataType::Wjets] = 8630776; // previous
    // numberOfMuonSkimEvents[DataType::Wjets] = 8587398; 	// previous

    // numberOfProcessedEvents[DataType::Wjets] = 49484941;
    // numberOfPattplSkimEvents[DataType::Wjets] = 49335978;
    // numberOfNtplSkimEvents[DataType::Wjets] = 15175745;
    // numberOfElectronSkimEvents[DataType::Wjets] = 7365443;
    // numberOfMuonSkimEvents[DataType::Wjets] = 7817858;
		// PATtuple V8 /WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/srappocc-ttbsm_v8_Summer11-PU_S4_START42_V11-v1-87037ef7c828ea57e128f1ace23a632e/USER 
    numberOfProcessedEvents[DataType::Wjets] = 77105816;
    numberOfPattplSkimEvents[DataType::Wjets] = 76978604;
    numberOfNtplSkimEvents[DataType::Wjets] = 24580667;
    numberOfElectronSkimEvents[DataType::Wjets] = 11449725;
    numberOfMuonSkimEvents[DataType::Wjets] = 13144761;
		// PATtuple V9 /WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/dstrom-prod_2011_10_05_17_14_11-bf57a985b107a689982b667a3f2f23c7/USER

    numberOfProcessedEvents[DataType::WToENu] = 5334220;	// Checked 21/07/11
    numberOfPattplSkimEvents[DataType::WToENu] = 5334220;
    numberOfNtplSkimEvents[DataType::WToENu] = 2662025;
    numberOfElectronSkimEvents[DataType::WToENu] = 2662025;
    numberOfMuonSkimEvents[DataType::WToENu] = 0;

    numberOfProcessedEvents[DataType::QCD_Flat_15to3000] = 10960800;
    // numberOfPattplSkimEvents[DataType::QCD_Flat_15to3000] = 10929635; // According to PATtple twiki
    numberOfPattplSkimEvents[DataType::QCD_Flat_15to3000] = 10751721;	// Actual received
    numberOfNtplSkimEvents[DataType::QCD_Flat_15to3000] = 163974;
    numberOfElectronSkimEvents[DataType::QCD_Flat_15to3000] = 163906;
    numberOfMuonSkimEvents[DataType::QCD_Flat_15to3000] = 70;
		// Previous values before LoosePFlow
    // numberOfNtplSkimEvents[DataType::QCD_Flat_15to3000] = 1010202;
    // numberOfElectronSkimEvents[DataType::QCD_Flat_15to3000] = 775300;
    // numberOfMuonSkimEvents[DataType::QCD_Flat_15to3000] = 256834;
		// PATtuple V8 /QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6/srappocc-ttbsm_v8_Summer11-PU_S3_-START42_V11-v2-d870fa9b0dd695e8eb649b7e725d070f/USER
    numberOfProcessedEvents[DataType::QCD_EMEnriched_Pt20to30] = 22529376;
    numberOfProcessedEvents[DataType::QCD_EMEnriched_Pt30to80] = 36409308;

    numberOfProcessedEvents[DataType::QCD_EMEnriched_Pt80to170] = 8150672;
    numberOfPattplSkimEvents[DataType::QCD_EMEnriched_Pt80to170] = 5078282;	// Actual received
    numberOfNtplSkimEvents[DataType::QCD_EMEnriched_Pt80to170] = 123843;
    numberOfElectronSkimEvents[DataType::QCD_EMEnriched_Pt80to170] = 123744;
    numberOfMuonSkimEvents[DataType::QCD_EMEnriched_Pt80to170] = 101;

    numberOfProcessedEvents[DataType::QCD_BCtoE_Pt20to30] = 2081560;	// Checked 21/07/11

    numberOfProcessedEvents[DataType::QCD_BCtoE_Pt30to80] = 2030033;	// Checked 19/08/11
    numberOfPattplSkimEvents[DataType::QCD_BCtoE_Pt30to80] = 2004505;	// Actual received
    numberOfNtplSkimEvents[DataType::QCD_BCtoE_Pt30to80] = 98332;
    numberOfElectronSkimEvents[DataType::QCD_BCtoE_Pt30to80] = 98318;
    numberOfMuonSkimEvents[DataType::QCD_BCtoE_Pt30to80] = 17;

    numberOfProcessedEvents[DataType::QCD_BCtoE_Pt80to170] = 1082691;	// Checked 19/08/11
    numberOfPattplSkimEvents[DataType::QCD_BCtoE_Pt80to170] = 362356;	// Actual received
    numberOfNtplSkimEvents[DataType::QCD_BCtoE_Pt80to170] = 63109;
    numberOfElectronSkimEvents[DataType::QCD_BCtoE_Pt80to170] = 63092;
    numberOfMuonSkimEvents[DataType::QCD_BCtoE_Pt80to170] = 22;

    numberOfProcessedEvents[DataType::QCD_MuEnriched_Pt20Pt15] = 24324525;	// Number missing, just use received
    numberOfPattplSkimEvents[DataType::QCD_MuEnriched_Pt20Pt15] = 24324525;	// Actual received
    numberOfNtplSkimEvents[DataType::QCD_MuEnriched_Pt20Pt15] = 2711630;
    numberOfElectronSkimEvents[DataType::QCD_MuEnriched_Pt20Pt15] = 80370;
    numberOfMuonSkimEvents[DataType::QCD_MuEnriched_Pt20Pt15] = 2650602;

    numberOfProcessedEvents[DataType::PhotonJets_Pt40to100] = 2217101;
    numberOfProcessedEvents[DataType::PhotonJets_Pt100to200] = 1065691;
    numberOfProcessedEvents[DataType::PhotonJets_Pt200toInf] = 1079950;

    numberOfProcessedEvents[DataType::WWtoAnything] = 2039440;
    numberOfProcessedEvents[DataType::WZtoAnything] = 2085696;
    numberOfProcessedEvents[DataType::ZZtoAnything] = 2108608;

    numberOfProcessedEvents[DataType::singleTop_And_W] = 489417;
    numberOfProcessedEvents[DataType::singleTopTChannel] = 484060;
    numberOfProcessedEvents[DataType::singleTopSChannel] = 494967;
    numberOfProcessedEvents[DataType::VQQ] = 720613;
    numberOfProcessedEvents[DataType::Zprime_M500GeV_W5GeV] = 227068;
    numberOfProcessedEvents[DataType::Zprime_M500GeV_W50GeV] = 238963;
    numberOfProcessedEvents[DataType::Zprime_M750GeV_W7500MeV] = 204819;
    numberOfProcessedEvents[DataType::Zprime_M1TeV_W10GeV] = 213384;
    numberOfProcessedEvents[DataType::Zprime_M1TeV_W100GeV] = 200387;
    numberOfProcessedEvents[DataType::Zprime_M1250GeV_W12500MeV] = 233361;
    numberOfProcessedEvents[DataType::Zprime_M1500GeV_W15GeV] = 193779;
    numberOfProcessedEvents[DataType::Zprime_M1500GeV_W150GeV] = 199121;
    numberOfProcessedEvents[DataType::Zprime_M2TeV_W20GeV] = 238752;
    numberOfProcessedEvents[DataType::Zprime_M2TeV_W200GeV] = 213363;
    numberOfProcessedEvents[DataType::Zprime_M3TeV_W30GeV] = 205270;
    numberOfProcessedEvents[DataType::Zprime_M3TeV_W300GeV] = 229034;
    numberOfProcessedEvents[DataType::Zprime_M4TeV_W40GeV] = 183920;
    numberOfProcessedEvents[DataType::Zprime_M4TeV_W400GeV] = 238142;

    numberOfProcessedEvents[DataType::WprimeTToTTD_M400] = 96990;	// V9
    numberOfPattplSkimEvents[DataType::WprimeTToTTD_M400] = 96698;
    numberOfNtplSkimEvents[DataType::WprimeTToTTD_M400] = 37674;
    numberOfElectronSkimEvents[DataType::WprimeTToTTD_M400] = 21097;
    numberOfMuonSkimEvents[DataType::WprimeTToTTD_M400] = 19681;

    // numberOfProcessedEvents[DataType::WprimeTToTTD_M600] = 96500;	// Muonless generation
    // numberOfPattplSkimEvents[DataType::WprimeTToTTD_M600] = 34552;	// After generator filter
    // numberOfNtplSkimEvents[DataType::WprimeTToTTD_M600] = 25535;
    // numberOfElectronSkimEvents[DataType::WprimeTToTTD_M600] = 22774;

    // numberOfMuonSkimEvents[DataType::WprimeTToTTD_M600] = 7012;//
    //low stats sample
    //    numberOfProcessedEvents[DataType::WprimeTToTTD_M600] = 86991;	// With muons
    // numberOfPattplSkimEvents[DataType::WprimeTToTTD_M600] = 86601;
    //numberOfNtplSkimEvents[DataType::WprimeTToTTD_M600] = 29439;
    // numberOfElectronSkimEvents[DataType::WprimeTToTTD_M600] = 19580;
    //numberOfMuonSkimEvents[DataType::WprimeTToTTD_M600] = 11656;
    //850k sample
    //    numberOfProcessedEvents[DataType::WprimeTToTTD_M600] =849920;	// With muons
    //numberOfPattplSkimEvents[DataType::WprimeTToTTD_M600] = 849920;
    //numberOfNtplSkimEvents[DataType::WprimeTToTTD_M600] = 340080;
    //numberOfElectronSkimEvents[DataType::WprimeTToTTD_M600] = 19580;
    //numberOfMuonSkimEvents[DataType::WprimeTToTTD_M600] = 11656;
    //500k sample
    numberOfProcessedEvents[DataType::WprimeTToTTD_M600] =499960;  // With muons                                                                                            
    numberOfPattplSkimEvents[DataType::WprimeTToTTD_M600] = 499960;
    numberOfNtplSkimEvents[DataType::WprimeTToTTD_M600] = 199893;
    numberOfElectronSkimEvents[DataType::WprimeTToTTD_M600] = 19580;
    numberOfNtplSkimEvents[DataType::WprimeTToTTD_M600] = 29439;	// V8 prescription
    numberOfMuonSkimEvents[DataType::WprimeTToTTD_M600] = 11656;  // V8
    // numberOfNtplSkimEvents[DataType::WprimeTToTTD_M600] = 35447;// V9
    // numberOfMuonSkimEvents[DataType::WprimeTToTTD_M600] = 19127; // V9



    //low stats
    //    numberOfProcessedEvents[DataType::WprimeTToTTD_M800] = 84993;
    //numberOfPattplSkimEvents[DataType::WprimeTToTTD_M800] = 84416;
    //numberOfNtplSkimEvents[DataType::WprimeTToTTD_M800] = 29269;
    //numberOfElectronSkimEvents[DataType::WprimeTToTTD_M800] = 19619;
    //numberOfMuonSkimEvents[DataType::WprimeTToTTD_M800] = 11363;
    //high stats 700k
    // numberOfProcessedEvents[DataType::WprimeTToTTD_M800] = 689950;
    //numberOfPattplSkimEvents[DataType::WprimeTToTTD_M800] = 689950;
    //numberOfNtplSkimEvents[DataType::WprimeTToTTD_M800] = 283455;
    //numberOfElectronSkimEvents[DataType::WprimeTToTTD_M800] = 19619;
    //numberOfMuonSkimEvents[DataType::WprimeTToTTD_M800] = 11363;
    //high stats 500k
    numberOfProcessedEvents[DataType::WprimeTToTTD_M800] = 489965;
    numberOfPattplSkimEvents[DataType::WprimeTToTTD_M800] = 489965;
    numberOfNtplSkimEvents[DataType::WprimeTToTTD_M800] = 201242;

    numberOfElectronSkimEvents[DataType::WprimeTToTTD_M800] = 19619;

    // numberOfProcessedEvents[DataType::WprimeTToTTD_M1000] = 93000;	// From Khristian's e-mail, muonless generation
    // numberOfPattplSkimEvents[DataType::WprimeTToTTD_M1000] = 35600;	// After generator filter
    // numberOfNtplSkimEvents[DataType::WprimeTToTTD_M1000] = 26563;
    // numberOfElectronSkimEvents[DataType::WprimeTToTTD_M1000] = 23025;
    // numberOfMuonSkimEvents[DataType::WprimeTToTTD_M1000] = 8849;
    numberOfProcessedEvents[DataType::WprimeTToTTD_M1000] = 94996; // With muons
    numberOfPattplSkimEvents[DataType::WprimeTToTTD_M1000] = 94145;
    numberOfElectronSkimEvents[DataType::WprimeTToTTD_M1000] = 21854;
    // numberOfNtplSkimEvents[DataType::WprimeTToTTD_M1000] = 40756;  // V9
    // numberOfMuonSkimEvents[DataType::WprimeTToTTD_M1000] = 22898;  // V9
    numberOfNtplSkimEvents[DataType::WprimeTToTTD_M1000] = 32039;  // V8
    numberOfMuonSkimEvents[DataType::WprimeTToTTD_M1000] = 12054;  // V8

    numberOfProcessedEvents[DataType::WprimeToTBbar_M1000] = 110000;	// Checked 21/07/11
}

EventWeightProvider::~EventWeightProvider() {

}

double EventWeightProvider::getWeight(DataType::value type) {
    if (type == DataType::DATA)
        return 1.;
    else
        return xsection[type] * lumiInInversePb / numberOfProcessedEvents[type];
}

float EventWeightProvider::reweightPileUp(unsigned int numberOfVertices){
    if(numberOfVertices >= pileUpWeights.size()){
        ++numberOfEventsWithTooHighPileUp;
        return 0.;
    }

    return pileUpWeights.at(numberOfVertices);
}

boost::shared_ptr<TH1D> EventWeightProvider::getPileUpHistogram(std::string pileUpEstimationFile){
    std::cout << "Using pile-up estimation file " << pileUpEstimationFile << std::endl;
    boost::scoped_ptr<TFile> file(TFile::Open(pileUpEstimationFile.c_str()));
    boost::shared_ptr<TH1D> pileUp((TH1D*) file->Get("pileup")->Clone());
    file->Close();
    return pileUp;
}

void EventWeightProvider::generate_flat10_weights(){
    // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
    //only use this for Wprime with true Mc compared to true data Cert_160404-180252_7TeV_PromptReco_Collisons11_JSON.pileup_v2.root
  //  const boost::array<double, 25> npu_probs = {{0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584 /* <-- 10*/,
  //       0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322,
  //       0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 /* <-- 24 */}};


  //USE THIS FOR TTBAR with Cert_160404_180252_7TeV_Collisions11_JSON.pileup_v2.root
  // Summer11 PU_S3 and PU_S4, distribution obtained by only looking at the in-time crossing.  This is the "spike+smear" distribution, RECOMMENDED FOR REWEIGHTING.
  const boost::array<double, 35> npu_probs = {{1.45346E-01,    6.42802E-02,    6.95255E-02,    6.96747E-02,    6.92955E-02,    6.84997E-02,    6.69528E-02,    6.45515E-02,    6.09865E-02,    5.63323E-02,    5.07322E-02,    4.44681E-02,    3.79205E-02,    3.15131E-02,    2.54220E-02,    2.00184E-02,    1.53776E-02,    1.15387E-02,    8.47608E-03,    6.08715E-03,    4.28255E-03,    2.97185E-03,    2.01918E-03,    1.34490E-03,    8.81587E-04,    5.69954E-04,    3.61493E-04,    2.28692E-04,    1.40791E-04,    8.44606E-05,    5.10204E-05,    3.07802E-05,    1.81401E-05,    1.00201E-05,    5.80004E-06  }};


    double s = 0.0;
    for (unsigned int npu = 0; npu < npu_probs.size(); ++npu) {
        double npu_estimated = estimatedPileUp->GetBinContent(estimatedPileUp->GetXaxis()->FindBin(npu));
        pileUpWeights[npu] = npu_estimated / npu_probs[npu];
	//	std::cout<<npu_estimated
	//std::cout<<npu<<std::endl;
        s += npu_estimated;
    }
    // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    for (unsigned int npu = 0; npu < pileUpWeights.size(); ++npu) {
        pileUpWeights[npu] /= s;
    }
}

unsigned long EventWeightProvider::getNumberOfEventsWithTooHighPileUp() const{
    return numberOfEventsWithTooHighPileUp;
}

} // namespace BAT
