#include "DQMServices/Diagnostic/test/HDQMInspector.h"
#include "DQM/TrackingHistoricInfoClient/test/HDQMInspectorConfigTracking.h"
#include <string>

void TrackingHDQMInspector (const std::string & tagName, std::string const& Password, int const NRuns) {
/////////////////////////////////////////////////////////////////
//
// Extraction of the summary information using 
// DQMServices/Diagnostic/test/HDQMInspector.
// The sqlite database should have been filled using the new
// TrackingHistoryDQMService.   
//
/////////////////////////////////////////////////////////////////


  //std::map<int, std::string> pixelTranslator = sipixelsummary::GetMap();

  //pixelTranslator Translator;

  //AutoLibraryLoader::enable();

  HDQMInspectorConfigTracking StripConfig;
  std::vector<std::string> ItemsForIntegration;
  ItemsForIntegration.push_back("Chi2overDoF_CosmicTk_entries");
  ItemsForIntegration.push_back("NumberOfTracks_CosmicTk_entries");
  StripConfig.computeIntegralList(ItemsForIntegration);
  HDQMInspector A(&StripConfig);
  //HDQMInspector A;
  //A.setDB("sqlite_file:dbfile.db",tagName,"cms_cond_strip","w3807dev","");
  A.setDB("oracle://cms_orcoff_prep/CMS_DQM_31X_OFFLINE",tagName,"cms_dqm_31x_offline", Password,"");


  A.setDebug(1);
  A.setDoStat(1);



  //A.setBlackList("68286");
  // 268435456
  A.createTrendLastRuns("268435456@NumberOfTracks_CosmicTk@entries", "NumberOfTracks_CosmicTk_entries.gif", 0, "268435456@NumberOfTracks_CosmicTk@entries > 0", NRuns);
  A.createTrendLastRuns("268435456@Chi2overDoF_CosmicTk@entries", "Chi2overDoF_CosmicTk_entries.gif", 0, "268435456@Chi2overDoF_CosmicTk@entries > 0", NRuns);
  A.createTrendLastRuns("268435456@NumberOfTracks_CKFTk@entries", "NumberOfTracks_CKFTk_entries.gif", 0, "268435456@NumberOfTracks_CKFTk@entries > 0", NRuns);
  A.createTrendLastRuns("268435456@NumberOfRecHitsPerTrack_CKFTk@entries", "NumberOfRecHitsPerTrack_CKFTk_entries.gif", 0, "268435456@NumberOfRecHitsPerTrack_CKFTk@entries > 0", NRuns);
  A.createTrendLastRuns("268435456@NumberOfTracks_CKFTk@mean", "NumberOfTracks_CKFTk_mean.gif", 0, "268435456@NumberOfTracks_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@NumberOfRecHitsPerTrack_CKFTk@entries", "NumberOfRecHitsPerTrack_CKFTk_entries.gif", 0, "268435456@NumberOfRecHitsPerTrack_CKFTk@entries > 0", NRuns);
  A.createTrendLastRuns("268435456@Chi2overDoF_CosmicTk@mean", "Chi2overDoF_CosmicTk.gif", 0, "268435456@Chi2overDoF_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@Chi2overDoF_CKFTk@mean", "Chi2overDoF_CKFTk.gif", 0, "268435456@Chi2overDoF_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@Chi2overDoF_RSTk@mean", "Chi2overDoF_RSTk.gif", 0, "268435456@Chi2overDoF_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@NumberOfTracks_CosmicTk@mean", "NumberOfTracks_CosmicTk.gif", 0, "268435456@NumberOfTracks_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@NumberOfRecHitsPerTrack_CosmicTk@mean", "NumberOfRecHitsPerTrack_CosmicTk.gif", 0, "268435456@NumberOfRecHitsPerTrack_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPt_CosmicTk@mean", "TrackPt_CosmicTk.gif", 0, "268435456@TrackPt_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPz_CosmicTk@mean", "TrackPz_CosmicTk.gif", 0, "268435456@TrackPz_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPx_CosmicTk@mean", "TrackPx_CosmicTk.gif", 0, "268435456@TrackPx_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPy_CosmicTk@mean", "TrackPy_CosmicTk.gif", 0, "268435456@TrackPy_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPhi_CosmicTk@mean", "TrackPhi_CosmicTk.gif", 0, "268435456@TrackPhi_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackEta_CosmicTk@mean", "TrackEta_CosmicTk.gif", 0, "268435456@TrackEta_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@DistanceOfClosestApproach_CosmicTk@mean", "DistanceOfClosestApproach_CosmicTk.gif", 0, "268435456@DistanceOfClosestApproach_CosmicTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPt_CKFTk@mean", "TrackPt_CKFTk.gif", 0, "268435456@TrackPt_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPz_CKFTk@mean", "TrackPz_CKFTk.gif", 0, "268435456@TrackPz_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPx_CKFTk@mean", "TrackPx_CKFTk.gif", 0, "268435456@TrackPx_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPy_CKFTk@mean", "TrackPy_CKFTk.gif", 0, "268435456@TrackPy_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPhi_CKFTk@mean", "TrackPhi_CKFTk.gif", 0, "268435456@TrackPhi_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackEta_CKFTk@mean", "TrackEta_CKFTk.gif", 0, "268435456@TrackEta_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@DistanceOfClosestApproach_CKFTk@mean", "DistanceOfClosestApproach_CKFTk.gif", 0, "268435456@DistanceOfClosestApproach_CKFTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@NumberOfTracks_RSTk@mean", "NumberOfTracks_RSTk.gif", 0, "268435456@NumberOfTracks_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@NumberOfRecHitsPerTrack_RSTk@mean", "NumberOfRecHitsPerTrack_RSTk.gif", 0, "268435456@NumberOfRecHitsPerTrack_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPt_RSTk@mean", "TrackPt_RSTk.gif", 0, "268435456@TrackPt_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPz_RSTk@mean", "TrackPz_RSTk.gif", 0, "268435456@TrackPz_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPx_RSTk@mean", "TrackPx_RSTk.gif", 0, "268435456@TrackPx_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPy_RSTk@mean", "TrackPy_RSTk.gif", 0, "268435456@TrackPy_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackPhi_RSTk@mean", "TrackPhi_RSTk.gif", 0, "268435456@TrackPhi_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@TrackEta_RSTk@mean", "TrackEta_RSTk.gif", 0, "268435456@TrackEta_RSTk@mean > 0", NRuns);
  A.createTrendLastRuns("268435456@DistanceOfClosestApproach_RSTk@mean", "DistanceOfClosestApproach_RSTk.gif", 0, "268435456@DistanceOfClosestApproach_RSTk@mean > 0", NRuns);

  A.closeFile();


  return;


}


void TrackingHDQMInspector (const std::string &tagName, std::string const& Password, int const Start, int const End) {
/////////////////////////////////////////////////////////////////
//
// Extraction of the summary information using 
// DQMServices/Diagnostic/test/HDQMInspector.
// The sqlite database should have been filled using the new
// TrackingHistoryDQMService.   
//
/////////////////////////////////////////////////////////////////


  //std::map<int, std::string> pixelTranslator = sipixelsummary::GetMap();

  //pixelTranslator Translator;

  //AutoLibraryLoader::enable();

  HDQMInspectorConfigTracking StripConfig;
  std::vector<std::string> ItemsForIntegration;
  ItemsForIntegration.push_back("Chi2overDoF_CosmicTk_entries");
  ItemsForIntegration.push_back("NumberOfTracks_CosmicTk_entries");
  StripConfig.computeIntegralList(ItemsForIntegration);
  HDQMInspector A(&StripConfig);
  //A.setDB("sqlite_file:dbfile.db",tagName,"cms_cond_strip","w3807dev","");
  A.setDB("oracle://cms_orcoff_prep/CMS_DQM_31X_OFFLINE",tagName,"cms_dqm_31x_offline", Password,"");


  A.setDebug(1);
  A.setDoStat(1);

  //A.setBlackList("68286");
  A.createTrend("268435456@NumberOfTracks_CosmicTk@entries", "NumberOfTracks_CosmicTk_entries.gif", 0, "268435456@NumberOfTracks_CosmicTk@entries > 0", Start, End);
  A.createTrend("268435456@Chi2overDoF_CosmicTk@entries", "Chi2overDoF_CosmicTk_entries.gif", 0, "268435456@Chi2overDoF_CosmicTk@entries > 0", Start, End);
  A.createTrend("268435456@NumberOfTracks_CKFTk@entries", "NumberOfTracks_CKFTk_entries.gif", 0, "268435456@NumberOfTracks_CKFTk@entries > 0", Start, End);
  A.createTrend("268435456@NumberOfRecHitsPerTrack_CKFTk@entries", "NumberOfRecHitsPerTrack_CKFTk_entries.gif", 0, "268435456@NumberOfRecHitsPerTrack_CKFTk@entries > 0", Start, End);
  A.createTrend("268435456@NumberOfTracks_CKFTk@mean", "NumberOfTracks_CKFTk_mean.gif", 0, "268435456@NumberOfTracks_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@NumberOfRecHitsPerTrack_CKFTk@entries", "NumberOfRecHitsPerTrack_CKFTk_entries.gif", 0, "268435456@NumberOfRecHitsPerTrack_CKFTk@entries > 0", Start, End);
  A.createTrend("268435456@Chi2overDoF_CosmicTk@mean", "Chi2overDoF_CosmicTk.gif", 0, "268435456@Chi2overDoF_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@Chi2overDoF_CKFTk@mean", "Chi2overDoF_CKFTk.gif", 0, "268435456@Chi2overDoF_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@Chi2overDoF_RSTk@mean", "Chi2overDoF_RSTk.gif", 0, "268435456@Chi2overDoF_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@NumberOfTracks_CosmicTk@mean", "NumberOfTracks_CosmicTk.gif", 0, "268435456@NumberOfTracks_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@NumberOfRecHitsPerTrack_CosmicTk@mean", "NumberOfRecHitsPerTrack_CosmicTk.gif", 0, "268435456@NumberOfRecHitsPerTrack_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPt_CosmicTk@mean", "TrackPt_CosmicTk.gif", 0, "268435456@TrackPt_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPz_CosmicTk@mean", "TrackPz_CosmicTk.gif", 0, "268435456@TrackPz_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPx_CosmicTk@mean", "TrackPx_CosmicTk.gif", 0, "268435456@TrackPx_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPy_CosmicTk@mean", "TrackPy_CosmicTk.gif", 0, "268435456@TrackPy_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPhi_CosmicTk@mean", "TrackPhi_CosmicTk.gif", 0, "268435456@TrackPhi_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackEta_CosmicTk@mean", "TrackEta_CosmicTk.gif", 0, "268435456@TrackEta_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@DistanceOfClosestApproach_CosmicTk@mean", "DistanceOfClosestApproach_CosmicTk.gif", 0, "268435456@DistanceOfClosestApproach_CosmicTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPt_CKFTk@mean", "TrackPt_CKFTk.gif", 0, "268435456@TrackPt_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPz_CKFTk@mean", "TrackPz_CKFTk.gif", 0, "268435456@TrackPz_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPx_CKFTk@mean", "TrackPx_CKFTk.gif", 0, "268435456@TrackPx_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPy_CKFTk@mean", "TrackPy_CKFTk.gif", 0, "268435456@TrackPy_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPhi_CKFTk@mean", "TrackPhi_CKFTk.gif", 0, "268435456@TrackPhi_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackEta_CKFTk@mean", "TrackEta_CKFTk.gif", 0, "268435456@TrackEta_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@DistanceOfClosestApproach_CKFTk@mean", "DistanceOfClosestApproach_CKFTk.gif", 0, "268435456@DistanceOfClosestApproach_CKFTk@mean > 0", Start, End);
  A.createTrend("268435456@NumberOfTracks_RSTk@mean", "NumberOfTracks_RSTk.gif", 0, "268435456@NumberOfTracks_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@NumberOfRecHitsPerTrack_RSTk@mean", "NumberOfRecHitsPerTrack_RSTk.gif", 0, "268435456@NumberOfRecHitsPerTrack_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPt_RSTk@mean", "TrackPt_RSTk.gif", 0, "268435456@TrackPt_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPz_RSTk@mean", "TrackPz_RSTk.gif", 0, "268435456@TrackPz_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPx_RSTk@mean", "TrackPx_RSTk.gif", 0, "268435456@TrackPx_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPy_RSTk@mean", "TrackPy_RSTk.gif", 0, "268435456@TrackPy_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackPhi_RSTk@mean", "TrackPhi_RSTk.gif", 0, "268435456@TrackPhi_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@TrackEta_RSTk@mean", "TrackEta_RSTk.gif", 0, "268435456@TrackEta_RSTk@mean > 0", Start, End);
  A.createTrend("268435456@DistanceOfClosestApproach_RSTk@mean", "DistanceOfClosestApproach_RSTk.gif", 0, "268435456@DistanceOfClosestApproach_RSTk@mean > 0", Start, End);




  A.closeFile();


  return;


}







int main (int argc, char* argv[])
{
  if (argc != 4 && argc != 5) {
    std::cerr << "Usage: " << argv[0] << " [TagName] [Password] [NRuns] " << std::endl;
    std::cerr << "Or:    " << argv[0] << " [TagName] [Password] [FirstRun] [LastRun] " << std::endl;
    return 1;
  }

  if (argc == 4) {
    std::cout << "Creating trends for NRuns = " << argv[3] << " for tag: " << argv[1] << std::endl;
    TrackingHDQMInspector( argv[1], argv[2], atoi(argv[3]) );
  } else if(argc == 5) {
    std::cout << "Creating trends for range:  " << argv[3] << " " << argv[4] << " for tag: " << argv[1] << std::endl;
    TrackingHDQMInspector( argv[1], argv[2], atoi(argv[3]), atoi(argv[4]) );
  }

  return 0;
}
