//fast eventdisplay
// after Chiohia
//
int display_event(int run_nr)
{

    gROOT->Reset();
    gStyle->SetPalette(1,0);
    gStyle->SetOptStat(0);

    // open file
    Char_t filename[80];
    sprintf(filename,"/home/l_tester/log/bt05r%06d/run_%06d.root",run_nr, run_nr);
    cout << filename << endl ;
    TFile *f = new TFile(filename);
    f->ls();

    TTree *t = (TTree*)f->Get("t1");

    Int_t PixN, PixCol[1000], PixRow[1000], PixAna[1000];
    Float_t PixVcal[1000];
  
    
    Int_t CluN, CluSize[1000];
    Float_t CluQ[1000], CluCol[1000], CluRow[1000];
    
    Int_t zeroValueAna = -400;
    Float_t zeroValueVcal = -40.;
    
    t->SetBranchAddress("PixN",&PixN);
    t->SetBranchAddress("PixCol",&PixCol[0]);
    t->SetBranchAddress("PixRow",&PixRow[0]); 
    t->SetBranchAddress("PixAna",&PixAna[0]);
    t->SetBranchAddress("PixVcal",&PixVcal[0]);            
 
    t->SetBranchAddress("CluN",&CluN);    
    t->SetBranchAddress("CluSize",&CluSize[0]);        
    t->SetBranchAddress("CluQ",&CluQ[0]);    
    t->SetBranchAddress("CluCol",&CluCol[0]);        
    t->SetBranchAddress("CluRow",&CluRow[0]);        
        
    TH2S *hPixelAna  = new TH2S("hPixelAna","ADC values",52,0,52,80,0,80);
    TH2S *hPixelCal  = new TH2S("hPixelCal","Vcal values",52,0,52,80,0,80);
    TH2F *hPixelClu  = new TH2F("hPixelClu","Clusters",52,0,52,80,0,80);
    
    // init the histograms
    for ( Int_t i=0; i<52; i++) {
      for (Int_t j=0; j<80; j++) {
        hPixelAna->SetBinContent(i+1,j+1,zeroValueAna); 
	hPixelCal->SetBinContent(i+1,j+1,zeroValueVcal);
	hPixelClu->SetBinContent(i+1,j+1,0);
      }
    }        
    int k,n;
    char weiter;

    n = t->GetEntries();
    TCanvas *lego = new TCanvas("lego2","Pixel Map",100,100,1000,600);
    lego->Divide(3,1);
    
    


    for (k=0; k<n; k++){

      t1->GetEntry(k);
      cout << "Event No " << k << " of " << n <<endl;
      
      if (PixN > 0){
      
        cout << "PixN "<< PixN << endl;
        // fill pixel hits
        for (Int_t pixCounter=0; pixCounter < PixN; pixCounter++){
      
	  cout << "Pix Nr. " << pixCounter ;
          cout << " PixCol "<<PixCol[pixCounter]<<" PixRow "<<PixRow[pixCounter];
	  cout << " PixAna "<<PixAna[pixCounter]<<" PixVcal "<<PixVcal[pixCounter]<<endl;
		
          hPixelAna->SetBinContent(PixCol[pixCounter]+1,PixRow[pixCounter]+1,PixAna[pixCounter]);
          hPixelCal->SetBinContent(PixCol[pixCounter]+1,PixRow[pixCounter]+1,PixVcal[pixCounter]);	  
        }      

        //fill the his of CLusters
	
        for (Int_t cluCounter=0; cluCounter < CluN; cluCounter++){
      
          cout << "CluSize "<<CluSize[cluCounter]<<" CluQ "<<CluQ[cluCounter]<<
	          " CluCol "<<CluCol[cluCounter]<<" CluRow "<<CluRow[cluCounter]<<endl;
		
          hPixelClu->SetBinContent(CluCol[cluCounter]+1,CluRow[cluCounter]+1,CluQ[cluCounter]);
        }
      
      
        // Plot histogram

        lego->cd(1);
        hPixelAna->Draw("colz");
        lego->cd(2);
        hPixelCal->Draw("colz");
        lego->cd(3);
        hPixelClu->Draw("colz");
        lego->Update();     
      }else{
        cout << "No hit in this event";
	continue;
      }
      
      //if (k==529) lego->Print("event_529.eps");
      //if (k==680) lego->Print("event_680.eps");

      cout << "Type in <cr> for next event, h <cr>: +100, t<cr> + 1000, b<cr> 1 back \n";
      weiter = getchar();
      if (weiter == 'q') {
        f->Close();
        return 0;
      };
      if (weiter == 'h') {
        k+=100;
      };
      if (weiter == 't') {
        k+=1000;
      };
      if (weiter == 'b') {
        k-=2;
      };
      
      // now delete the hits from the histograms for the next event
      
      if (PixN > 0){
      
        // fill pixel hits
        for (Int_t pixCounter=0; pixCounter < PixN; pixCounter++){
		
          hPixelAna->SetBinContent(PixCol[pixCounter]+1,PixRow[pixCounter]+1,zeroValueAna);
          hPixelCal->SetBinContent(PixCol[pixCounter]+1,PixRow[pixCounter]+1,zeroValueVcal);	  
        }      

        //fill the his of CLusters
	
        for (Int_t cluCounter=0; cluCounter < CluN; cluCounter++){
 
          hPixelClu->SetBinContent(CluCol[cluCounter]+1,CluRow[cluCounter]+1,0);
        }
      }      
      
      


    }// event loop

    f->Close();
    return 1;
}
