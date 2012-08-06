void testFit(int runNr){
  char fileName[100];

  sprintf(fileName,"/data/disk1/rohe/source_test_2008/rawdata/bt05r%06d/spectra_%06d.root\0",runNr, runNr);

  TFile *f;
  f = new TFile(fileName,"read");
  if (f->IsZombie()){
    cout<< fileName << " does not exist\n";
    delete(f);
    return;
  } 
  
  f->ls(); 
  
  // get the histograms
  
  char hName[50], fName[50];
  
  TH1D *hq[8];
  TF1  *fq[8];
  
  
  for (int i=0; i<8; i++){
    //create the names out of the runnumber
    sprintf(hName,"cQ%d_%d",i,runNr);
    cout<< "Histogram name "<< hName << endl;
    sprintf(fName,"Fitfcn_cQ%d_%d",i,runNr);
    cout<< "Function name "<< fName << endl;
        
    hq[i] = (TH1D*)f->Get(hName);
    fq[i] = (TF1*)f->Get(fName);
  }
  
  // get Fit parameters  
  TH2D *hFit, *hErr;
  double mpw[8], sMpw[8], gauss[8], sGauss[8];
  sprintf(hName,"hFitPar_%d",runNr);
  hFit = (TH2D*) f->Get(hName);
  sprintf(hName,"hFitErr_%d",runNr);
  hErr = (TH2D*) f->Get(hName);
  
  for (int i=0; i<8; i++){
    mpw[i] = hFit->GetBinContent(i+1,2);
    sMpw[i] = hErr->GetBinContent(i+1,2);
    gauss[i] = hFit->GetBinContent(i+1,4);
    sGauss[i] = hErr->GetBinContent(i+1,4);
  
  }  
  
  // Plot the stuff
  
 // Window 1 for the dirtributions
 TCanvas c1("c1","distribution",1);
 c1.cd();
 hq[0]->SetLineWidth(2);
 hq[0]->Draw();
 
 for (int i = 1; i<7; i++){
    hq[i]->SetLineColor(i);
    hq[i]->Draw("same");
 }
 
 char plotName[100];
 sprintf(plotName,"plots/cluster_run_%d.pdf",runNr);
 c1.Print(plotName);  
 
 // Window 2 for the Fits
 TCanvas c2("c2","Fits",0,0, 630,891);// aspect ratio as A4
 c2.Divide(1,4);
 
 TPaveText *pt[4];
 
 for (int i=0; i<4; i++){
    
   c2.cd(i+1);
   
   int j = (i==3)?7:i;
   
   hq[j]->SetLineWidth(1);
   hq[j]->SetLineColor(1);
   hq[j]->SetAxisRange(0.,1000.,"X");
   hq[j]->Draw();
   fq[j]->Draw("samel");
   
   char text1[100], text2[100];
   sprintf(text1,"MPW = %3.2f +/- %1.2f \n",mpw[j],sMpw[j]);
   sprintf(text2,"Gauss = %3.2f +/- %1.2f \0",gauss[j], sGauss[j]);
      
   double yMin, yMax;
   yMin = (hq[j]->GetMaximum() * 0.1);
   yMax = yMin * 4.;
   
   pt[i] = new TPaveText(600,yMin,1000,yMax);
   pt[i]->AddText(text1);
   pt[i]->AddText(text2);   
   pt[i]->Draw();
      
 }
 

 sprintf(plotName,"plots/fits_run_%d.pdf",runNr); 
 c2.Print(plotName);
 
 f->Close();
 
 return;
}