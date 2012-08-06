void plotClusterSize(){

  TFile *f[4];
  int runNr[4] = {3052,3602,4272,4457};
  double normFactor = 100000.;
  
  char entryText[4][50]={"Not irradiated, 150V","6.1E14, 600V","1.1E15, 600V","2.8E15, 800V" };
  
  TLegend *leg;
  leg = new TLegend(.5,.65,.85,.85);
  leg->SetFillColor(0);
  
  TH1D *his[4];
 
  TCanvas c;
  c.SetFillColor(0);
  
  TH2D axis("axis","",10,0,10,10,0,40000);
  axis.SetStats(0);
  axis.GetXaxis()->SetTitle("Cluster Size");
  //axis.GetYaxis()->SetTitle("Number of Entries"); 
  axis.Draw();
  
  for(int i=0; i<4; i++)
  {
  
    char filename[200];
    sprintf(filename,"../rawdata/bt05r%06d/spectra_%06d.root",runNr[i],runNr[i]);
    f[i] = new TFile(filename,"read");
    
    char hisname[50];
    sprintf(hisname,"cSize_%d",runNr[i]);
    his[i] = (TH1D*) f[i]->Get(hisname);
    double scaleFactor = normFactor/ his[i]->Integral();
    his[i]->Scale(scaleFactor);
    his[i]->SetLineWidth(3);
    if (i%2==0)  his[i]->SetLineColor(1); else his[i]->SetLineColor(14);
    if (i>1) his[i]->SetLineStyle(2);
    his[i]->SetStats(0);
    his[i]->Draw("same");
    
    leg->AddEntry(his[i],entryText[i],"l");
  };
  
  leg->Draw();
  
  c.Print("clusterSize.eps");
  c.Print("clusterSize.pdf");
  c.Print("clusterSize.jpg");
}
