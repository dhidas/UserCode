plotSpectrum(int rN){
  
  TCanvas c;
  c.SetFillColor(0);
  c.SetLogy(1);

  char filename[200];
  sprintf(filename,"../rawdata/bt05r%06d/spectra_%06d.root",rN,rN);
  TFile *f;
  f = new TFile(filename,"read");
  
  char hisname[50];
  sprintf(hisname,"cQ0_%d",rN);
  
  TH2D axis("axis","",10,0,2000,10,50,3000);
  axis.SetStats(0);
  axis.GetXaxis()->SetTitle("Signal Height [a.u.]");  
  axis.Draw(); 
  
  TLegend leg(.65,.5,.85,.85);
  leg.SetFillColor(0);
  
  char legText[100];
  
  TH1D *his[8];
  
  for (int i=0; i<7;i++){
    sprintf(hisname,"cQ%d_%d",i,rN);
    his[i] = (TH1D*) f->Get(hisname);
    
    // width
    if (i==0) his[i]->SetLineWidth(4);
    else if (i>3)  his[i]->SetLineWidth(2);
    else his[i]->SetLineWidth(3);
    // type

    if (i==2) his[i]->SetLineStyle(2);    
    else his[i]->SetLineStyle(1);
    
    if (i==3) his[i]->SetLineColor(15); 
    else his[i]->SetLineColor(1);
    
    int n=2;
    his[i]->Smooth(n);
    
    if (i<5) his[i]->Draw("same");
      
    if (i==0) sprintf(legText,"All Clusters");
    else if (i==1) sprintf(legText,"%d Pixel",i);
    else if (i<5) sprintf(legText,"%d Pixels",i);
       
    if (i<5) leg.AddEntry(his[i],legText,"l");
    
  }
  
  TH1D *sum;
  sum = new TH1D("sum","",501, -100, 2400);
  sum->Add(his[5],his[6],1,1);
  sum->SetLineColor(15);
  sum->SetLineStyle(1);
  sum->Draw("same");
  leg.AddEntry(sum,"> 4 Pixels","l");
    
  leg.Draw();
  c.Print("spectrum.eps");
  c.Print("spectrum.pdf");  
  c.Print("spectrum.jpg");  





}