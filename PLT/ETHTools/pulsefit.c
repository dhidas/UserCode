{

double params[6][4];
int run[6];

run[0] = 4930;
run[1] = 4936;
run[2] = 4941;
run[3] = 4913;
run[4] = 4946;
run[5] = 4951;

//Get parameters
for (int i = 0; i < 6; i++)
  {
    TFile f(Form("/home/l_tester/log/bt05r00%i/fitParameterDistributionRun_%i.root", run[i], run[i]));
    //f.Open();
    if (!f.IsOpen()){cout << "Error: could not open file for run " << run[i] << endl; break;}

    for (int j = 0; j < 4; j++)
      {
	TH1D *h = (TH1D*)f->Get(Form("hParOld%i",j));
	params[i][j] = h->GetMean(1);
	//params[j] = h->GetMaximumBin();
      }
    f.Close();
    
  }


for (int k = 0; k < 6; k++)
  {
    TF1* c = new TF1("c","[3]+[2]*tanh([0]*x-[1])",0,1600);
    c->SetParameters( params[k][0], params[k][1], params[k][2], params[k][3]);
    if (k > 2){c->SetLineColor(2);c->SetLineStyle(2);}
    if (k==0){c->Draw("l");}
    else{c->Draw("lsame");}
  }


}
