TH1* h1;
Double_t Npoints;
Double_t  callikelihood(Double_t x,Double_t *par){
  double val=0;
   val=+1.0/(par[1]*sqrt(2*3.14))*exp((-1./2.)*pow(((x-par[0])/par[1]),2)); 
  // val=+(par[2]/(par[1]*sqrt(2*3.14)))*exp((-1./2.)*pow(((x-par[0])/par[1]),2)); 
  return val;
}

Double_t  functest(Double_t *x,Double_t *par){
  double val=0;
  // val=+1.0/(par[1]*sqrt(2*3.14))*exp((-1./2.)*pow(((x-par[0])/par[1]),2)); 
   val=+par[2]*exp((-1./2.)*pow(((x[0]-par[0])/par[1]),2)); 
  return val;
}

void fcn(Int_t& npar, Double_t* gin, Double_t& f, Double_t *par, Int_t iflag){
 
  Double_t logl = 0;
  Double_t temp;
  Double_t energy;
  for (int i=0;i<h1->GetNbinsX();i++){
  energy=h1->GetBinCenter(i);
  double pdf=callikelihood(energy,par);
  if (pdf>0){
    
    double insidelog=pdf;
    
    cout<<insidelog<<endl;
    temp = log(insidelog);
    
    logl += temp; 
  } 
  }
  
  f=-2*logl;
  
}
void likelihood_norm(){
  TF1 *f1=new TF1("f1","gaus(0)",0,10);
  f1->FixParameter(0,100.0);
  f1->FixParameter(1,5);
  f1->FixParameter(2,2);

  f1->Draw();
  h1=f1->GetHistogram();
  Npoints=h1->GetNbinsX();
  cout<<Npoints<<"-------"<<endl;
  h1->Sumw2();
  h1->Draw();
  //h1->Rebin(2);
   TMinuit *gMinuit = new TMinuit(2);
   gMinuit->SetFCN(fcn);
   Double_t arglist[10];
   Int_t ierflg=0;
   // arglist[0]=0.5; //for maximum likelihood for some reason it is 0.5 for chisq =1.0
   // gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
   //cout<<ierflg<<"------------"<<endl;
  
const int npar = 3;
 double par[npar];   
 double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];
  
  par[0]=4;
  par[1]=2;
  par[2]=100;
  
   parName[0]="mean";
   parName[1]="choose";
   parName[2]="norm";
  

   minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.
   maxVal[0] = 0; 
   minVal[1] = 0;   // if min and max values = 0, parameter is unbounded.
   maxVal[1] = 0;
    minVal[2] =0;   // if min and max values = 0, parameter is unbounded.
  maxVal[2] = 0;
   

   stepSize[0] =0.001;
   stepSize[1]=0.001;
   stepSize[2]=0.001;

for (int i=0; i<npar; i++){
  gMinuit.DefineParameter(i, parName[i].c_str(), 
      par[i], stepSize[i], minVal[i], maxVal[i]);
  }
 gMinuit.Migrad();
   // stepSize[2]=0.001;
// gMinuit->mnparm(0,"mu",par[0],stepSize[0],0.,0.,ierflg);
// gMinuit->mnparm(1,"sigma",par[1],stepSize[1],0.,0.,ierflg);
   //gMinuit->mnparm(2,"norm",par[2],stepSize[2],0.,0.,ierflg);

// arglist[0] = 1 ;  
  
// gMinuit->mnexcm("SCAN",arglist,0,ierflg);
//  cout<<ierflg<<"------------"<<endl;

   // Now ready for minimization step
/* arglist[0] = 500;
   arglist[1] = 0.0;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   cout<<ierflg<<"------------"<<endl;

   gMinuit->mnexcm("HESSE", arglist ,1,ierflg);
   cout<<ierflg<<"------------"<<endl;
   // Print results
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   gMinuit->mnprin(3,amin);
   /*  
   Double_t arglist[10];
   Int_t ierflg = 0;
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

   arglist[0]=500;
   arglist[1] = 1.;
   
   for (int i=0; i<3; i++){
     gMinuit->DefineParameter(i, parName[i].c_str(), par[i], stepSize[i], minVal[i], maxVal[i]);   
   }


gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

   */
   //for (int i=0; i<3; i++){
   // gMinuit->DefineParameter(i, parName[i].c_str(), par[i], stepSize[i], minVal[i], maxVal[i]);   
   // }
   //gMinuit->Migrad();       // Minuit's best minimization algorithm
// gMinuit->Draw();
/*TF1* funcfit = new TF1("funcfit",callikelihood,0.0, 10, 3);
funcfit->SetParameter(0,5);
funcfit->SetParameter(1,2);
funcfit->SetParameter(2,100);
 funcfit->Draw();
 h1->Draw();
 funcfit->Draw("same");
 //h1->Draw("same");
 */

 TF1 *f2=new TF1("f2",functest,0,10,3);
  f2->FixParameter(2,100.0);
  f2->FixParameter(0,5);
  f2->FixParameter(1,2);
  f2->Draw("same");

}
