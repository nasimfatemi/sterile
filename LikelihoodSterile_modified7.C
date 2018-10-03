#include <TROOT.h>
#include <string>
#include <TFile.h>
#include <TApplication.h>
#include <TF1.h>
#include <TROOT.h>
#include <TFile.h>
#include "TF1.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TChain.h"
#include <vector>
#include <iostream>
#include <TCanvas.h>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <sstream>
#include <ostream>
#include "TH2F.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TFitResultPtr.h"
#include "RooAbsReal.h"
#include "RooGlobalFunc.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooExponential.h"
 #include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooNLLVar.h"
#include "RooHistPdf.h"
#include "RooAbsData.h"
#include "RooPrintable.h"
#include "RooDirItem.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooChebychev.h"
#include "TMinuit.h"
#include "RooRealSumPdf.h"
#include "RooNLLVar.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"
#include "RooTFnBinding.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
using namespace RooFit;
using namespace std;
#define PI 3.14159265
double Q=0.565;
double snu_mass=0.03;
TH1F *hmag;
TH1F *hmom;
TH1F *hfermi;
TH1F *hetot;
TH1F *pileup;
TH1F *hpileupsmear;
double sn_mass=0.1;

int chisqcount=0;
double Min=0.4;
double Max=0.6;
std::vector <TH1F*> histo(10);

TH1F *htest;
TH1F *hbetadecaystrileneutrino;
TH1F *hbetadecayredrawn;
TH1F *hmcarsmear;
TH1F *hfakedataressterile;

//TH1::SetDefaultSumw2();
TH1F *hfakedatabetadecay;
TH1F *hfakedatasterile;
TRandom3* Rand;


Double_t EnergyScale(Double_t E,Double_t *par){

  Double_t PE=0;
  PE=E;

  return PE;
}

Double_t GetSigma(Double_t E,Double_t *par){

  Double_t Sigma=0;
  //par[3]=par[0];
  Sigma=sqrt(EnergyScale(E,par)/7400);

  return Sigma;
}
Double_t GausFuncCoin(Double_t *x,Double_t *par){
  Double_t resultsmultiple = 0.;
  for (int i=1;i<pileup->GetNbinsX();i++){
 double tempx=hbetadecayredrawn->GetBinCenter(i);
    Double_t E=pileup->GetBinCenter(i);
    Double_t BinofCoin=pileup->FindBin(E);
    Double_t PE=EnergyScale(E,par);
    Double_t Sigma=sqrt(double(E)/7400);;

    resultsmultiple+=(pileup->GetBinContent(BinofCoin)*1.0*exp(-1.0/2.0*pow((x[0]-PE)/Sigma,2))/((sqrt(2*3.14))*Sigma));
  }
  return resultsmultiple;
}

Double_t GausFuncSingle(Double_t *x,Double_t *par){

  Double_t val=0.0;
  Double_t sumval=0.0;
  Double_t tempsum=0;
  
  Double_t resultssingle = 0.;
  Double_t resultsmultiple = 0.;
  Double_t resultbi214=0.;
  
  //cout<<x[0]<<endl;


for(int j=1; j<=hbetadecayredrawn->GetNbinsX(); j++)
    {

      double tempx=hbetadecayredrawn->GetBinCenter(j);
      double x1=tempx*Q/par[1];
      //double tempx=j/1000.;
      //double mom=sqrt((tempx*tempx)+(2*tempx*0.511));
      double mom=sqrt((x1*x1)+(2*x1*0.511));
      double etot=x1+0.511;
      double epsilon=(Q-x1)*(Q-x1);
      double betaenergycoef=mom*epsilon*etot;
      double mag=hmag->GetBinContent(hmag->FindBin(x1));
      double fermi=mag/betaenergycoef; //need to intorduce a way to give the magnitude

      mom=sqrt((tempx*tempx)+(2*tempx*0.511));
      etot=tempx+0.511;
      epsilon=(par[1]-tempx)*(par[1]-tempx);
      betaenergycoef=mom*epsilon*etot;

      double testsqrt= par[1]-tempx;
      if (testsqrt>snu_mass){

 
	val=mom*etot*fermi*(((testsqrt)*(1-(sin(par[2])*sin(par[2])))*sqrt((testsqrt)*(testsqrt)))+((testsqrt)*((sin(par[2])*sin(par[2]))*(sqrt(((testsqrt)*(testsqrt))-(snu_mass*snu_mass))))));
      
      }
      else val=mom*etot*fermi*(((testsqrt)*(1-(sin(par[2])*sin(par[2])))*sqrt((testsqrt)*(testsqrt))));
    Double_t E=hbetadecayredrawn->GetBinCenter(j);
   
    Double_t PE=EnergyScale(E,par);
    Double_t Sigma=sqrt(double(tempx)/7400);;
    resultssingle+=(val*1.0*exp(-1.0/2.0*pow((x[0]-PE)/Sigma,2))/((sqrt(2*3.14))*Sigma));
    // resultsmultiple+=(hcoinEdep->GetBinContent(BinofCoin)*1.0*exp(-1.0/2.0*pow((x[0]-PE)/Sigma,2))/((sqrt(2*3.14))*Sigma));
    //resultbi214+=(hBi214Edep->GetBinContent(i)*1.0*exp(-1.0/2.0*pow((x[0]-PE)/Sigma,2))/((sqrt(2*3.14))*Sigma));
    }
 
 return  resultssingle;
}

Double_t GausFunc(Double_t *x,Double_t *par){

 
  return (par[0]*GausFuncSingle(x,par)+par[3]*GausFuncCoin(x,par));

}



void PoissonData(std::vector<TH1F*>&   histo,TH1F* & hfakedata){
  Rand=new TRandom3(0);

  


  for (int ihist=0;ihist<histo.size();ihist++){

    int nbins=histo[ihist]->GetNbinsX();

    for (int i=0;i<nbins;i++){ 
      double     temp=Rand->Poisson(histo[ihist]->GetBinContent(i));
      hfakedata->Fill(histo[ihist]->GetBinCenter(i),temp);

    } 
  }

}

/*Double_t GausFuncCombined(Double_t *x,Double_t *par){

  Double_t val=0;
  val=par[0]*GausFuncSingle(x,par)+par[3]*GausFuncCoin(x,par);
  return val;

  }*/





void fcn(Int_t& npar, Double_t* gin, Double_t& f, Double_t *par, Int_t iflag){
 
  Double_t logl = 0;
  Double_t temp;
  Double_t energy[1];
  for (int i=0;i<hfakedatabetadecay->GetNbinsX();i++){
  energy[0]=hfakedatabetadecay->GetBinCenter(i);
  double scale=hfakedatabetadecay->Integral();
  //h1->Scale(1.0/scale);
  double y=hfakedatabetadecay->GetBinContent(i);
  Double_t pdf=GausFunc(energy,par);
  if ( y>0&& pdf>0 &&energy[0]<Max &&energy[0]>Min){
    
    double insidelog=pdf;
    
    // cout<<insidelog<<endl;
    // temp = y*(log(y)-log(insidelog));
    // temp=y*log(par[2]*insidelog)-par[2]*insidelog;
    temp=y*(log(insidelog))-insidelog;//-(((par[1]-0.565)*(par[1]-0.565))/2/0.001/0.001);
   
    // cout<<temp<<" "<<noconstrain<<endl;
    //temp=-log(insidelog);
    logl += temp; 
  } 
  }
  
  f=-logl;
  
}



void smearMC(TH1F *hmcsmear,TH1F *har,double smearvalue){

  
  double  minbin=hmcsmear->GetXaxis()->GetXmin();
  double  maxbin=hmcsmear->GetXaxis()->GetXmax(); 
  int Nbins=hmcsmear->GetNbinsX();

  double binwidth=hmcsmear->GetBinWidth(1);

  TH1F *tempa=new TH1F("TEMPa","temp",Nbins,minbin,maxbin);
	
  double dist=0;
  double ene=0;
  double sigma=0;
  double tempsum=0;
  for(int i=1; i<Nbins+1; i++)
    {
      ene = har->GetBinCenter(i);
	    
      sigma = sqrt(ene/smearvalue);
	   
      tempsum=0;
      for(int j=1; j<Nbins+1; j++)
	{
		  
	  dist = (double)(i-j)*binwidth;
		 
	  tempsum+=1./sqrt(2*3.1415)/sigma*exp(-1*dist*dist/sigma/sigma/2.);
	}
		
      for(int j=1; j<Nbins+1; j++)
	{
		 
	  dist = (double)(i-j)*binwidth;
		  
	  double wenergy=har->GetBinContent(i)*1./sqrt(2*3.1415)/sigma/tempsum*exp(-1*dist*dist/sigma/sigma/2.);
	  //cout<<wenergy<<endl;  
	  // if you dont want to smear
		 
	  tempa->Fill(har->GetBinCenter(j),wenergy);
		  
	}
    }
  for(int i=1; i<Nbins+1; i++)
    {
      hmcsmear->SetBinContent(i, tempa->GetBinContent(i));
    }
  //hbackpdf[nbkg]->Scale(1/hbackpdf[nbkg]->Integral());    

 
}

void  LikelihoodSterile_modified7(int ifile){
  //output histograms
  TH1F *hq=new TH1F("hq","Q par",50,0.56,0.57);
  TH1F *hqer=new TH1F("hqer","Q parerror",50,1.0E-4,1.0E-4);

  TH1F *hsintheta=new TH1F("hsintheta","sintheta",50,-8.0E-3,8.0E-3);
  TH1F *hsinthetaerror=new TH1F("hsinthetaerror","sintheta error",50,-1.0E-4,1.0E-3);



const int npar=4;
  //TH1::SetDefaultSumw2();
  cout<<"pased here"<<endl;
  TH1F *hbetadecay=new TH1F("hbetadecay","",1128,0.0,1.128);
  hbetadecayredrawn=new TH1F("hbetadecayredrawn","",1128,0.0,1.128);
  TH1F *hbetadecayneutrinomass=new TH1F("hbetadecayneutrinomass","",1128,0.0,1.128);
  hbetadecaystrileneutrino=new TH1F("hbetadecaystrileneutrino","sin2theta 0.2 and SN mass=10 keV",1128,0,1.128);
  TH1F *hfermi=new TH1F("hfermi","",564,0,0.2);
  TH1F *hbetaenergycoef=new TH1F("hbetaenergycoef","",50,0,30);
  //double Q=0.565;
  // double k=565.113;
  double k=1.0;
  double e_mass=0.511;
  double sin2theta=0.01;
  double cos2theta=1-sin2theta;
  double nu_mass=0.00000;
  

  std::vector<double> spec_E;
  std::vector<double> spec_mag;
  std::vector<double> spec_fermi;
  std::vector<double> spec_mom;
  std::vector<double> spec_etot;
  std::vector<double> spec_epsilon;
  std::vector<double> spec_epsilonneutrinomass;
  std::vector<double> spec_epsilonstrileneutrino;
  hmag=new TH1F("hmag","",565,0,0.565);
  // _lspec    = RAT::DB::Get()->GetLink("SPECTRUM","Ar39_beta");
  // spec_E=_lspec->GetDArray("spec_e");
  //spec_mag=_lspec->GetDArray("spec_mag");

ifstream spectrum_file;
  spectrum_file.open("/home/fatemigh/projects/rpp-jillings/fatemigh/sterile/Likelihood/39Ar_norm.dat");
  if (spectrum_file.fail()){
    cout << "Couldn't open file" << endl;
    exit(1);
  }

 int lines = 0;
  while ( !spectrum_file.eof() ){
    double E_temp,mag_temp;
    spectrum_file >> E_temp >> mag_temp;
    if ( !spectrum_file.eof() ) {
      lines++;
      //cout << E_temp << " " << mag_temp << endl;
      spec_E.push_back(E_temp);
      spec_mag.push_back(mag_temp);
    }
  }

  cout << lines << " lines read from file" << endl;
  spectrum_file.close();

  double rate=3.0*1.01*3300*24*3600*365;
  double minvalue=(snu_mass+0.001)*1000;

  cout<<minvalue<<endl;
  cout<<"spec_E"<<"  "<<spec_E.size()<<endl;
  for(unsigned int istep = 0; istep<spec_E.size()-minvalue; ++istep){
    spec_mag[istep]=spec_mag[istep]; 
    hbetadecay->SetBinContent(istep,spec_mag[istep]);
    double mom=sqrt((spec_E[istep]*spec_E[istep])+(2*spec_E[istep]*e_mass));
    double epsilon=(Q-spec_E[istep])*(Q-spec_E[istep]);
    double neutrinomassterm=sqrt((Q-spec_E[istep])*(Q-spec_E[istep])-(nu_mass*nu_mass));
    double epsilonneutrinomass=(Q-spec_E[istep])*neutrinomassterm;
    double epsilonstrileneutrino=(Q-spec_E[istep])*((cos2theta*neutrinomassterm)+sin2theta*(sqrt((Q-spec_E[istep])*(Q-spec_E[istep])-snu_mass*snu_mass)));
    double etot=spec_E[istep]+(0.511);
    hmag->SetBinContent(istep,spec_mag[istep]);
    
    double betaenergycoef=k*mom*epsilon*etot;
    double fermi=spec_mag[istep]/betaenergycoef;
 
    spec_fermi.push_back(fermi);
    spec_mom.push_back(mom);
    spec_etot.push_back(etot);
    spec_epsilon.push_back(epsilon);
    spec_epsilonneutrinomass.push_back(epsilonneutrinomass);
    spec_epsilonstrileneutrino.push_back(epsilonstrileneutrino);
     
    hfermi->Fill(fermi);
    hbetaenergycoef->Fill(betaenergycoef);
    
    //cout<<istep<<"  "<<epsilonneutrinomass<<"   "<<epsilonstrileneutrino<<endl;
    
  }
  
  for(unsigned int istep = spec_E.size()-minvalue; istep<spec_E.size(); ++istep){

    spec_mag[istep]=spec_mag[istep]; 
    //cout<<spec_mag[istep]<<"  "<<spec_E[istep]<<endl;
    hbetadecay->SetBinContent(istep,spec_mag[istep]);
    double mom=sqrt((spec_E[istep]*spec_E[istep])+(2*spec_E[istep]*e_mass));
    double epsilon=(Q-spec_E[istep])*(Q-spec_E[istep]);
    double neutrinomassterm=sqrt((Q-spec_E[istep])*(Q-spec_E[istep])-(nu_mass*nu_mass));
    double epsilonneutrinomass=(Q-spec_E[istep])*neutrinomassterm;
    double epsilonstrileneutrino=(Q-spec_E[istep])*(cos2theta*neutrinomassterm); //remove the heavy part from the equation
    double etot=spec_E[istep]+(0.511);
   
    hmag->SetBinContent(istep,spec_mag[istep]);
   
    double betaenergycoef=k*mom*epsilon*etot;
    double fermi=spec_mag[istep]/betaenergycoef;
  
    spec_fermi.push_back(fermi);
    spec_mom.push_back(mom);
    spec_etot.push_back(etot);
    spec_epsilon.push_back(epsilon);
    spec_epsilonneutrinomass.push_back(epsilonneutrinomass);
    spec_epsilonstrileneutrino.push_back(epsilonstrileneutrino);
    
    hfermi->Fill(fermi);
    hbetaenergycoef->Fill(betaenergycoef);    
  }

  for (int i=0;i<spec_fermi.size()-1;i++){
    //cout<<spec_fermi[i]<<"  "<<spec_mom[i]<<" "<<spec_etot[i]<<"  "<<spec_epsilonstrileneutrino[i]<<"  "<<spec_epsilonneutrinomass[i]<<endl;
    double dn=k*(spec_fermi[i]*spec_mom[i]*spec_etot[i]*spec_epsilon[i]);
    double dnneutrinomass=k*(spec_fermi[i]*spec_mom[i]*spec_etot[i]*spec_epsilonneutrinomass[i]);
    double dnstrileneutrino=k*(spec_fermi[i]*spec_mom[i]*spec_etot[i]*spec_epsilonstrileneutrino[i]);

    //cout<<dn<<"    "<<dnneutrinomass<<"  "<<dnstrileneutrino<<endl;
    //cout<<dnneutrinomass<<"  "<<dnstrileneutrino<<endl;
    hbetadecayredrawn->SetBinContent(i,dn);
    hbetadecayneutrinomass->SetBinContent(i,dnneutrinomass);
    hbetadecaystrileneutrino->SetBinContent(i,dnstrileneutrino);



  }
  hbetadecay->Scale(rate/hbetadecay->Integral());
  hbetadecayredrawn->Scale(rate/hbetadecayredrawn->Integral());
  hbetadecaystrileneutrino->Scale(rate/hbetadecaystrileneutrino->Integral());

  //-------------------------------


  //error bar calculation
  //TGraph *gr=new TGraph();
  for (int i=0;i<hbetadecayredrawn->GetNbinsX();i++){
    hbetadecayredrawn->SetBinError(i,sqrt(hbetadecayredrawn->GetBinContent(i)));
    hbetadecaystrileneutrino->SetBinError(i,sqrt(hbetadecaystrileneutrino->GetBinContent(i)));

  }

  /*

  //first histogram
  sprintf(temp,"histo%d",0);
  histo[0]=(TH1F*)hbetadecayredrawn->Clone(temp);

  //second histogram and all

  TH1F *hdatasmear=(TH1F*)histo[0]->Clone("hdatasmear");
  hdatasmear->Reset();

  smearMC(hdatasmear,histo[0],7800);
  histo[0]=hdatasmear;
  PoissonData(histo);
  */
  int  rebin=1;
  hbetadecay->Rebin(rebin);
  hbetadecayredrawn->Rebin(rebin);
  hbetadecaystrileneutrino->Rebin(rebin);
  hmag->Rebin(rebin);
  char temp[30];
  for (int i=0;i<10;i++){
    sprintf(temp,"histo%d",i);
    histo[i]=(TH1F*)hbetadecayredrawn->Clone(temp);
    histo[i]->Reset();
  }

  //first histogram
  sprintf(temp,"histo%d",0);
  histo[0]=(TH1F*)hbetadecayredrawn->Clone(temp);
  TH1F *hdatasmear=(TH1F*)histo[0]->Clone("hdatasmear"); // hmm,
 TH1F *hdatasmear2=(TH1F*)histo[0]->Clone("hdatasmear2"); 
  //pile up histogram

 TFile *finpileup=new TFile("/home/fatemigh/projects/rpp-jillings/fatemigh/sterile/Likelihood/pileup_theoritical.root","read");
 pileup=(TH1F*)finpileup->Get("final");
 pileup->Scale(1.0/pileup->Integral());
 pileup->Scale(0.01*hbetadecay->Integral()/pileup->Integral());
 hpileupsmear=(TH1F*)pileup->Clone("hpileupsmear");
 hpileupsmear->Reset();
 //histo[1]=(TH1F*)pileup->Clone("pileup1");

  smearMC(hdatasmear,histo[0],7400);
  //smearMC(pileup,histo[1],7400);
  histo[0]=hdatasmear;
  //histo[1]=pileup;
  double xmin=hbetadecayredrawn->GetXaxis()->GetXmin();
  double xmax=hbetadecayredrawn->GetXaxis()->GetXmax(); 

  hfakedatabetadecay=new TH1F("hfakedatabetadecay","",hbetadecayredrawn->GetNbinsX(),xmin,xmax);
  PoissonData(histo,hfakedatabetadecay);

 

  ///// sterile neutrino
  histo[0]=(TH1F*)hbetadecayredrawn->Clone(temp);

  hdatasmear=(TH1F*)histo[0]->Clone("hdatasmear");
  hdatasmear->Reset();
  
  hfakedatasterile=new TH1F("hfakedatasterile","",hbetadecayredrawn->GetNbinsX(),xmin,xmax);
  PoissonData(histo,hfakedatasterile);



  TH1F *hfakedatares=(TH1F*)hfakedatabetadecay->Clone("hfakedatares");
  hfakedatares->Reset();
  hfakedataressterile=(TH1F*)hfakedatasterile->Clone("hfakedataressterile");
  hfakedataressterile->Reset();
  TH1F *hfittedvalue=(TH1F*)hfakedatabetadecay->Clone("hfittedvalue");
  hfittedvalue->Reset();


  //hbetadecayredrawn->Integral(rate/hbetadecayredrawn->Integral());

  cout<<"Before entering fitting" <<endl;
  
 
  TF1 *func2=new TF1("func2",GausFunc,0.53,0.59,npar);

 TH1F *hpileuppar=new TH1F("hpileuppar","pile up par",200,0.5,1.5);
  
  TH1F *hnorm=new TH1F("hnorm","Norm",500,5.5e+05,5.7e+05);
  TH1F *hnormerror=new TH1F("hnormerror","Norm error",500,1.0E4,5.0E4);
  TH1F *hqerror=new TH1F("hqerror","Q Error",100,1.0E-3,0.06E-3);
  TH1F *hangle=new TH1F("hangle","Angle",500,-0.1E-3,1.0E-3);
  TH1F *hangleerror=new TH1F("hangleerror","Angle error",500,-0.1E-1,1.0E-1);
  TH2F *hangleerrorqpe2=new TH2F("hangleerrorqpe2","2D",1200,-1.0E-2,1E-2,400,0,1.13);
 

  TFile *fout=new TFile(Form("/home/fatemigh/projects/rpp-jillings/fatemigh/sterile/Likelihood/fitresultnoconstrain/fitresultswithpileupthetanoconstrain%i.root",ifile),"recreate");
 
  for (int ifit=0; ifit<10;ifit++){
 
 
  TF1 *func=new TF1("func",GausFunc,Min,Max,npar);
    for (int i=0;i<10;i++){
      sprintf(temp,"histo%d",i);
      histo[i]=(TH1F*)hbetadecayredrawn->Clone(temp);
      histo[i]->Reset();
    }

//first histogram
    sprintf(temp,"histo%d",0);
    // histo[0]=(TH1F*)hbetadecaystrileneutrino->Clone(temp);
    histo[0]=(TH1F*)hbetadecayredrawn->Clone(temp);
    hdatasmear=(TH1F*)histo[0]->Clone("hdatasmear");
    hdatasmear->Reset();
    hfakedatabetadecay->Reset();
    smearMC(hdatasmear,histo[0],7400);
    histo[0]=hdatasmear;
    //second histogram
    pileup=(TH1F*)finpileup->Get("final");
    pileup->Scale(1.0/pileup->Integral());
    pileup->Scale(0.01*hbetadecay->Integral()/pileup->Integral());
    hpileupsmear=(TH1F*)pileup->Clone("hpileupsmear");
    hpileupsmear->Reset();
    histo[1]=(TH1F*)pileup->Clone("pileup1");
    smearMC(hpileupsmear,histo[1],7400);
    histo[1]=hpileupsmear;
    hfakedatabetadecay=new TH1F("hfakedatabetadecay","",hbetadecayredrawn->GetNbinsX(),xmin,xmax);
    PoissonData(histo,hfakedatabetadecay);

    func=new TF1("func",GausFunc,Min,Max,npar);
    func->SetParameter(0,5.58e+05);
    func->SetParameter(1,0.565);
    func->SetParameter(2,0.0);
    func->SetParameter(3,9.99871E-4);
    //  func->FixParameter(3,0);
    // TFitResult r= (TFitResult)hfakedatabetadecay->Fit(func,"LVR");
     int fitattempt=0;
     /* while(r.Status()==4 &&fitattempt<4){
     
    r= (TFitResult)hfakedatabetadecay->Fit(func,"LVR");
    fitattempt++;
    
    }*/
     
 for (int i=0;i<hfakedatabetadecay->GetNbinsX();i++){
      hfakedatabetadecay->SetBinError(i,sqrt(histo[0]->GetBinContent(i)));
    }
    //  if (r.Status()!=0) continue;
    //if (func->GetParameter(1)>0.59) continue;
    //if (func->GetParameter(2)>1.0) continue;
  
 
 // cout<<"passed here----------------------"<<func->GetParameter(3)<<endl;
  cout<<"------------------------------------------"<<endl;
  cout<<"Enering ROOFIt-------"<<endl;
  
 
TMinuit *gMinuit = new TMinuit(npar);
   gMinuit->SetFCN(fcn);
   Double_t arglist[10];
   Int_t ierflg=0;
   
  
   fitattempt=0;
 double par[npar];   
 double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];
  
  par[0]=5.5923E5;
  par[1]=0.565;
  par[2]=0.0;
  par[3]=9.94E-4;
  
   parName[0]="norm";
   parName[1]="Q";
   parName[2]="theta";
   parName[3]="pileup";

   minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.
   maxVal[0] = 0; 
   minVal[1] = 0;   // if min and max values = 0, parameter is unbounded.
   maxVal[1] = 0;
   minVal[2] =0;   // if min and max values = 0, parameter is unbounded.
   maxVal[2] = 0;
    minVal[3] =0;   // if min and max values = 0, parameter is unbounded.
   maxVal[3] = 0;

   stepSize[0] =0.1;
   stepSize[1]=0.1;
   stepSize[2]=0.1;
   stepSize[3]=0.1;
   int fitstatus=4;
   TF1 *funcout=new TF1("funcout",GausFunc,Min,Max,npar);
   Double_t currentValue[npar];
   Double_t currentError[npar];
   while ((fitstatus==4 || currentError[2]>=0.5)&& fitattempt<2) {
   
for (int i=0; i<npar; i++){
  gMinuit->DefineParameter(i, parName[i].c_str(), 
      par[i], stepSize[i], minVal[i], maxVal[i]);
  }
 gMinuit->Migrad();
 fitstatus=gMinuit->Command("MIGRAD");
 fitattempt++;
 cout<<"status of the fit-----000------"<<fitstatus<<endl;
   // stepSize[2]=0.001;
 
 for (int ipar=0;ipar<npar;ipar++){
   Int_t test=gMinuit->GetParameter(ipar,currentValue[ipar],currentError[ipar]);
   //cout<<test<<"  "<<currentValue[ipar]<<endl;
   //funcout->SetParameter(ipar,currentValue[ipar]);
   //funcout->SetParError(ipar,currentError[ipar]);  
 }

 par[0]=currentValue[0];
 par[3]=currentValue[3];
 par[1]=currentValue[1];
 par[2]=currentValue[2];
 
  
   }

  if (fitstatus==0){
    for (int ipar=0;ipar<npar;ipar++){
    Int_t test=gMinuit->GetParameter(ipar,currentValue[ipar],currentError[ipar]);
    cout<<test<<"  "<<currentValue[ipar]<<endl;
    funcout->SetParameter(ipar,currentValue[ipar]);
    funcout->SetParError(ipar,currentError[ipar]);  
    }
    
    hfakedatabetadecay->Draw();
    funcout->SetLineColor(kBlue);
    funcout->Draw("same");
    funcout->Write();
    hfakedatabetadecay->Write();
  hq->Fill(currentValue[1]);
  hqerror->Fill(currentError[1]);

  hsintheta->Fill(currentValue[2]);
  hsinthetaerror->Fill(currentError[2]);
  }
  }

  fout->cd();
  hq->Write();
  hqerror->Write();
  hsintheta->Write();
  hsinthetaerror->Write();
  fout->Close();

}
