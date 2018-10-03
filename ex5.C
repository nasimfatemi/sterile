// ex5.C
//
// unbinned likelihood fit to data in NTuple
#include <cmath>
#include <iostream>

#include <TH1D.h>
#include <TF1.h>
#include <TMinuit.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TColor.h>

// tuple to which to fit
TNtuple* tuple;
// start and end event for fit (i.e. fit only to data in these entries)
unsigned long iStartEv, iEndEv;

// FCN to fit exp(-t/tau)/tau to the data in the tuple
void fcnTau(Int_t& npar, Double_t* gin, Double_t& f, 
	Double_t* par, Int_t iflag)
{
    if (1 == iflag) {
	// starting a new fit, initialize anything we need to have
	// (not needed here)
    }

    // calculate value for f
    double LogL = 0., dLogLdtau = 0.;
    for (unsigned long iEv = iStartEv; iEv < iEndEv; ++iEv) {
	tuple->GetEntry(iEv);
	const double t = tuple->GetArgs()[0];
	LogL += t/par[0] + std::log(par[0]);
	if (2 == iflag)
	    dLogLdtau += -t/(par[0]*par[0]) + 1./par[0];
    }
    f = LogL;

    if (2 == iflag) {
	// optionally, it is possible to calculate the gradient of FCN
	// wrt. the parameters and provide it to MINUIT in gin
	*gin = dLogLdtau;
    }
    if (3 == iflag) {
	// ending a fit, we free any resources we may hold, and do other
	// cleanup work
	// (not needed here)
    }
}

// FCN to fit lambda*exp(-lambda*t) to the data in the tuple
void fcnLambda(Int_t& npar, Double_t* gin, Double_t& f, 
	Double_t* par, Int_t iflag)
{
    if (1 == iflag) {
	// starting a new fit, initialize anything we need to have
	// (not needed here)
    }

    // calculate value for f
    double LogL = 0., dLogLdlambda = 0.;
    for (unsigned long iEv = iStartEv; iEv < iEndEv; ++iEv) {
	tuple->GetEntry(iEv);
	const double t = tuple->GetArgs()[0];
	LogL += t*par[0] - std::log(par[0]);
	if (2 == iflag)
	    dLogLdlambda += t - 1./par[0];
    }
    f = LogL;

    if (2 == iflag) {
	// optionally, it is possible to calculate the gradient of FCN
	// wrt. the parameters and provide it to MINUIT in gin
	*gin = dLogLdlambda;
    }
    if (3 == iflag) {
	// ending a fit, we free any resources we may hold, and do other
	// cleanup work
	// (not needed here)
    }
}

// do an unbinned fit to the data in tuple
// fitTau	true to fit tau, false to fit lambda
// par0		upon return, contains tau/lambda (depending on fitTau)
// par0err	upon return, contains uncertainty on fitted parameter
void unbinnedFit(bool fitTau, double& par0, double& par0err)
{
    // create a fitter object
    TMinuit* minu = new TMinuit(1);

    // set FCN
    if (fitTau)
	minu->SetFCN(fcnTau);
    else
	minu->SetFCN(fcnLambda);
    // set level of diagnostics output
    minu->SetPrintLevel(1);

    double arglist[10];	// used to pass parameters to Minuit
    int ierflg = 0;	// used by Minuit to indicate errors
    // set parameter name, starting value(2.), step size (0.1), no limits
    const char *parname = fitTau?"#tau":"#lambda";
    minu->mnparm(0, parname, 2., 0.1, 0., 0., ierflg);
    if (ierflg)
	std::cerr << "Error setting up parameters." << std::endl;

    // adjust the error definition for likelihood fits
    arglist[0] = 0.5;
    minu->mnexcm("SET ERR", arglist, 1, ierflg);
    if (ierflg)
	std::cerr << "Error setting error definition." << std::endl;
    // 500 steps, Simplex method for minimization to find good starting values
    arglist[0] = 500;
    minu->mnexcm("SIMPLEX", arglist, 1, ierflg);
    if (ierflg)
	std::cerr << "Error executing SIMPLEX commnad." << std::endl;
    // 500 steps, Migrad minimization
    minu->mnexcm("MIGRAD", arglist, 1, ierflg);
    if (ierflg)
	std::cerr << "Error executing MIGRAD commnad." << std::endl;
    // calculate Hesse Matrix (i.e. 2nd derivatives) for parameter uncertainties
    minu->mnexcm("HESSE", arglist, 1, ierflg);
    if (ierflg)
	std::cerr << "Error executing HESSE commnad." << std::endl;
    // calculate uncertainties with Minos (asymmetric errors)
    minu->mnexcm("MINOS", arglist, 1, ierflg);
    if (ierflg)
	std::cerr << "Error executing MINOS commnad." << std::endl;

    // get fit status
    double fmin, edm, errdef;
    int nvpar, nparx, icstat;
    minu->mnstat(fmin, edm, errdef, nvpar, nparx, icstat);
    // print Minuit results
    minu->mnprin(4, fmin);
    // get tau and error
    double dummy;
    TString dummystr;
    int dummyint;
    minu->mnpout(0, dummystr, par0, par0err, dummy, dummy, dummyint);

    // release Minuit object
    delete minu;
}

// fname	ROOT file name
// tuplename	name of tuple inside ROOT file
// fitTau	when true fit tau in ex. 5.2 otherwise, fit lambda = 1 / tau
// 		for ex. 5.3
void ex5(const char* fname, const char* tuplename, bool fitTau = true)
{
    // open ROOT file
    TFile infile(fname, "READ");
    // get tuple
    tuple = (TNtuple*)infile.Get(tuplename);
    if (0 == tuple) {
	std::cerr << "Cannot read tuple from file." << std::endl;
	return;
    }
    // get number of entries in tuple
    const unsigned long nentries = tuple->GetEntries();
    // set up variables and fit
    iStartEv = 0;
    iEndEv = nentries;
    double tau, tauerr;
    unbinnedFit(true, tau, tauerr);
   
    // make a histogram of the data in the tuple
    TH1D *htau = new TH1D("htau", "htau;t;entries", 50, 0., 5.);
    for (unsigned long iev = iStartEv; iev < iEndEv; ++iev) {
	tuple->GetEntry(iev);
	const double t = tuple->GetArgs()[0];
	htau->Fill(t);
    }
    // create a function according to the fitted parameters
    // (the tricky thing is the normalization; the bin width of the
    // histogram enters here)
    TF1* f = new TF1("expdecay", "[0]*exp(-x/[1])/[1]", 0., 5.);
    f->SetParameters(double(iEndEv - iStartEv) * htau->GetBinWidth(1), tau);
    // make a scan of the log likelihood function (we vary the fit parameter
    // and calculate the likelihood for each value of the fit parameter we try)
    TGraph *glh = new TGraph;
    glh->SetTitle("- Log Likelihood");
    // 100 values for the fit parameter
    for (unsigned k = 0; k < 100; ++k) {
	// values are between 0.5 and 5.0
	double mytau = .5 + (5. - .5) * double(k) / double(100);
	double y, dummy;
	int npar = 1;
	// work out - log likelihood
	fcnTau(npar, &dummy, y, &mytau, 4);
	// and add the pair (mytau, y) to the graph
	glh->SetPoint(k, mytau, y);
    }

    // plot distribution of fit parameters in many low statistics fits
    const unsigned nlowstat = 50;	// 50 events per fit
    TH1D* hpars = new TH1D("hpars",
	    "distribution of fitted parameter;param.;entries",
	    100, 0., 5.);
    TH1D* hparpull = new TH1D("hparpull",
	    "pull of fitted parameter;(param. - 1)/(param. uncertainty);entries",
	    100, -5., 5.);
    unsigned long istartev = iStartEv, iendev = iEndEv;
    for (unsigned long iev = istartev; iev < iendev; iev += nlowstat) {
	// set up fit
	iStartEv = iev;
	iEndEv = iev + nlowstat;
	if (iEndEv > iendev) iEndEv = iendev;
	// perform fit
	unbinnedFit(fitTau, tau, tauerr);
	// fill fit parameter
	hpars->Fill(tau);
	// calculate pull: (fitted - true) / uncertainty
	// why is this a good quantity to plot? which implementation
	// errors can be seen in a pull? (hint: for a fit (when implemented
	// correctly) with sufficient statistics, one expects a unit width
	// gaussian with mean zero)
	// we assume the true quantity to be 1.0 here, if you generated
	// something different, adapt below
	hparpull->Fill((tau - 1.) / tauerr);
    }

    // and then, we draw everything to a canvas
    TCanvas* c = new TCanvas;
    c->cd();
    c->Divide(2, 2);
    // start with the lifetime distribution and the fit function
    c->cd(1);
    htau->Draw("E");
    f->SetLineColor(kRed);
    f->SetLineWidth(1);
    f->Draw("SAME");
    // the -log likelihood scan is next
    c->cd(2);
    glh->Draw("AC");
    // last but not least, we have the distribution of fit parameters
    // in many low statistics fits (ex. 5.3)
    c->cd(3);
    hpars->Draw();
    // and its pull
    c->cd(4);
    hparpull->Draw();

    // and print it to an eps file (use "gv out.eps" to view the results)
    c->Update();
    c->Print("out.eps");
    delete c;

    infile.Close();
}
