// g++ -o likelihood_example likelihood_example.cpp -I$ROOTSYS/include -L$ROOTSYS/lib -lCore -lCint -lHist -lMatrix -lMinuit -lGraf -lm -ldl

#include <cstdio> 
#include <cmath> 
#include <TMinuit.h>

// global data
// better solution: read in from file
#define NDATA 18
  int r[NDATA]    = {1,1,  5, 4,2,0,3,2,  4,1,2,1,1,0,1,1,2,1};
  int rfak[NDATA] = {1,1,120,24,2,1,6,2, 24,1,2,1,1,1,1,1,2,1};

void results(TMinuit *minuit)
{
// helper function for nice formatting 

#define MAXPAR 50        

   double fmin,fedm;   
   double errdef;    
   int    nparv,nparx; 
   int    fstat;      

   TString pname;         
   double pvalue,perror; 
   double plbound,pubound;
   int    pvari;        
   double eplus,eminus;  
   double gcorr;         

   double emat[MAXPAR][MAXPAR]; 
   double kmat[MAXPAR][MAXPAR]; 

   int i,j;                                 
 
   minuit->mnstat(fmin,fedm,errdef,nparv,nparx,fstat);

   printf("\n\n");
   printf("Results of MINUIT minimisation\n");
   printf("-------------------------------------\n\n");
   
   printf(" Minimal function value:              %8.3lf  \n",fmin);
   printf(" Estimated difference to true minimum: %11.3le \n",fedm);
   printf(" Number of parameters:         %3i     \n",nparv);
   printf(" Error definition (Fmin + Delta):      %8.3lf  \n",errdef);
   if(fstat==3) {
    printf(" Exact covariance matrix.\n");
   }
   else {
    printf(" No/error with covariance matrix.\n");
    printf(" Error code: %i3\n",fstat);
   }
   printf("\n");


   printf("   Parameter     Value       Error    positive    negative    L_BND    U_BND\n");
   for(i=0;i<nparx;i++) 
    {
        minuit->mnpout(i,pname,pvalue,perror,plbound,pubound,pvari);
        if(pvari>0) // only variable parameters
            {                      
                minuit->mnerrs(i,eplus,eminus,perror,gcorr);
                printf("%2i %10s %10.3le %10.3le %+10.3le %10.3le %8.1le %8.1le\n",
                i,pname.Data(),pvalue,perror,eplus,eminus,plbound,pubound);
            }
   }

   minuit->mnemat(&emat[0][0],MAXPAR);

   for(i=0; i<nparv; i++) {
    for(j=0; j<nparv; j++){ 
     kmat[i][j] = sqrt(emat[i][i]*emat[j][j]); 
      if(kmat[i][j]>1E-80) {
       kmat[i][j] = emat[i][j]/kmat[i][j];
      }
      else kmat[i][j] = 0.0;
     }
   }

   printf("\n");
   printf("Covariance matrix: \n");
   for(i=0; i<nparv; i++) {
    for(j=0; j<nparv; j++){ printf(" %10.3le",emat[i][j]);}
    printf("\n");
   }
   printf("\n");

   printf("Correlation matrix: \n");
   for(i=0; i<nparv; i++) {
     for(j=0; j<nparv; j++){ printf(" %6.3lf",kmat[i][j]);}
     printf("\n");
   }
   printf("\n");
}

/*  negative log likelihood of Poisson distribution*/
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
   int i;
   double mu, lnL;

   mu = par[0];        //  current parameter value
   lnL = 0.0;          //  log likelihood sum
   for(i=0;i<NDATA; i++) {
     lnL += r[i]*log(mu) - mu - log(double(rfak[i]));
   }
   f = -lnL;            // remember: *negative* log likelihood
}

int main() 
{ 

   double arglist[10];    
   int ierflg = 0;        

   TMinuit minuit(1);  

   minuit.SetFCN(fcn);

   arglist[0] = 0.5;
   minuit.mnexcm("SET ERR",arglist,1,ierflg);

   minuit.mnparm(0,"Poisson mu",1.0,0.1,0,0,ierflg);

   minuit.mnexcm("MIGRAD",arglist,0,ierflg);

   minuit.mnexcm("MINOS",arglist,0,ierflg);

   results(&minuit);

   return 0;
}
