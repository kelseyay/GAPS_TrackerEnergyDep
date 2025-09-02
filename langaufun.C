#include <iostream>
#include <TMath.h>

Double_t langaufun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location


  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0];



  // Control constants
  Double_t np;      // number of convolution steps
  Double_t sc;      // convolution extends to +-sc Gaussian sigmas
  Double_t xlow,xupp;// Range of convolution integral
  Double_t step;

  np = 500.0;
  sc =   5.0;
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / (2.*np+1);




  // Convolution integral of Landau and Gaussian by sum
  fland = TMath::Landau(x[0],mpc,par[0]); // /par[0];
  sum += fland * TMath::Gaus(x[0],x[0],par[3]);

  for(i=1.0; i<=np; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]); // / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]); // / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[0]*par[2] * step * sum * invsq2pi/ (par[3]) );
};


Double_t langaufun_max(Double_t *par) {


    // // double mpv1 = sum2langaufun(par,par);
    // // double mpv2 = sum2langaufun(par+4,par);
    // // double start = (mpv1 > mpv2 ? par[0] : par[4]);
    // // start -= TMath::Max(par[3],par[7]);
    // double start = 0.;
    // double val  = 0.;
    // double val0 = 0.;
    // for(int i=0 ;i<100000; i++){
    // 	double x = start + i*par[3]/10.;
    // 	val = langaufun(&x ,par);
    // 	if(val<val0)return x-par[3]/10.;
    // 	val0=val;
    // }
    // return -1;

  double start = par[1];
  double val  = 0.;
  double val0 = langaufun(&start ,par);
  double step = TMath::Abs(par[3]/100.);
  cout << endl << " (step "<<step<<")";

  for(int i=0 ;i<10000; i++){
    double x = start + i*step;
    val = langaufun(&x ,par);
    if(val<val0)return x - step;
    val0=val;
  }
  return -1;

}




// TF1* FitHistoLandau(TH1F* h, char* label="He4", Double_t* peak=NULL, Double_t* sig=NULL, Double_t* epeak=NULL,int npeak=1, bool draw=true){



//   if(!h)return NULL;
//   if(h->Integral()==0)return NULL;

//   int nSearch = npeak;//2;
//   TSpectrum *s = new TSpectrum(nSearch);
//   npeak = s->Search(h,5,"nobackgroundnodraw");
//   Double_t *xpeaks = s->GetPositionX();
//   TF1* fPeaks[nSearch];
//   TF1* fPeaksL[nSearch];
//   Float_t wimin = 0.15;
//   Float_t wimax = 0.25;//0.1;
//   //    cout << " N.PEAK "<<npeak<<endl;

//   if(h->Integral()<100){
//     for (int p=0;p<npeak;p++) {
//       if(peak)peak[p]  = xpeaks[p];
//       if(epeak)epeak[p] = 1000;
//       if(sig)sig[p]   = 1000;
//     }
//     return NULL;
//   }



//   for (int p=0;p<npeak;p++) {

//     //     fPeaks[p]=NULL;
//     //     fPeaksL[p]=NULL;


//     Double_t xmin = xpeaks[p]*(1.-wimin);
//     Double_t xmax = xpeaks[p]*(1.+wimax);
//     fPeaks[p] = new TF1(Form("peak%i",p),"gaus",-1000,3000.);
//     fPeaks[p]->SetParameter(0,h->Integral());
//     fPeaks[p]->SetParameter(1,xpeaks[p]);
//     fPeaks[p]->SetParameter(2,xpeaks[p]*0.2);
//     h->Fit(fPeaks[p],"RQ0","0",xmin,xmax);
//     //
//     xmin = fPeaks[p]->GetParameter(1)-fPeaks[p]->GetParameter(2)*2.;//1.5;
//     xmax = fPeaks[p]->GetParameter(1)+fPeaks[p]->GetParameter(2)*1;

//     double par[]={0.5*fPeaks[p]->GetParameter(2),
//           fPeaks[p]->GetParameter(1),
// h->GetEntries()*h->GetBinWidth(1)/(TMath::Power(0.5*fPeaks[p]->GetParameter(2),2)),
//           0.5*fPeaks[p]->GetParameter(2)};
//     //    cout << par[0]<<" "<< par[1]<<" "<< par[2]<<" "<< par[3]<<endl;
//     fPeaksL[p] = new TF1(Form("FitFNC_%s",h->GetName()),langaufun,0.,200.,4);
// fPeaksL[p]->SetParNames("Width","MP","Area/Width^2","GSigma");
//     //     fPeaksL[p]->SetParLimits(3,0.,100.);
//     //     fPeaksL[p]->SetParLimits(1,0.,1000.);
//     //     fPeaksL[p]->SetParLimits(0,0.,100.);
//     fPeaksL[p]->SetParameters(par);
//     h->Fit(fPeaksL[p],"R0","0",xmin,xmax);
//     //    cout <<"secondo fatto"<<endl;
//     fPeaksL[p]->GetParameters(par);
//     double fit1;
//     double fit2;

//     /// per variare l'intervallo del fit devi modificare questi parametri:

//     fit1 = par[1]*0.5;
//     fit2 = par[1]*1.6;
//     h->Fit(fPeaksL[p],"R","",fit1,fit2);
//     fPeaksL[p]->SetRange(fit1,fit2);
//     fPeaksL[p]->GetParameters(par);

//     cout << endl << " MAX "<< langaufun_max(par);


//     if(peak)peak[p]  =fPeaksL[p]->GetParameter(1);
//     if(epeak)epeak[p] =fPeaksL[p]->GetParError(1);
//     if(sig)sig[p]   =fPeaksL[p]->GetParameter(3);
//   }

//   //    cout << "disegna "<<endl;

//   TF1 *f = NULL;
//   if(fPeaksL[0]){
//     f = new TF1();
//     fPeaksL[0]->Copy(*f);

//     if(draw){

//       gStyle->SetOptStat(1110);
//       gStyle->SetOptFit(1);
//       TCanvas *c = new TCanvas("c");
//       c->SetTicks(1,1);
//       h->GetXaxis()->SetTitle("Segnale (canali ADC)");
//       h->GetYaxis()->SetTitle("Conteggi");
//       h->Draw();
//       c->Print(Form("Fit-%s-%s.png",h->GetName(),label));
//       c->Print(Form("Fit-%s-%s.root",h->GetName(),label));
//     }
//   }
//   //    cout << "fatto "<<endl;
//   s->Delete();
//   return f;
// }
