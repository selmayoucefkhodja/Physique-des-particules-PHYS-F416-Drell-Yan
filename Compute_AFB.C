#include <iostream>
#include <math.h>
#include <string>
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"
#include <fstream>
#include <string>
#include <typeinfo>
#include <sstream>
#include "TCanvas.h"
#include "TPad.h"
#include <complex> 
// In this exercise, you should:
// 1) Assign the correct values to the physical constants
// 2) Implement a function that returns the real part of R
// 3) Make everything work together to give you the AFB for u quarks
// 3) Add the necessary parts in order to compute AFB for d quarks

// Computes axial-vector coupling
float Compute_gA(float I3){
   return -I3;
}

// Computes vector coupling
float Compute_gV(float I3, float Q, float sinWeakAngleSqr){
   return I3-(2*Q*sinWeakAngleSqr);
}

// Computes square module of R
float Compute_ModR2(float Q_l, float Q_q, float sinWeakAngleSqr, float MZ, float GammaZ, float sprime){

   float sin2WeakAngleSqr = pow(sin(2*asin(sqrt(sinWeakAngleSqr))),2);

   float RFactor = 1/(Q_l*Q_q*sin2WeakAngleSqr);
   float num = sprime;
   float ReDen = sprime-pow(MZ,2);
   float ImDen = sprime*(GammaZ/MZ);
   return pow(RFactor,2)*pow(num,2)/(pow(ReDen,2)+pow(ImDen,2));
}

// Computes real part of R
float Compute_RealPArt(float Q_l, float Q_q, float sinWeakAngleSqr, float MZ, float GammaZ, float sprime){

   float sin2WeakAngleSqr = pow(sin(2*asin(sqrt(sinWeakAngleSqr))),2);

   float RFactor = 1/(Q_l*Q_q*sin2WeakAngleSqr);
   float num = sprime*(sprime-pow(MZ,2));
   float ReDen = sprime-pow(MZ,2);
   float ImDen = sprime*(GammaZ/MZ);

   return RFactor*num/(pow(ReDen,2)+pow(ImDen,2));
}



// Draws a 1D histogram and saves the result in pdf and root format
void DrawHisto1D(TCanvas *c, TH1F *Histo, string name){

   gStyle->SetOptStat(0);

   c->cd();

   Histo->GetXaxis()->SetTitleOffset(1.2);
   Histo->Draw();

   std::string pdfname = name+".pdf";
   const char *pdfcharname = pdfname.c_str();
   c->SaveAs(pdfcharname);
}

void Compute_AFB(){

   // Define physical constants
   float Q_e = -1.0 ; float Q_u = 2.0/3.0 ; float Q_d = -1.0/3.0;
   float I3_e = -1.0/2.0; float I3_u = 1.0/2.0; float I3_d = -1.0/2.0;
   float sinWeakAngleSqr = 0.22290;
   float MZ = 91.1876 ; float GammaZ = 2.5 ;

   // Declare histograms
   TH1F *Histo_AFB_u = new TH1F("Histo_AFB_u","AFB_u;#sqrt{s'} [GeV];A_{FB}",10000,0,1000);
   TH1F *Histo_AFB_d = new TH1F("Histo_AFB_d","AFB_d;#sqrt{s'} [GeV];A_{FB}",10000,0,1000);


   // Compute couplings
   float gVe = Compute_gV(I3_e,Q_e,sinWeakAngleSqr);
   float gAe = Compute_gA(I3_e);
   float gVu = Compute_gV(I3_u,Q_u,sinWeakAngleSqr);
   float gAu = Compute_gA(I3_u);
   float gVd = Compute_gV(I3_d,Q_d,sinWeakAngleSqr);
   float gAd = Compute_gA(I3_d);

   // Loop over the bins in the histogramCompute_RealPArt(Q_l, Q_q,sinWeakAngleSqr, MZ, GammaZ, sprime)
   for(int bin=1; bin<=Histo_AFB_u->GetNbinsX(); bin++){

     // Turns each bin into a sprime value
     float sprime = pow(Histo_AFB_u->GetBinCenter(bin),2); 

     float ModR2_u = Compute_ModR2(Q_e,Q_u,sinWeakAngleSqr,MZ,GammaZ,sprime);
     float ModR2_d = Compute_ModR2(Q_e,Q_d,sinWeakAngleSqr,MZ,GammaZ,sprime);
     
     float R_u_real = Compute_RealPArt(Q_e, Q_u,sinWeakAngleSqr, MZ, GammaZ, sprime);
     float R_d_real = Compute_RealPArt(Q_e, Q_d,sinWeakAngleSqr, MZ, GammaZ, sprime);

     float c1_u = 1 + 2*R_u_real*gVe*gVu + ModR2_u*(pow(gVe,2)+pow(gAe,2))*(pow(gVu,2)+pow(gAu,2));
     float c1_d = 1 + 2*R_d_real*gVe*gVd + ModR2_d*(pow(gVe,2)+pow(gAe,2))*(pow(gVd,2)+pow(gAd,2));
     float c2_u = 4*R_u_real*gAe*gAu + 8*ModR2_u*gVe*gAe*gVu*gAu;
     float c2_d = 4*R_d_real*gAe*gAd + 8*ModR2_d*gVe*gAe*gVd*gAd;


     float AFB_u = (3*c2_u)/(8*c1_u);
     float AFB_d = (3*c2_d)/(8*c1_d);

     Histo_AFB_u->SetBinContent(bin,AFB_u);
     Histo_AFB_d->SetBinContent(bin,AFB_d);

   }

   // Draw histogram
   TCanvas *c_u = new TCanvas("c_u","c_u",900,800);
   TCanvas *c_d = new TCanvas("c_d","c_d",900,800);
   DrawHisto1D(c_u,Histo_AFB_u,"AFB_u");
   DrawHisto1D(c_d,Histo_AFB_d,"AFB_d");

   // Save histogram in a ROOT file
   TFile *file_AFB_u = new TFile("file_AFB_u.root","RECREATE");
   TFile *file_AFB_d = new TFile("file_AFB_d.root","RECREATE");
   Histo_AFB_u->Write();
   Histo_AFB_d->Write();
   file_AFB_u->Close();
   file_AFB_d->Close();

}

