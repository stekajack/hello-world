#include <TTree.h>
#include <TFile.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TF1.h>
#include <vector>
#include "TCanvas.h"
#include <cmath>
#include <TH2D.h>
#include <TDatabasePDG.h>
#include <cstring>
#include "TWignerd.h"
#include "TSmallWignerd.h"
#include <TGraph.h>
#include <sstream>


void wig_graph(Int_t l,Int_t m,Int_t n)
{
   // Use this macr to draw a wigner-d function. For more information
   // see the classes TWignerd and TSmallWignerd

   Int_t L = 1024;
   if (l > L) L = l;

   // The number of points on which the function will be sampled
   // is 2*L.

   const Int_t size = 2*L;

   TWignerd wig(l);
   wig.Advance(l);

   TSmallWignerd smallwd(L);
   Int_t nDim(0);
   Double_t *vecy = smallwd.Get(wig,m,n,nDim);

   Double_t vecx[size];
   for (Int_t j = 0; j < size; ++j) 
      vecx[j] = TMath::Pi()*(j + 1/2)/size;

   TGraph *gr = new TGraph(size,vecx,vecy);
   std::ostringstream oss;
   oss << "Wigner-d function" << l << m<< n;
   std::string var = oss.str();
   const char* foobar2 = var.c_str();

   gr->SetTitle(foobar2);
   gr->SetLineColor(kSpring+3);
   gr->SetFillColor(kSpring+3);
   gr->SetLineWidth(2);

   TCanvas *canvas = new TCanvas("canvas","Wigned d function",100,10,1000,500);
   gr->Draw("ac");
}

