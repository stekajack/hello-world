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
#include <TGraph2D.h>
#include <sstream>
#include "TStyle.h"
typedef std::vector<std::vector<std::complex<Double_t>>> Matrica;
typedef std::vector<std::vector<Matrica>> matricaMatrice;

int noRows(Matrica &m)
{
    return m.size();
}

int noColumns(Matrica &m)
{
    return m[0].size();
}
Matrica makeMatrix(int no_rows, int no_columns)
{
    return Matrica(no_rows, std::vector<std::complex<Double_t>>(no_columns));
}
Matrica multiplyMatrices_alt(Matrica& m1, Matrica &m2,int setting=0)
{
    auto m3(makeMatrix(noRows(m1), noColumns(m2)));
    for(int i = 0; i < noRows(m1); i++)
        for(int j = 0; j < noColumns(m2); j++) {
            std::complex<Double_t> suma(0,0);
            if(setting==1){
                std::cout<< m1[i][j] *m2[i][j]<<"  "<<m1[i][j]<<"  "<< m2[i][j]<<std::endl;
            }
            suma += m1[i][j] * m2[i][j];
            m3[i][j] = suma;
            
        }
    return m3;
}


void wigner_plotter(Int_t l,Int_t m,Int_t n)
{
    TCanvas *canvas = new TCanvas("canvas","Wigned d function",100,10,1000,500);
    canvas->Divide(2);
    
        canvas->cd(1);
        Int_t L = 256;
        if (l > L) L = l;

        const Int_t size = 2*L;

        TWignerd wig(l);
        wig.Advance(l);

        TSmallWignerd smallwd(L);
        Int_t nDim(0);
        Double_t *vecy = smallwd.Get(wig,m,n,nDim);
        Double_t *ptrMiddleMan= vecy;
        int counter(0);
        smallwd.setNoOfElements();

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
   gr->Draw("ac");


        Matrica containterTheta;
        std::cout<<smallwd.getNoOfElements()<<" length set to containerTheta "<<std::endl;
        containterTheta.resize(smallwd.getNoOfElements());
        for (auto &element:containterTheta)
        {
            element.resize(1);
        }
        counter=0;
        for (auto &element:containterTheta)
        {
            for (auto &subElement:element)
            {
                subElement=(*(ptrMiddleMan+counter));
                counter++;
            }
        }
        std::cout<<noRows(containterTheta)<<" number of rows for containerTheta after filling "<<std::endl;
    std::cout<<noColumns(containterTheta)<<" number of columns for containerTheta after filling "<<std::endl;

        Matrica containterPhi;
        containterPhi.resize(smallwd.getNoOfElements());
        for (auto &element:containterPhi) {
            element.resize(1);
        }
        
        Int_t dummyPhi(0);
        for (auto &element:containterPhi)
        {
            for (auto &subElement: element)
            {
                subElement= std::exp(std::complex<Double_t>(0,m*(2*(TMath::Pi())*(dummyPhi + 1/2)/size)));
                dummyPhi++;
            }
            
        }
std::cout<<noRows(containterPhi)<<" number of rows for containterPhi after filling "<<std::endl;
    std::cout<<noColumns(containterPhi)<<" number of columns for containterPhi after filling "<<std::endl;

        Matrica aElementMatrix;
        aElementMatrix.resize(smallwd.getNoOfElements());
        for (auto &element:aElementMatrix) {
            element.resize(smallwd.getNoOfElements());
        }
        
        for (int i(0); i<noRows(containterPhi); i++)
        {
            for (auto &element1:containterTheta[i])
            {
                for (int j(0); j<noRows(containterPhi); j++)
                {
                    for (auto &element2:containterPhi[j])
                    {
                        aElementMatrix[i][j]=element1*element2;
                    }
                }
            }
        }
        std::cout<<noRows(aElementMatrix)<<" number of rows for aElementMatrix after filling "<<std::endl;
        std::cout<<noColumns(aElementMatrix)<<" number of columns for aElementMatrix after filling "<<std::endl;


        Matrica aElementMatrix_conj;
        aElementMatrix_conj.resize(noRows(aElementMatrix));
        for (auto &element:aElementMatrix_conj) {
            element.resize(noColumns(aElementMatrix));
        }
        
        for (int i(0);i<noRows(aElementMatrix);i++) {
            for (int j(0); j<noColumns(aElementMatrix); j++) {
                aElementMatrix_conj[i][j]=std::conj(aElementMatrix[i][j]);
            }
        }
        
        Matrica allTogetherReadyForFill;
        allTogetherReadyForFill=multiplyMatrices_alt(aElementMatrix, aElementMatrix_conj);

        TGraph2D *g1= new TGraph2D(2*L*2*L);
        oss.str("");
        oss.clear();
        oss << "Real Angular Distribution for " << l <<" " << m <<" "<< n;
        std::string var1 = oss.str();
        const char* foobar1 = var1.c_str();
        g1->SetTitle(foobar1);
        
        int iter(0);
        Double_t iks, ips;
        for (int j(0); j<noRows(allTogetherReadyForFill); j++)
        {
            iks=(TMath::Pi()*(j + 1/2)/(2*L));
            for (int k(0); k<noColumns(allTogetherReadyForFill); k++)
            {
                ips=(2*(TMath::Pi())*(k + 1/2)/(2*L));
                g1->SetPoint(iter,iks, ips, real(allTogetherReadyForFill[j][k]));
                iter++;
            }
        }
        canvas->cd(2);
        gStyle->SetPalette(1);
        g1->Draw("surf");

        

        
        
}


