//
//  main.cpp
//  d_Wigner_calc
//
//  Created by Deniz Mostarac on 10/07/2017.
//  Copyright © 2017 Deniz Mostarac. All rights reserved.
//

#include <iostream>
#include <TFile.h>
#include <TList.h>
#include <TH1D.h>
#include <TF1.h>
#include <vector>
#include "TCanvas.h"
#include <cmath>
#include <TMath.h>
#include <TH2D.h>
#include <cstring>
#include "TWignerd.h"
#include "TSmallWignerd.h"
#include <TGraph.h>
#include <TGraph2D.h>
#include "TComplex.h"
#include <complex>
#include "TStyle.h"
typedef std::vector<std::vector<std::complex<Double_t>>> Matrica;
typedef std::vector<std::vector<Matrica>> matricaMatrice;

// matrix manipulation metods ////////////////////////////////////////////////////////////////////
int noRows(Matrica &m)
{
return m.size();
}

int noColumns(Matrica &m)
{
    return m[0].size();
}
bool ifAllowed(Matrica &m1, Matrica &m2)
{
    return noColumns(m1) == noRows(m2);
}
Matrica makeMatrix(int no_rows, int no_columns)
{
    return Matrica(no_rows, std::vector<std::complex<Double_t>>(no_columns));
}
Matrica multiplyMatrices(Matrica& m1, Matrica &m2)
{
    auto m3(makeMatrix(noRows(m1), noColumns(m2)));
    for(int i = 0; i < noRows(m1); i++)
        for(int j = 0; j < noColumns(m2); j++) {
            std::complex<Double_t> suma(0,0);
            for(int k = 0; k < noRows(m2); k++) suma += m1[i][k] * m2[k][j]; m3[i][j] = suma;
	    
        }
    return m3;
}
Matrica multiplyMatrices_alt(Matrica& m1, Matrica &m2,int setting=0)
{
    auto m3(makeMatrix(noRows(m1), noColumns(m2)));
    for(int i = 0; i < noRows(m1); i++)
        for(int j = 0; j < noColumns(m2); j++) {
            std::complex<Double_t> suma(0,0);
if(setting==1){
	    std::cout<< m1[i][j] *m2[i][j]<<"  "<<m1[i][j]<<"  "<< m2[i][j]<<endl;
}
            suma += m1[i][j] * m2[i][j];
            m3[i][j] = suma;
            
        }
    return m3;
}
void printMatrix(Matrica &m) {
    for(int i = 0; i < noRows(m); i++)
    {
        for(int j = 0; j < noColumns(m); j++)
        std::cout << std::setw(10) << m[i][j];
        std::cout << std::endl;
    } }
bool canBeSummed(Matrica &m1, Matrica &m2)
{
    return (noRows(m1) == noRows(m2)) && (noColumns(m1) == noColumns(m2));
}

Matrica addMatrices(Matrica& m1, Matrica& m2)
{
    auto m3(makeMatrix(noRows(m1), noColumns(m1))); for(int i = 0; i < noRows(m1); i++)
        for(int j = 0; j < noColumns(m1); j++) m3[i][j] = m1[i][j] + m2[i][j]; return m3;
}
///////////////////////////////////////////////////////////////////////////////////////////////

Matrica wig_graph(Int_t l,Int_t m,Int_t n,Int_t L = 1024)
{
    // Use this macr to draw a wigner-d function. For more information
    // see the classes TWignerd and TSmallWignerd
    
    
    if (l > L) L = l;
    
    // The number of points on which the function will be sampled
    // is 2*L.
    
    const Int_t size = 2*L;
    
    TWignerd wig(l);
    wig.Advance(l);
    
    TSmallWignerd smallwd(L);
    
    //printout section
    //don't delete "unused" variables, using return by reference!!!
    
    Int_t nDim(0);
    Double_t *vecy = smallwd.Get(wig,m,n,nDim);
    Double_t *ptrMiddleMan= vecy;
    int counter(0);
    smallwd.setNoOfElements();
    while (counter<smallwd.getNoOfElements()) {
        counter++;
    }
    
    // container generation and fill
    /*idea is to generate a square matrix that will have for each d function,a square matrix of values for 2L theta and phi values*/
    
    Matrica containterTheta;
    std::cout<<smallwd.getNoOfElements()<<" length set to containerTheta "<<std::endl;
    containterTheta.resize(smallwd.getNoOfElements());
    for (auto &element:containterTheta) {
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

    /* Here we perform matrix multiplication for one pair of d functions and the presumably complementary phi information.As a chech of the matrices have been properly generated and passed for multiplication, we employ a checker funtion!*/

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
        return aElementMatrix;
    
    Double_t vecx[size];
    for (Int_t j = 0; j < size; ++j)
        vecx[j] = TMath::Pi()*(j + 1/2)/size;
}

int sandwichGen(int J,int m,int n, int m_prim,Int_t L = 1024)
{
/////// Canvas preparation and generation, as well as histo declarations! ///////////////

    TCanvas* c1= new TCanvas("c1","myhisto",1000,800);
    c1->Divide(2);
    //TH2D* deadwood1 = new TH2D("deadwood1","a 2D hiato",2048,0,TMath::Pi(),\
    //                          2048,0,2*(TMath::Pi()));
    //TH2D* deadwood2 = new TH2D("deadwood2","a 2D hiatoo",2048,0,TMath::Pi(),\
    //                          2048,0,2*(TMath::Pi()));

////////// D_matrix generation and fill /////////////////////////////////////////////////

    matricaMatrice D_matrix;
    D_matrix.resize(2*J+1);
    for (auto &element:D_matrix) {
        element.resize(1);
    }
    int cnt(-J);
    for (auto &element:D_matrix) {
        element[0]=wig_graph(J, cnt,n,L);
        cnt++;
	std::cout<<"tadaaaa"<<std::endl;
    }
///////// transposed_D_matrix generation and fill ///////////////////////////////////////////

    matricaMatrice D_matrix_conj;
    D_matrix_conj.resize(2*J+1);
    for (auto &element:D_matrix_conj) {
        element.resize(1);
        element[0].resize(2*L);
    }
    
    for (auto &element1: D_matrix)
    {
        for (auto &element2: D_matrix_conj)
        {
            for (auto &subelement1:element1)
            {
                for (auto &subelement2:element2)
                {
                    subelement2.resize(subelement1.size());
                    for (auto &resizer1:subelement1) {
                        for (auto &resizer2:subelement2) {
                            resizer2.resize(resizer1.size());
                        }
                    }

                }
            }
        }
    }
    int selector(0);
    for (auto &element1: D_matrix)
    {	
        auto& element2=D_matrix_conj[selector];
	selector++;
        
            for (auto &subelement1:element1)
            {
                for (auto &subelement2:element2)
                {
                    for (int i(0); i<noRows(subelement1); i++)
                    {
                        for (int j(0); j<noColumns(subelement1); j++)
                        {
                            subelement2[i][j]=std::conj(subelement1[i][j]);
			   
                        }
                    }
                }
            }
    }
    
    
    
    
    matricaMatrice transposed_D_matrix;
    transposed_D_matrix.resize(1);
    for (auto &element:transposed_D_matrix) {
        element.resize(2*J+1);
    }
    int browser(0);
    for (auto &element:D_matrix_conj) {
        transposed_D_matrix[0][browser]=element[0];
        browser++;
        std::cout<<"lalalalala"<<std::endl;
    }

//////// generation of identity matrix of a matrix containing a matrix of all elements 1 ////

	std::cout<<"jedinice"<<std::endl;
    Matrica jedinice;
    jedinice.resize(2*L);
    for (auto &element:jedinice) {
        element.resize(2*L);
    }
    for (auto &element:jedinice)
    {
        for (auto &subElement: element)
        {
            subElement=std::complex<Double_t>(1,0);
        }
        
    }
	std::cout<<"id_matrix"<<std::endl;
    
    matricaMatrice id_matrix;
    id_matrix.resize(2*J+1);
    for (auto &element:id_matrix) {
        element.resize(2*J+1);
    }
    for (int i(0); i<(2*J+1); i++) {
        for (int j(0); j<(2*J+1); j++) {
            if(i==j)
            {
                id_matrix[i][j]=jedinice;
            }
        }
    }
    std::cout<<"D_matrix_id"<<std::endl;
//////////// doing the sandwich ///////////////////////////////////////////////////////

    matricaMatrice D_matrix_id;
    D_matrix_id.resize(m_prim);
    for (auto &element:D_matrix_id) {
        element.resize(1);
    }
    for (int i(0); i<m_prim; i++) {
        for (int j(0); j<m_prim;j++) {
	if(i==j)
        {
            D_matrix_id[i][0]=multiplyMatrices_alt(id_matrix[i][j], D_matrix[i][0]);
	}
        }
    }

        std::cout<<"allTogetherReadyForFill"<<std::endl;
    int convolutor(0);
    Matrica allTogetherReadyForFill;
    selector=0;
    for (auto &element:D_matrix_id)
    {
	
        auto &subelement=D_matrix_conj[selector];
        
	    std::cout<<"entered second loop "<<std::endl;
            if (convolutor==0)
            {
            	
                allTogetherReadyForFill=multiplyMatrices_alt(element[0], subelement[0]);
 		
                convolutor++;
                
            }
            else
            {
                Matrica privremeni=multiplyMatrices_alt(element[0], subelement[0]);
                if (canBeSummed(allTogetherReadyForFill, privremeni))
                {
                    allTogetherReadyForFill=addMatrices(allTogetherReadyForFill, privremeni);
                    std::cout<<convolutor<<std::endl;
                    convolutor++;
                }
                else
                {
                    std::cout<<"huston, we have a problem!!!"<<std::endl;
                    return 8;
                }
            }

         selector++;
        
    }

////////// print and plot !!! ////////////////////////////////////////////////////////

    
    //printMatrix(allTogetherReadyForFill);

    TGraph2D *g1= new TGraph2D(2*L*2*L);
    std::ostringstream oss;
    oss << "Real Angular Distribution for " << J <<" " << m<<" "<< n;
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
    c1->cd(1);
    gStyle->SetPalette(1);
    g1->Draw("surf");

    TGraph2D *g2= new TGraph2D(2*L*2*L);
    oss.str("");
    oss.clear();
    oss << "Imaginary Angular Distribution for " << J <<" " << m<<" "<< n;
    std::string var2 = oss.str();
    const char* foobar2 = var2.c_str();
    g2->SetTitle(foobar2);

    iter = 0;
    for (int j(0); j<noRows(allTogetherReadyForFill); j++)
    {
        iks=(std::cos(TMath::Pi()*(j + 1/2)/(2*L)));
        for (int k(0); k<noColumns(allTogetherReadyForFill); k++)
        {
 	    ips=(2*(TMath::Pi())*(k + 1/2)/(2*L));
            g2->SetPoint(iter,iks, ips, std::abs(imag(allTogetherReadyForFill[j][k])));
	    iter++;
        }
    }
    c1->cd(2);
    gStyle->SetPalette(1);
    g2->Draw("surf");
    



    return 0;
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
int main()
{
    wig_graph(2, -2, 0);
    return 0;
}
