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
#include <TGraph.h>
#include "TComplex.h"
#include <complex>
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
Matrica multiplyMatrices_alt(Matrica& m1, Matrica &m2)
{
    auto m3(makeMatrix(noRows(m1), noColumns(m2)));
    for(int i = 0; i < noRows(m1); i++)
        for(int j = 0; j < noColumns(m2); j++) {
            std::complex<Double_t> suma(0,0);
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
void program()
{
TCanvas* c1= new TCanvas("c1","myhisto",1000,800);
    c1->Divide(2);
    TH2D* deadwood1 = new TH2D("deadwood1","a 2D hiato",2,0,2,2,0,2);
Matrica teta;
    teta.resize(2);
    for (auto &element:teta)
    {
        element.resize(2);
    }
    int ja(1);
    for (auto &element:teta)
    {
        for (auto &subelement:element)
        {
            subelement=ja;
            ja++;
        }
    }
for (int i(0); i<noRows(teta); i++) {
        for (int j(0); j<noColumns(teta); j++) {
		
            deadwood1->Fill(i,j,real(teta[i][j]));
			
        }
c1->cd(1);
deadwood1->Draw("box");
}}
