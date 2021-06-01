// 1. start data analysis tool root with "root"
// 2. load the analysis macro Analyzer.C with ".L Analyzer.C"
// 3. run the macro with "Analyzer()"
// 4. after it finished exit root with ".q"
// 5. open the pdf file with "evince MyAnalysis.pdf &"


#include "TList.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

void Analyzer (TString partname)
{
  
  // some initialisations
  Int_t evn = 0;
  Int_t pdg = 0;
  Double_t px = 0;
  Double_t py = 0;
  Double_t pz = 0;
  

  // open the data file and prepare it for reading
  TFile *file = new TFile(partname+"_EvtGendecays.root", "READ");
  TTree *tree = (TTree*)file->Get("hidden treasure");
 
  TBranch *b1 = tree->GetBranch("EVN");
  b1->SetAddress(&evn);
  TBranch *b2 = tree->GetBranch("PDG");
  b2->SetAddress(&pdg);
  TBranch *b3 = tree->GetBranch("px");
  b3->SetAddress(&px);
  TBranch *b4 = tree->GetBranch("py");
  b4->SetAddress(&py);
  TBranch *b5 = tree->GetBranch("pz");
  b5->SetAddress(&pz);

  // loop over all events in the data file
  Int_t nev = b1->GetEntries();
  for (Int_t ii=0; ii<nev; ii++) {
    b1->GetEntry(ii);
    b2->GetEntry(ii);
    b3->GetEntry(ii);
    b4->GetEntry(ii);
    b5->GetEntry(ii);
    
    printf("event %i pdg %i p (%f, %f, %f)\n",evn,pdg,px,py,pz);
        
  }

  file->Close();
  
}
