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
#include "TLorentzVector_plus.h"

void treeRedefinition()
{
    // some initialisation
    Double_t E = 0;
    Int_t evn = 0;
    Int_t pdg = 0;
    Double_t px = 0;
    Double_t py = 0;
    Double_t pz = 0;
	Int_t nuki(0);
    Int_t nev;
    Int_t logic(1);
    // open the data file and prepare it for reading
    TFile *file = new TFile("f_2_EvtGendecays.root","READ");
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
    
    //generation of loader slots
    TList* list_slot=new TList();
    list_slot->SetOwner(kTRUE);
    TFile target("rafined_tree.root","recreate");
    TTree tree_target("tree_target","target tree");
    tree_target.Branch("Redefined_Tlist",list_slot);
    nev = b1->GetEntries();
    TDatabasePDG database;
	gSystem->Load("TLorentzVector_plus_h.so");
    //tree loader logic
    for (int ii(0); ii<nev; ii++)
    {
        getoverhere:
        b1->GetEntry(ii);
        b2->GetEntry(ii);
        b3->GetEntry(ii);
        b4->GetEntry(ii);
        b5->GetEntry(ii);
        E=sqrt(pow((database.GetParticle(pdg))->Mass(),2)+pow(px,2)+pow(py, 2)+pow(pz, 2));
        
        TLorentzVector_plus* slot= new TLorentzVector_plus;
        slot->SetPxPyPzE(px, py, pz, E);
	slot->SetPDGId(pdg);
        
        if (evn==logic)
        {
            list_slot->Add(slot);
	    std::cout<<slot->GetPDGId()<< " ";
        }
        else
        {
            logic+=1;
            tree_target.Fill();
            list_slot->Clear();
            goto getoverhere;
        }
        if(ii==((nev-1)))
        {
            tree_target.Fill();
            list_slot->Clear();
            target.Write();
            target.Close();
            delete list_slot;
        }
        if (ii%10==0) {
            std::cout<<(ii/10)<<" ";
        }
    }
    
    
}
void Printer ()
{
    Int_t nev, npart;
    
    TList *ll = new TList();
    TLorentzVector lv;
    
    TFile *file = new TFile("rafined_tree.root", "READ");
    TTree *tree = (TTree*)file->Get("tree_target");
    TBranch *branch = tree->GetBranch("Redefined_Tlist");
    branch->SetAddress(&ll);
    
    nev = branch->GetEntries();
    for (Int_t ii=0; ii<nev; ii++) {
        printf("Event number %i\n",ii+1);
        branch->GetEntry(ii);
        
        npart = ll->GetEntries();
        for (Int_t jj=0; jj<npart; jj++)
        {
            printf("part[%i] ",jj+1);
            ll->At(jj)->Print();
	    std::cout<<((TLorentzVector_plus*)(ll->At(jj)))->GetPDGId()<< " ";
        }
        
        printf("\n");
    }
    
    file->Close();    
}
void Plotter_H(Int_t PDGId)
{
    TCanvas* c1 = new TCanvas("c1","myhisto",1000, 800);
    c1->Divide(2,2);
    TH2D* htwoD = new TH2D("htwoD","a 2D hiato",100,-3,3,\
                           100,-3,3);
    TH1D *h3= new TH1D("h3","h3",100,-3,5);
    TH1D *h2 = new TH1D("h2","h2",100,-3,3);
    TH1D *h = new TH1D("h","h",100,-3,3);
    Int_t nev, npart;
    TList *ll = new TList();
    TFile *file = new TFile("transformed_tree.root", "update");
    TTree *tree = (TTree*)file->Get("tree_target");
    TBranch *branch = tree->GetBranch("Redefined_Tlist");
    branch->SetAddress(&ll);
    std::vector<TLorentzVector_plus*> container;
    container.reserve(4);
    std::vector<TLorentzVector_plus*> kanter;
    kanter.reserve(3);
    double beamAxis(0);
    TLorentzVector_plus motherMomenta;
    TLorentzVector_plus bAxis;
    
    
    nev = branch->GetEntries();
    for (Int_t ii=0; ii<nev; ii++)
    {
        
        branch->GetEntry(ii);
        
        npart = ll->GetEntries();
        TLorentzVector_plus lvsum;
	lvsum= TLorentzVector_plus(TLorentzVector(0,0,0,0));
        ////////////////////////////////////////////     LOADER!!!
        for (Int_t jj=0; jj<npart; jj++)
        {
            container.push_back((TLorentzVector_plus*)(ll->At(jj)));
            
        }
        /////////////////////////////////////////////     GENERATOR!!!
        for (Int_t i=0; i<static_cast<int>(container.size()); i++)
        {
	    std::cout<<lvsum.Pt()<<" ";
            lvsum+=(*(container[i]));
	    std::cout<<lvsum.Pt()<<" ";
   
        }
        h->Fill(lvsum.M());
	//std::cout<<"lvsum momenta are now beeing printed!!!!"<<endl;
        std::cout<<lvsum.Pt()<< " ";
        ////////////////////////////////////////////      OFLOADER!!!
        for (int k(0); k<npart; k++)
        {
            container.pop_back();
        }
        /////////////////////////////////////////////     printf("\n");
        
        //        if (ii%100==0)
        //        {
        //            std::cout<<(ii/100);
        //            std::cout<< " ";
        //        }
    }
    c1->cd(1);
    h->Draw();
    h->SetLineColor(kGreen);
    h->SetLineWidth(2);
    double fitmin=1;
    double fitmax=2;
    TF1* func = new TF1("Fit", "[0]*([1]/((x-[2])*(x-[2]) + [1]*[1]/4))", fitmin, fitmax);
    func->SetParameters(10000,1.335);
    func->SetLineColor(2);
    h->Fit(func);
    func->Draw("same");
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////      PHI!!!!!!!
    nev = branch->GetEntries();
    for (Int_t ii=0; ii<nev; ii++)
    {
        //    printf("Event number %i\n",ii+1);
        branch->GetEntry(ii);
        npart = ll->GetEntries();
	////////////////////////////////////////////     LOADER!!!
        for (Int_t jj=0; jj<npart; jj++)
        {
            container.push_back((TLorentzVector_plus*)(ll->At(jj)));
            
        }
        /////////////////////////////////////////////     GENERATOR!!!
        for (Int_t i=0; i<static_cast<int>(container.size()); i++)
        {
            h2->Fill(container[i]->Phi());
            
        }
	////////////////////////////////////////////      OFLOADER!!!
        for (int k(0); k<npart; k++)
        {
            container.pop_back();
        }
        /////////////////////////////////////////////     printf("\n");
        
        //        if (ii%100==0)
        //        {
        //            std::cout<<(ii/100);
        //            std::cout<< " ";
        //        }
    }
    c1->cd(2);
    h2->SetLineColor(kBlue);
    h2->SetLineWidth(2);
    h2->Draw();
////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////     ETA(THETA)!!!!!
    nev = branch->GetEntries();
    for (Int_t ii=0; ii<nev; ii++)
    {
        //    printf("Event number %i\n",ii+1);
        branch->GetEntry(ii);
        npart = ll->GetEntries();
	////////////////////////////////////////////     LOADER!!!
        for (Int_t jj=0; jj<npart; jj++)
        {
            container.push_back((TLorentzVector_plus*)(ll->At(jj)));
            
        }
        /////////////////////////////////////////////     GENERATOR!!!
        for (Int_t i=0; i<static_cast<int>(container.size()); i++)
        {
            h3->Fill(2*atan(exp(-(container[i]->Eta()))));
            
        }
	////////////////////////////////////////////      OFLOADER!!!
        for (int k(0); k<npart; k++)
        {
            container.pop_back();
        }
        /////////////////////////////////////////////     printf("\n");
        
        //        if (ii%100==0)
        //        {
        //            std::cout<<(ii/100);
        //            std::cout<< " ";
        //        }
    }
    c1->cd(3);
    h3->SetLineColor(kBlack);
    h3->SetLineWidth(3);
    h3->Draw(); 
	/////////////////
    /////////////////////////////////////////////////////////////////////HISTO!!!!!!
    
    nev = branch->GetEntries();
    for (Int_t ii=0; ii<nev; ii++)
    {
        //    printf("Event number %i\n",ii+1);
        branch->GetEntry(ii);
        
        npart = ll->GetEntries();
        for (Int_t jj=0; jj<npart; jj++)
        {
		
            container.push_back((TLorentzVector_plus*)(ll->At(jj)));
            kanter.push_back((TLorentzVector_plus*)(ll->At(jj)));
		
            
        }
        /////////////////////////////////////////////     GENERATOR!!!
        for (Int_t i=0; i<static_cast<int>(container.size()); i++)
        {
		if(container.size()==2)
		{
			std::cout<<(container[i]->GetPDGId())<<std::endl;
			 if((container[i]->GetPDGId())==PDGId)
				{
			htwoD->Fill(2*atan(exp(container[i]->Eta())),kanter[i]->Phi());		
				}
		}
            
        }
        for (int k(0); k<npart; k++)
        {
            container.pop_back();
            kanter.pop_back();
        }
        /////////////////////////////////////////////     printf("\n");
        
        //        if (ii%100==0)
        //        {
        //            std::cout<<(ii/100);
        //            std::cout<< " ";
        //        }
    }
    c1->cd(4);
    //htwoD->Rebin2D(2,2);
    htwoD->Draw("col");
    //    h->Write();
    //    h2->Write();
    //    file->Close();
}

void Transform_Calc()
{
    
    Int_t nev, npart;
    TList *ll = new TList();
    TFile *file = new TFile("rafined_tree.root", "update");
    TTree *tree = (TTree*)file->Get("tree_target");
    TBranch *branch = tree->GetBranch("Redefined_Tlist");
    branch->SetAddress(&ll);
    std::vector<TLorentzVector_plus*> container;
    container.reserve(4);
    double beamAxis(0);
    TLorentzVector_plus motherMomenta;
    TLorentzVector_plus bAxis;
    TLorentzVector_plus xAxis;
    TVector3 motherMomenta_TV3;
    
    
    
    nev = branch->GetEntries();
    int brojac(0);
    for (Int_t ii=0; ii<nev; ii++)
    {
        //                printf("Event number %i\n",ii+1);
        branch->GetEntry(ii);
        
        npart = ll->GetEntries();
        for (Int_t jj=0; jj<npart; jj++)
        {
            container.push_back((TLorentzVector_plus*)(ll->At(jj)));
        }
        
        ///////////////////////////////////////////     GENERATOR!!!
	//reset motherMomenta;
      	motherMomenta = TLorentzVector_plus(TLorentzVector(0,0,0,0));
        for (Int_t i=0; i<static_cast<int>(container.size()); i++)
        {
            motherMomenta+=(*(container[i]));
            
        }
        //            std::cout<<motherMomenta.Pz()<<"    ";
        
        //            if (container.size()==npart) {
        //                std::cout<<(brojac)<< " ";
        //                brojac++;
        //            }
        
        bAxis.SetPxPyPzE(0,0,1,0);
	xAxis.SetPxPyPzE(0,1,0,0);
	motherMomenta_TV3=motherMomenta.Vect();
	Int_t jj=0;
        for (int lll(0); lll<static_cast<int>(container.size()); lll++)
        {
	    container[lll]->Boost(-(motherMomenta.BoostVector()));
	    container[lll]->RotateX((Double_t)bAxis.Angle(motherMomenta.Vect()));
            container[lll]->RotateY((Double_t)xAxis.Angle(motherMomenta.Vect()));
	    std::cout<<container[lll]->Px()<<" "<<container[lll]->Py()<<" "<<container[lll]->Pz()<<" "<<container[lll]->Pt()<<std::endl;
	npart = ll->GetEntries();
       
		std::cout<<((TLorentzVector_plus*)(ll->At(jj)))->Px()<<" "<<((TLorentzVector_plus*)(ll->At(jj)))->Py()<<" "<<((TLorentzVector_plus*)(ll->At(jj)))->Pz()<<" "<<((TLorentzVector_plus*)(ll->At(jj)))->Pt()<<std::endl;
	    jj++;

        }
	std::cout<<"next dude"<<std::endl;
        for (int k(0); k<npart; k++)
        {
            container.pop_back();
            
        }
    }

}
void Tester_Ocm()
{
    
    Int_t nev, npart;
    TList *ll = new TList();
    TFile *file = new TFile("transformed_tree.root", "update");
    TTree *tree = (TTree*)file->Get("tree_target");
    TBranch *branch = tree->GetBranch("Redefined_Tlist");
    branch->SetAddress(&ll);
    std::vector<TLorentzVector_plus*> container;
    container.reserve(4);
    TLorentzVector_plus motherMomenta;
    nev = branch->GetEntries();
    int brojac(0);
    for (Int_t ii=0; ii<nev; ii++)
    {
        branch->GetEntry(ii);
        npart = ll->GetEntries();
        for (Int_t jj=0; jj<npart; jj++)
        {
            container.push_back((TLorentzVector_plus*)(ll->At(jj)));
        }
      	motherMomenta = TLorentzVector_plus(TLorentzVector(0,0,0,0));
        for (Int_t i=0; i<static_cast<int>(container.size()); i++)
        {
            motherMomenta+=(*(container[i]));
            
        }
        
        std::cout<< motherMomenta.Pz() <<" "<< motherMomenta.Px() <<" "<< motherMomenta.Py() <<std::endl;
        for (int k(0); k<npart; k++)
        {
            container.pop_back();
            
        }
    }

}


void Checker()
{
    
    Int_t nev, npart;
    TList *ll = new TList();
    TFile *file = new TFile("rafined_tree.root", "update");
    TTree *tree = (TTree*)file->Get("tree_target");
    TBranch *branch = tree->GetBranch("Redefined_Tlist");
    branch->SetAddress(&ll);
    std::vector<TLorentzVector_plus*> container;
    container.reserve(4);
    double beamAxis(0);
    TLorentzVector_plus motherMomenta;
    TLorentzVector_plus bAxis;
    TLorentzVector_plus yAxis;
    TVector3 motherMomenta_TV3;
    
    
    
    nev = branch->GetEntries();
    for (Int_t ii=0; ii<nev; ii++)
    {
        //                printf("Event number %i\n",ii+1);
        branch->GetEntry(ii);
        
        npart = ll->GetEntries();
        for (Int_t jj=0; jj<npart; jj++)
        {
            container.push_back((TLorentzVector_plus*)(ll->At(jj)));
        }
        
        ///////////////////////////////////////////     GENERATOR!!!
        motherMomenta = TLorentzVector_plus(TLorentzVector(0,0,0,0));
        for (Int_t i=0; i<static_cast<int>(container.size()); i++)
        {
            motherMomenta+=(*(container[i]));
            
        }
        //            std::cout<<motherMomenta.Pz()<<"    ";
        
        //            if (container.size()==npart) {
        //                std::cout<<(brojac)<< " ";
        //                brojac++;
        //            }
        
        bAxis.SetPxPyPzE(0,0,1,0);
        yAxis.SetPxPyPzE(0, 1, 0, 0);
        motherMomenta_TV3=motherMomenta.Vect();
        for (int lll(0); lll<static_cast<int>(container.size()); lll++)
        {
            
            std::cout<<motherMomenta.Px()<<" "<<motherMomenta.Py()<<" "<<motherMomenta.Pz()<<" "<<motherMomenta.Pt()<<std::endl;
            
            motherMomenta.Boost(-(motherMomenta.BoostVector()));
            
            
            std::cout<<motherMomenta.Px()<<" "<<motherMomenta.Py()<<" "<<motherMomenta.Pz()<<" "<<motherMomenta.Pt()<<std::endl;
        }
	std::cout<<"next dude"<<std::endl;
        for (int k(0); k<npart; k++)
        {
            container.pop_back();
            
        }
    }
    
}


void transformTreeRedefinition()
{
    Int_t nev, npart;
    TList *ll = new TList();
    TFile *file = new TFile("rafined_tree.root", "update");
    TTree *tree = (TTree*)file->Get("tree_target");
    TBranch *branch = tree->GetBranch("Redefined_Tlist");
    branch->SetAddress(&ll);
    std::vector<TLorentzVector_plus*> container;
    container.reserve(4);
    double beamAxis(0);
    TLorentzVector_plus motherMomenta;
    TLorentzVector_plus bAxis;
    TLorentzVector_plus xAxis;
    TVector3 motherMomenta_TV3;
    
    
    TList* list_slot=new TList();
    list_slot->SetOwner(kTRUE);
    TFile target("transformed_tree.root","recreate");
    TTree tree_target("tree_target","target tree");
    tree_target.Branch("Redefined_Tlist",list_slot);
    nev = branch->GetEntries();
    npart = ll->GetEntries();
    
    for (int ii(0); ii<nev; ii++)
    {
        branch->GetEntry(ii);
        npart = ll->GetEntries();
        for (Int_t jj=0; jj<npart; jj++)
        {
            container.push_back((TLorentzVector_plus*)(ll->At(jj)));
        }
        motherMomenta = TLorentzVector_plus(TLorentzVector(0,0,0,0));
        for (Int_t i=0; i<static_cast<int>(container.size()); i++)
        {
            motherMomenta+=(*(container[i]));
            
        }
        
        bAxis.SetPxPyPzE(0,0,1,0);
        xAxis.SetPxPyPzE(0,1,0,0);
        for (int lll(0); lll<static_cast<int>(container.size()); lll++)
        {
            container[lll]->Boost(-(motherMomenta.BoostVector()));
            container[lll]->RotateX((Double_t)bAxis.Angle(motherMomenta.Vect()));
            container[lll]->RotateY((Double_t)xAxis.Angle(motherMomenta.Vect()));
            TLorentzVector_plus* slot= new TLorentzVector_plus(*container[lll]);
            list_slot->Add(slot);
        }
        tree_target.Fill();
	list_slot->Clear();
	for (int k(0); k<npart; k++)
        {
            container.pop_back();
            
        }
        if (ii%10==0) {
            std::cout<<(ii/10)<< " ";
        }
        
    }
    
    target.Write();
    target.Close();
    
    delete list_slot;
}

#include <cstdlib> // for rand() and srand()
#include <ctime>

int getRandomNumber(int min, int max)
{
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);  // static used for efficiency, so we only calculate this value once
    // evenly distribute the random number across our range
    return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}






