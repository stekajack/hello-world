#include <TLorentzVector.h>

#ifndef TLorentzVector_plus_h
#define TLorentzVector_plus_h

class TLorentzVector_plus:public TLorentzVector
{
private:
    Int_t PDGId;
public:
    TLorentzVector_plus()
    {
        PDGId=0;
    }
    TLorentzVector_plus(const TLorentzVector & lorentzvector): TLorentzVector(lorentzvector)
    {
        PDGId=0;
    }
    void SetPDGId(Int_t x)
    {
        PDGId=x;
    }
    Double_t GetPDGId()
    {
        return PDGId;
    }
    TLorentzVector_plus & operator = (const TLorentzVector &);
};
TLorentzVector_plus &TLorentzVector_plus::operator = (const TLorentzVector & q)
{
    
    return *this;
}

#endif 
