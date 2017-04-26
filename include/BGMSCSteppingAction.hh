#ifndef BGMSCSTEPPINGACTION_H
#define BGMSCSTEPPINGACTION_H

#include "G4UserSteppingAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

class PrimaryGeneratorAction;

class BGMSCSteppingAction : public G4UserSteppingAction {

public:
    BGMSCSteppingAction()
    {
        i = 0;
        j = 0;
        Edep = 0;
        Edepout = 0;
        SelfAbsE = 0;
        SelfAbsG = 0;
    }

    ~BGMSCSteppingAction()
    {
        //G4cout << "Number of ionization is " << i << G4endl;
        G4cout << "Total energy deposition is " << Edep/MeV << G4endl;
        //G4cout << "Number of electrons escaping GNP is " << j << G4endl;
        G4cout << "Energy deposition by electrons proudced inside GNP is " << Edepout/MeV << G4endl;
        G4cout << "Self-absorption of electron is " << SelfAbsE/MeV << G4endl;
        G4cout << "Self-absorption of gamma is " << SelfAbsG/MeV << G4endl;
    }

    void UserSteppingAction(const G4Step* aStep);

private:
    G4int i;
    G4int j;
    G4double Edep;
    G4double Edepout;
    G4double SelfAbsE;
    G4double SelfAbsG;

    G4int InEnergy[200];
    G4int OutEnergy[200];
};

#endif
