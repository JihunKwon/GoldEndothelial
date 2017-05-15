#include "BGMSCPrimaryGeneratorAction.hh"

#include "G4Proton.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

#include <G4SingleParticleSource.hh>

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#define WORLD_SIDE (100. * um)

using namespace CLHEP;

BGMSCPrimaryGeneratorAction::BGMSCPrimaryGeneratorAction()
{

//    //Electron Source
//        CircleSource = new G4SingleParticleSource;
//        CircleSource->SetParticleDefinition(G4Electron::Definition());
//        this->setEnergyToBeta();

////        G4SPSEneDistribution* energy = CircleSource->GetEneDist();
////        energy->SetMonoEnergy(10*MeV);

//        G4SPSPosDistribution* position = CircleSource->GetPosDist();
//        position->SetPosDisType("Plane");
//        position->SetPosDisShape("Circle");
//        position->SetRadius(radius*nm);
//        position->SetCentreCoords(G4ThreeVector(0*m, 0*m,  -radius*nm));

//        G4SPSAngDistribution* angular = CircleSource->GetAngDist();
//        angular->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.0));
//        angular->DefineAngRefAxes("angref1", G4ThreeVector(-1.0,0.0,0.0));
//        CircleSource->SetNumberOfParticles(1);

        //Psudo Gamma
        CircleSourceG = new G4SingleParticleSource;
        CircleSourceG->SetParticleDefinition(G4Gamma::Definition());//G4Gamma G4Electron
        ;

        // Mono Energy
        G4SPSEneDistribution* energyG = CircleSourceG->GetEneDist();
        energyG->SetMonoEnergy(100*keV);
        // Energy Distribution
//        this->setEnergyToGamma()


//        /////X-Z////Vertical
//        G4SPSPosDistribution* positionG = CircleSourceG->GetPosDist();
//        positionG->SetPosDisType("Plane");
//        positionG->SetPosDisShape("Square");
//        positionG->SetHalfX(WORLD_SIDE/2);
//        positionG->SetHalfY(WORLD_SIDE/2);
//        positionG->SetHalfZ(WORLD_SIDE/2);
//        positionG->SetCentreCoords(G4ThreeVector(0*m, WORLD_SIDE/2,  0*um));
//        positionG->SetPosRot1(G4ThreeVector(1, 0, 0));
//        positionG->SetPosRot2(G4ThreeVector(0, 0, 1));
//        G4SPSAngDistribution* angularG = CircleSourceG->GetAngDist();
//        angularG->SetParticleMomentumDirection(G4ThreeVector(0.0,-1.0,0.0));
//        angularG->DefineAngRefAxes("angref1", G4ThreeVector(0.0,1.0,0.0));
//        CircleSourceG->SetNumberOfParticles(1);

        /////X-Y////Paralell
        G4SPSPosDistribution* positionG = CircleSourceG->GetPosDist();
        positionG->SetPosDisType("Plane");
        positionG->SetPosDisShape("Square");
        positionG->SetHalfX(WORLD_SIDE/2);
        positionG->SetHalfY(WORLD_SIDE/2);
        positionG->SetHalfZ(WORLD_SIDE/2);
        positionG->SetCentreCoords(G4ThreeVector(0*m, 0*um,  WORLD_SIDE/2));
        positionG->SetPosRot1(G4ThreeVector(1, 0, 0));
        positionG->SetPosRot2(G4ThreeVector(0, 1, 0));
        G4SPSAngDistribution* angularG = CircleSourceG->GetAngDist();
        angularG->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,-1.0));
        angularG->DefineAngRefAxes("angref1", G4ThreeVector(0.0,0.0,1.0));
        CircleSourceG->SetNumberOfParticles(1);



//        /////Y-Z////
//        G4SPSPosDistribution* positionG = CircleSourceG->GetPosDist();
//        positionG->SetPosDisType("Plane");
//        positionG->SetPosDisShape("Square");
//        positionG->SetHalfX(WORLD_SIDE/2);
//        positionG->SetHalfY(WORLD_SIDE/2);
//        positionG->SetHalfZ(WORLD_SIDE/2);
//        positionG->SetCentreCoords(G4ThreeVector(0*m, 0,  0*um));
//        positionG->SetPosRot1(G4ThreeVector(0, 1, 0));
//        positionG->SetPosRot2(G4ThreeVector(0, 0, 1));
//        G4SPSAngDistribution* angularG = CircleSourceG->GetAngDist();
//        angularG->SetParticleMomentumDirection(G4ThreeVector(1.0,0.0,0.0));
//        angularG->DefineAngRefAxes("angref1", G4ThreeVector(-1.0,0.0,0.0));
//        CircleSourceG->SetNumberOfParticles(1);

}

BGMSCPrimaryGeneratorAction::~BGMSCPrimaryGeneratorAction()
{
    //delete CircleSource;
    delete CircleSourceG;
    //delete fParticleGun;
}

void BGMSCPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    //CircleSource->SetNumberOfParticles(3);

        // CircleSource->GeneratePrimaryVertex(event);   //Changed
       CircleSourceG->GeneratePrimaryVertex(event);
        //fParticleGun->GeneratePrimaryVertex(event);
}


//G4SPSEneDistribution* BGMSCPrimaryGeneratorAction::setEnergyToBeta() {
//    G4SPSEneDistribution* eneDist = CircleSource->GetEneDist();
//    eneDist->SetEnergyDisType("Arb");
//    eneDist->ArbEnergyHistoFile("Electron_EnergyDistribution_20nm_100keVX-ray_allon_2115_10^6.txt");
//    eneDist->ArbInterpolate("Lin");
//    return eneDist;
//}

//G4SPSEneDistribution* BGMSCPrimaryGeneratorAction::setEnergyToGamma() {
//    G4SPSEneDistribution* eneDistG = CircleSourceG->GetEneDist();
//    eneDistG->SetEnergyDisType("Arb");
//    //eneDistG->ArbEnergyHistoFile("Gamma_EnergyDistribution_allon_10nm_100keVX-ray_2.115_10^9.txt");
//    eneDistG->ArbEnergyHistoFile("Electron_EnergyDistribution_10nm_100keVX-ray.txt");
//    eneDistG->ArbInterpolate("Lin");
//    return eneDistG;
//}




