#include "BGMSCDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4VSolid.hh"
#include "G4Sphere.hh"
#include "Randomize.hh"
#include "G4UserLimits.hh"

#include "math.h"
using namespace CLHEP;

/*
10mg 3106
20mg 6211
30mg 9317
40mg 12422
50mg 15528
*/
///////////////////////////////////////////////////////////////////
#define WORLD_SIDE (25. * um)
#define VESSEL_OUTER_DIAM (20.0 * um)
#define VESSEL_INNER_DIAM (16.0 * um)
#define VESSEL_HIGHT (10.0 * um)
#define GNP_DIAM (100 * nm) // diameter!!
#define GNP_COUNT 15528 // Number!
#define CELL_COUNT 60
//#define VESSEL_SPACE (5. * um)
//#define VESSEL_SHELL 5
//#define SHELL1_DIAM (5.0 * um)
//#define SHELL2_DIAM (10.0 * um)
//#define SHELL3_DIAM (15.0 * um)
//#define SHELL4_DIAM (20.0 * um)
//#define SHELL5_DIAM (25.0 * um)
//#define SHELL_FRC20NM_0_10 0.4332
//#define SHELL_FRC20NM_10_20 0.2594
//#define SHELL_FRC20NM_20_30 0.1431
//#define SHELL_FRC20NM_30_40 0.08949
//#define SHELL_FRC20NM_40_50 0.07479
///////////////////////////////////////////////////////////////////

BGMSCDetectorConstruction::BGMSCDetectorConstruction()
    :fStepLimit(NULL)
{
    m_dWorldSide = WORLD_SIDE;
    m_nGnpCount = GNP_COUNT;
    m_strDistribution = "Constrained";   // Constrained  Random None 8HPI
    m_nCellCount = CELL_COUNT;
}

G4VPhysicalVolume* BGMSCDetectorConstruction::Construct()
{
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    G4VisAttributes* pVisAttributesEndotherial = new G4VisAttributes;
    pVisAttributesEndotherial->SetForceWireframe(true);
    pVisAttributesEndotherial->SetForceAuxEdgeVisible(true);
    pVisAttributesEndotherial->SetForceSolid(false);
    pVisAttributesEndotherial->SetVisibility(true);
    pVisAttributesEndotherial->SetColor(1.0, 0.0, 0.0); // red

    // Build materials
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* Au = nistManager->FindOrBuildMaterial("G4_Au");
    G4Material* vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
    G4Material* water = nistManager->FindOrBuildMaterial("G4_WATER");

    // Define World
    G4Box* pWorldBox = new G4Box("WorldBox", WORLD_SIDE/2, WORLD_SIDE/2, WORLD_SIDE/2);
    G4LogicalVolume *pWorldLogic = new G4LogicalVolume(pWorldBox, water, "WorldLog");
    G4VPhysicalVolume *pWorldPhys = new G4PVPlacement(0, G4ThreeVector(), pWorldLogic, "WorldPhys", 0, false, 0);

//    // Define Endotherial cell
//    G4Tubs* pEndotherialTubs = new G4Tubs("EndotherialTubs", VESSEL_INNER_DIAM/2, VESSEL_OUTER_DIAM/2, VESSEL_HIGHT/2, 0*deg, 360*deg);
//    G4LogicalVolume *pEndotherialLogic = new G4LogicalVolume(pEndotherialTubs, water, "EndotherialLog");
//    G4VPhysicalVolume *pEndotherialPhys = new G4PVPlacement(0, G4ThreeVector(), pEndotherialLogic, "EndotherialPhys", pWorldLogic, 0, 0);
//    pEndotherialLogic->SetVisAttributes(pVisAttributesEndotherial);

    // Parameterize each Endotherial cell.
    for (int nCellIdx = 0; nCellIdx < m_nCellCount; nCellIdx++)
    {
        G4RotationMatrix* rotm = new G4RotationMatrix;
        rotm->rotateZ((nCellIdx+1) * 6*deg - 90*deg); //Set the first position 0 oclock., clock wise.

        G4Tubs* pEndotherialTubs = new G4Tubs("EndotherialTubs", VESSEL_INNER_DIAM/2, VESSEL_OUTER_DIAM/2, VESSEL_HIGHT/2, 0*deg, 6*deg); // Last deg is the size of cell!
        G4LogicalVolume *pEndotherialLogic = new G4LogicalVolume(pEndotherialTubs, water, "EndotherialLog");
        G4VPhysicalVolume *EndotherialPhys = new G4PVPlacement(rotm, G4ThreeVector(), pEndotherialLogic, "EndotherialPhys", pWorldLogic, 0, nCellIdx);
        pEndotherialLogic->SetVisAttributes(pVisAttributesEndotherial);
    }

    if (m_strDistribution == "Constrained")
    {
        DistributeGnpsSurface(pWorldLogic);
    }
    else if (m_strDistribution == "Random")
    {
        DistributeGnpsRandom(pWorldLogic);
    }
    else if (m_strDistribution == "8HPI")
    {
        //DistributeGnps8PHI(pWorldLogic);
    }
    else if (m_strDistribution == "None")
    {
    }
    else
    {
        G4Exception("Invalid distribution name", "001", G4ExceptionSeverity::FatalException, "Comment");
    }

    G4double maxStep = 1*nm;
    fStepLimit = new G4UserLimits(maxStep);

    return pWorldPhys;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
///////////////////////////CONSTRAINED////////////////////////////////////////
// Distribute GNP randomly on the inner surface of the blood vessel.
void BGMSCDetectorConstruction::DistributeGnpsSurface(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributes = new G4VisAttributes;
    pVisAttributes->SetForceWireframe(true);
    pVisAttributes->SetForceAuxEdgeVisible(true);
    pVisAttributes->SetForceSolid(false);
    pVisAttributes->SetVisibility(true);
    pVisAttributes->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...");
    G4int ary[360]={0};

    for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
    {
        retry:
        // Compute a random position for the GNP
        double dTheta = G4UniformRand() * 360.*deg;
        double dGnpX = ((VESSEL_INNER_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
        double dGnpY = ((VESSEL_INNER_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
        double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Special case for one GNP at the center.
                if (m_nGnpCount == 1)
                {
                    dGnpX = 0.0;
                    dGnpY = 0.0;
                    dGnpZ = 0.0;
                }

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry;
                }
                printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

                G4double r = VESSEL_INNER_DIAM/2.;
                G4double Cos = cos(dGnpX/r);
                G4double theta = acos(Cos)*180/M_PI;
//                if (dGnpX<0 && dGnpY>0) theta = theta+90; //2nd
//                else if(dGnpX<0 && dGnpY<0) theta = theta+180; //3rd
//                else if(dGnpX>0 && dGnpY<0) theta = theta+270; //4th

                for (int i=0; i<60; i++)
                {
                    if((i<=theta) && (theta<(i+1)))
                    {
                        ary[i] += 1;
                    }
                }

                int nCopyNumber = nGnpIdx; //Why?

                new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);

                pGnpLog->SetVisAttributes(pVisAttributes);

                aryGnpInfo[nGnpIdx].dPosX = dGnpX;
                aryGnpInfo[nGnpIdx].dPosY = dGnpY;
                aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
    }
    FILE* fp =fopen("position.txt", "wt");
    for (int i=0; i<60; i++)
    {
        printf("%d %d\n", i, ary[i]);
        fprintf(fp, "%d %d\n", i, ary[i]);
    }
    fclose(fp);
    delete [] aryGnpInfo;
}

///////////////////////////RANDOM////////////////////////////////////////
// Distribute GNP randomly in the blood vessel.
void BGMSCDetectorConstruction::DistributeGnpsRandom(G4LogicalVolume *pWorldLogic)
{
    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
    assert(pMaterialWater != NULL);
    assert(pMaterialGold != NULL);

    struct SGnpInfo
    {
        double dPosX;
        double dPosY;
        double dPosZ;
    };

    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

    // Set visible attributes
    G4VisAttributes* pVisAttributes = new G4VisAttributes;
    pVisAttributes->SetForceWireframe(true);
    pVisAttributes->SetForceAuxEdgeVisible(true);
    pVisAttributes->SetForceSolid(false);
    pVisAttributes->SetVisibility(true);
    pVisAttributes->SetColor(255. / 255., 215. / 255., 0.); // gold

    // Create the Sphere object
    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

    printf("Distributing GNPs randomly...");

    for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
    {
        retry:
        // Compute a random position for the GNP
        double dTheta = G4UniformRand() * 2*M_PI;
        double dRand = G4UniformRand();
        double dGnpX = (sqrt(dRand) * (VESSEL_INNER_DIAM/2.-GNP_DIAM/2) * cos(dTheta));
        double dGnpY = (sqrt(dRand) * (VESSEL_INNER_DIAM/2.-GNP_DIAM/2) * sin(dTheta));
        double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

                // Special case for one GNP at the center.
                if (m_nGnpCount == 1)
                {
                    dGnpX = 0.0;
                    dGnpY = 0.0;
                    dGnpZ = 0.0;
                }

                // Check if this GNP is over-lapping an existing GNP
                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
                {
                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

                    double dx = dExistingX - dGnpX;
                    double dy = dExistingY - dGnpY;
                    double dz = dExistingZ - dGnpZ;

                    // Don't bother computing sqrt if the distance is large on any axis.
                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
                    if (fabs(dx) > GNP_DIAM) continue;
                    if (fabs(dy) > GNP_DIAM) continue;
                    if (fabs(dz) > GNP_DIAM) continue;

                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dDist <= GNP_DIAM)
                        goto retry;
                }
                printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

                int nCopyNumber = nGnpIdx; //Why?

                new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);

                pGnpLog->SetVisAttributes(pVisAttributes);

                aryGnpInfo[nGnpIdx].dPosX = dGnpX;
                aryGnpInfo[nGnpIdx].dPosY = dGnpY;
                aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
    }
    delete [] aryGnpInfo;
}

///////////////////////////8HPI////////////////////////////////////////
// Distribute GNP with migration after 8 hour post injection (8HPI)
//void BGMSCDetectorConstruction::DistributeGnps8PHI(G4LogicalVolume *pWorldLogic)
//{
//    G4VisAttributes* pVisAttributesGold = new G4VisAttributes;
//    pVisAttributesGold->SetForceWireframe(true);
//    pVisAttributesGold->SetForceAuxEdgeVisible(true);
//    pVisAttributesGold->SetForceSolid(false);
//    pVisAttributesGold->SetVisibility(true);
//    pVisAttributesGold->SetColor(255. / 255., 215. / 255., 0.); // gold

//    G4Material *pMaterialWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
//    G4Material *pMaterialGold = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
//    assert(pMaterialWater != NULL);
//    assert(pMaterialGold != NULL);

////////This is just for visualization//////
////    G4VisAttributes* pVisAttributes = new G4VisAttributes;
////    pVisAttributes->SetForceWireframe(true);
////    pVisAttributes->SetForceAuxEdgeVisible(true);
////    pVisAttributes->SetForceSolid(false);
////    pVisAttributes->SetVisibility(true);
////    pVisAttributes->SetColor(0.0, 1.0, 1.0); // cyan

////    for (int nVesselShellIdx = 0; nVesselShellIdx < VESSEL_SHELL; nVesselShellIdx++)
////    {
////        G4Tubs* pVesselShellTubs = new G4Tubs("VesselShellTubs", VESSEL_SPACE*nVesselShellIdx, VESSEL_SPACE*(nVesselShellIdx+1), VESSEL_HIGHT/2., 0*deg, 360*deg);
////        G4LogicalVolume *pVesselShellLogic = new G4LogicalVolume(pVesselShellTubs, pMaterialWater, "VesselShellLog");
////        G4VPhysicalVolume *pVesselShellPhys = new G4PVPlacement(0, G4ThreeVector(), pVesselShellLogic, "VesselShellLog", pWorldLogic, 0, 0);
////        pVesselShellLogic->SetVisAttributes(pVisAttributes);
////    }
/////////////////////////////////////////////

//    G4VisAttributes* pVisAttributesEndotherial = new G4VisAttributes;
//    pVisAttributesEndotherial->SetForceWireframe(true);
//    pVisAttributesEndotherial->SetForceAuxEdgeVisible(true);
//    pVisAttributesEndotherial->SetForceSolid(false);
//    pVisAttributesEndotherial->SetVisibility(true);
//    pVisAttributesEndotherial->SetColor(1.0, 0.0, 0.0); // red

//    // Define Endotherial cell
//    G4Tubs* pEndotherialTubs = new G4Tubs("EndotherialTubs", VESSEL_INNER_DIAM / 2, VESSEL_OUTER_DIAM / 2, VESSEL_HIGHT / 2, 0*deg, 360*deg);
//    G4LogicalVolume *pEndotherialLogic = new G4LogicalVolume(pEndotherialTubs, pMaterialWater, "EndotherialLog");
//    G4VPhysicalVolume *pEndotherialPhys = new G4PVPlacement(0, G4ThreeVector(), pEndotherialLogic, "EndotherialPhys", pWorldLogic, 0, 0);
//    pEndotherialLogic->SetVisAttributes(pVisAttributesEndotherial);


//    struct SGnpInfo
//    {
//        double dPosX;
//        double dPosY;
//        double dPosZ;
//    };

//    SGnpInfo *aryGnpInfo = new SGnpInfo [m_nGnpCount];

//    // Create the Sphere object
//    G4Sphere* pGnpSphere = new G4Sphere("GNP", 0., GNP_DIAM / 2, 0*deg, 360*deg, 0*deg, 180*deg);
//    G4LogicalVolume *pGnpLog = new G4LogicalVolume(pGnpSphere, pMaterialGold, "GNPLogic");

//    printf("Distributing GNPs randomly...");

//    if (GNP_DIAM == 20 * nm)
//    {
//        for (int nGnpIdx = 0; nGnpIdx < m_nGnpCount; nGnpIdx++)
//        {
//            G4int GNP_COUNT_0_10 = GNP_COUNT*SHELL_FRC20NM_0_10;
//            if (0 <= nGnpIdx && nGnpIdx < GNP_COUNT_0_10) // Shell1
//            {
//                retry1:
//                // Compute a random position for the GNP
//                double dTheta = G4UniformRand() * 360.*deg;
//                double dGnpX = (SHELL1_DIAM/2. * G4UniformRand() * cos(dTheta));
//                double dGnpY = (SHELL1_DIAM/2. * G4UniformRand() * sin(dTheta));
//                double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

//                // Check if this GNP is over-lapping an existing GNP
//                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
//                {
//                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
//                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
//                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

//                    double dx = dExistingX - dGnpX;
//                    double dy = dExistingY - dGnpY;
//                    double dz = dExistingZ - dGnpZ;

//                    // Don't bother computing sqrt if the distance is large on any axis.
//                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
////                    if (fabs(dx) > GNP_DIAM) continue;
////                    if (fabs(dy) > GNP_DIAM) continue;
////                    if (fabs(dz) > GNP_DIAM) continue;

//                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
//                    if (dDist <= GNP_DIAM)
//                        goto retry1;

//                    printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

//                    int nCopyNumber = nGnpIdx;
//                    new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
//                    pGnpLog->SetVisAttributes(pVisAttributesGold);

//                    aryGnpInfo[nGnpIdx].dPosX = dGnpX;
//                    aryGnpInfo[nGnpIdx].dPosY = dGnpY;
//                    aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
//                }
//            }
//            else if (GNP_COUNT * SHELL_FRC20NM_0_10 <= nGnpIdx && GNP_COUNT*SHELL_FRC20NM_10_20) // Shell2
//            {
//                retry2:
//                // Compute a random position for the GNP
//                double dTheta = G4UniformRand() * 360.*deg;
//                double dGnpX = (SHELL2_DIAM/2. * G4UniformRand() * cos(dTheta));
//                double dGnpY = (SHELL2_DIAM/2. * G4UniformRand() * sin(dTheta));
//                double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

//                // Check if this GNP is in the inner shell
//                double dWhichShell = sqrt(dGnpX*dGnpX + dGnpY*dGnpY);
//                if (dWhichShell < SHELL1_DIAM)
//                    goto retry2;

//                // Check if this GNP is over-lapping an existing GNP
//                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
//                {
//                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
//                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
//                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

//                    double dx = dExistingX - dGnpX;
//                    double dy = dExistingY - dGnpY;
//                    double dz = dExistingZ - dGnpZ;

//                    // Don't bother computing sqrt if the distance is large on any axis.
//                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
//                    if (fabs(dx) > GNP_DIAM) continue;
//                    if (fabs(dy) > GNP_DIAM) continue;
//                    if (fabs(dz) > GNP_DIAM) continue;

//                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
//                    if (dDist <= GNP_DIAM)
//                        goto retry2;

//                    printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

//                    int nCopyNumber = nGnpIdx;
//                    new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
//                    pGnpLog->SetVisAttributes(pVisAttributesGold);

//                    aryGnpInfo[nGnpIdx].dPosX = dGnpX;
//                    aryGnpInfo[nGnpIdx].dPosY = dGnpY;
//                    aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
//                }
//            }
//            else if (GNP_COUNT * SHELL_FRC20NM_10_20 <= nGnpIdx && GNP_COUNT*SHELL_FRC20NM_20_30) // Shell3
//            {
//                retry3:
//                // Compute a random position for the GNP
//                double dTheta = G4UniformRand() * 360.*deg;
//                double dGnpX = (SHELL3_DIAM/2. * G4UniformRand() * cos(dTheta));
//                double dGnpY = (SHELL3_DIAM/2. * G4UniformRand() * sin(dTheta));
//                double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

//                // Check if this GNP is in the inner shell
//                double dWhichShell = sqrt(dGnpX*dGnpX + dGnpY*dGnpY);
//                if (dWhichShell < SHELL2_DIAM)
//                    goto retry3;

//                // Check if this GNP is over-lapping an existing GNP
//                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
//                {
//                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
//                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
//                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

//                    double dx = dExistingX - dGnpX;
//                    double dy = dExistingY - dGnpY;
//                    double dz = dExistingZ - dGnpZ;

//                    // Don't bother computing sqrt if the distance is large on any axis.
//                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
//                    if (fabs(dx) > GNP_DIAM) continue;
//                    if (fabs(dy) > GNP_DIAM) continue;
//                    if (fabs(dz) > GNP_DIAM) continue;

//                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
//                    if (dDist <= GNP_DIAM)
//                        goto retry3;

//                    printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

//                    int nCopyNumber = nGnpIdx;
//                    new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
//                    pGnpLog->SetVisAttributes(pVisAttributesGold);

//                    aryGnpInfo[nGnpIdx].dPosX = dGnpX;
//                    aryGnpInfo[nGnpIdx].dPosY = dGnpY;
//                    aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
//                }
//            }
//            else if (GNP_COUNT * SHELL_FRC20NM_20_30 <= nGnpIdx && GNP_COUNT*SHELL_FRC20NM_30_40) // Shell4
//            {
//                retry4:
//                // Compute a random position for the GNP
//                double dTheta = G4UniformRand() * 360.*deg;
//                double dGnpX = (SHELL4_DIAM/2. * G4UniformRand() * cos(dTheta));
//                double dGnpY = (SHELL4_DIAM/2. * G4UniformRand() * sin(dTheta));
//                double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

//                // Check if this GNP is in the inner shell
//                double dWhichShell = sqrt(dGnpX*dGnpX + dGnpY*dGnpY);
//                if (dWhichShell < SHELL3_DIAM)
//                    goto retry4;

//                // Check if this GNP is over-lapping an existing GNP
//                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
//                {
//                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
//                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
//                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

//                    double dx = dExistingX - dGnpX;
//                    double dy = dExistingY - dGnpY;
//                    double dz = dExistingZ - dGnpZ;

//                    // Don't bother computing sqrt if the distance is large on any axis.
//                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
//                    if (fabs(dx) > GNP_DIAM) continue;
//                    if (fabs(dy) > GNP_DIAM) continue;
//                    if (fabs(dz) > GNP_DIAM) continue;

//                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
//                    if (dDist <= GNP_DIAM)
//                        goto retry4;

//                    printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

//                    int nCopyNumber = nGnpIdx;
//                    new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
//                    pGnpLog->SetVisAttributes(pVisAttributesGold);

//                    aryGnpInfo[nGnpIdx].dPosX = dGnpX;
//                    aryGnpInfo[nGnpIdx].dPosY = dGnpY;
//                    aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
//                }
//            }
//            else if (GNP_COUNT * SHELL_FRC20NM_30_40 <= nGnpIdx && GNP_COUNT*SHELL_FRC20NM_40_50) // Shell5
//            {
//                retry5:
//                // Compute a random position for the GNP
//                double dTheta = G4UniformRand() * 360.*deg;
//                double dGnpX = (SHELL5_DIAM/2. * G4UniformRand() * cos(dTheta));
//                double dGnpY = (SHELL5_DIAM/2. * G4UniformRand() * sin(dTheta));
//                double dGnpZ = (VESSEL_HIGHT * G4UniformRand()) - VESSEL_HIGHT/2.0;

//                // Check if this GNP is in the inner shell
//                double dWhichShell = sqrt(dGnpX*dGnpX + dGnpY*dGnpY);
//                if (dWhichShell < SHELL4_DIAM)
//                    goto retry5;

//                // Check if this GNP is over-lapping an existing GNP
//                for (int nExistingGnpIdx = 0; nExistingGnpIdx < nGnpIdx; nExistingGnpIdx++)
//                {
//                    double dExistingX = aryGnpInfo[nExistingGnpIdx].dPosX;
//                    double dExistingY = aryGnpInfo[nExistingGnpIdx].dPosY;
//                    double dExistingZ = aryGnpInfo[nExistingGnpIdx].dPosZ;

//                    double dx = dExistingX - dGnpX;
//                    double dy = dExistingY - dGnpY;
//                    double dz = dExistingZ - dGnpZ;

//                    // Don't bother computing sqrt if the distance is large on any axis.
//                    //if ((fabs(dx) <= GNP_DIAM) && (fabs(dy) <= GNP_DIAM) && (fabs(dz) <= GNP_DIAM))
//                    if (fabs(dx) > GNP_DIAM) continue;
//                    if (fabs(dy) > GNP_DIAM) continue;
//                    if (fabs(dz) > GNP_DIAM) continue;

//                    double dDist = sqrt(dx*dx + dy*dy + dz*dz);
//                    if (dDist <= GNP_DIAM)
//                        goto retry5;

//                    printf("%d (%le, %le, %le)\n", nGnpIdx, dGnpX, dGnpY, dGnpZ);

//                    int nCopyNumber = nGnpIdx;
//                    new G4PVPlacement(0, G4ThreeVector(dGnpX, dGnpY, dGnpZ), pGnpLog, "GnpPhys", pWorldLogic, false, nCopyNumber);
//                    pGnpLog->SetVisAttributes(pVisAttributesGold);

//                    aryGnpInfo[nGnpIdx].dPosX = dGnpX;
//                    aryGnpInfo[nGnpIdx].dPosY = dGnpY;
//                    aryGnpInfo[nGnpIdx].dPosZ = dGnpZ;
//                }
//            }
//            else
//            {
//                G4cout << "Shell index error!!" << G4endl;
//            }
//        }

//    }
//    delete [] aryGnpInfo;
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BGMSCDetectorConstruction::SetMaxStep(G4double maxStep)
{
    if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}
