/**
 * @file Poisson 3D in hexahedra with shock problem
 */
#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "pzintel.h"


#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"


#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzmatmixedpoisson3d.h"
#include <tuple>
#include <memory>
#include "Tools.h"



bool gmeshreader=true;
bool isH1=true;
int matId=1;

bool datatest=false;

int bc0 = -1;
int bc1 = -2;
int bc2 = -3;
int bc3 = -4;
int bc4 = -5;
int bc5 = -6;

//#define DUALMat

#define UnitSquare

int dirichlet = 0;
int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif

// Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    
    TPZGeoMesh *gmesh = nullptr;
    
    ProblemConfig config;
    config.porder = 1;
   // config.problemname = "LaplaceOnSphere";
    config.dimension = 2;
    config.ndivisions = 0;
    TLaplaceExample1 example;
    config.exact.fExact = example.EConst;//ESinSinCircle;//EBubble;////ESinSin;//ESinSinDirNonHom;//;//ESinMark;//EArcTanSingular;//EArcTan;//
    
    int dim = 3;//somente para ler o arquivo
    
    if(gmeshreader){

      //  std::string meshfilename = "../Mysquare.msh";
       //  std::string meshfilename = "../Quad.msh";
        std::string meshfilename = "../Circle.msh";
        
          TPZGmshReader gmsh;
        
        if(dim==3)
        {
            meshfilename = "../esfera2.msh";
            //meshfilename = "../Mysphere.msh";
           // meshfilename = "../Cube.msh";
           // gmsh.GetDimNamePhysical()[2]["dirichlet"] =2;
            gmsh.GetDimNamePhysical()[1]["boundary"] =2;
            gmsh.GetDimNamePhysical()[2]["domain"] = 1;
            

            
        }
        else{
            
            gmsh.GetDimNamePhysical()[1]["dirichlet"] =2;
            gmsh.GetDimNamePhysical()[2]["domain"] = 1;
        }
      
     
        config.materialids.insert(1);
        config.bcmaterialids.insert(2);
        
        
        gmsh.SetFormatVersion("4.1");
        gmesh = gmsh.GeometricGmshMesh(meshfilename);
        gmsh.PrintPartitionSummary(std::cout);
        gmesh->SetDimension(config.dimension);
        config.gmesh = gmesh;
        
        
    }
    else{
       
        
#ifdef UnitSquare
        gmesh= CreateGeoMesh(2);
#else
        gmesh= GMeshSphere2();
#endif
    
    
        config.gmesh = gmesh;
        config.materialids.insert(1);
        config.bcmaterialids.insert(bc1);
        config.bcmaterialids.insert(bc2);
         config.bcmaterialids.insert(bc3);
         config.bcmaterialids.insert(bc4);
    
    }
    
    UniformRefinement(config.ndivisions, gmesh);
    

    {
        std::ofstream out("Sgmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out,true);
        ofstream arg1("GeoMesh.txt");
        gmesh->Print(arg1);
    }

    

    if (isH1) {
        example.fSignConvention = -1;
        TPZCompMesh *cmeshH1 = CMeshH1(config);
        {
        ofstream arg1("CompMeshH1.txt");
            cmeshH1->Print(arg1);
        }
        
        SolveH1Problem(cmeshH1,config);
        
       // return EXIT_SUCCESS;
    }
    
    //REsolvendo problema misto

    
 
    TPZMultiphysicsCompMesh *mixedmesh=CreateHDivMesh(config);
  {

        ofstream file3("MultPhysics.txt");
        mixedmesh->Print(file3);

    }


    SolveMixedProblem(mixedmesh,config);
    
    //Calculo de erro
    
    
//  std::ofstream Errorfile("Error.txt");
//    ErrorHDiv(fluxmesh, Errorfile, config.porder, ndiv);
//    ErrorL2(pressuremesh, Errorfile, config.porder, ndiv);
    
    

        
        
    return 0;
}





void SolveMixSystem(TPZCompMesh *fCmesh)
{

    TPZAnalysis an(fCmesh);
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(fCmesh);
    strmat.SetNumThreads(2);
    //        strmat.SetDecomposeType(ELDLt);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(fCmesh);
    strmat.SetNumThreads(2);

#endif
    
    int dim=fCmesh->Reference()->Dimension();
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
//    delete direct;
//    direct = 0;
    an.Assemble();
    an.Solve();
//
//    fCmesh->LoadSolution(an.Solution());
//    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, fCmesh);
    
//    meshvector[0]->Solution().Print("p = ",std::cout,EMathematicaInput);
    
    //Pos processamento
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    scalnames.Push("ExactPressure");
    vecnames.Push("ExactFlux");
    scalnames.Push("Divergence");
    scalnames.Push("ExactDiv");
    an.DefineGraphMesh(dim, scalnames, vecnames, "MixedSolution.vtk");
    // Post processing
    an.PostProcess(2,dim);

}



TPZGeoMesh *CreateGeoMesh() {
    
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    int matID=1;
    
    
    // Creates matrix with quadrilateral node coordinates.
    const int quadNodeNumber = 4;
    REAL coordinates[quadNodeNumber][3] = {
        {0., 0., 1.},
        {1., 0., 1.},
        {1., 1., 1.},
        {0., 1., 1.}
    };
    
    // Inserts coordinates in the TPZGeoMesh object.
    for(int i = 0; i < quadNodeNumber; i++) {
        int64_t nodeID = gmesh->NodeVec().AllocateNewElement();
        
        TPZVec<REAL> nodeCoord(3);
        nodeCoord[0] = coordinates[i][0];
        nodeCoord[1] = coordinates[i][1];
        nodeCoord[2] = coordinates[i][2];
        
        gmesh->NodeVec()[nodeID] = TPZGeoNode(i, nodeCoord, *gmesh);
    }
    
    gmesh->SetDimension(2);
    // Creates quadrilateral element.
    int64_t index;
    TPZManVector<int64_t> nodeIDs(quadNodeNumber);
    
    for(int n = 0; n < quadNodeNumber; n++) {
        nodeIDs[n] = n;
    }
    gmesh->CreateGeoElement(EQuadrilateral, nodeIDs, matID, index);
    
    // Creates line elements where boundary conditions will be inserted.
    nodeIDs.Resize(2);
    for (int i = 0; i < 4; i++) {
        
        nodeIDs[0] = i % 4;
        nodeIDs[1] = (i + 1) % 4;
        
        gmesh->CreateGeoElement(EOned, nodeIDs, -1, index);
    }
    
    gmesh->BuildConnectivity();
    return gmesh;
}

