//
//  Tools.h
//  ErrorEstimation
//
//  Created by Denise De Siqueira on 28/03/19.
//

#ifndef Tools_h
#define Tools_h

#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "ProblemConfig.h"

#include "mixedpoisson.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"

#include "pzintel.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZHybridizeHDiv.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZMultiphysicsCompMesh.h"


class TPZMultiphysicsCompMesh;

#include <tuple>
#include <memory>

#include <stdio.h>

#endif /* Tools_hpp */
TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem);
TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem);
TPZMultiphysicsCompMesh *CreateHDivMesh(const ProblemConfig &problem);

void UniformRefinement(int nDiv, TPZGeoMesh *gmesh);
TPZGeoMesh *CreateGeoMesh(int nelem);

void RandomRefine(ProblemConfig &config,int numelrefine);


void PrintSolAndDerivate(const ProblemConfig config);
void FunctionTest();

void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder);



void SolveHybridProblem(TPZCompMesh *Hybridmesh,int InterfaceMatId,const ProblemConfig &problem);
void SolveMixedProblem(TPZCompMesh *cmesh_HDiv,const ProblemConfig &config);
void PlotLagrangreMultiplier(TPZCompMesh *cmesh, const ProblemConfig &problem);

TPZCompMesh * CMeshH1(const ProblemConfig &problem);
void SolveH1Problem(TPZCompMesh *cmeshH1,const ProblemConfig &config);
void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv);
void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv);



