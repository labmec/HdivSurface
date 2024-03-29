//
//  ProblemConfig.h
//  ErrorEstimation
//
//  Created by Philippe Devloo on 09/05/18.
//

#ifndef ProblemConfig_h
#define ProblemConfig_h

#include <set>
#include "TPZAnalyticSolution.h"

/// class to guide the error estimator
struct ProblemConfig
{
    
    /// geometric mesh on which the computational meshes are based
    TPZGeoMesh *gmesh = 0;
    /// polynomial order of the original mesh
    int porder = 1;
    /// increment in internal order of flux and pressure
    int hdivmais = 1;
    /// option to compute the error based on continuous pressures or not
    bool makepressurecontinuous = 0;
    
    /// number of uniform refinements applied to the mesh
    int ndivisions=1;
    
    bool prefine=false;
    bool steklovexample=false;
    bool GalvisExample;
    STATE alpha=1;
    /// directory where the files will be stored
    std::string dir_name;
    /// name identifying the problem
    std::string problemname;
    /// set of materialids in the mesh
    std::set<int> materialids;
    /// set of boundary condition material ids
    std::set<int> bcmaterialids;
    /// exact solution
    TLaplaceExample1 exact;
    TLaplaceExample1 forcingCte;
    int dimension=2;
    
    
    ProblemConfig(){};
    ProblemConfig(const ProblemConfig &cp) : gmesh(cp.gmesh), porder(cp.porder), hdivmais(cp.hdivmais),
    makepressurecontinuous(cp.makepressurecontinuous),
    problemname(cp.problemname),
    materialids(cp.materialids), bcmaterialids(cp.bcmaterialids),exact(cp.exact),forcingCte(cp.forcingCte),dimension(cp.dimension),ndivisions(cp.ndivisions)

    {
    }
    
    ProblemConfig &operator=(const ProblemConfig &cp)
    {
        gmesh = cp.gmesh;
        porder = cp.porder;
        hdivmais = cp.hdivmais;
        makepressurecontinuous = cp.makepressurecontinuous;
        problemname = cp.problemname;
        materialids = cp.materialids;
        bcmaterialids = cp.bcmaterialids;
        exact = cp.exact;

        forcingCte = cp.forcingCte;
        dimension = cp.dimension;
        ndivisions = cp.ndivisions;
        prefine = cp.prefine;
        alpha = cp.alpha;
        dir_name = cp.dir_name;
        steklovexample=cp.steklovexample;
        GalvisExample=cp.GalvisExample;
        

        return *this;
    }
};

#endif /* ProblemConfig_h */
