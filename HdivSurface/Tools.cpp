//
//  Tools.cpp
//  ErrorEstimation
//
//  Created by Denise De Siqueira on 28/03/19.
//

#include "Tools.h"
#include "pzgengrid.h"

TPZCompMesh *CMeshH1(const ProblemConfig &problem) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = 0;
    for (auto matid : problem.materialids) {
        TPZMatPoisson3d *mix = new TPZMatPoisson3d(matid, cmesh->Dimension());
        mix->SetForcingFunctionExact(problem.exact.Exact());
        mix->SetForcingFunction(problem.exact.ForcingFunction());
        
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
        
    }
    
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        val2.Zero();
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        
        cmesh->InsertMaterialObject(bc);
    }
    
    cmesh->SetDefaultOrder(problem.porder);//ordem
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    
    

    
    
    return cmesh;
}

TPZCompMesh *CreatePressureMesh(const ProblemConfig &problem) {
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = 0;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    cmesh->SetDefaultOrder(problem.porder+problem.hdivmais);//ordem + hdivmais
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    int64_t n_connects = cmesh->NConnects();
    for (int64_t i = 0; i < n_connects; ++i) {
        cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
    }
    
   
    
    return cmesh;
}

TPZCompMesh *CreateFluxHDivMesh(const ProblemConfig &problem) {
    int dim = problem.gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    problem.gmesh->ResetReference();
    for (auto matid : problem.materialids) {
        TPZVecL2 *mix = new TPZVecL2(matid);
        mix->SetDimension(dim);
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 1.);
        TPZBndCond *bc = mat->CreateBC(mat, matid, 0, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());
        cmesh->InsertMaterialObject(bc);
    }
    cmesh->SetDefaultOrder(problem.porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();
    if (problem.hdivmais) {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == dim) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, problem.porder + problem.hdivmais);//seta ordem +hdivmais
                intel->SetPreferredOrder(problem.porder+problem.hdivmais);
            }
        }
 
    }
   
    
    cmesh->InitializeBlock();
    return cmesh;
    
}

TPZMultiphysicsCompMesh *CreateHDivMesh(const ProblemConfig &problem) {
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(problem.gmesh);
    TPZMaterial *mat = NULL;
    for (auto matid : problem.materialids) {
        TPZMixedPoisson *mix = new TPZMixedPoisson(matid, cmesh->Dimension());
        mix->SetForcingFunction(problem.exact.ForcingFunction());
        mix->SetForcingFunctionExact(problem.exact.Exact());
        
        if (!mat) mat = mix;
        cmesh->InsertMaterialObject(mix);
    }
    for (auto matid : problem.bcmaterialids) {
        TPZFNMatrix<1, REAL> val1(1, 1, 0.), val2(1, 1, 0.);
        int bctype = 0;
        val2.Zero();
        TPZBndCond *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
        bc->TPZMaterial::SetForcingFunction(problem.exact.Exact());

        cmesh->InsertMaterialObject(bc);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<int> active(2,1);
    TPZManVector<TPZCompMesh *> meshvector(2,0);
    
    meshvector[0] = CreateFluxHDivMesh(problem);
    meshvector[1] = CreatePressureMesh(problem);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    bool keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true, keepmatrix);
    
    return cmesh;
}


/// Increase the approximation orders of the sides of the flux elements


void UniformRefinement(int nDiv, TPZGeoMesh *gmesh) {
    
    TPZManVector<TPZGeoEl*> children;
    for(int division = 0; division < nDiv; division++) {
        
        int64_t nels = gmesh->NElements();
        
        for(int64_t elem = 0; elem < nels; elem++) {
            
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            
            if(!gel || gel->HasSubElement()) continue;
            if(gel->Dimension() == 0) continue;
            gel->Divide(children);
        }
    }
}


TPZGeoMesh *CreateGeoMesh(int nel) {
    
    TPZManVector<int> nx(2,nel);
    TPZManVector<REAL> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    TPZGenGrid gen(nx,x0,x1);

    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -1);
    gen.SetBC(gmesh, 6, -1);
    gen.SetBC(gmesh, 7, -1);

    
    

    return gmesh;
}


void MultiPhysicsCompel(const ProblemConfig &config){
    
    TPZManVector<TPZCompMesh *,2> MeshesHDiv(2);
    TPZMultiphysicsCompMesh * mixed_cmesh = CreateHDivMesh(config);
    MeshesHDiv = mixed_cmesh->MeshVector();
    
    TPZMultiphysicsCompMesh *mphysicCompMesh = new TPZMultiphysicsCompMesh(config.gmesh);
    std::ofstream outgeo("geometria.txt");
    mphysicCompMesh->Reference()->Print(outgeo);
    
    
    //Have to include the materials. Here we just did a copy of previous materials
    TPZCompMesh * cmesh =  dynamic_cast<TPZCompMesh *>(mphysicCompMesh);
    mixed_cmesh->CopyMaterials(*cmesh);
    
    TPZManVector<TPZCompMesh *,3> mp_meshes_vec(3);
    mp_meshes_vec[0] = mixed_cmesh;
    mp_meshes_vec[1] = MeshesHDiv[0];
    mp_meshes_vec[2] = MeshesHDiv[1];
    
    mphysicCompMesh->SetDimModel(2);
    TPZManVector<int,5>  active_approx_spaces(3,1);//teste usando todos os espaços
    mphysicCompMesh->BuildMultiphysicsSpace( active_approx_spaces, mp_meshes_vec);
    
    {
        std::ofstream out("mixed.txt");
        mphysicCompMesh->MeshVector()[0]->Print(out);
        
        std::ofstream out2("hdiv.txt");
        mphysicCompMesh->MeshVector()[1]->Print(out2);
        
        std::ofstream out3("L2.txt");
        mphysicCompMesh->MeshVector()[2]->Print(out3);
        
    }
    
    
    
}


void RandomRefine(ProblemConfig &config,int numelrefine){
    
    int64_t nel = config.gmesh->NElements();
    if (numelrefine > nel/2) {
        numelrefine = 1;
    }
    int count = 0;
    while(count < numelrefine) {
        int64_t elindex = random()%nel;
        TPZGeoEl *gel = config.gmesh->Element(elindex);
        if(gel && gel->Dimension() == config.gmesh->Dimension() && !gel->Father())
        {
            TPZStack<TPZGeoEl *> subels;
            gel->Divide(subels);
            count++;
        }
    }
    nel = config.gmesh->NElements();
    bool changed = true;
    while(changed)
    {
        changed = false;
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = config.gmesh->Element(el);
            if(gel && gel->Dimension() < config.gmesh->Dimension())
            {
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->HasSubElement() != 0 && !gel->HasSubElement()) {
                        TPZStack<TPZGeoEl *> subels;
                        gel->Divide(subels);
                        changed = true;
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
        }
    }
    
}

 
void Print(const FADREAL &a, std::ostream &out)
{
    out << " val " << a.val() << std::endl;
    for (int i=0; i< a.dx().size(); i++) {
        out << a.d(i) << " ";
    }
    out << std::endl;
}
void Print(const FADFADREAL &a, std::ostream &out)
{
    out << "Value ";
    Print(a.val(),out);
    out << "Derivatives\n";
    for (int i=0; i< a.dx().size(); i++) {
        Print(a.d(i),out);
    }
    out << "End\n";
    
}

void PrintSolAndDerivate(const ProblemConfig config){
    
    TPZManVector<REAL,3> x(3,0.25);
    
    TPZManVector<Fad<REAL>,3> xfad(x.size()), graduxy(x.size());
    TPZManVector<FADFADREAL,3> xfadfad(x.size()), uxyfadfad(1);
    for(int i=0; i<3; i++)
    {
        xfad[i] = Fad<REAL>(3,i,x[i]);
        xfadfad[i] = FADFADREAL(3,i,xfad[i]);
        for(int j=0; j<3; j++)
        {
            xfadfad[i].fastAccessDx(j) = Fad<REAL>(3,xfadfad[i].val().dx(j));
        }
    }
    std::cout << "xfadfad = \n";
    for(int i=0; i<3; i++)
    {
        Print(xfadfad[i],std::cout);
    }
    std::cout << std::endl;
    config.exact.graduxy(xfad, graduxy);
    config.exact.uxy(xfadfad, uxyfadfad);
    for(int i=0; i<3; i++)
    {
        std::cout << "xfad = ";
        Print(xfad[i],std::cout);
        std::cout << std::endl;
    }
    std::cout << "graduxy = \n";
    for(int i=0; i<3; i++)
    {
        Print(graduxy[i],std::cout);
    }
    std::cout << std::endl;
    std::cout << "uxyfadfad = \n";
    for(int i=0; i<uxyfadfad.size(); i++)
    {
        Print(uxyfadfad[i],std::cout);
    }
    REAL laplace = uxyfadfad[0].dx(0).dx(0)+uxyfadfad[0].dx(1).dx(1)+uxyfadfad[0].dx(2).dx(2);
    std::cout << "Laplacian " << laplace << std::endl;
    }


void FunctionTest(){
        TLaplaceExample1 Denise;
    Denise.fExact = TLaplaceExample1::ESinSinDirNonHom;//TLaplaceExample1::ESinMark;//
        TPZVec<FADFADREAL> x(3);
        FADFADREAL x0 = (FADFADREAL) 0.001;
        FADFADREAL x1 = (FADFADREAL) 0.5;
        FADFADREAL x2 = (FADFADREAL) 0;
        x[0]= x0;
        x[1]= x1;
        x[2]= x2;
        TPZVec<FADFADREAL> disp(1);
        Denise.uxy(x, disp);
    std::cout<< "Pto x[0] "<<x[0]<<std::endl;
    std::cout<< "Pto x[1] "<<x[1]<<std::endl;
    std::cout<< "Pto x[2] "<<x[2]<<std::endl;
    
    std::cout<<"valor de ur0 "<<disp[0]<<std::endl;
    
        TPZVec<REAL> x_r(3);
        x_r[0] = x[0].val().val();
        x_r[1] = x[1].val().val();
        x_r[2] = x[2].val().val();
        TPZManVector<REAL,3> grad(3);
        Denise.graduxy(x_r, grad);
    
        REAL force;
        Denise.DivSigma(x_r, force);
    
}


void Prefinamento(TPZCompMesh * cmesh, int ndiv, int porder){
    if(ndiv<1) return;
    int nel = cmesh->NElements();
    for(int iel = 0; iel < nel; iel++){
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if(!cel) continue;
        
        TPZInterpolationSpace *sp = dynamic_cast<TPZInterpolationSpace *>(cel);
        if(!sp) continue;
        int level = sp->Reference()->Level();
        TPZGeoEl * gel = sp->Reference();
        if((gel->Dimension()==2) && (iel % 2==0)){
            int ordem= 0;
            ordem=porder + (ndiv-1 ) + (level);
            std::cout<<"level "<< level<<" ordem "<<ordem<<std::endl;
            sp->PRefine(ordem);
        }
    }
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
        std::stringstream sout;
        sout<<"malha computacional apos pRefinamento\n";
        cmesh->Print(sout);

}


void SolveMixedProblem(TPZCompMesh *cmesh_HDiv,const ProblemConfig &config)
{
    
    TPZAnalysis an(cmesh_HDiv);
    
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh_HDiv);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmesh_HDiv);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema misto ate aqui
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");
    
    int dim = config.gmesh->Dimension();
    
    std::stringstream sout;
    
    sout << "MixedSolution_Order_"<<config.problemname<<"Order"<< config.porder<<"Nref_"<<config.ndivisions<<".vtk";
    
    //an.DefineGraphMesh(2, scalnames, vecnames, "Original_Misto.vtk");
    an.DefineGraphMesh(dim, scalnames, vecnames, sout.str());
    int resolution=2;
    an.PostProcess(resolution,dim);
}

void SolveH1Problem(TPZCompMesh *cmeshH1,const ProblemConfig &config){
    
    TPZAnalysis an(cmeshH1);
    
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        strmat.SetDecomposeType(ELDLt);
    an.SetStructuralMatrix(strmat);
#else
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    an.Assemble();
    an.Solve();//resolve o problema misto ate aqui
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Solution");
    vecnames.Push("Derivative");
    scalnames.Push("ExactSolution");
    
    int dim = cmeshH1->Reference()->Dimension();
    
    std::stringstream sout;
    
    sout << "H1Solution_"<<config.problemname<<"Order"<< config.porder<<"Nref_"<<config.ndivisions<<".vtk";

    an.DefineGraphMesh(dim, scalnames, vecnames, sout.str());
    int resolution=2;
    an.PostProcess(resolution,dim);
    
    
    
}


//void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out, int p, int ndiv)
//{
//    int64_t nel = hdivmesh->NElements();
//    int dim = hdivmesh->Dimension();
//    TPZManVector<REAL,10> globalerrors(10,0.);
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = hdivmesh->ElementVec()[el];
//        if(cel->Reference()->Dimension()!=dim) continue;
//        TPZManVector<REAL,10> elerror(10,0.);
//        elerror.Fill(0.);
//        cel->EvaluateError(SolExata, elerror, NULL);
//        int nerr = elerror.size();
//        for (int i=0; i<nerr; i++) {
//            globalerrors[i] += elerror[i]*elerror[i];
//        }
//
//    }
//    out << "Errors associated with HDiv space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
//    //    out << "L2 Norm for flux - "<< endl;
//    //    out<< "L2 Norm for divergence - Hdiv Norm for flux " << endl;
//    //    out <<  setw(16) << sqrt(globalerrors[1]) <<endl;
//    //    out<< setw(25)  << sqrt(globalerrors[2]) << setw(21)  << sqrt(globalerrors[3]) << endl;
//
//    out << "L2 Norm for flux = "    << sqrt(globalerrors[1]) << endl;
//    out << "L2 Norm for divergence = "    << sqrt(globalerrors[2])  <<endl;
//    out << "Hdiv Norm for flux = "    << sqrt(globalerrors[3])  <<endl;
//
//}
//
//void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out, int p, int ndiv)
//{
//    int64_t nel = l2mesh->NElements();
//    //int dim = l2mesh->Dimension();
//    TPZManVector<REAL,10> globalerrors(10,0.);
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = l2mesh->ElementVec()[el];
//        TPZManVector<REAL,10> elerror(10,0.);
//        cel->EvaluateError(SolExata, elerror, NULL);
//        int nerr = elerror.size();
//        globalerrors.resize(nerr);
//        for (int i=0; i<nerr; i++) {
//            globalerrors[i] += elerror[i]*elerror[i];
//        }
//
//    }
//    out << "Errors associated with L2 space - ordem polinomial = " << p << "- divisoes = " << ndiv << endl;
//    out << "L2 Norm for pressure = "    << sqrt(globalerrors[1]) << endl;
//}
