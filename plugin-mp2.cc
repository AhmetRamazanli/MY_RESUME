/*
 * @BEGIN LICENSE
 *
 * plugin_test by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/lib3index/df_helper.h"
#include "psi4/lib3index/dftensor.h"
#include "tensors.h"

namespace psi{ namespace plugin_mp2 {


extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "PLUGIN_BASICS"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
	/*- Basis set -*/
	options.add_str("DF_BASIS_MP2", "");
	/*- REF -*/
	options.add_str("REFERENCE","RHF", "RHF UHF");
    }

    return true;
}

extern "C"
SharedWavefunction plugin_mp2(SharedWavefunction ref_wfn, Options& options)
{
    // Get options
    int print = options.get_int("PRINT");
    std::string reference = options.get_str("REFERENCE");
    //options.print();
    // Grab the global (default) PSIO object, for file I/O
    std::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Defs
    SharedTensor2d CmoA, FockA,FockB, Hso, Tso, Vso, Sso, bQso, bQij, bQia, bQab, bQiaA, bQiaB;
    SharedTensor1d eps_orbA, eps_orbB;

    int nQ;

    //=============== RHF =========================
    if(reference == "RHF"){

    /* Your code goes here */
    int nmo = ref_wfn->nmo(); // MO sayısı
    int nso = ref_wfn->nso(); // AO sayısı
    int noccA = ref_wfn->doccpi()[0]; // Toplam OCC sayısı
    int nfrzc = ref_wfn->frzcpi()[0]; // Frozen core sayısı (inaktif OCC) 
    int nfrzv = ref_wfn->frzvpi()[0];  // Frozen virtuals
    int nalpha = ref_wfn->nalphapi()[0]; // alpha e sayısı
    int nbeta = ref_wfn->nbetapi()[0];   // beta e sayısı
    int natom = ref_wfn->molecule()->natom(); // natom

    // figure out number of MO spaces
    int nvirA = nmo - noccA;         // Number of virtual orbitals
    int naoccA = noccA - nfrzc;       // Number of active occupied orbitals
    int navirA = nvirA - nfrzv;       // Number of active virtual orbitals
    int ntri_so = 0.5*nso*(nso+1);

    // Read in nuclear repulsion energy
    double Enuc = ref_wfn->molecule()->nuclear_repulsion_energy();

    // Read SCF energy
    double Escf = ref_wfn->reference_energy();
    double Eref = Escf;
    double Eelec = Escf - Enuc;

    outfile->Printf("\n");
    outfile->Printf("\tNMO  : %2d \n", nmo);

    // Read MO energies
    SharedVector epsilon_a_ = ref_wfn->epsilon_a();
    //epsilon_a_->print();

    eps_orbA = std::shared_ptr<Tensor1d>(new Tensor1d("epsilon <P|Q>", nmo));
    for (int p = 0; p < nmo; ++p) eps_orbA->set(p, epsilon_a_->get(0, p));

    // Build Initial fock matrix
    FockA = SharedTensor2d(new Tensor2d("MO-basis alpha Fock matrix", nmo, nmo));
    for (int i = 0; i < noccA; ++i) FockA->set(i, i, epsilon_a_->get(0, i));
    for (int a = 0; a < nvirA; ++a) FockA->set(a + noccA, a + noccA, epsilon_a_->get(0, a + noccA));
    //FockA->print();

    // MO Coefficient matrix
    SharedMatrix Ca = SharedMatrix(ref_wfn->Ca());
    CmoA = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso, nmo));
    CmoA->set(Ca);
    Ca.reset();
    //CmoA->print();

    // Build Hso
    std::shared_ptr<Matrix> Hso_ = std::shared_ptr<Matrix>(new Matrix("SO-basis One-electron Ints", nso, nso));
    std::shared_ptr<Matrix> Tso_ = std::shared_ptr<Matrix>(new Matrix("SO-basis Kinetic Energy Ints", nso, nso));
    std::shared_ptr<Matrix> Vso_ = std::shared_ptr<Matrix>(new Matrix("SO-basis Potential Energy Ints", nso, nso));
    std::shared_ptr<Matrix> Sso_ = std::shared_ptr<Matrix>(new Matrix("SO-basis Overlap Ints", nso, nso));
    Hso_->zero();
    Tso_->zero();
    Vso_->zero();
    Sso_->zero();

    // Read SO-basis one-electron integrals
    MintsHelper mints(MintsHelper(ref_wfn->basisset(), options, 0));
    Sso_ = mints.ao_overlap();
    Tso_ = mints.ao_kinetic();
    Vso_ = mints.ao_potential();
    Hso_->copy(Tso_);
    Hso_->add(Vso_);
    Tso_.reset();
    Vso_.reset();
    Hso = SharedTensor2d(new Tensor2d("SO-basis One-electron Ints", nso, nso));
    Hso->set(Hso_);
    Hso_.reset();
    Sso = SharedTensor2d(new Tensor2d("SO-basis Overlap Ints", nso, nso));
    Sso->set(Sso_);
    Sso_.reset();
    //Hso->print();

    // Read DF-CC Integrals
    //std::shared_ptr<Molecule> molecule = ref_wfn->molecule();
    //std::shared_ptr<BasisSet> primary = ref_wfn->basisset();
    //outfile->Printf("\tI am here \n");
    std::shared_ptr<BasisSet> auxiliary = ref_wfn->get_basisset("DF_BASIS_CC");
    std::shared_ptr<DFTensor> DF (new DFTensor(ref_wfn->basisset(),auxiliary,ref_wfn->Ca(),noccA,nvirA,naoccA,navirA,options));
    nQ = auxiliary->nbf();

    // read AO basis B tensor
    std::shared_ptr<Matrix> tmp = DF->Qso();
    double ** Qso = tmp->pointer();
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso, nso));
    bQso->set(Qso);
    bQso->write(psio, PSIF_DFOCC_INTS, true, true);
    //bQso->print();
    bQso.reset();
    tmp.reset();

    // Read MO basis B tensors
    // OO Block
    /*
    std::shared_ptr<Matrix> tmpOO = DF->Qoo();
    double ** Qoo = tmpOO->pointer();
    bQij = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQij->set(Qoo);
    bQij->write(psio, PSIF_DFOCC_INTS);
    //bQij->print();
    bQij.reset();
    tmpOO.reset();
    */

    // OV Block
    std::shared_ptr<Matrix> tmpOV = DF->Qov();
    double ** Qov = tmpOV->pointer();
    bQia = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQia->set(Qov);
    bQia->write(psio, PSIF_DFOCC_INTS);
    //bQia->print();
    tmpOV.reset();

    // VV Block
    /*
    std::shared_ptr<Matrix> tmpVV = DF->Qvv();
    double ** Qvv = tmpVV->pointer();
    bQab = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQab->set(Qvv);
    bQab->write(psio, PSIF_DFOCC_INTS, true, true);
    //bQab->print();
    bQab.reset();
    tmpVV.reset();
    */

    // ahmet_df_mp2_rhf

    tstart();

    SharedTensor2d J,T,U;
    J = SharedTensor2d(new Tensor2d("J (IA|JB)", naoccA, navirA, naoccA, navirA));
    T = SharedTensor2d(new Tensor2d("T (IA|JB)", naoccA, navirA, naoccA, navirA)); 
    U = SharedTensor2d(new Tensor2d("U (IA|JB)", naoccA, navirA, naoccA, navirA));

    // Forming J
    J->gemm(true,false,bQia,bQia,1.0,0.0);
    J->print();
    bQia.reset();

    double Ecorr  = 0.0;
    double value  = 0.0;
    double denom  = 0.0;
   
    //computing T
    T->zero();
    T->copy(J);
    T->apply_denom_chem(nfrzc,noccA,FockA);
    T->print();
    
    //forming U
    U->zero();
    U->axpy(T,2.0);
    U->sort(3214,T,-1.0,1.0);

    //computing Ecorr
    Ecorr = U->vector_dot(J);
 
    J.reset();
    eps_orbA.reset();
    T.reset();
    U.reset(); 

    double Emp2 = Escf + Ecorr;

    outfile->Printf("\tSCF Energy        : %20.14f \n", Escf);
    outfile->Printf("\tCorrelation Energy: %20.14f \n", Ecorr);
    outfile->Printf("\tMP2 Energy        : %20.14f \n", Emp2);
    tstop();
}// end rhf

 else if (reference == "UHF") {

    // Ahmet_df_mp2_uhf

        int nmo = ref_wfn->nmo(); // MO sayısı
        int nso = ref_wfn->nso(); // AO sayısı
        int nfrzc = ref_wfn->frzcpi()[0]; // Frozen core sayısı (inaktif OCC) 
        int nfrzv = ref_wfn->frzvpi()[0];  // Frozen virtuals

        int noccA = ref_wfn->doccpi()[0] + ref_wfn->soccpi()[0];  // Toplam OCC sayısı-Alfa
	int noccB = ref_wfn->doccpi()[0];         // Toplam OCC sayısı-Beta
        int nvirA = nmo - noccA;         // Number of virtual orbitals-Alfa
        int nvirB = nmo - noccB;         // Number of virtual orbitals-Beta

        int naoccA = noccA - nfrzc;       // Number of active occupied orbitals-Alfa
        int naoccB = noccB - nfrzc;       // Number of active occupied orbitals-Beta
        int navirA = nvirA - nfrzv;       // Number of active virtual orbitals-Alfa
        int navirB = nvirB - nfrzv;       // Number of active virtual orbitals-Beta

        // Read in nuclear repulsion energy
        double Enuc = ref_wfn->molecule()->nuclear_repulsion_energy();

        // Read SCF energy
        double Escf = ref_wfn->reference_energy();
        double Eref = Escf;
        double Eelec = Escf - Enuc;

        outfile->Printf("\tMO spaces... \n\n");
        outfile->Printf("\t FC   AOCC   BOCC  AVIR   BVIR   FV \n");
        outfile->Printf("\t------------------------------------------\n");
        outfile->Printf("\t%3d   %3d   %3d   %3d    %3d   %3d\n", nfrzc, naoccA, naoccB, navirA, navirB, nfrzv);

        // Read MO energies
        SharedVector epsilon_a_ = ref_wfn->epsilon_a();
        //epsilon_a_->print();

        SharedVector epsilon_b_ = ref_wfn->epsilon_b();
        //epsilon_b_->print();

        eps_orbA = std::shared_ptr<Tensor1d>(new Tensor1d("epsilon <P|Q>", nmo));
        for (int p = 0; p < nmo; ++p) eps_orbA->set(p, epsilon_a_->get(0, p));
        //eps_orbA->print();

        eps_orbB = std::shared_ptr<Tensor1d>(new Tensor1d("epsilon <p|q>", nmo));
        for (int p = 0; p < nmo; ++p) eps_orbB->set(p, epsilon_b_->get(0, p));

	// Build Initial fock matrix
        FockA = SharedTensor2d(new Tensor2d("MO-basis alpha Fock matrix", nmo, nmo));
        for (int i = 0; i < noccA; ++i) FockA->set(i, i, epsilon_a_->get(0, i));
        for (int a = 0; a < nvirA; ++a) FockA->set(a + noccA, a + noccA, epsilon_a_->get(0, a + noccA));
        //FockA->print();

        // Build Initial fock matrix
        FockB = SharedTensor2d(new Tensor2d("MO-basis alpha Fock matrix", nmo, nmo));
        for (int i = 0; i < noccB; ++i) FockB->set(i, i, epsilon_b_->get(0, i));
        for (int a = 0; a < nvirB; ++a) FockB->set(a + noccB, a + noccB, epsilon_b_->get(0, a + noccB));
        //FockB->print();

        // Read DF-CC Integrals
        std::shared_ptr<BasisSet> auxiliary = ref_wfn->get_basisset("DF_BASIS_CC");
        std::shared_ptr<DFTensor> DFA (new DFTensor(ref_wfn->basisset(),auxiliary,ref_wfn->Ca(),noccA,nvirA,naoccA,navirA,options));
        nQ = auxiliary->nbf();

       // Read MO basis B tensors
        std::shared_ptr<Matrix> tmpOVA = DFA->Qov();
        double ** QovA = tmpOVA->pointer();
        bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
        bQiaA->set(QovA);
        //bQiaA->write(psio, PSIF_DFOCC_INTS);
        //bQiaA->print();
        tmpOVA.reset();
        DFA.reset();

        // Beta part
        std::shared_ptr<DFTensor> DFB (new DFTensor(ref_wfn->basisset(),auxiliary,ref_wfn->Cb(),noccB,nvirB,naoccB,navirB,options));
        std::shared_ptr<Matrix> tmpOVB = DFB->Qov();
        double ** QovB = tmpOVB->pointer();
        bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB));
        bQiaB->set(QovB);
        //bQiaB->write(psio, PSIF_DFOCC_INTS);
        //bQiaB->print();
        tmpOVB.reset();
        DFB.reset();

        /*
        // MO Coefficient matrix
        SharedMatrix Ca = SharedMatrix(ref_wfn->Ca());
        CmoA = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso, nmo));
        CmoA->set(Ca);
        Ca.reset();
        CmoA->print();

        SharedMatrix Cb = SharedMatrix(ref_wfn->Cb());
        CmoB = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso, nmo));
        CmoB->set(Cb);
        Cb.reset();
        CmoA->print();
        */
    
          
    tstart();
    SharedTensor2d J,I,Taa,AI,K,N,Tab,L,M,Tbb,AN;

    J  = SharedTensor2d(new Tensor2d("J (IA|JB)", naoccA, navirA, naoccA, navirA));
    I  = SharedTensor2d(new Tensor2d("I (IJ|AB)", naoccA, naoccA, navirA, navirA));
    Taa  = SharedTensor2d(new Tensor2d("T2 (IJ|AB)", naoccA, naoccA, navirA, navirA));
    AI = SharedTensor2d(new Tensor2d("AI (IJ||AB)",naoccA, naoccA, navirA, navirA));

    K  = SharedTensor2d(new Tensor2d("J (IA|jb)", naoccA, navirA, naoccB, navirB));
    M  = SharedTensor2d(new Tensor2d("I (Ij|Ab)", naoccA, naoccB, navirA, navirB));
    Tab  = SharedTensor2d(new Tensor2d("T2 (Ij|Ab)", naoccA, naoccB, navirA, navirB));

    L  = SharedTensor2d(new Tensor2d("J (ia|jb)", naoccB, navirB, naoccB, navirB));
    N  = SharedTensor2d(new Tensor2d("I (ij|ab)", naoccB, naoccB, navirB, navirB));
    Tbb  = SharedTensor2d(new Tensor2d("T2 (ij|ab)", naoccB, naoccB, navirB, navirB));
    AN = SharedTensor2d(new Tensor2d("AN (ij||ab)",naoccB, naoccB, navirB, navirB));

    double Ecorr = 0.0;
    double Eaa  =  0.0;
    double Eab  =  0.0;
    double Ebb  =  0.0;

    // AA part
    J->gemm(true,false,bQiaA,bQiaA,1.0,0.0);
    J->print();

    //from J to I (chemist to physicist)
    I->sort(1324,J,1.0,0.0);
    I->print();

    //forming Antisymmetrized 
    AI->copy(I);
    AI->sort(1243,I,-1.0,1.0);
    AI->print();

    //forming T2
    Taa->copy(AI);
    Taa->apply_denom(nfrzc,noccA,FockA);
    Taa->print();

    //E_AA
    Eaa = 0.25*(Taa->vector_dot(AI));

    J.reset();
    I.reset();
    AI.reset();
    Taa.reset();

    // AB part

    //forming J
    K->gemm(true,false,bQiaA,bQiaB,1.0,0.0);
    K->print();
    bQiaA.reset();
   

    //from J to I (chemist to physicist)
    M->sort(1324,K,1.0,0.0);
    M->print();

    //forming T2
    Tab->copy(M);
    Tab->apply_denom_os(nfrzc,noccA,noccB,FockA,FockB);
    Tab->print();

    // E_AB
    Eab = Tab->vector_dot(M); 

    K.reset();
    M.reset();
    Tab.reset();
    eps_orbA.reset();

    // BB part
    
    //forming J
    L->gemm(true,false,bQiaB,bQiaB,1.0,0.0);
    L->print();
    bQiaB.reset();

    // from J to I (chemist to physicist)
    N->sort(1324,L,1.0,0.0);
    N->print();

    //forming Antisymmetrized AI
    AN->copy(N);
    AN->sort(1243,N,-1.0,1.0);

    // forming T2
    Tbb->copy(AN);
    Tbb->apply_denom(nfrzc,noccB,FockB);
    Tbb->print();

    // E_BB
    Ebb = 0.25*(Tbb->vector_dot(AN));

    L.reset();
    N.reset();
    Tbb.reset();
    AN.reset();

    eps_orbB.reset();
   
    Ecorr = Eaa + Eab + Ebb;

    double Emp2 = Escf + Ecorr;
    outfile->Printf("\tSCF Energy        : %20.14f \n", Escf);
    outfile->Printf("\tAA Energy         : %20.14f \n", Eaa);
    outfile->Printf("\tAB Energy         : %20.14f \n", Eab);
    outfile->Printf("\tBB Energy         : %20.14f \n", Ebb);
    outfile->Printf("\tCorrelation Energy: %20.14f \n", Ecorr);
    outfile->Printf("\tMP2 Energy        : %20.14f \n", Emp2);
    tstop();
}
    // Typically you would build a new wavefunction and populate it with data
    return ref_wfn;
 }//end uhf

}} // End namespaces

