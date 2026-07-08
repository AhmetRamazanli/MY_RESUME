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

namespace psi{ namespace plugin_cis {


#define BIGNUM 1E100
#define maxiter 1000

// Defs
SharedTensor2d CmoA, FockA, Hso, Tso, Vso, Sso, bQso, bQij, bQia, bQab, bQabA, bQabB, bQiaA, bQiaB, bQijA, bQijB;
SharedTensor2d FockB, cQia, cQiaA, cQiaB, sigma0A, sigma0B, sigmaA, sigmaB, cia, ciaA, ciaB, Cvec;
SharedTensor1d eps_orbA, eps_orbB, cQ, HII, sigma_temp, diag;
SharedTensor1i HNII,list;

int nQ, naoccA, navirA, noccA, nvirA, nmo, nso, nfrzc, nfrzv, nalpha, nbeta, natom, ntri_so, noccB, nvirB, naoccB, navirB;
int Enuc, Eref, Escf, Eelec, nroot;
double reg_sigma;

//void sigma_ov();
SharedTensor1d sigma_ov(SharedTensor2d cia);

void sigma_ovUHF();
//extern int davidson(SharedTensor2d &A, int M, SharedTensor1d &eps, SharedTensor2d &v, double cutoff, int print);
int davidson2(SharedTensor1d &diagA, int M, SharedTensor1d &eps, SharedTensor2d &v, double cutoff, int print);

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
        /*- delta used for sigma -*/
        options.add_double("REG_SIGMA", 1e-4);
    }

    return true;
}

extern "C"
SharedWavefunction plugin_cis(SharedWavefunction ref_wfn, Options& options)
{
    // Get options
    int print = options.get_int("PRINT");
    std::string reference = options.get_str("REFERENCE");
    //options.print();
    // Grab the global (default) PSIO object, for file I/O
    std::shared_ptr<PSIO> psio(_default_psio_lib_);
    reg_sigma = options.get_double("REG_SIGMA");
  
    // Read in nuclear repulsion energy
    Enuc = ref_wfn->molecule()->nuclear_repulsion_energy();

    // Read SCF energy
    Escf = ref_wfn->reference_energy();
    Eref = Escf;
    Eelec = Escf - Enuc;

    nmo = ref_wfn->nmo(); // MO sayısı
    nso = ref_wfn->nso(); // AO sayısı
    nfrzc = ref_wfn->frzcpi()[0]; // Frozen core sayısı (inaktif OCC) 
    nfrzv = ref_wfn->frzvpi()[0];  // Frozen virtuals

//============================================
//======== RHF ===============================
//============================================
if(reference == "RHF"){

    
    /* Your code goes here */
    noccA = ref_wfn->doccpi()[0]; // Toplam OCC sayısı 
    nalpha = ref_wfn->nalphapi()[0]; // alpha e sayısı
    nbeta = ref_wfn->nbetapi()[0];   // beta e sayısı
    natom = ref_wfn->molecule()->natom(); // natom

    // figure out number of MO spaces
    nvirA = nmo - noccA;         // Number of virtual orbitals
    naoccA = noccA - nfrzc;       // Number of active occupied orbitals
    navirA = nvirA - nfrzv;       // Number of active virtual orbitals
    ntri_so = 0.5*nso*(nso+1);

    outfile->Printf("\n");
    outfile->Printf("\tNMO  : %2d \n", nmo);

    // Read MO energies
    SharedVector epsilon_a_ = ref_wfn->epsilon_a();
    //epsilon_a_->print();

    eps_orbA = std::shared_ptr<Tensor1d>(new Tensor1d("epsilon <P|Q>", nmo));
    for (int p = 0; p < nmo; ++p) eps_orbA->set(p, epsilon_a_->get(0, p));
    //eps_orbA->print();

    // Build Initial fock matrix
    FockA = SharedTensor2d(new Tensor2d("MO-basis alpha Fock matrix", nmo, nmo));
    for (int i = 0; i < noccA; ++i) FockA->set(i, i, epsilon_a_->get(0, i));
    for (int a = 0; a < nvirA; ++a) FockA->set(a + noccA, a + noccA, epsilon_a_->get(0, a + noccA));
    //FockA->print();

    // MO Coefficient matrix
    SharedMatrix Ca = SharedMatrix(ref_wfn->Ca());
    CmoA = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso, nmo));
    CmoA->set(Ca);
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
    std::shared_ptr<Matrix> tmpOO = DF->Qoo();
    double ** Qoo = tmpOO->pointer();
    bQij = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQij->set(Qoo);
    //bQij->write(psio, PSIF_DFOCC_INTS);
    //bQij->print();
    //bQij.reset();
    tmpOO.reset();

    // OV Block
    std::shared_ptr<Matrix> tmpOV = DF->Qov();
    double ** Qov = tmpOV->pointer();
    bQia = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQia->set(Qov);
    //bQia->write(psio, PSIF_DFOCC_INTS);
    //bQia->print();
    tmpOV.reset();

    // VV Block
    std::shared_ptr<Matrix> tmpVV = DF->Qvv();
    double ** Qvv = tmpVV->pointer();
    bQab = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQab->set(Qvv);
    //bQab->write(psio, PSIF_DFOCC_INTS, true, true);
    //bQab->print();
    //bQab.reset();
    tmpVV.reset();

    // ahmet_cis_rhf

    tstart();

    nroot = 5;
    // malloc
    cQia     = SharedTensor2d(new Tensor2d("C (Q|IA)", nQ, naoccA, navirA));
    ciaA     = SharedTensor2d(new Tensor2d("C (I|A)", naoccA, navirA));
    sigmaA   = SharedTensor2d(new Tensor2d("Sigma (I|A)", naoccA, navirA));
    cQ = std::shared_ptr<Tensor1d>(new Tensor1d("c <Q>", nQ));
    HII      = SharedTensor1d(new Tensor1d("DiagonalH (I|A)", naoccA*navirA));
    HNII     = SharedTensor1i(new Tensor1i("DiagonalH (I|A)", naoccA*navirA));
    cia      = SharedTensor2d(new Tensor2d ("C (I|A)", naoccA, navirA));
    list     = SharedTensor1i(new Tensor1i ("DiagonalH root (I|A)",nroot));
    Cvec     = SharedTensor2d(new Tensor2d("Cia", naoccA*navirA, nroot));

    // HII
    for(int i=0; i<naoccA; i++){
        int ii = i*naoccA + i;
	for(int a=0; a<navirA; a++){
            int ia = i*navirA + a;
	    HII->set(ia, eps_orbA->get(a+noccA) - eps_orbA->get(i+nfrzc));
            int aa = a*navirA + a;
            double sum = 0.0;
	    for(int Q=0; Q<nQ; Q++){
                sum += bQia->get(Q,ia)*bQia->get(Q,ia);
                sum -= bQij->get(Q,ii)*bQab->get(Q,aa);
            }
            HII->add(ia,sum);
	    HNII->set(ia,ia+1);
	}
    }
    HII->print();
    HNII->print();

    double swap = 0.0;
    for(int ia = 0; ia<naoccA*navirA; ia++){
	for(int j = ia+1; j<naoccA*navirA; j++){
	    if(HII->get(ia) > HII->get(j)){
	       swap = HII->get(ia);
	       HII->set(ia, HII->get(j));
	       HII->set(j, swap);
	       swap = HNII->get(ia);
	       HNII->set(ia, HNII->get(j));
	       HNII->set(j, swap);
	    }
	}
    }
    
    HII->print();
    HNII->print();

    for(int ia = 0; ia<nroot; ia++){
	list->set(ia, HNII->get(ia));
    }
    for(int ia = 0; ia<nroot; ia++){
	for(int j=ia+1; j<nroot; j++){
	    if(list->get(ia) > list->get(j)){
	       swap = list->get(ia);
	       list->set(ia, list->get(j));
	       list->set(j, swap);
	    }
	}
    }

    list->print();
    
    //form the Cvec(sigma initial)
    Cvec->zero();
    for(int root = 0; root<nroot; root++){
	int ia = list->get(root) - 1;
	Cvec->set(ia,root,1);
    }

    Cvec->print();
 
    // call davidson
    SharedTensor1d diag = std::shared_ptr<Tensor1d>(new Tensor1d("eigenvalues", naoccA*navirA));
    
    // initial Cvec
    
    davidson2(HII, nroot, diag, Cvec, 1e-10, 0);
   
    diag.reset();
    

    //resetting tensors and vectors    
    sigmaA.reset();
    cQia.reset();
    bQia.reset();
    bQij.reset();
    bQab.reset();
    cQ->zero();
    eps_orbA.reset();
    HII.reset();
    HNII.reset();
    list.reset();

    tstop();

}// end rhf

//============================================
//======== UHF ===============================
//============================================
else if (reference == "UHF") {
    

    /* Your code goes here */
    noccA = ref_wfn->doccpi()[0] + ref_wfn->soccpi()[0];  // Toplam OCC sayısı-Alfa
    noccB = ref_wfn->doccpi()[0];         // Toplam OCC sayısı-Beta
    nvirA = nmo - noccA;         // Number of virtual orbitals-Alfa
    nvirB = nmo - noccB;         // Number of virtual orbitals-Beta

    naoccA = noccA - nfrzc;       // Number of active occupied orbitals-Alfa
    naoccB = noccB - nfrzc;       // Number of active occupied orbitals-Beta
    navirA = nvirA - nfrzv;       // Number of active virtual orbitals-Alfa
    navirB = nvirB - nfrzv;       // Number of active virtual orbitals-Beta


    // ahmet_cis_uhf
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
    //eps_orbB->print();

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

    //Read MO basis B tensors
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

    
    // MO Coefficient matrix
    SharedMatrix Ca = SharedMatrix(ref_wfn->Ca());
    //CmoA = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso, nmo));
    //CmoA->set(Ca);
    //Ca.reset();
    //CmoA->print();

    SharedMatrix Cb = SharedMatrix(ref_wfn->Cb());
    //CmoB = SharedTensor2d(new Tensor2d("Alpha MO Coefficients", nso, nmo));
    //CmoB->set(Cb);
    //Cb.reset();
    //CmoA->print();
    
    std::shared_ptr<DFTensor> DFooA (new DFTensor(ref_wfn->basisset(),auxiliary,ref_wfn->Ca(),noccA,nvirA,naoccA,navirA,options));
    std::shared_ptr<Matrix> tmpOOA = DFooA->Qoo();
    double ** QooA = tmpOOA->pointer();
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQijA->set(QooA);
    //bQijA->write(psio, PSIF_DFOCC_INTS);
    //bQijA->print();
    //bQijA.reset();
    tmpOOA.reset();
    DFooA.reset();

    // OO Block
    std::shared_ptr<DFTensor> DFooB (new DFTensor(ref_wfn->basisset(),auxiliary,ref_wfn->Cb(),noccB,nvirB,naoccB,navirB,options));
    std::shared_ptr<Matrix> tmpOOB = DFooB->Qoo();
    double ** QooB = tmpOOB->pointer();
    bQijB = SharedTensor2d(new Tensor2d("DFOO_BASIS_CC B (Q|IJ)", nQ, naoccB, naoccB));
    bQijB->set(QooB);
    //bQijB->write(psio, PSIF_DFOCC_INTS);
    //bQijB->print();
    //bQijB.reset();
    tmpOOB.reset();
    DFooB.reset();

    // VV Block alpha
    std::shared_ptr<DFTensor> DFvvA (new DFTensor(ref_wfn->basisset(),auxiliary,ref_wfn->Ca(),noccA,nvirA,naoccA,navirA,options));
    std::shared_ptr<Matrix> tmpVVA = DFvvA->Qvv();
    double ** QvvA = tmpVVA->pointer();
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->set(QvvA);
    //bQabA->write(psio, PSIF_DFOCC_INTS, true, true);
    //bQabA->print();
    //bQabA.reset();
    tmpVVA.reset();   
    DFvvA.reset();
  
    // VV Block beta
    std::shared_ptr<DFTensor> DFvvB (new DFTensor(ref_wfn->basisset(),auxiliary,ref_wfn->Cb(),noccB,nvirB,naoccB,navirB,options));
    std::shared_ptr<Matrix> tmpVVB = DFvvB->Qvv();
    double ** Qvv = tmpVVB->pointer();
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirB, navirB));
    bQabB->set(Qvv);
    //bQabB->write(psio, PSIF_DFOCC_INTS, true, true);
    //bQabB->print();
    //bQabB.reset();
    tmpVVB.reset();
    DFvvB.reset();

    tstart();

    // malloc
    cQiaA = SharedTensor2d(new Tensor2d("C (Q|IA)", nQ, naoccA, navirA));
    cQiaB = SharedTensor2d(new Tensor2d("C (Q|ia)", nQ, naoccB, navirB));
    ciaA = SharedTensor2d(new Tensor2d("C (I|A)", naoccA, navirA));
    ciaB = SharedTensor2d(new Tensor2d("C (i|a)", naoccB, navirB));
    sigmaA = SharedTensor2d(new Tensor2d("Sigma (I|A)", naoccA, navirA));
    sigmaB = SharedTensor2d(new Tensor2d("Sigma (i|a)", naoccB, navirB));
    cQ = std::shared_ptr<Tensor1d>(new Tensor1d("c <Q>", nQ));
    sigma0A = SharedTensor2d(new Tensor2d("Sigma initial for alpha (I|A)",naoccA,navirA));
    sigma0B = SharedTensor2d(new Tensor2d("Sigma initial for beta (i|a)",naoccB,navirB));

    sigma0A->zero();    
    // compute initial sigma for alpha
    for(int i=0; i<naoccA; i++){
        for(int a=0; a<navirA; a++){
            sigma0A->set(i,a,reg_sigma*(eps_orbA->get(a+noccA) - eps_orbA->get(i+nfrzc)));
	}
    }
    sigma0A->print();

    // compute initial sigma for beta
    for(int i=0; i<naoccB; i++){
        for(int a=0; a<navirB; a++){
            sigma0B->set(i,a,reg_sigma*(eps_orbB->get(a+noccB) - eps_orbB->get(i+nfrzc)));
        }
    }
    sigma0B->print();

    // call sigma
    sigma_ovUHF();

    //call davidson

    //print sigmaA
    sigmaA->print();
    sigmaB->print();

    //resetting tensors and vectors
    bQijA.reset();
    bQijB.reset();
    sigmaA.reset();
    sigmaB.reset();
    sigma0A.reset();
    sigma0B.reset();
    cQ.reset();
    cQiaA.reset();
    cQiaB.reset();
    ciaA.reset();
    ciaB.reset();
    eps_orbA.reset();
    eps_orbB.reset();
    Ca.reset();
    Cb.reset();

    tstop();

}// end uhf


    // Typically you would build a new wavefunction and populate it with data
    return ref_wfn;

}

//=========================================
//======= Sigma_ia & Sigma initial ========
//=========================================

SharedTensor1d sigma_ov(SharedTensor2d cia) //RHF
{   
    cQ->gemv(false,bQia,cia,2.0,0.0);
    cQia->gemv(false,bQab,cia,1.0,0.0);

    //computing sigmaA 
    for(int i=0; i<naoccA; i++){
        for(int a=0; a<navirA; a++){
            sigmaA->set(i,a,cia->get(i,a)*(eps_orbA->get(a+noccA) - eps_orbA->get(i+nfrzc)));
        }
    }
    sigmaA->gemv(true,bQia,cQ,1.0,1.0);
    sigmaA->contract332(true,false,naoccA,bQij,cQia,-1.0,1.0);
  // sigmaA->print();
 
    // form sigma vec for a root
    sigma_temp = SharedTensor1d(new Tensor1d("DiagonaltempH (I|A)", naoccA*navirA));
    for(int i=0; i<naoccA; i++){
        for(int a=0; a<navirA; a++){
            int ia = i*navirA + a;
            sigma_temp->set(ia,sigmaA->get(i,a));
        }
    }    
    return sigma_temp;

    
}// end sigma_ov

void sigma_ovUHF()  //UHF
{   
    ciaA->zero(); 
    cQ->gemv(false,bQiaA,ciaA,1.0,0.0);
    cQ->gemv(false,bQiaB,ciaB,1.0,1.0);
    cQiaA->gemv(false,bQabA,ciaA,1.0,0.0);
    cQiaB->gemv(false,bQabB,ciaB,1.0,0.0);

    sigmaA->zero();
    //computing sigmaAA 
    for(int i=0; i<naoccA; i++){
        for(int a=0; a<navirA; a++){
            sigmaA->set(i,a,ciaA->get(i,a)*(eps_orbA->get(a+noccA) - eps_orbA->get(i+nfrzc)));
        }
    }
    sigmaA->gemv(true,bQiaA,cQ,1.0,1.0);
    sigmaA->contract332(true,false,naoccA,bQijA,cQiaA,-1.0,1.0);

    //computing sigmaAB 
    for(int i=0; i<naoccB; i++){
        for(int a=0; a<navirB; a++){
            sigmaB->set(i,a,ciaB->get(i,a)*(eps_orbB->get(a+noccB) - eps_orbB->get(i+nfrzc)));
        }
    }
    sigmaB->gemv(true,bQiaB,cQ,1.0,1.0);
    sigmaB->contract332(true,false,naoccB,bQijB,cQiaB,-1.0,0.0);
}

/*!
** davidson(): Computes the lowest few eigenvalues and eigenvectors of a
** symmetric matrix, A, using the Davidson-Liu algorithm.  
**
** The matrix must be small enough to fit entirely in core.  This algorithm 
** is useful if one is interested in only a few roots of the matrix
** rather than the whole spectrum.
**
** NB: This implementation will keep up to eight guess vectors for each
** root desired before collapsing to one vector per root.  In
** addition, if smart_guess=1 (the default), guess vectors are
** constructed by diagonalization of a sub-matrix of A; otherwise,
** unit vectors are used.
** UB, November 2017
**
** \param diagA      = diagonals of A
** \param N      = dimension of A
** \param M      = number of roots desired
** \param eps    = eigenvalues
** \param v      = eigenvectors
** \param cutoff = tolerance for convergence of eigenvalues
** \param print  = Boolean for printing additional information
**
** Returns: number of converged roots
*/

int davidson2(SharedTensor1d &diagA, int M, SharedTensor1d &eps, SharedTensor2d &v, double cutoff, int print) {
  int i, j, k, L, I;
  double minimum;
  int min_pos, numf, iter, converged, maxdim, skip_check;
  int init_dim;
  //int smart_guess = 0;
  double norm, denom, diff;

  SharedTensor2d B, Bnew, sigma, sigmaia, G, alpha, Res;
  SharedTensor1d diag,lambda, lambda_old;
  SharedTensor1i small2big, conv;

  // get dimension of matrix A
  int N = diagA->dim1();

  // define maximum dimension
  maxdim = 8 * M;

  // Current set of guess vectors, stored by row
  B = SharedTensor2d(new Tensor2d("B matrix", maxdim, N));
  // Guess vectors formed from old vectors, stored by row
  Bnew = SharedTensor2d(new Tensor2d("New B matrix", M, N));
  // Sigma vectors, stored by column
  sigma = SharedTensor2d(new Tensor2d("sigma matrix", N, maxdim));
  // Davidson mini-Hamitonian
  G = SharedTensor2d(new Tensor2d("G matrix", maxdim, maxdim));
  // Residual eigenvectors, stored by row
  Res = SharedTensor2d(new Tensor2d("f matrix", maxdim, N));
  // Eigenvectors of G
  alpha = SharedTensor2d(new Tensor2d("alpha matrix", maxdim, maxdim));
  // Eigenvalues of G
  lambda = SharedTensor1d(new Tensor1d("lambda vector", maxdim));
  // Approximate roots from previous iteration
  lambda_old = SharedTensor1d(new Tensor1d("lambda vector", maxdim));

  sigmaia = SharedTensor2d(new Tensor2d("cia", naoccA, navirA));
    // init dim
    if (N > 7*M) init_dim = 7*M;
    else init_dim = M;

    // Diagonals of A
    diag = SharedTensor1d(new Tensor1d("A diag vector", N));
    // Addressing
    small2big = SharedTensor1i(new Tensor1i("A diag vector", 7*M));

    // diag of A
    for (i=0; i < N; i++) diag->set(i,diagA->get(i)); 

    // Form addressing array 
    for (i=0; i < init_dim; i++) {
         minimum = diag->get(0);
         min_pos = 0;
         for (j=1; j < N; j++) {
	      if(diag->get(j) < minimum) {	
                 minimum = diag->get(j);
	         min_pos = j; 
		 small2big->set(i,j);
	      } 
         }
         B->set(i, min_pos, 1.0);
         diag->set(min_pos, BIGNUM);
         lambda_old->set(i, minimum);
    }
  diag.reset();
  //===================================
  //======== Initial guess ============
  //===================================
  // Use eigenvectors of a sub-matrix as initial guesses
  /*
  if (smart_guess) { 

    // init dim
    if (N > 7*M) init_dim = 7*M;
    else init_dim = M;

    // Diagonals of A
    diag = SharedTensor1d(new Tensor1d("A diag vector", N));
    // Addressing
    small2big = SharedTensor1i(new Tensor1i("A diag vector", 7*M));

    // diag of A
    for (i=0; i < N; i++) diag->set(i, A->get(i,i)); 

    // Form addressing array 
    for (i=0; i < init_dim; i++) {
         minimum = diag->get(0);
         min_pos = 0;
         for (j=1; j < N; j++) {
	      if(diag->get(j) < minimum) {	
                 minimum = diag->get(j);
	         min_pos = j; 
		 small2big->set(i,j);
	      } 
         }
	 diag->set(min_pos, BIGNUM);
	 lambda_old->set(i, minimum);
    }

    // Form G
    for (i=0; i < init_dim; i++) {
	 int ii = small2big->get(i);
         for (j=0; j < init_dim; j++) {
	      int jj = small2big->get(j);
	      G->set(i, j, A->get(ii,jj));

         }
    }

    // Diagonalize G
    G->diagonalize(init_dim, alpha, lambda, 1e-12, true);

    // Form B 
    for (i=0; i < init_dim; i++) {
         for (j=0; j < init_dim; j++) {
	      int jj = small2big->get(j);
	      B->set(i, jj, alpha->get(j,i));
	 }
    }

    diag.reset();
    small2big.reset();
  }// end if(smart_guess)
  

  // Use unit vectors as initial guesses
  //else { 
    // Diagonals of A
    diag = SharedTensor1d(new Tensor1d("H diag vector", N));
    for (i=0; i < N; i++) diag->set(i, diagA->get(i));
	 
    // Form B
    for (i=0; i < M; i++) {
         minimum = diag->get(0);
         min_pos = 0;
         for (j=1; j < N; j++) {
	      if (diag->get(j) < minimum) { 
		  minimum = diag->get(j); 
	          min_pos = j; 
	      }
         }

	 B->set(i, min_pos, 1.0);
	 diag->set(min_pos, BIGNUM);
	 lambda_old->set(i, minimum);
    }
    diag.reset();
  //} // end else
*/

  //===================================
  //======== Loop =====================
  //===================================
  // start
  L = init_dim;
  iter = 0;
  converged = 0;

  // boolean array for convergence of each root
  conv = SharedTensor1i(new Tensor1i("conv vector", M));

  // Head of Loop
  while(converged < M && iter < maxiter) {

    skip_check = 0;
    if(print) printf("\niter = %d\n", iter); 

     // Form G
     // sigma = A*B'
     // sigma->gemm(false, true, A, B, 1.0, 0.0);
     // call sigma_ov;
     int K,g;
     for (k=0; k < maxdim; k++){
          for (i=0; i < naoccA; i++ ) {
	       for (int a=0; a < navirA; a++ ) {
                    int ia = i*navirA + a;
        	    if(k < M){
                       sigma->set(ia,k,v->get(ia,k));
        	    }
		    else{
		       K=k;
		       k=k-M;
		       sigma->set(ia,K,sigma_ov(sigma)->get(ia));
		       k=k+M;
		    }   
	       }
	  }
     }
     outfile->Printf("converged %d ",converged);
     outfile->Printf(" iter %d", iter);
     sigma->print();
     

    // G = B*sigma
    G->gemm(false, false, B, sigma, 1.0, 0.0);

    // Diagonalize G
    G->diagonalize(L, alpha, lambda, 1e-12, true);

    // Form preconditioned residue vectors
    for (k=0; k < M; k++) {
         for (I=0; I < N; I++) {
	      Res->set(k,I, 0.0);
	      double value = 0.0;
	      for(i=0; i < L; i++) {
		  value += alpha->get(i,k) * (sigma->get(I,i) - lambda->get(k) * B->get(i,I));
	      }
              Res->add(k,I,value);
	      denom = lambda->get(k) - diagA->get(I);
	      if (fabs(denom) > 1e-6) Res->set(k,I, Res->get(k,I)/denom);
	      else Res->set(k,I,0.0);
         }
    }

    // Normalize each residual
    for (k=0; k < M; k++) {
         norm = 0.0;
         for (I=0; I < N; I++) {
	     norm += Res->get(k,I) * Res->get(k,I);
         }
      
	 norm = std::sqrt(norm);
         for (I=0; I < N; I++) {
	      if (norm > 1e-6) Res->set(k,I, Res->get(k,I)/norm);
	      else Res->set(k,I,0.0);
         }
    }

    // Schmidt orthogonalize the Res[k] against the set of B[i] and add new vectors
    for (k=0,numf=0; k < M; k++) {
         SharedTensor1d Rk = SharedTensor1d(new Tensor1d("Res[k] vector", N));
	 Rk->row_vector(Res,k);
         if (B->gs_add(L,Rk)) { 
             L++; 
	     numf++; 
	 }
    }

    // If L is close to maxdim, collapse to one guess per root
    if (maxdim - L < M) {
        if (print) {
	    printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
	    printf("Collapsing eigenvectors.\n");
        }

	// form Bnew
        for (i=0; i < M; i++) {
	     for (k=0; k < N; k++) {
		  double sum = 0.0;
	          for (j=0; j < L; j++) {
	               sum += alpha->get(j,i) * B->get(j,k);
	          } 
		  Bnew->set(i,k,sum);
	     }
        }

        // Copy new vectors into place
	B->copy(Bnew);
        skip_check = 1;
        L = M;
    } // end if

    // check convergence on all roots
    if (!skip_check) {
        converged = 0;
	conv->zero();
        if (print) {
	    printf("Root      Eigenvalue       Delta  Converged?\n");
	    printf("---- -------------------- ------- ----------\n");
        }  

        for (k=0; k < M; k++) {
	     diff = std::fabs(lambda->get(k) - lambda_old->get(k));
	     if (diff < cutoff) {
		 conv->set(k,1);
	         converged++;
	     }
	     lambda_old->set(k,lambda->get(k));
	     if (print) {
	         printf("%3d  %20.14f %4.3e    %1s\n", k, lambda->get(k), diff, conv->get(k) == 1 ? "Y" : "N");
	     }
        }
    }// end if (!skip_check)
      
      //alpha->print();
      //B->print();
      
      
    // get eigen vecs
      for (i=0; i < M; i++) {
           eps->set(i,lambda->get(i));
	   for (I=0; I < N; I++) {
	        double sum = 0.0;
                for (j=0; j < L; j++) {
		     sum += alpha->get(j,i) * B->get(j,I);
	        }  
	        v->add(I,i,sum);
           }
      }
    iter++;
  }// end of loop 
      outfile->Printf("iter %d", iter); 
     // get sigma    
   
  //v->print();
  // Generate final eigenvalues and eigenvectors
  if (converged == M) {
      for (i=0; i < M; i++) {
           eps->set(i,lambda->get(i));
	   for (I=0; I < N; I++) {
	        double sum = 0.0;
                for (j=0; j < L; j++) {
		     sum += alpha->get(j,i) * B->get(j,I);
	        }  
	        v->add(I,i,sum);
           }
      } 
      if (print) printf("Davidson algorithm converged in %d iterations.\n", iter);
  }

 
  conv.reset();
  B.reset();
  Bnew.reset();
  sigma.reset();
  G.reset();
  Res.reset();
  alpha.reset();
  lambda.reset();
  lambda_old.reset();

  return converged;
}

}} //namespace
