//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file disk.cpp
//  \brief Initializes stratified Keplerian accretion disk in both cylindrical and
//  spherical polar coordinates.  Initial conditions are in vertical hydrostatic eqm.

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min
#include <cstdlib>    // srand
#include <cfloat>     // FLT_MIN

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
static Real DenProfileCyl(const Real rad, const Real phi, const Real z);
static Real PoverR(const Real rad, const Real phi, const Real z);
static void VelProfileCyl(const Real rad, const Real phi, const Real z,
  Real &v1, Real &v2, Real &v3,const Real rad2);
static void Rotation(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
static void VaryingViscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, int is, int ie, int js, int je, int ks, int ke);
// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

// problem parameters which are useful to make global to this file
static Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
static Real dfloor;
static Real companionmass, companiondistance;
//static Real Disk_Truncate,TruncateWidth;
//static Real InnerLimbo;
static Real RInnerBoundary, ROuterBoundary;
//static Real SmearMode;
//static Real Truncate_Scale;
static Real Perturb_Type;
static Real ViscHalfRecess;
static Real kr,kz,amp,phi0,theta,rref,phiref,zref;
//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);
  companiondistance = pin->GetReal("problem","companiondistance");
  companionmass = pin->GetReal("problem","companionmass");
  RInnerBoundary = pin->GetReal("problem","RInnerBoundary");
  ROuterBoundary = pin->GetReal("problem","ROuterBoundary");
  ViscHalfRecess = pin->GetReal("problem","ViscHalfRecess");
  kr = pin->GetReal("problem","kr");
  kz = pin->GetReal("problem","kz");
  amp = pin->GetReal("problem","amp");
  theta = std::acos(kz/sqrt(kz*kz+kr*kr));
  rref = pin->GetReal("problem","rref");
  phiref = pin->GetReal("problem","phiref");
  zref = pin->GetReal("problem","zref");
  phi0 = pin->GetReal("problem","phi0");

  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  //EnrollUserExplicitSourceFunction(Rotation);
  //EnrollViscosityCoefficient(VaryingViscosity);
  // enroll user-defined boundary condition
  if (mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X1, DiskInnerX1);
  }
  if (mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X1, DiskOuterX1);
  }
  if (mesh_bcs[INNER_X2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X2, DiskInnerX2);
  }
  if (mesh_bcs[OUTER_X2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X2, DiskOuterX2);
  }
  if (mesh_bcs[INNER_X3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X3, DiskInnerX3);
  }
  if (mesh_bcs[OUTER_X3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X3, DiskOuterX3);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad, phi, z;
  Real v1, v2, v3;
  Real rad2=0,phi2,z2;
  Real xsheet,ysheet,zsheet;
  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
      // compute initial conditions in cylindrical coordinates
      xsheet = rad - rref;
      zsheet = z;
      Real Ux,Uy,Uz,drho,dp,Al,kappa,Bl,Omega,omegak;
      Al = - 0.25 * (3*gm0/pow(rref,3.0) - (dslope+pslope)*(pslope-2)*p0_over_r0*pow(rref,pslope-2)) / sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
      Omega = sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
      kappa = sqrt(4*Omega*(Omega+Al));
      Bl = Omega + Al;
      omegak = kappa * cos(theta);
      Ux = sqrt(amp)*cos(theta)*cos(phi0)*cos(kr*xsheet+kz*zsheet);
      Uy = -2*Bl/kappa * sqrt(amp)*sin(phi0)*cos(kr*xsheet+kz*zsheet);
      Uz = -kr/kz * Ux;
      drho = 0;//omegak / kz * rho0 / p0_over_r0 * pow(rref,dslope - pslope) * Uz;
      dp = 0;//drho * p0_over_r0 * pow(rad,pslope);
      phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z) + drho;
      VelProfileCyl(rad,phi,z,v1,v2,v3,rad2);
  
      phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*(v1+Ux);
      phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*(v2+Uy);
      phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*(v3+Uz);
      if (NON_BAROTROPIC_EOS) {
        Real p_over_r = PoverR(rad,phi,z);
        phydro->u(IEN,k,j,i) = (p_over_r*phydro->u(IDN,k,j,i) + dp)/(gamma_gas - 1.0);
        phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                   + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);

        
      }
    }
  }
}
  return;
}

//----------------------------------------------------------------------------------------
//~\f to conduct the cooling process. It would incur singularity near inner boundaries of meshblocks 
//if done in the explcit source funciton rather than here, unknown why.
/*void MeshBlock::UserWorkInLoop(){
  Real rad, phi, z;
  Real v1, v2, v3;
  Real dt=pmy_mesh->dt;
  Real KeplerianOmega;
  KeplerianOmega = sqrt((gm0+companionmass)/pow(companiondistance,3.0));
  Real g = gamma_gas;
  Real tau = 1.0 / KeplerianOmega;
  int flag3 = int(block_size.nx3>1);
  int flag2 = int(block_size.nx2>1);
  //  Initialize density and momenta
  for (int k=ks-NGHOST*flag3; k<=ke+NGHOST*flag3; ++k) {
    for (int j=js-NGHOST*flag2; j<=je+NGHOST*flag2; ++j) {
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        //if(i<is+block_size.nx3/16.0){
        //  phydro->u(IM3,k,j,i)=0;
        //  phydro->w(IVZ,k,j,i)=0;
        //}
        Real temp;
        Real poverr = PoverR(rad,phi,z); 
        temp = phydro->w(IEN,k,j,i)/phydro->w(IDN,k,j,i);
        phydro->u(IEN,k,j,i) -= dt*phydro->w(IDN,k,j,i)*(temp - poverr)/tau/(g-1.0);
        phydro->w(IEN,k,j,i) -= dt*phydro->w(IDN,k,j,i)*(temp - poverr)/tau;
        }
      }
    }
  return;
}*/
//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate

static void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (COORDINATE_SYSTEM == "cylindrical") {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    rad=fabs(pco->x1v(i)*sin(pco->x2v(j)));
    phi=pco->x3v(i);
    z=pco->x1v(i)*cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes density in cylindrical coordinates

static Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*pow(rad/r0,dslope);
  Real dentem = denmid*exp(gm0/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = denmid;//dentem;
  /*if(rad>Disk_Truncate){
     //Real den0=DenProfileCyl(Disk_Truncate,phi,z);
     //den=den0*exp(-(rad - Disk_Truncate)*(rad - Disk_Truncate)/(TruncateWidth*TruncateWidth));
    den *= exp(-(rad - Disk_Truncate)/Truncate_Scale); 
  }*/
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! \f  computes pressure/density in cylindrical coordinates
static Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  //poverr = p0_over_r0*pow(DenProfileCyl(rad,phi,z), gamma_gas-1);
  poverr = p0_over_r0*pow(rad/r0, pslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates

static void VelProfileCyl(const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3, const Real rad2) {
  Real p_over_r = PoverR(rad, phi, z);
  Real p_over_rref = PoverR(rref,phi,z);
  Real vel,vel0;
  /*if(rad > Disk_Truncate){
     Real PressGradient=(PoverR(rad,phi,z)*DenProfileCyl(rad,phi,z)-PoverR(rad2,phi,z)*DenProfileCyl(rad2,phi,z))/(rad-rad2);
     vel = 1.0 + rad/DenProfileCyl(rad,phi,z)*PressGradient/(gm0/rad)+pslope - pslope*rad/std::sqrt(rad*rad+z*z);
     //std::cout<<rad<<" "<<PressGradient<<" "<<vel<<std::endl;
  }*/
  vel = (dslope+pslope)*p_over_r/(gm0/rad) + 1.0;//(1.0+pslope) - pslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(gm0/rad)*std::sqrt(vel);//
  vel0 = (dslope+pslope)*p_over_rref/(gm0/rref) + 1.0;
  vel0 = std::sqrt(gm0/rref)*std::sqrt(vel0);
  //vel-= vel0/rref*rad;//vel/(vel+KeplerianOmega*rad)*;
  if (COORDINATE_SYSTEM == "cylindrical") {
    v1=0.0;
    v2=vel;
    v3=0.0;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    v1=0.0;
    v2=0.0;
    v3=vel;
  }
  return;
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: sets solution in ghost zones to initial values
//

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  Real rad,phi,z,rad2;
  Real v1, v2, v3;
  Real xsheet,ysheet,zsheet;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
        GetCylCoord(pco,rad,phi,z,is-i,j,k);
        GetCylCoord(pco,rad2,phi,z,is-i+1,j,k);
        xsheet = rad - rref;
        zsheet = z;
        Real Ux,Uy,Uz,drho,dp,Al,kappa,Bl,Omega,omegak;
        Al = - 0.25 * (3*gm0/pow(rref,3.0) - (dslope+pslope)*(pslope-2)*p0_over_r0*pow(rref,pslope-2)) / sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        Omega = sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        kappa = sqrt(4*Omega*(Omega+Al));
        Bl = Omega + Al;
        omegak = kappa * cos(theta);
        Ux = sqrt(amp)*cos(theta)*cos(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uy = -2*Bl/kappa * sqrt(amp)*sin(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uz = -kr/kz * Ux;
        drho = 0;//omegak / kz * rho0 / p0_over_r0 * pow(rref,dslope - pslope) * Uz;
        dp = 0;//drho * p0_over_r0 * pow(rad,pslope);
        VelProfileCyl(rad,phi,z,v1,v2,v3,rad2);
        prim(IDN,k,j,is-i) = DenProfileCyl(rad,phi,z) + drho;
        prim(IM1,k,j,is-i) = v1 + Ux;
        prim(IM2,k,j,is-i) = v2 + Uy;
        prim(IM3,k,j,is-i) = v3 + Uz;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,is-i) = PoverR(rad, phi, z)*prim(IDN,k,j,is-i) + dp;
      }
    }
  }
}

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  Real rad,phi,z,rad2;
  Real v1, v2, v3;
  Real xsheet,ysheet,zsheet;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
        GetCylCoord(pco,rad,phi,z,ie+i,j,k);
        GetCylCoord(pco,rad2,phi,z,ie+i-1,j,k);
        xsheet = rad - rref;
        zsheet = z;
        Real Ux,Uy,Uz,drho,dp,Al,kappa,Bl,Omega,omegak;
        Al = - 0.25 * (3*gm0/pow(rref,3.0) - (dslope+pslope)*(pslope-2)*p0_over_r0*pow(rref,pslope-2)) / sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        Omega = sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        kappa = sqrt(4*Omega*(Omega+Al));
        Bl = Omega + Al;
        omegak = kappa * cos(theta);
        Ux = sqrt(amp)*cos(theta)*cos(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uy = -2*Bl/kappa * sqrt(amp)*sin(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uz = -kr/kz * Ux;
        drho = 0;//omegak / kz * rho0 / p0_over_r0 * pow(rref,dslope - pslope) * Uz;
        dp = 0;//drho * p0_over_r0 * pow(rad,pslope);
        VelProfileCyl(rad,phi,z,v1,v2,v3,rad2);
        prim(IDN,k,j,ie+i) = DenProfileCyl(rad,phi,z) + drho;
        prim(IM1,k,j,ie+i) = v1 + Ux;
        prim(IM2,k,j,ie+i) = v2 + Uy;
        prim(IM3,k,j,ie+i) = v3 + Uz;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,ie+i) = PoverR(rad, phi, z)*prim(IDN,k,j,ie+i) + dp;
      }
    }
  }
}

void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  Real rad,phi,z;
  Real v1, v2, v3;
  Real xsheet,ysheet,zsheet;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pco,rad,phi,z,i,js-j,k);
        xsheet = rad - rref;
        zsheet = z;
        Real Ux,Uy,Uz,drho,dp,Al,kappa,Bl,Omega,omegak;
        Al = - 0.25 * (3*gm0/pow(rref,3.0) - (dslope+pslope)*(pslope-2)*p0_over_r0*pow(rref,pslope-2)) / sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        Omega = sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        kappa = sqrt(4*Omega*(Omega+Al));
        Bl = Omega + Al;
        omegak = kappa * cos(theta);
        Ux = sqrt(amp)*cos(theta)*cos(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uy = -2*Bl/kappa * sqrt(amp)*sin(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uz = -kr/kz * Ux;
        drho = 0;//omegak / kz * rho0 / p0_over_r0 * pow(rref,dslope - pslope) * Uz;
        dp = 0;//drho * p0_over_r0 * pow(rad,pslope);
        VelProfileCyl(rad,phi,z,v1,v2,v3,0);
        prim(IDN,k,js-j,i) = DenProfileCyl(rad,phi,z) + drho;
        prim(IM1,k,js-j,i) = v1 + Ux;
        prim(IM2,k,js-j,i) = v2 + Uy;
        prim(IM3,k,js-j,i) = v3 + Uz;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,js-j,i) = PoverR(rad, phi, z)*prim(IDN,k,js-j,i) + dp;
      }
    }
  }
}

void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  Real rad,phi,z;
  Real v1, v2, v3;
  Real xsheet,ysheet,zsheet;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pco,rad,phi,z,i,je+j,k);
        xsheet = rad - rref;
        zsheet = z;
        Real Ux,Uy,Uz,drho,dp,Al,kappa,Bl,Omega,omegak;
        Al = - 0.25 * (3*gm0/pow(rref,3.0) - (dslope+pslope)*(pslope-2)*p0_over_r0*pow(rref,pslope-2)) / sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        Omega = sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        kappa = sqrt(4*Omega*(Omega+Al));
        Bl = Omega + Al;
        omegak = kappa * cos(theta);
        Ux = sqrt(amp)*cos(theta)*cos(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uy = -2*Bl/kappa * sqrt(amp)*sin(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uz = -kr/kz * Ux;
        drho = 0;//omegak / kz * rho0 / p0_over_r0 * pow(rref,dslope - pslope) * Uz;
        dp = 0;//drho * p0_over_r0 * pow(rad,pslope);
        VelProfileCyl(rad,phi,z,v1,v2,v3,0);

        prim(IDN,k,je+j,i) = DenProfileCyl(rad,phi,z) + drho;
        prim(IM1,k,je+j,i) = v1 + Ux;
        prim(IM2,k,je+j,i) = v2 + Uy;
        prim(IM3,k,je+j,i) = v3 + Uz;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,je+j,i) = PoverR(rad, phi, z)*prim(IDN,k,je+j,i) + dp;
      }
    }
  }
}

void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  Real rad,phi,z;
  Real v1, v2, v3;
  Real xsheet,ysheet,zsheet;
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        //GetCylCoord(pco,rad,phi,z,i,j,ks-k);
        GetCylCoord(pco,rad,phi,z,i,j,ks);
        xsheet = rad - rref;
        zsheet = z;
        Real Ux,Uy,Uz,drho,dp,Al,kappa,Bl,Omega,omegak;
        Al = - 0.25 * (3*gm0/pow(rref,3.0) - (dslope+pslope)*(pslope-2)*p0_over_r0*pow(rref,pslope-2)) / sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        Omega = sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        kappa = sqrt(4*Omega*(Omega+Al));
        Bl = Omega + Al;
        omegak = kappa * cos(theta);
        Ux = sqrt(amp)*cos(theta)*cos(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uy = -2*Bl/kappa * sqrt(amp)*sin(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uz = -kr/kz * Ux;
        drho = 0;//omegak / kz * rho0 / p0_over_r0 * pow(rref,dslope - pslope) * Uz;
        dp = 0;//drho * p0_over_r0 * pow(rad,pslope);
        VelProfileCyl(rad,phi,z,v1,v2,v3,0);
        prim(IDN,ks-k,j,i) = DenProfileCyl(rad,phi,z) + drho;
        prim(IM1,ks-k,j,i) = v1 + Ux;
        prim(IM2,ks-k,j,i) = v2 + Uy;
        prim(IM3,ks-k,j,i) = v3 + Uz;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,ks-k,j,i) = PoverR(rad, phi, z)*prim(IDN,ks-k,j,i) + dp;
      }
    }
  }
}

void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  Real rad,phi,z;
  Real v1, v2, v3;
  Real xsheet,ysheet,zsheet;
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        //GetCylCoord(pco,rad,phi,z,i,j,ke+k);
        GetCylCoord(pco,rad,phi,z,i,j,ke);
        xsheet = rad - rref;
        zsheet = z;
        Real Ux,Uy,Uz,drho,dp,Al,kappa,Bl,Omega,omegak;
        Al = - 0.25 * (3*gm0/pow(rref,3.0) - (dslope+pslope)*(pslope-2)*p0_over_r0*pow(rref,pslope-2)) / sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        Omega = sqrt(gm0/pow(rref,3.0)+(dslope+pslope)*p0_over_r0*pow(rref,pslope-2));
        kappa = sqrt(4*Omega*(Omega+Al));
        Bl = Omega + Al;
        omegak = kappa * cos(theta);
        Ux = sqrt(amp)*cos(theta)*cos(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uy = -2*Bl/kappa * sqrt(amp)*sin(omegak*time+phi0)*cos(kr*xsheet+kz*zsheet);
        Uz = -kr/kz * Ux;
        drho = 0;//omegak / kz * rho0 / p0_over_r0 * pow(rref,dslope - pslope) * Uz;
        dp = 0;//drho * p0_over_r0 * pow(rad,pslope);
        VelProfileCyl(rad,phi,z,v1,v2,v3,0);
        prim(IDN,ke+k,j,i) = DenProfileCyl(rad,phi,z) + drho;
        prim(IM1,ke+k,j,i) = v1 + Ux;
        prim(IM2,ke+k,j,i) = v2 + Uy;
        prim(IM3,ke+k,j,i) = v3 + Uz;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,ke+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ke+k,j,i) + dp;
      }
    }
  }
}


void Rotation(MeshBlock *pmb,const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons){
  // MeshBlock *pmb = pmy_hydro_->pmy_block;
  // pmb is already assigned before calling this source term function
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real den = prim(IDN,k,j,i);
        Real rad, phi, z;
        Real KeplerianOmega;
        Real temp, tau;
        Real g = pmb->peos->GetGamma();
        Real poverr = PoverR(rad,phi,z);
        GetCylCoord(pmb->pcoord,rad,phi,z,i,j,k);
        // distance to the perturber
        // acceleration in r and phi
        Real p_over_rref = PoverR(rref,phi,z);
        Real vel0;
        vel0 = (dslope+pslope)*p_over_rref/(gm0/rref) + 1.0;
        vel0 = std::sqrt(gm0/rref)*std::sqrt(vel0);
        KeplerianOmega = vel0/rref;
        tau = 10 / KeplerianOmega;
        Real acc1,acc2;
        acc1 =  2 * KeplerianOmega * prim(IVY,k,j,i) + KeplerianOmega*KeplerianOmega*rad; //- KeplerianOmega*KeplerianOmega*companiondistance*(companionmass/(gm0+companionmass))*cos(phi) + companionmass*(companiondistance*cos(phi)-rad)/pow(dist,3.0);
        acc2 =  - 2 * KeplerianOmega * prim(IVX,k,j,i);// - companionmass*companiondistance*sin(phi)/pow(dist,3.0) + KeplerianOmega*KeplerianOmega*companiondistance*(companionmass/(gm0+companionmass))*sin(phi);

        Real src1 = dt*den*acc1;
        Real src2 = dt*den*acc2;
        // update
        cons(IM1,k,j,i) += src1;
        cons(IM2,k,j,i) += src2;
        // source = rho u. (-grad Phi_p) = rho ur g_r + rho us gs
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) +=
          src1*prim(IVX,k,j,i) + src2*prim(IVY,k,j,i);
      }
    }
  }
  return;
}





