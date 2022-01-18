//
//  main.cpp
//  simple test
//
//  Created by Sajjad Sadeghizadeh on 7/26/19.
//  Copyright © 2019 Sajjad Sadeghizadeh. All rights reserved.
//


#include <iostream>
#include <string>
#include <vector>
#include "OpenMM.h"

using namespace std ;

const double step_fs = 10;
const double report_fs = 40 ;
const double totaltime_ps = 50 ;
static const bool wantenergy = true ;

// if k1>>k2 --> Kelvin-Voigt with stiffness k2 parallel to dashpot
const double k1 = 20;
const double k2 = 30;
const double eta = 100;
double k_eff = (k1*k2)/(k1+k2);
double free_length = 3;
double eta1=1000;
double eta2=10;
double tau1 = eta1/k1;
double tau2 = eta2/k1;



struct Atomtype
{
    double mass , charge , vdwradius_nm , vdwenergy_KJ;
};

Atomtype atomtypes[] = { 40 , 0 , 0.01 , 0.004 ,
    1000000000000 , 0 , 0.01 , 0.004
};


struct Bondtype
{
    double length_nm , k_KJperNm2 ;
};

Bondtype bondtypes[] = {free_length , k1 } ;

struct Atominfo
{
    int type ; const char* pdb ; double initpos_nm[3] ; double pos_nm[3]; double vel_nmpers[3];
};

const int endlist = -1;

Atominfo atoms[] = {0 , " A " , {0.5,0,0} , {0,0,0} ,{0,0,0},
    1 , " B " , {0 , 0 , 0} , {0,0,0} ,{0,0,0},
    endlist
} ;

double dist[3] = {-1,-1,-1} ;
double last_f_length;


static struct {int type ; int bond_atoms[2];} bonds[] = { 0 , 0 , 1 , endlist} ;

static void writepdb(int framenum , double time_ps , double energy_KJ , double potential, const Atominfo atoms[])
{
    /* printf("MODEL     %d\n", framenum);
     printf("REMARK 250 time=%.3f ps; energy=%.3f kJ/mole; potential=%.3f kJ/mole\n",
     time_ps, energy_KJ , potential);
     for (int n=0; atoms[n].type != endlist; ++n)
     printf("ATOM  %5d %4s ETH     1    %8.3f%8.3f%8.3f  1.00  0.00\n",
     n+1, atoms[n].pdb,
     atoms[n].pos_nm[0] *OpenMM::AngstromsPerNm, atoms[n].pos_nm[1] * OpenMM::AngstromsPerNm, atoms[n].pos_nm[2] * OpenMM::AngstromsPerNm);
     printf("ENDMDL\n");*/
    std::cout << time_ps << ' ' << atoms[0].pos_nm[0] * OpenMM::AngstromsPerNm << '\n';
    
}

OpenMM::HarmonicBondForce& bondstretch = *new OpenMM::HarmonicBondForce();
OpenMM::CustomBondForce& mybond = *new OpenMM::CustomBondForce("step(4 -r) *0.5*k*(r-r0)^2");
OpenMM::NonbondedForce& nonbond = *new OpenMM::NonbondedForce();


struct openmm_data ;

static openmm_data* initializer(const Atominfo atoms[] , double stepsize_fs , std::string& platformName) ;
static void step(openmm_data* , int stepnum , Atominfo atoms[]) ;
static void getstate(openmm_data* , bool wantEnergy , double& time , double& energy, double& potential , Atominfo atoms[]);
static void cheap_getstate(openmm_data* ,Atominfo atoms[]);
static void terminateopenmm(openmm_data*) ;
static void Harmonic_update(openmm_data*  , OpenMM::HarmonicBondForce& bond  , double distance[3], int Update_Method) ;

bool flag = true ;
int main(){

    try {
        std::string   platformName;
        openmm_data* omm = initializer( atoms , step_fs , platformName) ;
        
        
        //printf("REMARK  Using OpenMM platform %s\n", platformName.c_str());
        
        const int NumSilentSteps = (int)(report_fs / step_fs + 0.5);
        //const int NumSteps = (int)(totaltime_ps / step_fs + 0.5);
        for (int frame=1; ; ++frame) {
            double time, energy , potential;
            getstate(omm, wantenergy, time, energy, potential, atoms);
            writepdb(frame, time, energy,potential, atoms);
            
            
            
            if (time >= totaltime_ps)
                break;
            
            step(omm, NumSilentSteps , atoms);
        }
        
        terminateopenmm(omm);
        
        return 0 ;
    }
    
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
    
}




using OpenMM::Vec3;


struct openmm_data {
    openmm_data() : system(0), context(0), integrator(0) {}
    ~openmm_data() {delete context; delete integrator; delete system;}
    OpenMM::System*         system;
    OpenMM::Integrator*     integrator;
    OpenMM::Context*  context;
};


static openmm_data* initializer( const Atominfo atoms[] , double stepsize_fs , std::string& platformName)
{
    OpenMM::Platform::loadPluginsFromDirectory
    (OpenMM::Platform::getDefaultPluginsDirectory());
    
    openmm_data* omm = new openmm_data();
    
    OpenMM::System& system = *(omm->system = new OpenMM::System()) ;
    
    //system.addForce(&bondstretch);
    
    //system.addForce(&nonbond);
    
    
    std::vector<Vec3> initpos_nm;
    for (int n=0; atoms[n].type != endlist; ++n) {
        const Atomtype& atype = atomtypes[atoms[n].type];
        system.addParticle(atype.mass);
        
        // Convert the initial position to nm and append to the array.
        const Vec3 pos_nm(atoms[n].initpos_nm[0] ,
                          atoms[n].initpos_nm[1] ,
                          atoms[n].initpos_nm[2] );
        initpos_nm.push_back(pos_nm);
        
        
        //nonbond.addParticle(atype.charge, atype.vdwradius_nm* OpenMM::SigmaPerVdwRadius, atype.vdwenergy_KJ);
        
    }
    
    mybond.addGlobalParameter("k", 10);
    mybond.addGlobalParameter("r0", 2);
    
    
    std::vector< std::pair<int,int> > bondPairs;
    for (int i=0; bonds[i].type != endlist; ++i) {
        const int*      atom = bonds[i].bond_atoms;
        const Bondtype& bond = bondtypes[bonds[i].type];
        
        bondstretch.addBond(atom[0], atom[1], bond.length_nm , bond.k_KJperNm2);
        mybond.addBond(atom[0], atom[1]);
        
        
        bondPairs.push_back(std::make_pair(atom[0], atom[1]));
    }
    
    system.addForce(&mybond);
    
    
    omm->integrator = new OpenMM::VerletIntegrator(stepsize_fs * OpenMM::PsPerFs);
    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator);
    omm->context->setPositions(initpos_nm);
    
    platformName = omm->context->getPlatform().getName();
    return omm;
}





static void
getstate(openmm_data* omm, bool wantEnergy,
         double& timeInPs, double& energy_KJ, double& potential_KJ ,
         Atominfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    if (wantEnergy) {
        infoMask += OpenMM::State::Velocities;
        infoMask += OpenMM::State::Energy;
    }
    
    
    const OpenMM::State state = omm->context->getState(infoMask);
    timeInPs = state.getTime(); // OpenMM time is in ps already
    
    
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    for (int i=0; i < (int)positionsInNm.size(); ++i)
        for (int j=0; j < 3; ++j)
            atoms[i].pos_nm[j] = positionsInNm[i][j] ;
    
    const std::vector<Vec3>& velInNmpers = state.getVelocities();
    for (int i=0; i < (int)velInNmpers.size(); ++i)
        for (int j=0; j < 3; ++j)
            atoms[i].vel_nmpers[j] = velInNmpers[i][j] ;
    
    
    energy_KJ = 0;
    potential_KJ = 0 ;
    if (wantEnergy)
        potential_KJ = state.getPotentialEnergy() ;
    energy_KJ = (state.getPotentialEnergy() + state.getKineticEnergy());
}


static void
cheap_getstate(openmm_data* omm,
               Atominfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    
    
    const OpenMM::State state = omm->context->getState(infoMask);
    //timeInPs = state.getTime(); // OpenMM time is in ps already
    
    
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    for (int i=0; i < (int)positionsInNm.size(); ++i)
        for (int j=0; j < 3; ++j)
            atoms[i].pos_nm[j] = positionsInNm[i][j] ;
}



static void Harmonic_update(openmm_data* omm , OpenMM::HarmonicBondForce& bond , double distance[3], int Update_Method )
{
    
    const int a = bond.getNumBonds();
    int atom1, atom2 ;
    double length, stiffness, damp_ps , viscosity;
    //damp_ps = 4 ;
    //viscosity = 10 ;
    for(int i=0 ; i<a ; ++i)
    {
        
        bond.getBondParameters(i, atom1, atom2, length, stiffness);
        
        last_f_length = length;
        
        switch (Update_Method) {
                /** No viscoelasticity */
            case 0:
                break;
                /** Kelvin-Voigt */
            case 1:
                length = free_length -  ((distance[1] - distance[0])) * ((eta/k_eff)/(step_fs*OpenMM::PsPerFs)) ;
                break;
                /** Standard linear model */
            case 2:
                length = free_length*(step_fs*OpenMM::PsPerFs / (step_fs*OpenMM::PsPerFs + eta/(k1+k2))) + length*((eta/(k1+k2)) / (step_fs*OpenMM::PsPerFs + eta/(k1+k2))) - ((distance[1] - distance[0]))*((eta/(k1+k2)) / (step_fs*OpenMM::PsPerFs + eta/(k1+k2))) * (k1/k2)  ;
                break;
            case 3:
                length = free_length*(step_fs*OpenMM::PsPerFs/(tau1+tau2+step_fs*OpenMM::PsPerFs)) + length* ( (tau1+tau2)/(tau1+tau2+step_fs*OpenMM::PsPerFs)) + distance[2] * ( (tau1+ step_fs*OpenMM::PsPerFs - (tau1*tau2/step_fs*OpenMM::PsPerFs))/ (tau1+tau2+step_fs*OpenMM::PsPerFs) ) + distance[1] * ( (-tau1 + (2*tau1*tau2/step_fs*OpenMM::PsPerFs))/ (tau1+tau2+step_fs*OpenMM::PsPerFs)) + distance[0] * ( (-(tau1*tau2/step_fs*OpenMM::PsPerFs))/ (tau1+tau2+step_fs*OpenMM::PsPerFs) ) ;
                break;
            case 4:
                length = free_length*(step_fs*OpenMM::PsPerFs/(tau1+step_fs*OpenMM::PsPerFs))+ length* ( (tau1)/(tau1+step_fs*OpenMM::PsPerFs)) - ((distance[2] - distance[1])) * ( (tau2+ (tau1*tau2/step_fs*OpenMM::PsPerFs))/ (tau1+step_fs*OpenMM::PsPerFs) ) + ((distance[1] - distance[0])) * ( ( (tau1*tau2/step_fs*OpenMM::PsPerFs))/ (tau1+step_fs*OpenMM::PsPerFs) ) ;
                
                
                
                
            default:
                break;
        }
        
        
        
        bond.setBondParameters(i, atom1, atom2, length , stiffness);
    }
    bond.updateParametersInContext(*omm->context);
    
    
}







static void
step(openmm_data* omm, int numSteps , Atominfo atoms[]) {
    double bond_distance;
    
//    for(int i=0 ; i<numSteps ; i++)
//    {
//        cheap_getstate(omm, atoms);
//
//        bond_distance= sqrt( (atoms[0].pos_nm[0] - atoms[1].pos_nm[0]) * (atoms[0].pos_nm[0] - atoms[1].pos_nm[0]) + (atoms[0].pos_nm[1] - atoms[1].pos_nm[1])*(atoms[0].pos_nm[1] - atoms[1].pos_nm[1]) + (atoms[0].pos_nm[2] - atoms[1].pos_nm[2])*(atoms[0].pos_nm[2] - atoms[1].pos_nm[2])  ) ;
//
//        dist[0] = dist[1];
//        dist[1] = dist[2];
//        dist[2] = bond_distance;
//
//        if(dist[0] != -1)
//        {
//            Harmonic_update(omm , bondstretch , dist , 0) ;
//        }
//
//        //dist.push_back(distance);
//        //Harmonic_update(omm , bondstretch , dist , 1 , frame) ;
//        omm->integrator->step(1);
//
//    }
    
    omm->integrator->step(numSteps);
}




static void
terminateopenmm(openmm_data* omm) {
    delete omm;
}

