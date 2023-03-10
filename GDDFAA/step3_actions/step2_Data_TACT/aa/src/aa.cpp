#include "aa.h"
#include "spherical_aa.h"
#include <limits>

VecDoub planar_sphr_actions(const VecDoub &x,Potential_JS* Pot){
    PlanarAxisymPotential PP(Pot);
    Actions_Spherical AS(&PP);
    return AS.actions(x);
}
VecDoub planar_sphr_angles(const VecDoub &x,Potential_JS* Pot){
    PlanarAxisymPotential PP(Pot);
    Actions_Spherical AS(&PP);
    return AS.angles(x);
}

/*!
    Checks input for actions
*/
int action_check(const VecDoub &x, VecDoub &acts, Potential_JS *Pot){

	// printf("action_check() called.\n"); //gjy add
    // printf("x[0] = %f, acts[0] = %f\n", x[0], acts[0]); //gjy add
    // printf("here1\n"); //gjy add
    // Check length of x
    if(x.size()<6){
        std::cerr<<"action_check(): Must pass 6D vector to action routine\n";
        for(unsigned i=0;i<3;++i) acts[i]=std::numeric_limits<double>::infinity();
        return 1;
    }
    // printf("here2\n"); //gjy add
    // Check if unbound
    //can not MultipoleExpansion //?? //gjy note
    if(Pot->H(x)>Pot->Phi({1.e7*(x[0]==0.?1.:x[0]),1.e7*(x[1]==0.?1.:x[1]),1.e7*(x[2]==0.?1.:x[2])})){
        std::cerr<<"action_check(): Orbit passed to action routine is unbound\n";
        for(unsigned i=0;i<3;++i)
            acts[i]=std::numeric_limits<double>::infinity();
        return 1;
    }

    VecDoub Polar = conv::CartesianToPolar(x);
    // std::cout<<"[aa.cpp: action_check() Polar = "<<Polar[0]<<" "<<Polar[1]<<" "<<Polar[2]<<"  "<<Polar[3]<<" "<<Polar[4]<<" "<<Polar[5]<<"]\n"; //gjy add
    if(fabs(Polar[0])<SMALL){ //SMALL=1e-5
        std::cout<<"fabs(Polar[0])<SMALL"<<"\n"; //gjy add
        acts[0]=0.;acts[1]=0.;acts[2]=0.;
        return 1;
    }

    acts[1]=Polar[0]*Polar[4];
    if(fabs(x[2])<SMALL and fabs(x[5])<SMALL){
        acts[2]=0.;
        std::cout<<"fabs(x[2])<SMALL and fabs(x[5])<SMALL"<<"\n"; //gjy add
        // std::cout<<"actions_check acts = "<<acts[0]<<" "<<acts[1]<<" "<<acts[2]<<"  "<<acts[3]<<"\n"; //gjy add
        if(fabs(Polar[3])<TINY and fabs(Polar[0])<0.01*SMALL){
            acts[0]=0.;
            return 1;
        }
        else{
            acts=planar_sphr_actions(x,Pot);
            acts[2]=0.;
            return 1;
        }
    }
    return 0;
}
/*!
    Checks input for angles
*/
int angle_check(const VecDoub &x, VecDoub &angs, Potential_JS *Pot){
    // Check length of x
    if(x.size()<6){
        std::cerr<<"angle_check(): Must pass 6D vector to angle routine\n";
        for(unsigned i=0;i<6;++i) angs[i]=std::numeric_limits<double>::infinity();
        return 1;
    }
    // Check if unbound
    if(Pot->H(x)>Pot->Phi({1.e7*(x[0]==0.?1.:x[0]),1.e7*(x[1]==0.?1.:x[1]),1.e7*(x[2]==0.?1.:x[2])})){
        std::cerr<<"angle_check(): Orbit passed to angle routine is unbound\n";
        for(unsigned i=0;i<6;++i)
            angs[i]=std::numeric_limits<double>::infinity();
        return 1;
    }
    VecDoub Polar = conv::CartesianToPolar(x);
    if(fabs(Polar[0])<SMALL){
        angs[0]=0.;angs[1]=0.;angs[2]=0.;
        angs[3]=0.;angs[4]=0.;angs[5]=0.;
        return 1;
    }
    if(fabs(x[2])<SMALL and fabs(x[5])<SMALL){
        angs[2]=0.;
        angs[5]=0.;
        if(fabs(Polar[3])<SMALL){
            angs[0]=0.;
            angs[3]=0.;
            angs[1]=Polar[1];
            angs[4]=Polar[4]/Polar[0];
            return 1;
        }
        else{
            angs=planar_sphr_angles(x,Pot);
            angs[2]=0.;
            angs[5]=0.;
            return 1;
        }
    }
    return 0;
    }
// ============================================================================
VecDoub integrate_a_bit(VecDoub x, Potential_JS *Pot){
    Orbit Orb(Pot);
    double torb = Pot->torb(x);
    return Orb.integrate(x,torb/1000.,torb/1000.,false);
}
