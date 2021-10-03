// This code is written by Muralidhar Nalabothula 
// Please conside citing https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.045416 if you use the code.
// contact me in case of bug reports : muralidharrsvm7@gmail.com 
#include "purcell.h"

int main(){
    Print_details();
    double dipole_distance = 1.6*1e-9;
    System elements = { Interface(1.0,  1.85, 0.0, 0.0, 0.0, 0.0), Dielectric_slab(1.85, 5*1e-9) ,Interface(1.85,  1.5, 4.915721315948482e-07+3.596760413126592e-05*I, 0.0, 0.0,  1.0021936438284495e-07+8.464171939906124e-06*I) };
    //                  Air-dielectric interface                     dielectric slab               dielectric, 2D material and substate interface
    auto ppp = Purcell(1145476.66561009,elements, dipole_distance ,"xx");//}
    std::cout<<"                                                                                     "<<std::endl;
    std::cout<<"-------------------------------------------------------------------------------------"<<std::endl;
    std::cout<< "Total decay rate (Normalized to free decay rate) =                    " <<ppp.total  <<std::endl;
    std::cout<< "Total radiative decay rate (Normalized to free decay rate) =          " <<ppp.rad  <<std::endl;
    std::cout<< "Total non-radiative rate (Normalized to free decay rate) =            " <<ppp.non_rad  <<std::endl;
    std::cout<< "Relative numerical integration error for radiative calculation =       " <<ppp.err_r  <<std::endl;
    std::cout<< "Relative numerical integration error for non-radiative calculation =   " <<ppp.err_nr  <<std::endl;
    std::cout<<"-------------------------------------------------------------------------------------"<<std::endl;
    return 0;
    return 0;
}
