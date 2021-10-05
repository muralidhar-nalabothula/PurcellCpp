// This code was written by Muralidhar Nalabothula 
// Please conside citing https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.045416 if you use the code.
// contact me in case of bug reports : muralidharrsvm7@gmail.com 
#pragma once
//Std headers
#include <iostream>
#include <complex>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <string>
#include <vector>

void Print_details();
// Defining Constants 
const std::complex<double> I(0.0f,1.0f);
const std::complex<double> c_zero(0.0f,0.0f);
const std::complex<double> c_one(1.0f,0.0f);
const double pi = 3.14159265358979323846 ;
const int precision =4;
const double eta_0 = 376.991118430775;

struct struct_elements { std::string element_type; /* Structure to store input data of the complete system*/
                         std::complex<double> nl; std::complex<double> nr; 
                         std::complex<double> sxx_p; std::complex<double> sxy_p; 
                         std::complex<double> syx_p; std::complex<double> syy_p;
                         std::complex<double> nd; double thickness; 
                         virtual ~struct_elements(){}};

struct Interface : struct_elements { /* Structure to store sinlge inteface input data of the system*/
                                     Interface(std::complex<double> nl_, std::complex<double> nr_, 
                                     std::complex<double> sxx_p_, std::complex<double> sxy_p_, 
                                     std::complex<double> syx_p_, std::complex<double> syy_p_)
                                     {
                                        element_type = "interface";
                                        nl = nl_;
                                        nr = nr_;
                                        sxx_p = sxx_p_ * eta_0;
                                        sxy_p = sxy_p_ * eta_0;
                                        syx_p = syx_p_ * eta_0;
                                        syy_p = syy_p_ * eta_0;
                                     }

                                      ~Interface(){}
                                     
                                     };

struct Dielectric_slab : struct_elements { /* Structure to store dielectric input data of the system*/
                                           Dielectric_slab(std::complex<double> nd_, double thickness_)
                                           {
                                               element_type = "Dielectric";
                                               nd = nd_;
                                               thickness = thickness_;
                                           } 
                                           ~Dielectric_slab(){}
                                           };

typedef std::vector<struct_elements> System;

/* Output Structure to store computed decay rates (and their relative errors) of the system
It contains transface matrix and propagation angle*/
struct Decay_rates {double total; double rad; double non_rad; double err_r; double err_nr ;};

Decay_rates Purcell(double k0, /* Computing the rad and non rad decay rates by perfoming the kp integral*/
                    std::vector<struct_elements> elements,
                    double z0 , 
                    std::string dir, 
                    size_t limit = 128,size_t m_deg = 10);

