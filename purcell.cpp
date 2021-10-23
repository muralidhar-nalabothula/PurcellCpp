// This code was written by Muralidhar Nalabothula 
// Please conside citing https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.045416 if you use the code.
// contact me in case of bug reports : muralidharrsvm7@gmail.com 
//Std headers
#include <iostream>
#include <complex>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <string>
//Eigen headers
#include "eigen/unsupported/Eigen/Polynomials"
//Integration headers
#include "quadpack/workspace.hpp"
//Symengine headers
#include <symengine/basic.h>
#include <symengine/add.h>
#include <symengine/complex_double.h>
#include <symengine/symbol.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/series.h>
#include <symengine/matrix.h>
#include <symengine/solve.h>
#include <symengine/eval_double.h>
#include <symengine/derivative.h>
#include "purcell.h"
// Defining Constants 

void Print_details(){
std::cout<<"---------------------- PurcellCpp ------------------------"<<std::endl;
std::cout<<"|                                                         |"<<std::endl;
std::cout<<"------------ Author: Muralidhar Nalabothula --------------"<<std::endl;
std::cout<<"----- Please cite the following work if you use the code -"<<std::endl;
std::cout<<"|                                                         |"<<std::endl;
std::cout<<"---- https://doi.org/10.1103/PhysRevB.102.045416 ---------"<<std::endl;
std::cout<<"|                                                         |"<<std::endl;
std::cout<<"--------- Contact: muralidharrsvm7@gmail.com -------------"<<std::endl;
std::cout<<"|                                                         |"<<std::endl;
std::cout<<"----------------------------------------------------------"<<std::endl;
}

const auto sym_zero = SymEngine::complex_double(c_zero);
const auto sym_iden = SymEngine::complex_double(c_one);
const auto sym_I = SymEngine::complex_double(I);
const auto sym_minus_one = SymEngine::integer(-1); 
const SymEngine::DenseMatrix sym_I_mat_4= SymEngine::DenseMatrix(4, 4,
                      { sym_iden,   sym_zero,   sym_zero,  sym_zero,
                        sym_zero,   sym_iden,   sym_zero,  sym_zero,
                        sym_zero,   sym_zero,   sym_iden,  sym_zero,
                        sym_zero,   sym_zero,   sym_zero,  sym_iden });



int cart_sym_2_num(char cart_string){ /* This Function converts string to number (x,y,z)--> (0,1,2) */
    if (cart_string == 'x'){
        return 0;
    }
    else if (cart_string == 'y'){
        return 1;
    }
    else if (cart_string == 'z'){
        return 2;
    }
    else {
        throw std::invalid_argument( "Wrong direction ! direction can only be in x,y,z" );
    }
     
}


/* Output Structure to store the interface function of the system
It contains the transfer matrix and the propagation angle*/
struct M_struct {SymEngine::DenseMatrix tranfser_mat; std::complex<double> prop_angle;}; 

/* Function to compute the transfer matrix for a interface*/
M_struct M_sub(std::complex<double> vl, struct_elements *interface, SymEngine::DenseMatrix &R_mat){ 
    /* Unrotated surface conductivity at the interface*/
    SymEngine::DenseMatrix sig_p = SymEngine::DenseMatrix(2, 2, 
        {SymEngine::complex_double(interface->sxx_p), SymEngine::complex_double(interface->sxy_p), 
         SymEngine::complex_double(interface->syx_p), SymEngine::complex_double(interface->syy_p)});
    SymEngine::DenseMatrix sig_pt(2,2), sig_pt_temp(2,2),R_mat_T(2,2);

    /* Rotating surface conductivity at the interface*/
    mul_dense_dense(sig_p, R_mat, sig_pt_temp); // sigma_p * R
    R_mat.transpose(R_mat_T); 
    mul_dense_dense(R_mat_T,sig_pt_temp , sig_pt); // sig_pt = R^T * sigma_p *R 

    /* sig_pt is the rotated conductivity tensor*/
    auto sxx = sig_pt.get(0,0);
    auto sxy = sig_pt.get(0,1);
    auto syx = sig_pt.get(1,0);
    auto syy = sig_pt.get(1,1);
    // nl is the refractive index (R.I) of incoming medium and nr is R.I of outgoing medium
    std::complex<double> nl = interface->nl;
    std::complex<double> nr = interface->nr;
    std::complex<double> c_vl = std::cos(std::complex<double>(vl));
    std::complex<double> s_vl = std::sin(std::complex<double>(vl));
    std::complex<double> vr = std::asin(std::complex<double>((nl*s_vl)*c_one/nr));
    std::complex<double> c_vr = std::cos(std::complex<double>(vr));
    std::complex<double> s_vr = std::sin(std::complex<double>(vr));

    // vl is the incoming angle and vr is the outgoing angle at interface
    auto f_p = SymEngine::complex_double(c_one + ((nl*c_vl)*c_one/(nr*c_vr)));
    auto f_n = SymEngine::complex_double(c_one - ((nl*c_vl)*c_one/(nr*c_vr)));
    auto g_p = SymEngine::complex_double((nl*c_one/nr) + (c_vl*c_one/c_vr));
    auto g_n = SymEngine::complex_double((nl*c_one/nr) - (c_vl*c_one/c_vr));
    auto b1 = mul(syy,SymEngine::complex_double(c_one/(nr*c_vr)));
    auto b2 = mul(syx,SymEngine::complex_double(c_vl*c_one/(nr*c_vr)));
    auto b3 = mul(sxy,SymEngine::complex_double(c_one/nr));
    auto b4 = mul(sxx,SymEngine::complex_double(c_vl*c_one/(nr)));
    SymEngine::DenseMatrix tot
        = SymEngine::DenseMatrix(4, 4,
                      {      sub(f_p,b1),    sub(sym_zero,b2),          sub(f_n,b1),                  b2,
                        sub(sym_zero,b3),        sub(g_p,b4) ,     sub(sym_zero,b3),         add(g_n,b4),
                            add(f_n ,b1),                  b2,         add(f_p ,b1),    sub(sym_zero,b2),
                        sub(sym_zero,b3),         sub(g_n,b4),     sub(sym_zero,b3),      add(g_p,b4) });
    tot.mul_scalar(SymEngine::complex_double(0.5*c_one),tot);
    return M_struct{tot,vr};
    }
/* Function for computing transfer matrix matrix for dielectric */
SymEngine::DenseMatrix M_di(std::complex<double> v_d, double k0, struct_elements *dielectric_ )
{   // vd is the propagation angle, k0 is the wave vector
    double d = dielectric_->thickness ; 
    std::complex<double> n_d  = dielectric_->nd;
    std::complex<double> k_dielc = k0*n_d;

    SymEngine::RCP<const SymEngine::ComplexDouble> ele = 
    SymEngine::complex_double(std::exp(I*std::complex<double>(4.0)*pi*k_dielc*d*std::cos(std::complex<double>(v_d))));

    SymEngine::DenseMatrix M_free = SymEngine::DenseMatrix(4, 4,
                      { ele,        sym_zero,   sym_zero,  sym_zero,
                        sym_zero,   ele,        sym_zero,  sym_zero,
                        sym_zero,   sym_zero,   sym_iden,  sym_zero,
                        sym_zero,   sym_zero,   sym_zero,  sym_iden });
    return M_free;
    }

/* Integrand of Green's function */
double green_integrand (double kp, /* kp is the normalized to k0 !!*/ \
double k0, /* Wave vector*/
std::vector<struct_elements> elements, /* System configuation*/
double z0 , /* dipole distance*/
std::string dir /* Direction of the dipole matrix*/
){  
    SymEngine::RCP<const SymEngine::Symbol> z = SymEngine::symbol("z"); // temp complex variable
    std::complex<double> kz = std::sqrt(std::complex<double>(1-kp*kp)); //kz

    auto KP = SymEngine::complex_double(kp+c_zero); //kp in symbolic rep
    auto KZ = SymEngine::complex_double(kz); // kz in symbolic rep

    SymEngine::RCP<const SymEngine::Basic> kx, ky, sphi, cphi;

    cphi = mul(add(z, pow(z, sym_minus_one)), SymEngine::complex_double(0.5 + 0.0*I)); //cos(phi) = kx/kp
    sphi = mul(sub(z, pow(z, sym_minus_one)), SymEngine::complex_double(-0.5*I));  // , sin(phi) = ky/kp

    kx = mul(cphi, SymEngine::real_double(kp));
    ky = mul(sphi, SymEngine::real_double(kp));

    SymEngine::DenseMatrix R_mat = SymEngine::DenseMatrix(2, 2,   /* Rotation matrix*/
                                    {cphi,      mul(sphi, sym_minus_one),
                                     sphi,      cphi});

    SymEngine::DenseMatrix G(3,3), Mss(4,4), Msp(4,4), Mps(4,4), Mpp(4,4);

    std::complex<double> vl = std::acos(kz); //incoming angle

    SymEngine::DenseMatrix M_tot =  sym_I_mat_4; // Starting the transfer matrx with a identity matrix

    // Computing the total transfer matrix of the system from the user input
    for (auto i : elements)
    { 
        if (i.element_type =="interface"){
            auto i_temp  = M_sub(vl+c_zero ,&i, R_mat);
            mul_dense_dense(i_temp.tranfser_mat,M_tot, M_tot);
            vl = i_temp.prop_angle;
        }

        if (i.element_type =="Dielectric"){
            auto i_temp  = M_di(vl+c_zero ,k0, &i);
            mul_dense_dense(i_temp,M_tot, M_tot);
        }

    }

    SymEngine::DenseMatrix M_21(2,2), M_22(2,2), g(2,2);

    M_tot.submatrix(M_21, 2, 0, 4,2); // M21 submatrix
    M_tot.submatrix(M_22, 2, 2, 4,4); // M_22 submatrix

    SymEngine::DenseMatrix M22_inv = SymEngine::DenseMatrix(2, 2,
                    {M_22.get(1,1) ,           mul(M_22.get(0,1), sym_minus_one) , 
                     mul(M_22.get(1,0),         sym_minus_one), M_22.get(0,0)        });  
                     // M22^-1 (Note: This is not the true inverse it is |M_22| * M_22^-1)
                     // The det will be divided at the end 

    SymEngine::RCP<const SymEngine::Basic> det_M22 = mul(M_22.det(), sym_minus_one); // det of M22

    det_M22 = expand(mul(z,det_M22)); //det(M22) * z , this z is multiplied, as we do transformation of variable exp(I*phi) = z 

    M22_inv.mul_matrix(M_21, g ) ; // M22^-1 * M21

    Mss = SymEngine::DenseMatrix(3, 3,
                      {mul(ky,ky),mul(mul(kx,ky), sym_minus_one),sym_zero,
                      mul(mul(kx,ky), sym_minus_one),mul(kx,kx),sym_zero,
                      sym_zero,sym_zero,sym_zero });

    Mpp = SymEngine::DenseMatrix(3, 3,
                      { mul(mul(kx,kx), sym_minus_one), mul(mul(kx,ky), sym_minus_one), mul(mul(kx,sym_minus_one),mul(mul(KP,KP),pow(KZ, sym_minus_one))),
                       mul(mul(kx,ky), sym_minus_one), mul(mul(ky,ky), sym_minus_one), mul(mul(ky,sym_minus_one),mul(mul(KP,KP) ,pow(KZ, sym_minus_one))),
                       mul(kx,mul(mul(KP,KP) , pow(KZ, sym_minus_one))) ,mul(ky,mul(mul(KP,KP) , pow(KZ, sym_minus_one))), mul(pow(KP,SymEngine::integer(4)) , (pow(mul(KZ,KZ), sym_minus_one))) });

    Msp = SymEngine::DenseMatrix(3, 3,
                      { mul(mul(kx,ky), sym_minus_one), mul(mul(ky,ky), sym_minus_one), mul(mul(ky,sym_minus_one),mul(mul(KP,KP),pow(KZ, sym_minus_one))), 
                      mul(kx,kx), mul(kx,ky), mul(kx,mul(mul(KP,KP),pow(KZ, sym_minus_one))), 
                      sym_zero,sym_zero,sym_zero});
    Mps = SymEngine::DenseMatrix(3, 3,
                      { mul(kx,ky) , mul(mul(kx,kx), sym_minus_one) ,sym_zero, 
                      mul(ky,ky), mul(mul(kx,ky), sym_minus_one) ,sym_zero,
                       mul(mul(ky,sym_minus_one),mul(mul(KP,KP),pow(KZ, sym_minus_one))), mul(kx,mul(mul(KP,KP),pow(KZ, sym_minus_one))) ,sym_zero});

    Mss.mul_scalar(SymEngine::complex_double(c_one/(kz*kp*kp) ) ,Mss ) ;
    Mpp.mul_scalar(SymEngine::complex_double(kz * c_one/(kp*kp) ), Mpp )  ;
    Msp.mul_scalar(SymEngine::complex_double(c_one/(kp*kp) ) , Msp)  ;
    Mps.mul_scalar(SymEngine::complex_double(c_one/(kp*kp) ) , Mps)  ;

    Mss.mul_scalar(g.get(0,0),Mss);
    Msp.mul_scalar(g.get(0,1),Msp);
    Mps.mul_scalar(g.get(1,0),Mps);
    Mpp.mul_scalar(g.get(1,1),Mpp);
    Mss.add_matrix(Msp,G);
    Mps.add_matrix(G,G);
    Mpp.add_matrix(G,G);

    //std::cout<<*SymEngine::rcp_static_cast<const SymEngine::Add>(add(det_M22, SymEngine::complex_double(c_one*pow(10,-100)) ))<<std::endl;
    auto det_M22_dict = SymEngine::rcp_static_cast<const SymEngine::Add>(add(det_M22, SymEngine::complex_double(c_one*1e-50) )) ->get_dict(); // getting coeff of the det polynomial
    SymEngine::RCP<const SymEngine::Basic> gg_num = G.get(cart_sym_2_num(dir[0]),cart_sym_2_num(dir[1])); // get the G(x,y) element
    std::vector<signed long int> det_terms;
    std::vector<std::complex<double>> det_coeffs;
    det_M22 = mul(sym_zero,z);
    for(auto const& imap: det_M22_dict){ /* Filtering out co-eff >1E-8 and set others to 0. This is done due to machine prescision */
        std::complex<double> coeffs_temp = SymEngine::rcp_static_cast<const SymEngine::ComplexDouble>(add(imap.second,sym_zero))->i ;
        if (std::abs(coeffs_temp)>std::pow(10,-8)){
            det_terms.push_back((SymEngine::rcp_static_cast<const SymEngine::Integer>(SymEngine::rcp_static_cast<const SymEngine::Pow>(pow(imap.first,SymEngine::integer(2)))->get_exp())->as_int())/2);
            det_coeffs.push_back(coeffs_temp);     
            det_M22 = add(det_M22,mul(imap.second,imap.first));
        }
    }
    //std::cout<<*det_M22<<std::endl;
    auto minmax = std::minmax_element(det_terms.begin(), det_terms.end()); // finding min and max degrees of the det polynomial
    std::vector<std::complex<double>> poly_coeff(*minmax.second-*minmax.first+1, c_zero);
    // OPENMP?
    for (int i=0; i<det_terms.size(); ++i){
        poly_coeff[det_terms[i]-*minmax.first] = det_coeffs[i]; // storing all the elements in the vector
    }
    Eigen::Map<Eigen::VectorXcd> coeff(poly_coeff.data(),poly_coeff.size()); //mapping to eigen vector for finding the roots of the polynomial
    std::vector<std::complex<double>> roots;
    if (poly_coeff.size() > 1){
        Eigen::PolynomialSolver<std::complex<double>, Eigen::Dynamic> solver; // initialization of solver
        solver.compute(coeff); // solving for roots
        const Eigen::PolynomialSolver<std::complex<double>, Eigen::Dynamic>::RootsType &root = solver.roots(); // get the roots
        roots = std::vector<std::complex<double>>(root.data(), root.data() + root.size());
        if (root.size() != roots.size()){
            std::cout<<"Warning, some roots are missing"<<std::endl;
        }
    }
    
    auto integrand = mul(gg_num,pow(det_M22,sym_minus_one)); //Numerical evaulation of G(x,y) integrand 
    SymEngine::RCP<const SymEngine::Basic> der_det_M22 = diff(det_M22,z); // computing the derivative for finding the residue.
    // do a series expansion to find number of terms required
    SymEngine::RCP<const SymEngine::SeriesCoeffInterface> residue_series = SymEngine::series(mul(z,integrand), z,64); // series expanded to 64 terms
    auto unordered_dict = residue_series->as_dict(); // get coeff
    std::map<int, SymEngine::RCP<const SymEngine::Basic>> ordered_dict(unordered_dict.begin(), unordered_dict.end()); // order the coedd according to degree
    int series_order = ordered_dict.rbegin()->first; // sorting 
    if (series_order<0){ 
        residue_series = SymEngine::series(mul(z,integrand), z,200);
        unordered_dict = residue_series->as_dict();
        std::map<int, SymEngine::RCP<const SymEngine::Basic>> ordered_dict(unordered_dict.begin(), unordered_dict.end());
        series_order = ordered_dict.rbegin()->first;
    }
    SymEngine::RCP<const SymEngine::Basic> residue0 =  residue_series ->get_coeff(0);  // not get_coeff(-1) because we multiplied series with z and take the constant
    std::complex<double> residue = SymEngine::rcp_static_cast<const SymEngine::ComplexDouble>(add(residue0,sym_zero))->i; // z= 0 is a always pole
    for (auto i : roots){
         if (std::abs(i) < 1 && std::abs(i) > std::pow(10,-precision) ){
            auto pole = SymEngine::complex_double(i+c_zero); // only |z|>0 and |z|<1 residues contribute to integral
            std::complex<double> diff_at_pole = SymEngine::rcp_static_cast<const SymEngine::ComplexDouble>(expand(der_det_M22->subs({{z, pole}})))->i ;
            int counter = 1; // counts how many times det is differentiated with z
            while (std::abs(diff_at_pole)<std::pow(10,-precision)){ // if the pole  is higher order // be careful with prec, it can lead to infite loop 
               counter = counter+1 ;
               der_det_M22 = diff(der_det_M22,z);
               diff_at_pole = (SymEngine::rcp_static_cast<const SymEngine::ComplexDouble>(expand(der_det_M22->subs({{z, pole}})))->i)/(std::tgamma(counter+1)) ;
            }
            residue = residue + (SymEngine::rcp_static_cast<const SymEngine::ComplexDouble>(expand(gg_num->subs({{z, pole}})))->i)/diff_at_pole ;
         }
    // }
    }
    //auto integrand_p = integrand->subs({{z, add(z,p)}}); 
    return (residue*kp*std::exp(I*std::complex<double>(2.0)*kz*z0*std::complex<double>(2.0)*pi*k0)).real() ; //integrand of G(x,y) without phi integral
    ;
}


struct G_integrand_pars { double k0 ; // for performing the integrals 
    std::vector<struct_elements> elements;
    double z0 ;
    std::string dir; };

double integrand_rad (double kp, void * G_params) // for performing the radiative integrals 
  {
    struct G_integrand_pars * params = (struct G_integrand_pars *)G_params;
    double k0 = params->k0;
    std::vector<struct_elements> elements = params->elements;
    double z0 = params->z0;
    std::string dir= params->dir;

    return  green_integrand(kp,k0, elements, z0,dir);
  }

double integrand_non_rad (double t, void * G_params) // for performing the non-radiative integrals 
  { double kp = 1.0f/t;
    struct G_integrand_pars * params = (struct G_integrand_pars *)G_params;
    double k0 = params->k0;
    std::vector<struct_elements> elements = params->elements;
    double z0 = params->z0;
    std::string dir= params->dir;

    return  kp*kp*green_integrand(kp,k0, elements, z0,dir);
  }

Decay_rates Purcell(double k0, /* Computing the rad and non-rad decay rates by perfoming the kp integral*/
                    std::vector<struct_elements> elements,
                    double z0 , 
                    std::string dir, 
                    size_t limit,size_t m_deg){  // sets (2m+1)-point Gauss-Kronrod

    std::cout<<"                                                             "<<std::endl;
    std::cout<<"********  Computing " + dir + " component decay rate ********"<<std::endl;
    double err_r, err_nr, rad, non_rad;
    Workspace<double> Work(limit, m_deg);
    struct G_integrand_pars params = { k0, elements, z0,dir };    
	Function<double, void> F_r(integrand_rad, &params);
    Function<double, void> F_nr(integrand_non_rad, &params);
	Work.qag(F_r, 0.01, 1.001, 1e-6, 1e-8, rad, err_r);
    Work.qag(F_nr, 1e-10, 1.0/1.001, 1e-6, 1e-8, non_rad, err_nr);
    //std::cout<<rad<<"  "<<err_nr/non_rad<<std::endl;
    return Decay_rates { 1+ 1.5f*(rad + non_rad), 1+ 1.5f*rad, 1.5f * non_rad, err_r/rad, err_nr/non_rad}; // returns , total decay rate, rad decay, non rad decay, numerical integration erros for rad and non rad
	
}



