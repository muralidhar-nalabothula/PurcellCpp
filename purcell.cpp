//Std headers
#include <iostream>
#include <complex>
#include <cmath>
#include <chrono>
#include <stdexcept>
#include <string>
//Eigen headers
#include <unsupported/Eigen/Polynomials> 
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

// Defining Constants (!! Donot use these varibles any where)
const std::complex<double> I(0.0f,1.0f);
const std::complex<double> c_zero(0.0f,0.0f);
const std::complex<double> c_one(1.0f,0.0f);
const double pi = 3.14159265358979323846 ;
const int precision =4;
const auto sym_zero = SymEngine::complex_double(c_zero);
const auto sym_iden = SymEngine::complex_double(c_one);
const auto sym_I = SymEngine::complex_double(I);
const auto sym_minus_one = SymEngine::integer(-1); 

int cart_sym_2_num(char cart_string){
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

struct M_struct {SymEngine::DenseMatrix tranfser_mat; std::complex<double> prop_angle;};
struct Decay_rates {double total; double rad; double non_rad; double err_r; double err_nr ;};
M_struct M_sub(std::complex<double> nr, std::complex<double> nl, std::complex<double> vl, SymEngine::DenseMatrix &sig_un){
    auto sxx = sig_un.get(0,0);
    auto sxy = sig_un.get(0,1);
    auto syx = sig_un.get(1,0);
    auto syy = sig_un.get(1,1);
    std::complex<double> c_vl = std::cos(std::complex<double>(vl));
    std::complex<double> s_vl = std::sin(std::complex<double>(vl));
    std::complex<double> vr = std::asin(std::complex<double>((nl*s_vl)*c_one/nr));
    std::complex<double> c_vr = std::cos(std::complex<double>(vr));
    std::complex<double> s_vr = std::sin(std::complex<double>(vr));
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

SymEngine::DenseMatrix M_di(std::complex<double> v_d, double k0, double d, std::complex<double> n_d ){
    std::complex<double> k_dielc = k0*n_d;
    SymEngine::RCP<const SymEngine::ComplexDouble> ele = SymEngine::complex_double(std::exp(I*std::complex<double>(4.0)*pi*k_dielc*d*std::cos(std::complex<double>(v_d))));
    SymEngine::DenseMatrix M_free ;
    SymEngine::vec_basic grid_ele{ele, ele, sym_iden,sym_iden};
    diag(M_free, grid_ele, 4);
    return M_free;
    }

/* Integrand of Green's function */
double green_integrand (double kp, /* kp is normalized to k0 !!*/ \
double k0, \
std::complex<double> sxx_p,  \
std::complex<double> sxy_p, \
std::complex<double> syx_p, \
std::complex<double> syy_p, \
double z0 , \
std::string dir
){  
    SymEngine::RCP<const SymEngine::Symbol> z = SymEngine::symbol("z"); // temp complex variable
    std::complex<double> kz = std::sqrt(std::complex<double>(1-kp*kp));
    auto KP = SymEngine::complex_double(kp+c_zero);
    auto KZ = SymEngine::complex_double(kz);
    SymEngine::DenseMatrix sig_p = SymEngine::DenseMatrix(2, 2, \
        {SymEngine::complex_double(sxx_p), SymEngine::complex_double(sxy_p), \
         SymEngine::complex_double(syx_p), SymEngine::complex_double(syy_p)});
    SymEngine::RCP<const SymEngine::Basic> kx, ky, sphi, cphi;
    cphi = mul(add(z, pow(z, sym_minus_one)), SymEngine::complex_double(0.5 + 0.0*I));
    sphi = mul(sub(z, pow(z, sym_minus_one)), SymEngine::complex_double(-0.5*I));
    kx = mul(cphi, SymEngine::real_double(kp));
    ky = mul(sphi, SymEngine::real_double(kp));
    SymEngine::DenseMatrix R_mat = SymEngine::DenseMatrix(2, 2, \
        {cphi,mul(sphi, sym_minus_one),sphi,cphi});
    SymEngine::DenseMatrix sig_pt(2,2), sig_pt_temp(2,2),R_mat_T(2,2), G(3,3), \
    Mss(4,4), Msp(4,4), Mps(4,4), Mpp(4,4);
    mul_dense_dense(sig_p, R_mat, sig_pt_temp);
    R_mat.transpose(R_mat_T);
    mul_dense_dense(R_mat_T,sig_pt_temp , sig_pt);
    std::complex<double> in_ang = std::acos(kz);
    SymEngine::DenseMatrix M_tot = M_sub(1.5 + c_zero,c_one,in_ang,sig_pt).tranfser_mat;
    SymEngine::DenseMatrix M_21(2,2), M_22(2,2), g(2,2);
    M_tot.submatrix(M_21, 2, 0, 4,2);
    M_tot.submatrix(M_22, 2, 2, 4,4);
    SymEngine::DenseMatrix M22_inv = SymEngine::DenseMatrix(2, 2,
                    {M_22.get(1,1) , mul(M_22.get(0,1), sym_minus_one) , 
                    mul(M_22.get(1,0), sym_minus_one), M_22.get(0,0)});
    SymEngine::RCP<const SymEngine::Basic> det_M22 = mul(M_22.det(), sym_minus_one);
    det_M22 = expand(mul(z,det_M22));
    M22_inv.mul_matrix(M_21, g ) ;
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

    auto det_M22_dict = SymEngine::rcp_static_cast<const SymEngine::Add>(det_M22) ->get_dict();
    SymEngine::RCP<const SymEngine::Basic> gg_num = G.get(cart_sym_2_num(dir[0]),cart_sym_2_num(dir[1]));
    std::vector<signed long int> det_terms;
    std::vector<std::complex<double>> det_coeffs;
    det_M22 = mul(sym_zero,z);
    for(auto const& imap: det_M22_dict){
        std::complex<double> coeffs_temp = SymEngine::rcp_static_cast<const SymEngine::ComplexDouble>(add(imap.second,sym_zero))->i ;
        if (std::abs(coeffs_temp)>std::pow(10,-8)){
            det_terms.push_back((SymEngine::rcp_static_cast<const SymEngine::Integer>(SymEngine::rcp_static_cast<const SymEngine::Pow>(pow(imap.first,SymEngine::integer(2)))->get_exp())->as_int())/2);
            det_coeffs.push_back(coeffs_temp);     
            det_M22 = add(det_M22,mul(imap.second,imap.first));
        }
    }
    auto minmax = std::minmax_element(det_terms.begin(), det_terms.end());
    std::vector<std::complex<double>> poly_coeff(*minmax.second-*minmax.first+1, c_zero);
    // OPENMP?
    for (int i=0; i<det_terms.size(); ++i){
        poly_coeff[det_terms[i]-*minmax.first] = det_coeffs[i];
    }
    Eigen::Map<Eigen::VectorXcd> coeff(poly_coeff.data(),poly_coeff.size()); 
    Eigen::PolynomialSolver<std::complex<double>, Eigen::Dynamic> solver;
    solver.compute(coeff);
    const Eigen::PolynomialSolver<std::complex<double>, Eigen::Dynamic>::RootsType &roots = solver.roots();
    auto integrand = mul(gg_num,pow(det_M22,sym_minus_one));
    SymEngine::RCP<const SymEngine::Basic> der_det_M22 = diff(det_M22,z);
    // do a dummy series expansion to find number of terms required
    SymEngine::RCP<const SymEngine::SeriesCoeffInterface> residue_series = SymEngine::series(mul(z,integrand), z,64);
    auto unordered_dict = residue_series->as_dict();
    std::map<int, SymEngine::RCP<const SymEngine::Basic>> ordered_dict(unordered_dict.begin(), unordered_dict.end());
    int series_order = ordered_dict.rbegin()->first;
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
            auto pole = SymEngine::complex_double(i+c_zero);
            std::complex<double> diff_at_pole = SymEngine::rcp_static_cast<const SymEngine::ComplexDouble>(expand(der_det_M22->subs({{z, pole}})))->i ;
            if (std::abs(diff_at_pole)<std::pow(10,-precision)){
               throw std::invalid_argument( "Higher order poles detected within the set precession. The program only works for single order poles" );
            }
            residue = residue + (SymEngine::rcp_static_cast<const SymEngine::ComplexDouble>(expand(gg_num->subs({{z, pole}})))->i)/diff_at_pole ;
         }
    // }
    }
    //auto integrand_p = integrand->subs({{z, add(z,p)}}); 
    return (residue*kp*std::exp(I*std::complex<double>(2.0)*kz*z0*std::complex<double>(2.0)*pi*k0)).real() ;
    ;
}


struct G_integrand_pars { double k0 ;
    std::complex<double> sxx_p ;
    std::complex<double> sxy_p ;
    std::complex<double> syx_p ;
    std::complex<double> syy_p ;
    double z0 ;
    std::string dir; };

double integrand_rad (double kp, void * G_params)
  {
    struct G_integrand_pars * params = (struct G_integrand_pars *)G_params;
    double k0 = params->k0;
    std::complex<double> sxx_p = params->sxx_p;
    std::complex<double> sxy_p = params->sxy_p;
    std::complex<double> syx_p = params->syx_p;
    std::complex<double> syy_p = params->syy_p;
    double z0 = params->z0;
    std::string dir= params->dir;

    return  green_integrand(kp,k0, sxx_p, sxy_p, syx_p, syy_p, z0,dir);
  }

double integrand_non_rad (double t, void * G_params)
  { double kp = 1.0f/t;
    struct G_integrand_pars * params = (struct G_integrand_pars *)G_params;
    double k0 = params->k0;
    std::complex<double> sxx_p = params->sxx_p;
    std::complex<double> sxy_p = params->sxy_p;
    std::complex<double> syx_p = params->syx_p;
    std::complex<double> syy_p = params->syy_p;
    double z0 = params->z0;
    std::string dir= params->dir;

    return  kp*kp*green_integrand(kp,k0, sxx_p, sxy_p, syx_p, syy_p, z0,dir);
  }

Decay_rates Purcell(double k0, \
                    std::complex<double> sxx_p,  \
                    std::complex<double> sxy_p, \
                    std::complex<double> syx_p, \
                    std::complex<double> syy_p, \
                    double z0 , \
                    std::string dir){
    
    size_t limit = 128;
	size_t m_deg = 10; // sets (2m+1)-point Gauss-Kronrod
    double err_r, err_nr, rad, non_rad;
    Workspace<double> Work(limit, m_deg);
    struct G_integrand_pars params = { k0, sxx_p, sxy_p, syx_p, syy_p, z0,dir };    
	Function<double, void> F_r(integrand_rad, &params);
    Function<double, void> F_nr(integrand_non_rad, &params);
	Work.qag(F_r, 0.01, 1.001, 1e-6, 1e-8, rad, err_r);
    Work.qag(F_nr, 1e-10, 1.0/1.001, 1e-6, 1e-8, non_rad, err_nr);
    std::cout<<rad<<"  "<<err_nr/non_rad<<std::endl;
    return Decay_rates { 1+ 1.5f*(rad + non_rad), rad, non_rad, err_r/rad, err_nr/non_rad};
	
}



int main(){

    //double ppp;
    //double kp_in;
    std::string dir_str;
    //std::cin>>kp_in;
    std::cin>>dir_str;
    auto t1 = std::chrono::high_resolution_clock::now();
    //for (int i=0; i<10; ++i){
    //auto ppp = green_integrand(kp_in,1145476.66561009,std::complex<double>(0.0001853183276793419,0.0135594673087213),c_zero,c_zero,std::complex<double>(3.778181026711009*std::pow(10,-5),0.003190917646215592),1.6*std::pow(10,-9),dir_str);
    auto ppp = Purcell(1145476.66561009,std::complex<double>(0.0001853183276793419,0.0135594673087213),c_zero,c_zero,std::complex<double>(3.778181026711009*1e-5,0.003190917646215592),1.6*1e-9,"zz");//}
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Time difference (s) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count())*1e-6 <<std::endl;
    std::cout<<ppp.total<<std::endl;
    //std::cout<<ppp<<std::endl;
    return 0;
}
