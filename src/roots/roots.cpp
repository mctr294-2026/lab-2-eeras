#include "roots.hpp"
#include <cmath>
#include <functional>

bool bisection(std::function<double(double)> f, double a, double b, double *root) {
    
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0)
        return false;

    const double tolerance = 1e-6;
    const int max_iterations = 1e6;

    for (int i=0 ; i < max_iterations ; ++i) {
         double c = (a+b)/2.0;
        double fc = f(c);

        if (std::abs(fc) < tolerance){
        *root = c;
        return true;
        }

        if ((fa > 0 && fc >0) || (fa < 0 && fc < 0)){
            a = c;
            fa = fc;
        }
        else {
            b = c;
            fb = fc;
        }
    }
return false;

}


bool regula_falsi(std::function<double(double)> f, double a, double b, double *root) {
    
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0)
        return false;

    const double tolerance = 1e-6;
    const double max_iterations = 1e6;

    for (int i = 0 ; i < max_iterations ; ++i){
        double c = a - (f(a)*(b-a))/(f(b)-f(a));
        double fc = f(c);
    
    if (std::abs(fc) < tolerance){
        *root = c;
        return true;
    }
    if ((fa > 0 && fc > 0) || (fa < 0 && fc < 0)) {
            a = c;
            fa = fc;
        } else {
            b = c;
            fb = fc;
        }
    }
return false;
}

bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c_n, double *root) {
    
    const double tolerance = 1e-6;
    const double max_iterations = 1e6;

    for(int i = 0 ; i < max_iterations ; ++i){

        if (std::abs(f(c_n)) < tolerance){
            *root = c_n;
            return true;
        }
        
        double c_nplus1 = c_n - (f(c_n) / g(c_n));
        c_n = c_nplus1;

    }
return false;
}

bool secant(std::function<double(double)> f, double a, double b, double c, double *root) {

    const double tolerance = 1e-6;
    const int max_iterations = 1e6;

    double x_n = a;
    double x_nminus1 = b;
    double x_nplus1;

    for (int i =0 ; i < max_iterations ; ++i){
        x_nplus1 = x_n - f(x_n)*((x_n - x_nminus1)/(f(x_n)-f(x_nminus1)));
        
        if (std::abs(x_nplus1 - x_n) < tolerance){
            *root = x_nplus1;
            return true;
        }

        x_nminus1 = x_n;
        x_n = x_nplus1;
    }
return false;
}

