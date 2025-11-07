#include <iostream>
#include <cmath>
#include <functional>

using namespace std;

double simpson(function<double(double)> f, double a, double b){
    double c = (a + b) / 2.0;
    return (b-a)/6.0 * (f(a) + 4.0*f(c)+f(b));
}

double recurseSimp(function<double(double)> f, double a, double b, double eps, double s, double fa, double fb) {
    double c = (a+b)/2.0;
    double fc = f(c);
    double sleft = (b - a) / 12.0 * (fa + 4.0*f((a+c)/2.0) + fc);
    double sright = (b - a) / 12.0 * (fc + 4.0*f((c+b)/2.0) + fb);
    double s2 = sleft + sright;

    if (fabs(s2 - s) <= 15*eps)
        return s2+(s2-s)/15.0;
    return recurseSimp(f,a,c,eps/2.0,sleft,fa,fc) + recurseSimp(f, c, b, eps / 2.0, sright, fc, fb);
}
double adaptiveSimpson(function<double(double)>f, double a, double b, double eps) {
    double fa = f(a);
    double fb = f(b);
    double S = simpson(f,a,b);
    return recurseSimp(f,a,b,eps,S,fa,fb);
}

int main() {
    /*
     * == 2
    auto f = [](double x) {
        return std::sin(x);
    };
     == 107
     auto f = [](double x) {
        return std::sin(x) * exp(2*x);
    };
     == 1.1781
     auto f = [](double x) {
        return pow(std::sin(x), 4);
    };
     */

    double a = 0.0;
    double b = M_PI;
    double eps = 1e-6;
    double result = adaptiveSimpson(f,a,b,eps);
    cout.precision(10);
    cout << "Вычисленный интеграл:  " << result << endl;
    return 0;
}