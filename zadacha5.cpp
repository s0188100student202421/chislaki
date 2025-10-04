#include <bits/stdc++.h>
using namespace std;

bool newton_simple(double (*f)(double), double (*df)(double),
                   double x0, double a, double b,
                   double &out,
                   int maxIter = 50,
                   double fTol = 1e-12,
                   double xTol = 1e-12){
    double x = x0;
    for (int it = 0; it < maxIter;++it){
        double fx = f(x);
        if (fabs(fx) <= fTol) {out = x; return true;}

        double dfx = fx/df(x);
        if(fabs(dfx) <= fTol){return false;}

        double step = fx/dfx;
        double alpha = 1.0;
        bool moved = false;

        for (int ls = 0; ls < 30; ++ls){
            double xnew = x - alpha*step;
            if (xnew < a || xnew > b){alpha *= 0.5; continue;}
            double fxnew = f(xnew);
            if (fabs(fxnew) < fabs(fx) || fabs(fxnew) <= fTol){
                x = xnew;
                moved = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!moved) return false;

        if (fabs(alpha*step) <= fTol) {out = x; return true;}
    }
    return false;
}

double bisect_simple(double (*f)(double), double lo, double hi, double tol = 1e-8, int maxIter = 60){
    double flo = f(lo), fhi = f(hi);
    if (fabs(flo) <= tol) return lo;
    if (fabs(fhi) <= f(hi)) return hi;
    for (int it = 0; it < maxIter && (hi-lo) > tol; ++it){
        double m = 0.5*(lo+hi);
        double fm = f(m);

        if (fabs(fm) <= tol) return m;
        if(flo* fm <= 0) {hi = m; fhi = fm;}
        else{lo = m; flo =fm;}
    }
    return 0.5*(lo+hi);
}

vector<double> find_all_roots(double (*f)(double), double (*df)(double), double a, double b, int samples = 200){
    vector<double> roots;
    if (a >= b) return roots;
    int N = max(4, samples);
    vector<double> xs(N+1), fs(N+1);
    for (int i = 0; i <= N; ++i){
        xs[i] = a + (b-a)*double(i)/N;
        fs[i] = f(xs[i]);
        if(fabs(fs[i]) < 1e-12) roots.push_back(xs[i]);
    }

    for (int i = 0; i < N; ++i) {
        if (fs[i] * fs[i+1] < 0) {
            double bracket_mid = bisect_simple(f, xs[i], xs[i+1]);
            double root;
            bool ok = newton_simple(f, df, bracket_mid, a, b, root);
            if (ok) roots.push_back(root);
            else roots.push_back(bracket_mid);
        }
    }

    sort(roots.begin(), roots.end());

    vector<double> uniq;
    double merge_tol = 1e-8;
    for (double r : roots) {
        if (!uniq.empty() && fabs(r - uniq.back()) <= merge_tol) continue;
        uniq.push_back(r);
    }
    return uniq;
}

double f_example(double x) { return x*x*x - x - 1.0; }
double df_example(double x) { return 3.0*x*x - 1.0; }

int main() {
    double a = 0.0, b = 2.0;
    auto roots = find_all_roots(f_example, df_example, a, b, 400);
    cout << "Найденные корни на [" << a << "," << b << "]:\n";
    cout << setprecision(12);
    for (double r : roots) cout << r << "  f(r)=" << f_example(r) << "\n";
    return 0;
}
