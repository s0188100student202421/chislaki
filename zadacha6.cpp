#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
using namespace std;

class PicardSystem {
private:
    int n;
    double x0, X, h, eps;
    int max_iter;
    vector<double> y0;
    vector<function<double(double, const vector<double>&)>> f;

public:
    PicardSystem(int n_eq, double x0_, double X_, double h_, double eps_, int max_iter_)
            : n(n_eq), x0(x0_), X(X_), h(h_), eps(eps_), max_iter(max_iter_)
    {
        y0.resize(n, 0.0);
        f.resize(n);
    }

    void set_initial_conditions(const vector<double>& y0_) {
        y0 = y0_;
    }

    void set_function(int i, function<double(double, const vector<double>&)> func) {
        f[i] = func;
    }

    double integrate(function<double(double, const vector<double>&)> func,
                     const vector<double>& x,
                     const vector<vector<double>>& Y,
                     int i, double h, double y0)
    {
        double integral = 0.0;
        for (int j = 0; j < i; j++) {
            integral += func(x[j], Y[j]) * h;
        }
        return y0 + integral;
    }
    void solve() {
        int m = (int)((X - x0) / h);
        vector<double> x(m + 1);
        for (int i = 0; i <= m; i++)
            x[i] = x0 + i * h;

        vector<vector<double>> Y_prev(m + 1, vector<double>(n));
        vector<vector<double>> Y_curr(m + 1, vector<double>(n));

        // начальные условия
        for (int k = 0; k < n; k++)
            Y_prev[0][k] = y0[k];

        cout << "Метод Пикара для системы из " << n << " уравнений\n";

        for (int iter = 0; iter < max_iter; iter++) {
            for (int k = 0; k < n; k++)
                Y_curr[0][k] = y0[k];
            for (int i = 1; i <= m; i++) {
                for (int k = 0; k < n; k++) {
                    Y_curr[i][k] = integrate(f[k], x, Y_prev, i, h, y0[k]);
                }
            }

            double total_diff = 0.0;
            for (int i = 0; i <= m; i++)
                for (int k = 0; k < n; k++)
                    total_diff += fabs(Y_curr[i][k] - Y_prev[i][k]);

            cout << "Итерация " << iter + 1 << ": ошибка = " << total_diff << endl;

            if (total_diff < eps) {
                cout << "Метод сошелся после " << iter + 1 << " итераций.\n";
                break;
            }

            Y_prev = Y_curr;
        }
        cout << "\nРезультаты:\n";
        for (int i = 0; i <= m; i++) {
            cout << "x=" << x[i];
            for (int k = 0; k < n; k++)
                cout << "\ty" << (k + 1) << "=" << Y_curr[i][k];
            cout << endl;
        }
    }
};

int main() {
    int n;
    cout << "Введите количество уравнений в системе: ";
    cin >> n;

    // параметры задачи
    double x0 = 0.0, X = 1.0, h = 0.2, eps = 1e-4;
    int max_iter = 20;

    PicardSystem system(n, x0, X, h, eps, max_iter);

    // начальные условия
    vector<double> y0(n);
    cout << "Введите нач условия y0[i]:\n";
    for (int i = 0; i < n; i++) {
        cout << "y" << i + 1 << "(x0) = ";
        cin >> y0[i];
    }
    system.set_initial_conditions(y0);
    if (n == 2) {
        system.set_function(0, [](double x, const vector<double>& y) {
            return y[1] + x * x; // y1' = y2 + x^2
        });
        system.set_function(1, [](double x, const vector<double>& y) {
            return -y[0] * y[0] + 1; // y2' = -y1^2 + 1
        });
    }
    else if (n == 3) {
        system.set_function(0, [](double x, const vector<double>& y) {
            return y[1] + x; // y1' = y2 + x
        });
        system.set_function(1, [](double x, const vector<double>& y) {
            return -y[0] * y[2] + 1; // y2' = -y1*y3 + 1
        });
        system.set_function(2, [](double x, const vector<double>& y) {
            return sin(y[0]) - y[1]; // y3' = sin(y1) - y2
        });
    }
    else {
        cout << "Добавьте сюда свои функции для n=" << n << "!\n";
        return 0;
    }
    system.solve();

    return 0;
}
