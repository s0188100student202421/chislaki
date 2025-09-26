#include <bits/stdc++.h>

using namespace std;
const double EPS = 1e-15;

double x_to_t(double x, double x_min, double x_max){
    double length = x_max-x_min;
    if (fabs(length) < EPS) return 0.0;
    double x_mid = 0.5*(x_max+x_min);
    return 2.0*M_PI*(x-x_mid)/length;
}
vector<double> GaussRow(vector<vector<double>> A, vector<double> b){
    int n = A.size();
    vector<double> x(n);
    vector<double> colOrder(n); // поскольку для строк мы меняем места х-сов!!!

    for (int i = 0; i < n ; i++) colOrder[i] = i;

    for (int k = 0; k < n;k++){
        int max_col = k;
        for (int j = k+1; j < n; j++){
            if (abs(A[k][j]) > abs(A[k][max_col])) max_col = j;
        }
        if (A[k][max_col] == 0) throw runtime_error("Вырожденная система");

        for (int i = 0; i < n; i++){
            swap(A[i][k], A[i][max_col]);
        }
        swap(colOrder[k], colOrder[max_col]);

        for (int i = k+1; i<n;i++){
            double del = A[i][k] / A[k][k];
            for (int j = k; j< n;j++) A[i][j] -= del*A[k][j];
            b[i] -= del * b[k];
        }
    }

    for (int i = n - 1; i >= 0; i--){
        double sum = 0;
        for (int j = i+1; j < n; j++) sum += A[i][j] * x[j];
        x[i] = (b[i] - sum)/ A[i][i];
    }

    // надо восстановить порядок перменных
    vector<double> x_final(n);
    for(int i = 0; i < n; i++)
        x_final[colOrder[i]] = x[i];
    return x_final;
}
void printVector(const vector<double> &v, const string &label, int precision = 4) {
    cout << label << "\n";
    cout << fixed << setprecision(precision);
    for (double x : v) cout << x << "\t";
    cout << "\n";
    cout.unsetf(ios::fixed);
}
vector<vector<double>> buildMatrixPhi(const vector<double>& xs, int n, double x_min, double x_max){
    int m = (int)xs.size();
    int p = n+1;
    vector<vector<double>> Phi(m, vector<double>(p,0.0));
    for(int j = 0; j < m; j++){
        double t = x_to_t(xs[j], x_min, x_max);
        Phi[j][0] = 1.0;
        for (int k = 1; k <= n; k++) Phi[j][k] = cos(k*t);
    }
    return Phi;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout << "Введите количество гармоник" << endl;
    int n;
    if(!(cin >> n) || n < 0){cerr << "Неверное n" << endl;return 1;}
    int p = n+1;
    double x_min, x_max;
    cout << "Введите границы отрезка для нормирования в t" << endl;
    cin >> x_min >> x_max;
    if (x_min >= x_max) {cerr << "Неверные границы";}

    cout << "Введите узлы( пары х у)";
    vector<double> x(p), y(p);
    for (int i = 0; i < p; i++){
        cin >> x[i] >> y[i];
    }

    auto Phi = buildMatrixPhi(x,n,x_min, x_max);
    vector<vector<double>> A_system(p, vector<double>(p, 0.0));
    vector<double> b_system(p,0.0);
    for (int i = 0; i < p; i++){
        for(int k = 0; k <p;k++)
            A_system[i][k] = Phi[i][k];
        b_system[i] = y[i];
    }

    vector<double> result;
    try{
        result = GaussRow(A_system, b_system);
    }
    catch (const exception& ex){
        cerr << "Ошибка при решении системы методом Гаусса"<< ex.what() << endl;
        return 1;
    }

    cout << "\nКоэффициенты result_k:\n";
    for (int k = 0; k <= n; ++k)
        cout << "a[" << k << "] = " << fixed << setprecision(12) << result[k] << "\n";
    // Проверка в узлах
    for (int j = 0 ; j < p; j++){
        double t = x_to_t(x[j], x_min,x_max);
        double value = result[0];
        for (int k = 1;k <=n; k++)
            value += result[k]* cos(k*t);
        cout << "x[" << j << "]=" << x[j] << "  p(x)=" << value << "  y=" << y[j] << "  diff=" << value - y[j] << "\n";
    }
    return 0;
}



