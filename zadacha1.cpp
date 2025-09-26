#include "zadacha1.h"
using namespace std;

vector<double> GaussColumn(vector<vector<double>> A,  vector<double> b){
    int n = A.size();
    vector<double> x(n);

    for (int k = 0; k < n; k++){
        int max_row = k;
        for (int i = k+1; i < n; i++){
            if (abs(A[i][k]) > abs(A[max_row][k]))
                max_row = i;
        }

        if (A[max_row][k] == 0) throw runtime_error("Вырожденная система");

        swap(A[k], A[max_row]);
        swap(b[k], b[max_row]);

        for (int i = k+1; i < n; i++){
            double del = A[i][k] / A[k][k];
            for (int j = k; j < n; j++){
                A[i][j] -= del * A[k][j];
            }
            b[i] -= del * b[k];
        }
    }

    for (int i = n -1 ; i >= 0;i--){
        double sum = 0;
        for (int j = i + 1; j < n; ++j)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
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

vector<double> GaussMatrix(vector<vector<double>> A, vector<double> b){
    int n = A.size();
    vector<double> x(n);
    vector<int> col_order(n);
    for (int i = 0; i < n; i++) col_order[i] = i;

    for (int k = 0;k < n; k++){
        double max_val = 0;
        int iMax = k, jMax = k;
        for (int i = k; i < n; i++){
            for (int j = k;j < n;j++){
                if (abs(A[i][j]) > max_val){
                    max_val = abs(A[i][j]);
                    iMax = i;jMax = j;
                }
            }
        }

        if (max_val == 0) throw runtime_error("Вырожденная система");

        swap(A[k], A[iMax]);
        swap(b[k], b[iMax]);

        for (int i = 0; i < n; i++) swap(A[i][k], A[i][jMax]);
        swap(col_order[k], col_order[jMax]);

        for (int i = k+1; i< n; i++){
            double del = A[i][k]/ A[k][k];
            for (int j = k; j < n;j++)
                A[i][j] -= del*A[k][j];
            b[i] -= del*b[k];
        }
    }

    for (int i = n-1; i >= 0; i--){
        double sum = 0;
        for (int j = i+1; j<n;j++)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum)/A[i][i];
    }

    vector<double> x_final(n);
    for (int i = 0; i < n; i++) x_final[col_order[i]] = x[i];
    return x_final;
}

vector<vector<double>> readMatrix(int n) {
    vector<vector<double>> A(n, vector<double>(n));
    cout << "Введите матрицу A (по строкам):\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cin >> A[i][j];
    return A;
}

vector<double> readVector(int n) {
    vector<double> b(n);
    cout << "Введите вектор b:\n";
    for (int i = 0; i < n; i++)
        cin >> b[i];
    return b;
}

void printVector(const vector<double> &v, const string &label, int precision = 4) {
    cout << label << "\n";
    cout << fixed << setprecision(precision);
    for (double x : v) cout << x << "\t";
    cout << "\n";
    cout.unsetf(ios::fixed);
}


int main() {
    int precision = 5;
    int n;
    cout << "Введите размерность системы n:\n";
    cin >> n;

    auto A = readMatrix(n);
    auto b = readVector(n);

    cout << "\nМатрица A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << A[i][j] << " ";
        cout << "\n";
    }
    printVector(b, "Вектор b:",precision);

    auto x_col = GaussColumn(A, b);
    auto x_row = GaussRow(A, b);
    auto x_mat = GaussMatrix(A, b);

    printVector(x_col, "Решение (выбор по столбцу):",precision);
    printVector(x_row, "Решение (выбор по строке):",precision);
    printVector(x_mat, "Решение (полный выбор):",precision);

    return 0;
}

