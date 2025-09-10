#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
/*
 * 1 2 3
0 1 4
5 6 0

3 4
1 2

1 2 3
4 5 6
7 8 9(не сущ)
 */
const double INT_MIN = 1e-12;

vector<vector<double>> multiply(const vector<vector<double>>& A, const vector<vector<double>>& B){
    int n = A.size();
    vector<vector<double>> C(n, vector<double>(n,0.0));
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

void print_matrix(const vector<vector<double>>& M){
    int n = M.size();
    for (int i = 0; i < n; i++){
        for (int j = 0; j < (int)M[i].size();j++)
            cout << M[i][j] << "\t";
        cout << endl;
    }
}
bool inverMatrix(const vector<vector<double>>& A, vector<vector<double>>& A_inverted){
    int n = A.size();
    if (n==0) return false;
    for (auto &row:A)if ((int)row.size() != n) return false;

    vector<vector<double>> extended_matrix(n, vector<double>(2*n, 0.0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            extended_matrix[i][j] =  A[i][j];
        extended_matrix[i][n+i] = 1.0;
    }

    //приводим к единичной матрице
    for (int col = 0; col < n; col++){
        int pivot_row = col;
        double max_val = fabs(extended_matrix[col][col]);
        for (int r = col+1; r<n; r++){
            double val = fabs(extended_matrix[r][col]);
            if (val > max_val){
                max_val = val;
                pivot_row = r;
            }
        }

        if (max_val < INT_MIN) return false;

        if (pivot_row != col) swap(extended_matrix[pivot_row], extended_matrix[col]);

        //НОРМАЛИЗАЦИЯ
        double diag = extended_matrix[col][col];
        for (int j = 0; j < 2*n; j++)
            extended_matrix[col][j] /= diag;

        for (int r = 0; r<n; r++){
            if(r==col) continue;
            double del = extended_matrix[r][col];
            if (fabs(del) > 0.0){
                for (int j = 0; j<2*n;j++)
                    extended_matrix[r][j] -= del* extended_matrix[col][j];
            }
        }
    }
    A_inverted.assign(n, vector<double>(n,0.0));
    for (int i = 0; i < n;i++)
        for (int j = 0; j < n; j++)
            A_inverted[i][j] = extended_matrix[i][n+j];
    return true;
}
vector<vector<double>> roundToIdentity(const vector<vector<double>> &M) {
    int n = M.size();
    vector<vector<double>> R = M;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // если близко к 1 (в диагонали), ставим 1
            if (fabs(R[i][j] - 1.0) < 0.1)
                R[i][j] = 1.0;
                // если близко к 0, ставим 0
            else if (fabs(R[i][j]) < 0.1)
                R[i][j] = 0.0;
                // иначе можно оставить как есть или округлить до 1 знака
            else
                R[i][j] = round(R[i][j] * 10.0) / 10.0;
        }
    }
    return R;
}
int main() {
    int n;
    cout << "Введите размерность:\n";
    cin >> n;
    vector<vector<double>> A(n,vector<double>(n));
    cout << "Введите матрицу А:";
    for (int i = 0; i < n; i++)
        for(int j = 0;j < n;j++)
            cin >> A[i][j];

    vector<vector<double>> A_inverted;
    bool flag = inverMatrix(A, A_inverted);
    if (!flag){
        cout << "Обратной матрицы для введенной не существует";
        return 0;
    }
    cout << "Исходная матрица:\n";
    print_matrix(A);
    cout << "Инвертированная матрица:\n";
    print_matrix(A_inverted);
    cout << "Проверка умножением матриц A и A_inverted:\n";
    auto Iapprox = multiply(A, A_inverted);
    print_matrix(roundToIdentity(Iapprox));
}