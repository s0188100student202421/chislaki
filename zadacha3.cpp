#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
/*
 * 4 -5 7
1 -4 9
-4 0 5 (13)

3 1 -1 2
-1 2 1 -1
5 -1 2 1
1 3 -1 -1 (-67)

1 -2 3 5
10 2 -1 3
8 6 -7 -7(неквадратная)
 */

const double INT_MIN = 1e-12;

double determinantGaussian(vector<vector<double>> A){
    int n = (int)A.size();
    if (n==0) return 1.0;

    for (auto &row: A)
        if ((int)row.size() != n)
            throw invalid_argument("Матрица не квадратная");

    int swaps  = 0;

    for(int col = 0; col < n;col++){
        int pivot_row = col;
        double max_abs = fabs(A[col][col]);
        for (int r = col + 1; r < n; r++){
            double current = fabs(A[r][col]);
            if (current > max_abs){
                max_abs = current;
                pivot_row = r;
            }
        }

        if (max_abs < INT_MIN) return 0.0;

        if (pivot_row != col){
            swap(A[pivot_row],A[col]);
            ++swaps;
        }

        double pivot = A[col][col];
        for (int r = col + 1; r < n; r++){
            double factor = A[r][col]/pivot;
            for (int c = col; c < n; c++)
                A[r][c] -= factor * A[col][c];
        }
    }

    double det = 1.0;
    for (int i = 0; i < n; i++) det *= A[i][i];
    if (swaps%2) det  = -det;
    return det;
}

int main(){
    int n;
    cout << "Введите размерность n: ";
    if (!(cin >> n) || n <= 0) {
        cerr << "Неверная размерность\n";
        return 1;
    }

    vector<vector<double>> A(n, vector<double>(n));
    cout << "Введите матрицу A:\n";
    for (int i = 0; i < n; ++i) {
        int a; cin >> a;
        A[i].push_back(a);
    }
    cout << endl;
    double det  = determinantGaussian(A);
    cout << "Посчитанный определитель с помощью метода Гаусса:" << det;
}

