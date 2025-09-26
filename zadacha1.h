#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>
using namespace std;

vector<double> GaussColumn(vector<vector<double>> A,  vector<double> b);
vector<double> GaussRow(vector<vector<double>> A, vector<double> b);
vector<double> GaussMatrix(vector<vector<double>> A, vector<double> b);
vector<vector<double>> readMatrix(int n);
vector<double> readVector(int n);
void printVector(const vector<double> &v, const string &label, int precision);

#ifndef CHISLAKI_ZADACHA1_H
#define CHISLAKI_ZADACHA1_H

#endif //CHISLAKI_ZADACHA1_H
