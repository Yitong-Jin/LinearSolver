#pragma once
#include <cstdio>
#include <iostream>
#include "Matrix.h"
#include "CSRMatrix.h"

using namespace std;

class Interface{

public:
    Interface();
    ~Interface();

    void start();

private:
    Matrix* A;
    CSRMatrix* A_sparse;

    Matrix* b;
    Matrix* x;

    bool csr;

    void invalidInput(char c);
    void printLinearSystem(Matrix& A, Matrix& x, Matrix& b, bool unx);

    void readData();
    void generateMatrix();
    void inputMatrix();
    void selectSolver();

    void runDenseSolver(bool (*denseSolver)(Matrix&, Matrix&, Matrix&));
    void runSparseSolver(bool (*sparseSolver)(CSRMatrix&, Matrix&, Matrix&));

    void runTest(string fileName);
    void chooseTest();

};