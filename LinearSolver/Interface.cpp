#include <iostream>
#include <iomanip> 
#include "Interface.h"
#include "Matrix.h"
#include "Solver.h"
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <limits>

using namespace std;

Interface::Interface(){
    csr = false;
}

Interface::~Interface(){

}

void Interface::start(){
    system("CLS");
	cout << endl;
    cout << "-----------------------------------------------" << endl;
    cout << "|                                             |" << endl;
	cout << "|          Solvers for Linear System          |" << endl;
    cout << "|               Team Group Chat               |" << endl;
    cout << "|                                             |" << endl;
	cout << "-----------------------------------------------" << endl;
	cout << endl;


    cout << "-----------------------------------------------" << endl;
    cout << "Functions:" << endl << endl;
    cout << " 1: Browse the implemented solvers." << endl;
    cout << " 2: Try our solvers right now." << endl;
    cout << " q: Exit." << endl << endl;
    cout << ">> ";

    char input;
    cin >> input;

    switch (input) {
        case '1': chooseTest(); break;
        case '2': readData(); break;
        case 'q': exit(0);
        case 'Q': exit(0);
        default: invalidInput(input); start();
    }

}

void Interface::invalidInput(char c) {
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    system("CLS");
    cout << endl << " Invalid Selection : " << c << endl << endl;
    cout << "------------------------------" << endl;
    cout << endl;
    cout << " Press enter to return..." << endl;
    cin.get();
}

void Interface::readData() {
    system("CLS");
    cout << endl;
    cout << "Select the type of data:" << endl << endl;
    cout << " 1: Randomly generate a linear system." << endl;
    cout << " 2: Type your matrix through standard input." << endl;
    cout << " b: back" << endl;
    cout << " q: Exit." << endl << endl;
    cout << ">> ";

    char input;
    cin >> input;

    switch (input) {
    case '1': generateMatrix(); break;
    case '2': inputMatrix(); break;
    case 'b': start(); return;
    case 'B': start(); return;
    case 'q': return;
    case 'Q': return;
    default: invalidInput(input); readData(); return;
    }

    selectSolver();
}

void Interface::generateMatrix() {
    system("CLS");
    cout << endl;
    cout << "Randomly generate a linear system." << endl << endl;
    cout << "------------------------------" << endl;
    cout << "Use Dense Matrix or Sparse Matrix?" << endl << endl;
    cout << " 1: Dense Matrix" << endl;
    cout << " 2: Sparse Matrix" << endl;
    cout << " q: Exit." << endl << endl;
    cout << ">> ";

    char input;
    cin >> input;

    switch (input) {
    case '1': csr = false; break;
    case '2': csr = true; break;
    case 'q': exit(0);
    case 'Q': exit(0);
    default: invalidInput(input); generateMatrix(); return;
    }

    cout << "Please enter the size of A (n, integer): ";
    int rows = 0;
    while (rows <= 0 || rows > 1000) {
        string tmp;
        cin >> tmp;
        rows = atoi(tmp.c_str());
        if (rows <= 0 || rows > 1000) {
            cout << "The size MUST be greater than 0 and less than 1001" << endl;
            cout << "Please enter the size of A (n, integer): ";
        }
    }
    cout << "------------------------------" << endl;

    double* Avalues = new double[rows * rows];
    double* bvalues = new double[rows];
    double* xvalues = new double[rows];
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < i; j++) {
            if (csr) {
                if (rand() % 1000 <= 750) {
                    Avalues[i * rows + j] = 0;
                } else {
                    Avalues[i * rows + j] = rand() % 1000 + 0.2323 * (rand() % 5);
                }
            } else {
                Avalues[i * rows + j] = rand() % 1000 + 0.2323 * (rand() % 5);
            }
            Avalues[j * rows + i] = Avalues[i * rows + j];
        }
    }
    for (int i = 0; i < rows; i++) {
        bvalues[i] = rand() % 1000 + 0.234 * (rand() % 5);
        Avalues[i * rows + i] = (double) (rand() % 1000) + 9000;
    }
    this->A = new Matrix(rows, rows, Avalues);
    this->b = new Matrix(rows, 1, bvalues);
    printLinearSystem(*A, *x, *b, true);
}

void Interface::inputMatrix() {
    system("CLS");
    cout << endl;
    cout << "Type your matrix through standard input." << endl << endl;
    cout << "------------------------------" << endl;
    cout << "For convenience, we only accept dense matrix ipnut here." << endl;
    cout << "Please follow the given instructions." << endl;
    cout << "------------------------------" << endl;

    cout << "The linear system is defined as" << endl;
    cout << "Ax = b" << endl;
    cout << "where A is a nxn matrix and x, b are both nx1 vector." << endl;
    cout << "------------------------------" << endl;

    cout << "Please enter the size of A (n, integer, less than 30 is preferred): ";
    int rows = 0;
    while(rows <= 0 || rows > 1000) {
        string tmp;
        cin >> tmp;
        rows = atoi(tmp.c_str());
        if(rows <= 0 || rows > 1000){
            cout << "The size MUST be greater than 0 and less than 1001" << endl;
            cout << "Please enter the size of A (n, integer): ";
        }
    }
    cout << "------------------------------" << endl << endl;
    double* Avalues = new double[rows * rows];

    cout << "Please enter the " << rows * rows << " values of A (double)" << endl << endl;

    cout << "**********************************************" << endl;
    cout << "*** INVALID VALUES WILL BE TREATED AS ZERO ***" << endl;
    cout << "**********************************************" << endl << endl;

    for(int i = 0; i < rows; i++){
        cout << "Please enter the " << rows << " values on row " << i + 1  << " one by one seperated by empty space" << endl;
        for(int j = 0; j < rows; j++){
            string tmp;
            cin >> tmp;
            Avalues[i * rows + j] = atof(tmp.c_str());
        }
    }
    cout << "------------------------------" << endl;
    double* bvalues = new double[rows];
    cout << "Please enter the " << rows << " values of b (double)" << endl << endl;
    for(int i = 0; i < rows; i++){
        string tmp;
        cin >> tmp;
        bvalues[i] = atof(tmp.c_str());
    }

    this->A = new Matrix(rows, rows, Avalues);
    this->b = new Matrix(rows, 1, bvalues);
    this->x = new Matrix(rows, 1);

    cout << "Use Dense Matrix or Sparse Matrix to store your matrix?" << endl << endl;
    cout << " 1: Dense Matrix" << endl;
    cout << " 2: Sparse Matrix" << endl;
    cout << " q: Exit." << endl << endl;
    cout << ">> ";

    char input;
    cin >> input;

    switch (input) {
    case '1': csr = false; break;
    case '2': csr = true; break;
    case 'q': exit(0);
    case 'Q': exit(0);
    default: invalidInput(input); inputMatrix(); return;
    }

    cout << endl;

    printLinearSystem(*A, *x, *b, true);
}

void Interface::runDenseSolver(bool (*denseSolver)(Matrix&, Matrix&, Matrix&)) {
    system("CLS");
    cout << endl << "Running Solver on dense matrix" << endl;
    x = new Matrix(A->getRows(), 1);
    Matrix* A_tmp = new Matrix(A->getRows(), A->getCols());
    for (int i = 0; i < A->getRows(); i++) {
        for (int j = 0; j < A->getCols(); j++) {
            A_tmp->setValueAt(i, j, A->getValueAt(i, j));
        }
    }
    clock_t start = clock();
    bool re = denseSolver(*A_tmp, *b, *x);
    double t = (double)(clock() - start);
    printLinearSystem(*A, *x, *b, !re);
    cout << endl << "Time Elapsed: " << 1000 * t / CLOCKS_PER_SEC << " ms" << endl;
    delete A_tmp;
}

void Interface::runSparseSolver(bool (*sparseSolver)(CSRMatrix&, Matrix&, Matrix&)) {
    system("CLS");
    cout << endl << "Running Solver on sparse matrix" << endl;
    x = new Matrix(A->getRows(), 1);
    A_sparse = new CSRMatrix(*A);
    CSRMatrix* A_tmp = new CSRMatrix(A_sparse->getRows(), A_sparse->getCols());
    for (int i = 0; i < A_sparse->getRows(); i++) {
        for (int j = 0; j < A_sparse->getCols(); j++) {
            A_tmp->setValueAt(i, j, A_sparse->getValueAt(i, j));
        }
    }
    clock_t start = clock();
    bool re = sparseSolver(*A_tmp, *b, *x);
    double t = (double)(clock() - start);
    printLinearSystem(*A_sparse, *x, *b, !re);
    cout << endl << "Time Elapsed: " << 1000 * t / CLOCKS_PER_SEC << " ms" << endl;
    delete A_tmp;
}

void Interface::selectSolver() {
    cout << endl;
    cout << "Available Solvers:" << endl << endl;
    if(!csr) cout << " 1: Gaussian Elimination" << endl;
    cout << " 2: Jacobi" << endl;
    cout << " 3: Gauss-Seidel" << endl;
    cout << " 4: LU Fractorization" << endl;
    cout << " 5: LU with partial pivoting" << endl;
    cout << " 6: Cholesky Fractorization" << endl;
    cout << " b: Back" << endl;
    cout << " q: Exit" << endl;
    cout << ">> ";

    char input;
    cin >> input;
    
    if (!csr) {
        switch (input) {
        case '1': runDenseSolver(Solver::GaussElimination); break;
        case '2': runDenseSolver(Solver::Jacobi); break;
        case '3': runDenseSolver(Solver::Gauss_Seidel); break;
        case '4': runDenseSolver(Solver::LU); break;
        case '5': runDenseSolver(Solver::LU_pp); break;
        case '6': runDenseSolver(Solver::Cholesky); break;
        case 'b': readData(); return;
        case 'q': return;
        case 'Q': return;
        default: invalidInput(input); selectSolver(); return;
        }
    } else {
        switch (input) {
        //case '1': runDenseSolver(Solver::GaussElimination); break;
        case '2': runSparseSolver(Solver::Jacobi); break;
        case '3': runSparseSolver(Solver::Gauss_Seidel); break;
        case '4': runSparseSolver(Solver::LU); break;
        case '5': runSparseSolver(Solver::LU_pp); break;
        case '6': runSparseSolver(Solver::Cholesky); break;
        case 'b': readData(); return;
        case 'q': return;
        case 'Q': return;
        default: invalidInput(input); selectSolver(); return;
        }
    }

    cout << endl;
    cout << "------------------------------" << endl;
    cout << "Options:" << endl << endl;
    cout << " 1: Try another Solver." << endl;
    cout << " h: Home." << endl;
    cout << " q: Exit." << endl << endl;
    cout << ">> ";

    cin >> input;

    switch (input) {
    case '1': system("CLS"); selectSolver(); break;
    case 'h': start(); break;
    case 'H': start(); break;
    case 'q': return;
    case 'Q': return;
    default: invalidInput(input); cout << "return to Home" << endl; start(); return;
    }

}

void Interface::printLinearSystem(Matrix& A, Matrix& x, Matrix& b, bool unx){
    string emptySpace = "       ";
    int n_size = A.getRows();
    int cspace = n_size / 10 + 1;
    cout << " __ ";
    for(int i = 0; i < n_size; i++) cout << emptySpace;
    if(!unx) cout << "__   __ " << emptySpace << "__     __ " << emptySpace << "__ " << endl;
    else{
        cout << "__   __ " << " ";
        for(int i = 0; i < cspace; i++) cout << " ";
        cout << "__     __ " << emptySpace << "__ " << endl;
    }
    for(int i = 0; i < n_size; i++){
        if (i == n_size - 1) cout << "|__ ";
        else cout << "|   ";
        for(int j = 0; j < n_size; j++){
            double tmp = A.getValueAt(i, j);
            if(tmp >= 0) cout << " ";
            cout << setw(6) << setfill(' ') << setprecision(4) << round(tmp * 10000) / 10000; 
        }
        if (i == n_size - 1) cout << "__| ";
        else cout << "  | ";
        if (i == n_size - 1) cout << "|__ ";
        else cout << "|   ";
        if(!unx){
            double tmp = x.getValueAt(i, 0);
            if(tmp >= 0) cout << " ";
            else {
                cout << "-";
                tmp = -tmp;
            }
            cout << setw(6) << setfill('0') << setprecision(4) << round(tmp * 10000) / 10000;
        } else {
            cout << "x";
            cout << setw(cspace) << setfill('0') << setprecision(4) << i; 
        }
        if (i == n_size - 1) cout << "__| ";
        else cout << "  | ";
        if(i == n_size / 2) cout << "=";
        else cout << " ";
        if (i == n_size - 1) cout << " |__ ";
        else cout << " |   ";
        double tmp = b.getValueAt(i, 0);
        if(tmp >= 0) cout << " ";
        cout << setw(6) << setfill(' ') << setprecision(4) << round(tmp * 10000) / 10000;
        if (i == n_size - 1) cout << "__|" << endl;
        else cout << "  |" << endl;
    }
}

void Interface::chooseTest() {
    cout << endl;
    cout << "Available Test Files:" << endl << endl;
    cout << " 1: General" << endl;
    cout << " 2: Diagonal Dominant" << endl;
    cout << " 3: Symmetry" << endl;
    cout << " b: Back" << endl;
    cout << " q: Exit" << endl;
    cout << ">> ";

    char input;
    cin >> input;

    switch (input) {
    case '1': runTest("test_data/test_general.txt"); break;
    case '2': runTest("test_data/test_diagonal.txt"); break;
    case '3': runTest("test_data/test_symmetry.txt"); break;
    case 'b': start(); return;
    case 'B': start(); return;
    case 'q': return;
    case 'Q': return;
    default: invalidInput(input); chooseTest(); return;
    }

    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << endl << "------------------------------" << endl;
    cout << endl;
    cout << " Press enter to return" << endl;
    cin.get();

    start();
}

void Interface::runTest(string fileName) {
    system("CLS");
    cout << endl << endl;
    fstream fin, fout;
    fin.open(fileName, fstream::in);
    if (!fin) {
        cout << "failed to open test file" << endl;
        cout << "------------------------------" << endl;
        chooseTest();
        return;
    }
    int rows, cols;
    double* A_values, * b_values, * x_values;
    int count = 0;
    double prec = 1e-8;
    cout << "Testing: " << fileName << " with tolerance: " << prec << endl << endl;
    cout << "------------------------------" << endl << endl;
    cout << "GE: Gaussian Elimination" << endl;
    cout << "LU: LU Decomposition" << endl;
    cout << "LUP: LU Decomposition with Partial Pivoting" << endl;
    cout << "J: Jacobi - Iterative Method" << endl;
    cout << "GS: Gauss Seidel - Iterative Method" << endl;
    cout << "CH: Cholesky Decomposition" << endl << endl;
    cout << "------------------------------" << endl << endl;
    cout << "D: Dense Matirx S: Sparse Matirx" << endl << endl;
    cout << "------------------------------" << endl << endl;
    cout << "\t" << "GE-D\tGE-S\tLU-D\tLU-S\tLUP-D\tLUP-S\tJ-D\tJ-S\tGS-D\tGS-S\tCH-D\tCH-S" << endl;
    while (fin.good()) {
        count++;
        //cout << "Test" << count << ": ";
        fin >> rows >> cols;
        A_values = new double[rows * cols];
        b_values = new double[rows];
        x_values = new double[rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                fin >> A_values[i * cols + j];
            }
        }
        for (int i = 0; i < rows; i++) {
            fin >> b_values[i];
        }
        for (int i = 0; i < rows; i++) {
            fin >> x_values[i];
        }
        cout << rows << "x" << rows << "\t";
        Matrix* A1 = new Matrix(rows, cols, A_values);
        CSRMatrix* A1_s = new CSRMatrix(*A1);
        Matrix* b1 = new Matrix(rows, 1, b_values);
        Matrix* x_test = new Matrix(rows, 1);
        
        // Gaussian Elimination
        Solver::GaussElimination(*A1, *b1, *x_test);
        bool correct = true;
        for (int i = 0; i < rows; i++) {
            correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
        }
        cout << (correct ? "Pass" : "Fail") << "\t";

        b1 = new Matrix(rows, 1, b_values);
        x_test = new Matrix(rows, 1);
        Solver::GaussElimination(*A1_s, *b1, *x_test);
        correct = true;
        for (int i = 0; i < rows; i++) {
            correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
        }
        cout << (correct ? "Pass" : "Fail") << "\t";

        // LU
        A1 = new Matrix(rows, cols, A_values);
        A1_s = new CSRMatrix(*A1);
        b1 = new Matrix(rows, 1, b_values);
        x_test = new Matrix(rows, 1);

        Solver::LU(*A1, *b1, *x_test);
        correct = true;
        for (int i = 0; i < rows; i++) {
            correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
        }
        cout << (correct ? "Pass" : "Fail") << "\t";

        b1 = new Matrix(rows, 1, b_values);
        x_test = new Matrix(rows, 1);
        Solver::LU(*A1_s, *b1, *x_test);
        correct = true;
        for (int i = 0; i < rows; i++) {
            correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
        }
        cout << (correct ? "Pass" : "Fail") << "\t";

        // LUPP
        A1 = new Matrix(rows, cols, A_values);
        A1_s = new CSRMatrix(*A1);
        b1 = new Matrix(rows, 1, b_values);
        x_test = new Matrix(rows, 1);

        Solver::LU_pp(*A1, *b1, *x_test);
        correct = true;
        for (int i = 0; i < rows; i++) {
            correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
        }
        cout << (correct ? "Pass" : "Fail") << "\t";

        b1 = new Matrix(rows, 1, b_values);
        x_test = new Matrix(rows, 1);
        Solver::LU_pp(*A1_s, *b1, *x_test);
        correct = true;
        for (int i = 0; i < rows; i++) {
            correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
        }
        cout << (correct ? "Pass" : "Fail") << "\t";

        if (fileName.find("diagonal") != -1 || fileName.find("symmetry") != -1) {
            // Jacobi
            A1 = new Matrix(rows, cols, A_values);
            A1_s = new CSRMatrix(*A1);
            b1 = new Matrix(rows, 1, b_values);
            x_test = new Matrix(rows, 1);

            Solver::Jacobi(*A1, *b1, *x_test);
            correct = true;
            for (int i = 0; i < rows; i++) {
                correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
            }
            cout << (correct ? "Pass" : "Fail") << "\t";

            b1 = new Matrix(rows, 1, b_values);
            x_test = new Matrix(rows, 1);
            Solver::Jacobi(*A1_s, *b1, *x_test);
            correct = true;
            for (int i = 0; i < rows; i++) {
                correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
            }
            cout << (correct ? "Pass" : "Fail") << "\t";

            // Gauss_Seidel
            A1 = new Matrix(rows, cols, A_values);
            A1_s = new CSRMatrix(*A1);
            b1 = new Matrix(rows, 1, b_values);
            x_test = new Matrix(rows, 1);

            Solver::Gauss_Seidel(*A1, *b1, *x_test);
            correct = true;
            for (int i = 0; i < rows; i++) {
                correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
            }
            cout << (correct ? "Pass" : "Fail") << "\t";

            b1 = new Matrix(rows, 1, b_values);
            x_test = new Matrix(rows, 1);
            Solver::Gauss_Seidel(*A1_s, *b1, *x_test);
            correct = true;
            for (int i = 0; i < rows; i++) {
                correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
            }
            cout << (correct ? "Pass" : "Fail") << "\t";
        } else {
            cout << "/\t/\t/\t/\t";
        }

        if (fileName.find("symmetry") != -1) {
            // Cholesky
            A1 = new Matrix(rows, cols, A_values);
            A1_s = new CSRMatrix(*A1);
            b1 = new Matrix(rows, 1, b_values);
            x_test = new Matrix(rows, 1);

            Solver::Cholesky(*A1, *b1, *x_test);
            correct = true;
            for (int i = 0; i < rows; i++) {
                correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
            }
            cout << (correct ? "Pass" : "Fail") << "\t";

            b1 = new Matrix(rows, 1, b_values);
            x_test = new Matrix(rows, 1);
            Solver::Cholesky(*A1_s, *b1, *x_test);
            correct = true;
            for (int i = 0; i < rows; i++) {
                correct &= abs(x_test->getValueAt(i, 0) - x_values[i]) < prec;
            }
            cout << (correct ? "Pass" : "Fail") << "\t";
        } else {
            cout << "/\t/\t";
        }

        cout << endl;
        delete A1;
        delete A1_s;
        delete b1;
        delete x_test;
    }
}