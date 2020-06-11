#include <iostream>
#include <math.h>
#include "Solver.h"

#define MAX_ITER 100000
#define ACC_PREC 0.0000000001

bool Solver::GaussElimination(Matrix& A, Matrix& b, Matrix& x) {
    auto result = Matrix(A.getRows(), A.getCols());
    if(!A.inverse(result)){
		std::cout << "Gauss Elimination: Determination = 0. Unable to solve the linear system." << std::endl;
		return false;
	}
    result.matMul(b, x);
	return true;
}

bool Solver::LU(Matrix& A, Matrix& b, Matrix& x) {
	Matrix* L = new Matrix(A.getRows(), A.getCols());
	Matrix* y = new Matrix(A.getRows(), 1);
	if (!Solver::LU_decomposition(A, *L)) {
		std::cout << "LU: Unable to solve." << std::endl;
		return false;
	}
	if (!Solver::forward_substitution(*L, b, *y)) {
		std::cout << "LU: Unable to solve." << std::endl;
		return false;
	}
	if (!Solver::backward_substitution(A, *y, x)) {
		std::cout << "LU: Unable to solve." << std::endl;
		return false;
	}
	delete L;
	delete y;
	return true;
}

bool Solver::LU(CSRMatrix& A, Matrix& b, Matrix& x) {
	Matrix* L = new Matrix(A.getRows(), A.getCols());
	Matrix* y = new Matrix(A.getRows(), 1);
	if (!Solver::LU_decomposition(A, *L)) {
		std::cout << "LU: Unable to solve." << std::endl;
		return false;
	}
	CSRMatrix* L_Sparse = new CSRMatrix(*L);
	if (!Solver::forward_substitution(*L_Sparse, b, *y)) {
		std::cout << "LU: Unable to solve." << std::endl;
		return false;
	}
	if (!Solver::backward_substitution(A, *y, x)) {
		std::cout << "LU: Unable to solve." << std::endl;
		return false;
	}
	delete L;
	delete y;
	delete L_Sparse;
	return true;
}

bool Solver::LU_decomposition(Matrix& A, Matrix& L) {
	for (int i = 0; i < A.getRows(); i++) {
		for (int j = 0; j < A.getCols(); j++) {
			L.setValueAt(i, j, 0);
		}
	}
	for (int k = 0; k < A.getRows() - 1; k++) {
		for (int i = k + 1; i < A.getRows(); i++) {
			if (A.getValueAt(k, k) == 0) {
				return false;
			}
			double s = (A.getValueAt(i, k) / A.getValueAt(k, k));
			A.rocOp(k, -s, i, 1, false, k+1);
			L.setValueAt(i, k, s);
		}
	}
	for (int i = 0; i < A.getRows(); i++) {
		L.setValueAt(i, i, 1);
    }
	return true;
}

bool Solver::LU_decomposition(CSRMatrix& A, Matrix& L) {
	for (int k = 0; k < A.getRows() - 1; k++) {
		double diag = A.getValueAt(k, k);
		if (diag == 0) {
			return false;
		}
		for (int i = k + 1; i < A.getRows(); i++) {
			double s = (A.getValueAt(i, k) / diag);
			A.rocOp(k, -s, i, 1, false, k+1);
			ColNode* rowHead = A.getRowHead(i).head;
			for (; rowHead != nullptr && rowHead->col <= k; rowHead = rowHead->next);
			A.setRowHead(i, rowHead);
			L.setValueAt(i, k, s);
		}
	}
	for (int i = 0; i < A.getRows(); i++) {
		L.setValueAt(i, i, 1);
	}
	return true;
}
	
bool Solver::backward_substitution(Matrix& A, Matrix& b, Matrix& x) {
    for (int i = 0; i < A.getRows(); i++) {
		x.setValueAt(i, 0, 0);
	}
	for (int k = A.getRows() - 1; k > -1; k--) {
		double s = 0;
		for (int j = k + 1; j < A.getRows(); j++) {
			s += A.getValueAt(k, j) * x.getValueAt(j, 0);
        }
		if (A.getValueAt(k, k) == 0) {
			return false;
		}
        double tmp = (b.getValueAt(k, 0) - s) / A.getValueAt(k, k);
        x.setValueAt(k, 0, tmp);
	}
	return true;
}

bool Solver::backward_substitution(CSRMatrix& A, Matrix& b, Matrix& x) {
	for(int k = A.getRows() - 1; k >= 0; k--){
		double s = 0, diag = 0.0;
		ColNode* rowHead = A.getRowHead(k).head;
		for(; rowHead != nullptr && rowHead->col < k; rowHead = rowHead->next);
		if(rowHead != nullptr){
			diag = rowHead->value;
			rowHead = rowHead->next;
		}
		for(int j = k + 1; j < A.getRows(); j++){
			if(rowHead != nullptr && rowHead->col == j){
				s += rowHead->value * x.getValueAt(j, 0);
				rowHead = rowHead->next;
			}
		}
		if (diag == 0) {
			return false;
		}
		x.setValueAt(k, 0, (b.getValueAt(k, 0) - s) / diag);
	}
	return true;
}

bool Solver::forward_substitution(Matrix& A, Matrix& b, Matrix& x) {
    for (int i = 0; i < A.getRows(); i++) {
		x.setValueAt(i, 0, 0);
	}
	for (int k = 0; k < A.getRows(); k++) {
		double s = 0;
		for (int j = 0; j < k; j++) {
			s += A.getValueAt(k, j) * x.getValueAt(j, 0);
        }
		if (A.getValueAt(k, k) == 0) {
			return false;
		}
		double tmp = (b.getValueAt(k, 0) - s) / A.getValueAt(k, k);
		x.setValueAt(k, 0, tmp);
	}
	return true;
}

bool Solver::forward_substitution(CSRMatrix& A, Matrix& b, Matrix& x) {
	for(int k = 0; k < A.getRows(); k++){
		double s = 0, diag = 0.0;
		ColNode* rowHead = A.getRowHead(k).head;
		for(int j = 0; j <= k; j++){
			if(rowHead != nullptr && rowHead->col == j){
				if(j != k){
					s += rowHead->value * x.getValueAt(j, 0);
					rowHead = rowHead->next;
				} else {
					diag = rowHead->value;
					if (diag == 0) {
						return false;
					}
				}
			}
		}
		x.setValueAt(k, 0, (b.getValueAt(k, 0) - s) / diag);
	}
	return true;
}

bool Solver::LU_pp(Matrix& A, Matrix& b, Matrix& x) {
    Matrix* L = new Matrix(A.getRows(), A.getCols());
	Matrix* P_ = new Matrix(A.getRows(), A.getCols());
	Matrix* Pinvb = new Matrix(A.getRows(), 1);
	Matrix* y = new Matrix(A.getRows(), 1);
	if (!Solver::LU_decomposition_pp(A, *L, *P_)) {
		std::cout << "LU_pp: Unable to solve." << std::endl;
		return false;
	}
	P_->matMul(b, *Pinvb);
	if (!Solver::forward_substitution(*L, *Pinvb, *y)) {
		std::cout << "LU_pp: Unable to solve." << std::endl;
			return false;
	}
	if (!Solver::backward_substitution(A, *y, x)) {
		std::cout << "LU_pp: Unable to solve." << std::endl;
			return false;
	}
	delete L;
	delete P_;
	delete Pinvb;
	delete y;
	return true;
}

bool Solver::LU_pp(CSRMatrix& A, Matrix& b, Matrix& x) {
	Matrix* L = new Matrix(A.getRows(), A.getCols());
	Matrix* P_ = new Matrix(A.getRows(), A.getCols());
	CSRMatrix* Pinvb = new CSRMatrix(A.getRows(), 1);
	Matrix* y = new Matrix(A.getRows(), 1);
	if (!Solver::LU_decomposition_pp(A, *L, *P_)) {
		std::cout << "LU_pp: Unable to solve." << std::endl;
		return false;
	}
	CSRMatrix* L_sparse = new CSRMatrix(*L);
	CSRMatrix* P_sparse = new CSRMatrix(*P_);
	P_sparse->matMul(b, *Pinvb);
	if (!Solver::forward_substitution(*L_sparse, *Pinvb, *y)) {
		std::cout << "LU_pp: Unable to solve." << std::endl;
		return false;
	}
	if (!Solver::backward_substitution(A, *y, x)) {
		std::cout << "LU_pp: Unable to solve." << std::endl;
		return false;
	}
	delete L;
	delete P_;
	delete Pinvb;
	delete y;
	delete L_sparse;
	delete P_sparse;
	return true;
}

bool Solver::LU_decomposition_pp(Matrix& A, Matrix& L, Matrix& P_) {
	for (int i = 0; i < A.getRows(); i++) {
		for (int j = 0; j < A.getCols(); j++) {
			L.setValueAt(i, j, 0);
		}
	}
	for (int i = 0; i < A.getRows(); i++) {
		P_.setValueAt(i, i, 1);
    }
    for (int k = 0; k < A.getRows() - 1; k++) {
        double max = 0.0;
        int j = 0;
        for ( int l = k; l < A.getRows(); l++) {
            if (max < abs(A.getValueAt(l, k))) {
                max = abs(A.getValueAt(l, k));
                j = l;
            }
        }
        A.swap(j, k, false);
        P_.swap(j, k, false);
        L.swap(j, k, false);
        for (int i = k + 1; i < A.getRows(); i++) {
			if (A.getValueAt(k, k) == 0) {
				return false;
			}
            double s = A.getValueAt(i, k) / A.getValueAt(k, k);
			A.rocOp(k, -s, i, 1, false, k+1);
            L.setValueAt(i, k, s);
        }
    }
    for (int i = 0; i < A.getRows(); i++) {
        L.setValueAt(i, i, 1);
    }
	return true;
}

bool Solver::LU_decomposition_pp(CSRMatrix& A, Matrix& L, Matrix& P_) {
	for (int i = 0; i < A.getRows(); i++) {
		P_.setValueAt(i, i, 1);
	}
	for (int k = 0; k < A.getRows() - 1; k++) {
		double max = 0.0;
		int j = 0;
		for (int l = k; l < A.getRows(); l++) {
			if (max < abs(A.getValueAt(l, k))) {
				max = abs(A.getValueAt(l, k));
				j = l;
			}
		}
		A.swap(j, k, false);
		P_.swap(j, k, false);
		L.swap(j, k, false);
		double diag = A.getValueAt(k, k);
		if (diag == 0) {
			return false;
		}
		for (int i = k + 1; i < A.getRows(); i++) {
			double s = (A.getValueAt(i, k) / diag);
			A.rocOp(k, -s, i, 1, false, k+1);
			ColNode* rowHead = A.getRowHead(i).head;
			for (; rowHead != nullptr && rowHead->col <= k; rowHead = rowHead->next);
			A.setRowHead(i, rowHead);
			L.setValueAt(i, k, s);
		}
	}
	for (int i = 0; i < A.getRows(); i++) {
		L.setValueAt(i, i, 1);
	}
	return true;
}

bool Solver::Jacobi(Matrix& A, Matrix& b, Matrix& x) {
	double* next_x = new double[x.getRows()];
	for (int i = 0; i < A.getRows(); i++) {
		x.setValueAt(i, 0, 0);
	}
	//x.setValueAt(0, 0, 7);
	//x.setValueAt(1, 0, -3);
	bool converged = false;
	for (int k = 0; k < 10000 && !converged; k++) {
		converged = true;
		for (int i = 0; i < A.getRows(); i++) {
			double tmp = 0.0;
			for (int j = 0; j < A.getCols(); j++) {
				
				if (j != i)
					tmp += A.getValueAt(i, j) * x.getValueAt(j, 0);
			}
			if (A.getValueAt(i, i) == 0) {
				std::cout << "Jacobi: Unable to solve. Found zero diagonal value." << std::endl;
				return false;
			}
			next_x[i] = (b.getValueAt(i, 0) - tmp) / A.getValueAt(i, i);
			
		}
		for(int i = 0; i < x.getRows(); i++){
			converged &= abs(x.getValueAt(i, 0) - next_x[i]) <= ACC_PREC;
			x.setValueAt(i, 0, next_x[i]);
		}
		if (converged){
			//std::cout << " The number of iterations is:" << " " << k;
			delete[] next_x;
			return true;
		}
	}
	std::cout << "Jacobi: Unable to converge. Cannot solve this linear system." << std::endl;
	delete[] next_x;
	return false;
}

bool Solver::Jacobi(CSRMatrix& A, Matrix& b, Matrix& x) {
	int row_size = x.getRows();
	double* cur_x = new double[row_size];
	double* next_x = new double[row_size];
	for (int i = 0; i < row_size; i++) {
		cur_x[i] = 0;
	}
	int k = 0;
	for (bool converged = false; k < MAX_ITER && !converged; k++) {
		converged = true;
		for (int i = 0; i < row_size; i++) {
			double tmp = b.getValueAt(i, 0), diag = 0.0;
			for (ColNode* rowHead = A.getRowHead(i).head; rowHead != nullptr; rowHead = rowHead->next) {
				if (rowHead->col != i) {
					tmp -= rowHead->value * cur_x[rowHead->col];
				} else {
					diag = rowHead->value;
					if (diag == 0) {
						std::cout << "Jacobi: Unable to solve. Found zero diagonal value." << std::endl;
						return false;
					}
				}
			}
			tmp /= diag;
			converged &= (abs(cur_x[i] - tmp) <= ACC_PREC);
			next_x[i] = tmp;
		}
		for (int i = 0; i < row_size; i++) {
			cur_x[i] = next_x[i];
		}
	}
	for (int i = 0; i < row_size; i++) {
		x.setValueAt(i, 0, next_x[i]);
	}
	delete[] cur_x;
	delete[] next_x;
	if (k == MAX_ITER) {
		std::cout << "Jacobi: Unable to converge. Cannot solve this linear system." << std::endl;
		return false;
	}
	return true;
}

bool Solver::Gauss_Seidel(Matrix& A, Matrix& b, Matrix& x) {
	double* next_x = new double[x.getRows()];
	for (int i = 0; i < A.getRows(); i++) {
		x.setValueAt(i, 0, 0);
	}
	//x.setValueAt(0, 0, 7);
	//x.setValueAt(1, 0, -3);
	bool converged = false;
	for (int k = 0; k < MAX_ITER && !converged; k++) {
		converged = true;
		for (int i = 0; i < A.getRows(); i++) {
			double tmp1 = 0.0;
			double tmp2 = 0.0;
			
			for (int j = 0; j < i; j++) {
					
					tmp1 += A.getValueAt(i, j) * next_x[j];
			}

			for (int j = i + 1; j < A.getCols(); j++) {


					tmp2 += A.getValueAt(i, j) * x.getValueAt(j, 0);
			}
			if (A.getValueAt(i, i) == 0) {
				std::cout << "Jacobi: Unable to solve. Found zero diagonal value." << std::endl;
				return false;
			}
			next_x[i] = (b.getValueAt(i, 0) - tmp1 - tmp2) / A.getValueAt(i, i);
				
		}
		for (int i = 0; i < x.getRows(); i++) {
			converged &= abs(x.getValueAt(i, 0) - next_x[i]) <= ACC_PREC;
			x.setValueAt(i, 0, next_x[i]);
		}
		if (converged) {
			//std::cout << " The number of iterations is:" << " " << k;
			delete[] next_x;
			return true;
		}
	}
	std::cout << "Gauss Seidel: Unable to converge. Cannot solve this linear system." << std::endl;
	delete[] next_x;
	return false;
}

bool Solver::Gauss_Seidel(CSRMatrix& A, Matrix& b, Matrix& x) {
	int row_size = x.getRows();
	double* cur_x = new double[row_size];
	for (int i = 0; i < row_size; i++) {
		cur_x[i] = 0;
	}
	int k = 0;
	for (bool converged = false; k < MAX_ITER && !converged; k++) {
		converged = true;
		for (int i = 0; i < row_size; i++) {
			double tmp = b.getValueAt(i, 0), diag = 0.0;
			for (ColNode* rowHead = A.getRowHead(i).head; rowHead != nullptr; rowHead = rowHead->next) {
				int tmp_col = rowHead->col;
				if (tmp_col != i) {
					tmp -= rowHead->value * cur_x[tmp_col];
				} else {
					diag = rowHead->value;
					if (diag == 0) {
						std::cout << "Jacobi: Unable to solve. Found zero diagonal value." << std::endl;
						return false;
					}
				}
			}
			tmp /= diag;
			converged &= (abs(cur_x[i] - tmp) <= ACC_PREC);
			cur_x[i] = tmp;
		}
	}
	for (int i = 0; i < row_size; i++) {
		x.setValueAt(i, 0, cur_x[i]);
	}
	delete[] cur_x;
	if (k == MAX_ITER) {
		std::cout << "Gauss Seidel: Unable to converge. Cannot solve this linear system." << std::endl;
		return false;
	}
	return true;
}

bool Solver::Cholesky(Matrix& A, Matrix& b, Matrix& x) {
    Matrix* L = new Matrix(A.getRows(), A.getCols());
	Matrix* y = new Matrix(A.getRows(), 1);
	Matrix* L_T = new Matrix(L->getCols(), L->getRows());
	if (!Solver::CholeskyDecomposition(A, *L)) {
		std::cout << "Cholesky: Unable to solve." << std::endl;
		return false;
	}
	if (!Solver::forward_substitution(*L, b, *y)) {
		std::cout << "Cholesky: Unable to solve." << std::endl;
		return false;
	}
	L->Transpose(*L_T);
	if (!Solver::backward_substitution(*L_T, *y, x)) {
		std::cout << "Cholesky: Unable to solve." << std::endl;
		return false;
	}
	delete L;
	delete y;
	delete L_T;
	return true;
}

bool Solver::Cholesky(CSRMatrix& A, Matrix& b, Matrix& x) {
    Matrix* L = new Matrix(A.getRows(), A.getCols());
	Matrix* y = new Matrix(A.getRows(), 1);
	CSRMatrix* L_T = new CSRMatrix(L->getCols(), L->getRows());
	if (!Solver::CholeskyDecomposition(A, *L)) {
		std::cout << "Cholesky: Unable to solve." << std::endl;
		return false;
	}
	if (!Solver::forward_substitution(*L, b, *y)) {
		std::cout << "Cholesky: Unable to solve." << std::endl;
		return false;
	}
	CSRMatrix* L_sparse = new CSRMatrix(*L);
	L_sparse->Transpose(*L_T);
	if (!Solver::backward_substitution(*L_T, *y, x)) {
		std::cout << "Cholesky: Unable to solve." << std::endl;
		return false;
	}
	delete L;
	delete y;
	delete L_T;
	return true;
}


bool Solver::CholeskyDecomposition(Matrix& A, Matrix& L) {
    int n = A.getRows();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double tmp = A.getValueAt(i, j);
            for (int m = 0; m < j; m++) {
                tmp -= L.getValueAt(i, m) * L.getValueAt(j, m);
            }
            if (i != j) {
				if (L.getValueAt(j, j) == 0) {
					return false;
				}
                L.setValueAt(i, j, tmp / L.getValueAt(j, j));
                L.setValueAt(j, i, 0.0);
            } else {
                L.setValueAt(i, j, sqrt(tmp));
            }
        }
    }
	return true;
}

bool Solver::CholeskyDecomposition(CSRMatrix& A, Matrix& L) {
    int n = A.getRows();
    for (int i = 0; i < n; i++) {
		ColNode* rowHead = A.getRowHead(i).head;
        for (int j = 0; j <= i; j++) {
			double tmp = 0;
			if(rowHead != nullptr && rowHead->col == j){
				tmp = rowHead->value;
				rowHead = rowHead->next;
			}
            for (int m = 0; m < j; m++) {
                tmp -= L.getValueAt(i, m) * L.getValueAt(j, m);
            }
            if (i != j) {
				if (L.getValueAt(j, j) == 0) {
					return false;
				}
                L.setValueAt(i, j, tmp / L.getValueAt(j, j));
                L.setValueAt(j, i, 0.0);
            } else {
                L.setValueAt(i, j, sqrt(tmp));
            }
        }
    }
	return true;
}
