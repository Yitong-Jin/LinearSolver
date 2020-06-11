#include <iostream>
#include "Matrix.h"

Matrix::Matrix(int rows, int cols) : rows(rows), cols(cols), preallocate(true){
	values_ptr = new double[rows * cols];
	for (int i = 0; i < rows * cols; i++) {
		values_ptr[i] = 0.0;
	}
}

Matrix::Matrix(int rows, int cols, double* values_ptr) : rows(rows), cols(cols), preallocate(true){
	this->values_ptr = new double[rows * cols];
	for (int i = 0; i < rows * cols; i++) {
		this->values_ptr[i] = values_ptr[i];
	}
}

Matrix::Matrix(int rows, int cols, int* values_ptr) : rows(rows), cols(cols), preallocate(true) {
	this->values_ptr = new double[rows * cols];
	for (int i = 0; i < rows * cols; i++) {
		this->values_ptr[i] = (double) values_ptr[i];
	}
}

Matrix::~Matrix() {
	if (preallocate) {
		delete[] values_ptr;
	}
}

int Matrix::getRows() const { return rows; }

int Matrix::getCols() const { return cols; }

void Matrix::setValueAt(int row, int col, double value) {
	values_ptr[row * cols + col] = value;
}

double Matrix::getValueAt(int row, int col) const {
	return values_ptr[row * cols + col];
}

void Matrix::matAdd(const Matrix& rightMatrix, Matrix& result) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			result.setValueAt(i, j, this->getValueAt(i, j) + rightMatrix.getValueAt(i, j));
		}
	}
}

void Matrix::matMinus(const Matrix& rightMatrix, Matrix& result) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			result.setValueAt(i, j, this->getValueAt(i, j) - rightMatrix.getValueAt(i, j));
		}
	}
}

void Matrix::matMul(const Matrix& rightMatrix, Matrix& result) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rightMatrix.cols; j++) {
			double tmp = 0.0;
			for (int k = 0; k < cols; k++) {
				tmp += this->getValueAt(i, k) * rightMatrix.getValueAt(k, j);
			}
			result.setValueAt(i, j, tmp);
		}
	}
}

bool Matrix::inverse(Matrix& result) {
	auto elimated = Matrix(rows, 2 * cols);
	auto eyes = Matrix(rows, rows);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++) {
			if (i == j) {
				eyes.values_ptr[i * rows + j] = 1;
			}
			else {
				eyes.values_ptr[i * rows + j] = 0;
			}
		}
	}
	if(!GauEli(eyes, elimated)){
		std::cout << "Cannot inverse the matirx" << std::endl;
		return false;
	}
	for (int i = rows - 1; i >= 0; i--) {
		for (int j = 0; j < cols; j++) {
			double sumleft = 0.0;
			for (int k = i + 1; k < rows; k++) {
				sumleft += elimated.getValueAt(i, k) * result.getValueAt(k, j);
			}
			if (elimated.getValueAt(i, i) == 0) {
				return false;
			}
			result.setValueAt(i, j, (elimated.getValueAt(i, j + cols) - sumleft) / elimated.getValueAt(i, i));
		}
	}
	return true;
}

bool Matrix::GauEli(Matrix& RHS, Matrix& result) {
	int new_rows = rows, new_cols = cols + RHS.cols;
	result.values_ptr = new double[new_rows * new_cols];
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			result.setValueAt(i, j, this->getValueAt(i, j));
		}
		for (int j = 0; j < RHS.cols; j++) {
			result.setValueAt(i, cols + j, RHS.getValueAt(i, j));
		}
	}
	for (int i = 0; i < rows; i++) {
		if (result.getValueAt(i, i) == 0) {
			int j = i;
			do {
				for (int k = 0; k < j; k++) {
					if (result.getValueAt(j, k) != 0) {
						result.rocOp(k, result.getValueAt(j, k), j, -1 * result.getValueAt(k, k), false, 0);
					}
				}
				j++;
			} while (j < rows && result.getValueAt(j, i) == 0);
			// if j == rows then cannot solve
			if(j == rows){
				return false;
			}
			result.swap(i, j, false);
		}
		for (int j = 0; j < i; j++) {
			if (result.getValueAt(i, j) != 0) {
				if (result.getValueAt(j, j) == 0) {
					return false;
				}
				result.rocOp(j, -1 * result.getValueAt(i, j) / result.getValueAt(j, j), i, 1, false, 0);
			}
		}
	}
	return true;
}

void Matrix::rocOp(int left, double leftCoef, int target, double rightCoef, bool col, int start) {
	int size = col ? rows : cols;
	for (int i = start; i < size; i++) {
		if (!col) {
			values_ptr[target * size + i] = values_ptr[target * size + i] * rightCoef + values_ptr[left * size + i] * leftCoef;
		}
		else {
			values_ptr[i * size + target] = values_ptr[i * size + target] * rightCoef + values_ptr[i * size + left] * leftCoef;
		}
	}
}

void Matrix::swap(int index1, int index2, bool col) {
	int size = col ? rows : cols;
	if (!col) {
		double tmp = 0;
		for (int i = 0; i < cols; i++) {
			tmp = values_ptr[index1 * cols + i];
			values_ptr[index1 * cols + i] = values_ptr[index2 * cols + i];
			values_ptr[index2 * cols + i] = tmp;
		}
	}
	else {
		double tmp = 0;
		for (int i = 0; i < rows; i++) {
			tmp = values_ptr[i * cols + index1];
			values_ptr[i * cols + index1] = values_ptr[i * cols + index2];
			values_ptr[i * cols + index2] = tmp;
		}
	}
}

void Matrix::printMatrix() const {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			std::cout << values_ptr[i * cols + j] << " ";
		}
		std::cout << std::endl;
	}
}

void Matrix::Transpose(Matrix& result_T) {
	result_T.rows = cols;
	result_T.cols = rows;
	for (int j = 0; j < result_T.rows; j++) {
		for (int i = 0; i < result_T.cols; i++) {
			result_T.values_ptr[j * result_T.cols + i] = values_ptr[i * cols + j];
		}
	}
}
