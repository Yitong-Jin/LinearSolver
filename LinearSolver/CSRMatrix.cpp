#include <iostream>
#include <assert.h>
#include "CSRMatrix.h"

CSRMatrix::CSRMatrix(int rows, int cols) : Matrix(rows, cols) {
	rowNodes = new RowNode[rows];
	for (int i = 0; i < rows; i++) {
		rowNodes[i].head = nullptr;
	}
}

CSRMatrix::CSRMatrix(int rows, int cols, double* values, int* row_position, int* col_Index) : Matrix(rows, cols) {
	rowNodes = new RowNode[rows];
	for (int i = 0; i < rows; i++) {
		rowNodes[i].head = nullptr;
		ColNode* tmp = rowNodes[i].head;
		for (int j = row_position[i]; j < row_position[i + 1]; j++) {
			ColNode* new_node = new ColNode(col_Index[j], values[j]);
			if (tmp == nullptr) {
				rowNodes[i].head = new_node;
				tmp = rowNodes[i].head;
			} else {
				tmp->next = new_node;
				tmp = tmp->next;
			}
		}
	}
}

CSRMatrix::CSRMatrix(int rows, int cols, int* values, int* row_position, int* col_Index) : Matrix(rows, cols) {
	rowNodes = new RowNode[rows];
	for (int i = 0; i < rows; i++) {
		rowNodes[i].head = nullptr;
		ColNode* tmp = rowNodes[i].head;
		for (int j = row_position[i]; j < row_position[i + 1]; j++) {
			ColNode* new_node = new ColNode(col_Index[j], (double) values[j]);
			if (tmp == nullptr) {
				rowNodes[i].head = new_node;
				tmp = rowNodes[i].head;
			}
			else {
				tmp->next = new_node;
				tmp = tmp->next;
			}
		}
	}
}

CSRMatrix::CSRMatrix(Matrix& matrix) : Matrix(matrix.getRows(), matrix.getCols()) {
	this->rowNodes = new RowNode[matrix.getRows()];
	for (int i = 0; i < matrix.getRows(); i++) {
		rowNodes[i].head = nullptr;
		ColNode* tmp = rowNodes[i].head;
		for (int j = 0; j < matrix.getCols(); j++) {
			double tmpv = matrix.getValueAt(i, j);
			if (tmpv != 0) {
				ColNode* new_node = new ColNode(j, tmpv);
				if (tmp == nullptr) {
					rowNodes[i].head = new_node;
					tmp = rowNodes[i].head;
				}
				else {
					tmp->next = new_node;
					tmp = tmp->next;
				}
			}
		}
	}
}

CSRMatrix::~CSRMatrix() {
	for (int i = 0; i < rows; i++) {
		ColNode* tmp = rowNodes[i].head;
		ColNode* next = nullptr;
		while (tmp != nullptr) {
			next = tmp->next;
			delete tmp;
			tmp = next;
		}
	}
	delete[] rowNodes;
}

void CSRMatrix::setValueAt(int row, int col, double value) {
	ColNode* tmp = rowNodes[row].head;
	ColNode* prev = nullptr;
	while (tmp != nullptr && col > tmp->col) {
		prev = tmp;
		tmp = tmp->next;
	}
	if (tmp != nullptr && col == tmp->col) {
		tmp->value = value;
	} else {
		ColNode* newNode = new ColNode(col, value);
		if (prev == nullptr) {
			newNode->next = nullptr;
			rowNodes[row].head = newNode;
		}
		else {
			newNode->next = prev->next;
			prev->next = newNode;
		}
	}
}

double CSRMatrix::getValueAt(int row, int col) const {
	ColNode* tmp = rowNodes[row].head;
	for (; tmp != nullptr && col > tmp->col; tmp = tmp->next) {}
	if (tmp != nullptr && col == tmp->col) {
		return tmp->value;
	} else {
		return 0.0;
	}
}

void CSRMatrix::printMatrix() const {
	for (int i = 0; i < rows; i++) {
		ColNode* cur = rowNodes[i].head;
		for (int j = 0; j < cols; j++) {
			if (cur == nullptr || cur->col > j) {
				std::cout << 0 << " ";
			} else {
				std::cout << cur->value << " ";
				cur = cur->next;
			}
		}
		std::cout << std::endl;
	}
}

void CSRMatrix::rocOp(int left, double leftCoef, int target, double targetCoef, bool col, int start) {
	assert(!col); // Only support row operations
	ColNode* leftCol = rowNodes[left].head;
	ColNode* targetCol = rowNodes[target].head;
	ColNode* prevTarget = nullptr;
	for (; leftCol != nullptr && leftCol->col < start; leftCol = leftCol->next);
	for (; targetCol != nullptr && targetCol->col < start; targetCol = targetCol->next);
	while(!(leftCol == nullptr && targetCol == nullptr)){
		if(leftCol != nullptr && targetCol != nullptr && leftCol->col == targetCol->col){
			targetCol->value = leftCoef * leftCol->value + targetCoef * targetCol->value;
			prevTarget = targetCol;
			targetCol = targetCol->next;
			leftCol = leftCol->next;
		} else if(leftCol != nullptr && targetCol != nullptr){
			if(leftCol->col > targetCol->col){
				targetCol->value = targetCoef * targetCol->value;
				prevTarget = targetCol;
				targetCol = targetCol->next;
			} else {
				ColNode* tmp = new ColNode(leftCol->col, leftCoef * leftCol->value);
				tmp->next = targetCol;
				if(prevTarget == nullptr){
					rowNodes[target].head = tmp;
					prevTarget = rowNodes[target].head;
				} else {
					prevTarget->next = tmp;
					prevTarget = prevTarget->next;
				}
				leftCol = leftCol->next;
			}
		} else if(targetCol != nullptr){
			targetCol->value = targetCoef * targetCol->value;
			prevTarget = targetCol;
			targetCol = targetCol->next;
		} else {
			while(leftCol != nullptr){
				ColNode* tmp = new ColNode(leftCol->col, leftCoef * leftCol->value);
				if(prevTarget == nullptr){
					rowNodes[target].head = tmp;
					prevTarget = rowNodes[target].head;
				} else {
					prevTarget->next = tmp;
				}
				tmp->next = nullptr;
				prevTarget = prevTarget->next;
				leftCol = leftCol->next;
			}
		}
	}
}

void CSRMatrix::swap(int index1, int index2, bool col) {
	assert(!col); // Only support row swap
	ColNode* tmp = rowNodes[index1].head;
	rowNodes[index1].head = rowNodes[index2].head;
	rowNodes[index2].head = tmp;
}

void CSRMatrix::matAdd(const CSRMatrix& rightMatrix, CSRMatrix& result) {
	for (int i = 0; i < this->rows; i++) {
		ColNode* leftCol = rowNodes[i].head;
		ColNode* rightCol = rightMatrix.rowNodes[i].head;
		ColNode* resultCol = result.rowNodes[i].head;
		while (!(leftCol == nullptr && rightCol == nullptr)) {
			int c_col = 0;
			double c_v = 0;
			if (leftCol != nullptr && rightCol != nullptr && leftCol->col == rightCol->col) {
				c_v = leftCol->value + rightCol->value;
				c_col = leftCol->col;
				leftCol = leftCol->next;
				rightCol = rightCol->next;
			}
			else if ((leftCol != nullptr && rightCol != nullptr && leftCol->col > rightCol->col) || (leftCol == nullptr)) {
				c_v = rightCol->value;
				c_col = rightCol->col;
				rightCol = rightCol->next;
			}else {
				c_v = leftCol->value;
				c_col = leftCol->col;
				leftCol = leftCol->next;
			}
			ColNode* tmp = new ColNode(c_col, c_v);
			if (resultCol == nullptr) {
				result.rowNodes[i].head = tmp;
				resultCol = result.rowNodes[i].head;
			} else {
				resultCol->next = tmp;
				resultCol = resultCol->next;
			}
		}
	}
}

void CSRMatrix::matMinus(const CSRMatrix& rightMatrix, CSRMatrix& result) {
	for (int i = 0; i < this->rows; i++) {
		ColNode* leftCol = rowNodes[i].head;
		ColNode* rightCol = rightMatrix.rowNodes[i].head;
		ColNode* resultCol = result.rowNodes[i].head;
		while (!(leftCol == nullptr && rightCol == nullptr)) {
			int c_col = 0;
			double c_v = 0;
			if (leftCol != nullptr && rightCol != nullptr && leftCol->col == rightCol->col) {
				c_v = leftCol->value - rightCol->value;
				c_col = leftCol->col;
				leftCol = leftCol->next;
				rightCol = rightCol->next;
			}
			else if ((leftCol != nullptr && rightCol != nullptr && leftCol->col > rightCol->col) || (leftCol == nullptr)) {
				c_v = -1 * rightCol->value;
				c_col = rightCol->col;
				rightCol = rightCol->next;
			}
			else {
				c_v = leftCol->value;
				c_col = leftCol->col;
				leftCol = leftCol->next;
			}
			ColNode* tmp = new ColNode(c_col, c_v);
			if (resultCol == nullptr) {
				result.rowNodes[i].head = tmp;
				resultCol = result.rowNodes[i].head;
			}
			else {
				resultCol->next = tmp;
				resultCol = resultCol->next;
			}
		}
	}
}

void CSRMatrix::matMul(const CSRMatrix& rightMatrix, CSRMatrix& result) {
	ColNode** LNodes = new ColNode * [this->rows];
	for (int i = 0; i < this->rows; i++) {
		LNodes[i] = this->rowNodes[i].head;
	}
	for (int i = 0; i < this->cols; i++) {
		double* tmp = new double[rightMatrix.cols];
		for (int j = 0; j < rightMatrix.cols; j++) {
			tmp[j] = rightMatrix.getValueAt(i, j);
		}
		for (int j = 0; j < this->rows; j++) {
			if (LNodes[j] != nullptr && LNodes[j]->col == i) {
				double rv = LNodes[j]->value;
				ColNode* prev = nullptr;
				ColNode* CNodes = result.rowNodes[j].head;
				for (int k = 0; k < rightMatrix.cols; k++) {
					if (tmp[k] != 0) {
						for (; CNodes != nullptr && CNodes->col < k; prev = CNodes, CNodes = CNodes->next);
						if (CNodes == nullptr || CNodes->col > k) {
							ColNode* new_node = new ColNode(k, tmp[k] * rv);
							if (prev == nullptr) {
								result.rowNodes[j].head = new_node;
								new_node->next = CNodes;
							} else {
								prev->next = new_node;
								new_node->next = CNodes;
							}
							prev = new_node;
							CNodes = prev->next;
						} else {
							CNodes->value += tmp[k] * rv;
							prev = CNodes;
							CNodes = CNodes->next;
						}
					}
				}
				LNodes[j] = LNodes[j]->next;
			}
		}
	}
}

void CSRMatrix::Transpose(CSRMatrix& result){
	ColNode** LNodes = new ColNode * [result.rows];
	for (int i = 0; i < this->rows; i++) {
		LNodes[i] = result.rowNodes[i].head;
	}
	for(int i = 0; i < this->rows; i++){
		for(ColNode* rowHead = this->rowNodes[i].head; rowHead != nullptr; rowHead = rowHead->next){
			ColNode* new_node = new ColNode(i, rowHead->value);
			int ori_col = rowHead->col;
			if(LNodes[ori_col] == nullptr){
				result.rowNodes[ori_col].head = new_node;
				LNodes[ori_col] = result.rowNodes[ori_col].head;
			} else {
				LNodes[ori_col]->next = new_node;
				LNodes[ori_col] = LNodes[ori_col]->next;
			}
		}
	}
}

RowNode CSRMatrix::getRowHead(int rowId) const {
	return this->rowNodes[rowId];
}

void CSRMatrix::setRowHead(int rowId, ColNode* rowHead) {
	this->rowNodes[rowId].head = rowHead;
}