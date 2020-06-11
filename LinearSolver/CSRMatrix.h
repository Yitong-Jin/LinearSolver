#pragma once
#include "Matrix.h"

/**
	Column Node of CSRMatrix.
	The next points to the next column node on the same row.
*/
struct ColNode
{
	int col;	///< index of the column
	double value;	///< value of the matrix
	ColNode *next;	///< pointer to the next column node on the same row
	ColNode(int c, double v) {
		col = c;
		value = v;
		next = nullptr;
	}
};

/**
	Row Node of CSRMatrix.
	The head points to the first column node.
*/
struct RowNode
{
	ColNode *head;	///< pointer to the first column node on the row
	RowNode() {
		head = nullptr;
	}
};

/**
	The class CSRMatrix is designed for a sparse matrix which aims at
	storing the value of a Matrix and providing some basic operations of matrices.
	The CSRMatrix is improved using linked lists to store values on every row,
	which performs better than the one using basic arraies.
*/
class CSRMatrix : public Matrix {
public:
	/** Constructor of CSRMatrix

		Create a Sparse Matrix object with empty values.

		@param rows number of rows in the matrix
		@param cols number of columns in the matrix
	*/
	CSRMatrix(int rows, int cols);

	/** Constructor of CSRMatrix

		Create a Sparse Matrix object with given values.

		@param rows			number of rows in the matrix
		@param cols			number of columns in the matrix
		@param values		the 1-d array of values
		@param row_position	the 1-d array of row position
		@param col_Index	the 1-d array of corresponding column index
	*/
	CSRMatrix(int rows, int cols, double* values, int* row_position, int* col_Index);

	/** Constructor of CSRMatrix

		Create a Sparse Matrix object with given integer values.

		@param rows			number of rows in the matrix
		@param cols			number of columns in the matrix
		@param values		the 1-d array of values
		@param row_position	the 1-d array of row position
		@param col_Index	the 1-d array of corresponding column index
	*/
	CSRMatrix(int rows, int cols, int* values, int* row_position, int* col_Index);

	/** Constructor of CSRMatrix

		Create a Sparse Matrix object from a dense matrix.

		@param matrix	the basic matrix object
	*/
	CSRMatrix(Matrix& matrix);

	/** Deconstructor of CSRMatrix
	*/
	~CSRMatrix() override;

	/** Set a single value at a position

		@param row		the row index of the position
		@param col		the column index of the position
		@param value	the value to be set
	*/
	void setValueAt(int row, int col, double value) override;

	/** Retreive a value at a position

		@param row		the row index of the position
		@param col		the column index of the position
		@return			return the value at position row row, col column
	*/
	double getValueAt(int row, int col) const override;

	/** Conduct one row/column operation in sparse matrix

		Each value of the target row/column will be updated by
		Line[target] = Line[left] * leftCoef + Line[target] * targetCoef
		where the `Line` is row or column controlled by parameter col.

		@param left			the row/column index of the line
		@param leftCoef		the coefficient applied to the so-called left line
		@param target		the row/column index of the target line
		@param targetCoef	the coefficient applied to the target line
		@param col			a boolean value indicating whther to use column as the line
	*/
	void rocOp(int left, double leftCoef, int target, double targetCoef, bool col, int start) override;

	/** Swap two rows/columns in sparse matrix

		Swap each value of the two rows/columns.

		@param index1		the row/column index of the one line
		@param index2		the row/column index of the other line
		@param col			a boolean value indicating whther to use column as the line
		@param start		a integer value indicating the start position of operation
	*/
	void swap(int index1, int index2, bool col) override;

	/** Sparse Matrix Addition

		@param rightMatrix	the matrix to add
		@param result		the result of the addition
	*/
	void matAdd(const CSRMatrix& rightMatrix, CSRMatrix& result);

	/** Sparse Matrix Substraction

		@param rightMatrix	the matrix to minus
		@param result		the result of the substraction
	*/
	void matMinus(const CSRMatrix& rightMatrix, CSRMatrix& result);

	/** Sparse Matrix Multiplication

		@param rightMatrix	the matrix to mupltiply
		@param result		the result of the multiplication
	*/
	void matMul(const CSRMatrix& rightMatrix, CSRMatrix& result);

	/** Transpose the sparse matrix

		@param result		the transpose of the matrix
	*/
	void Transpose(CSRMatrix& result);

	/** Output the matrix

		Print each value of the matrix by the order of index.
	*/
	void printMatrix() const override;

	/** Retreive the row node at the certain row.
	*/
	RowNode getRowHead(int rowId) const;

	/** Set the row head node of certain row.
	*/
	void setRowHead(int rowId, ColNode* rowHead);

private:
	RowNode* rowNodes; ///< the basic array storing the head of the linked list of every row.
};

