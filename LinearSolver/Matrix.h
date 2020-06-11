#pragma once

/**
	The class Matrix is designed for a dense matrix which aims at 
	storing the value of a Matrix and providing some basic
	operations of matrices.
*/
class Matrix {
public:
	/** Constructor of Matrix
		
		Create a Matrix object with empty values.

		@param rows number of rows in the matrix
		@param cols number of columns in the matrix
	*/
	Matrix(int rows, int cols);

	/** Constructor of Matrix

		Create a Matrix object with given values.

		@param rows number of rows in the matrix
		@param cols number of columns in the matrix
		@param values_ptr the 1-d array of values
	*/
	Matrix(int rows, int cols, double* values_ptr);

	/** Constructor of Matrix

		Create a Matrix object with given integer values.

		@param rows number of rows in the matrix
		@param cols number of columns in the matrix
		@param values_ptr the 1-d array of values
	*/
	Matrix(int rows, int cols, int* values_ptr);

	/** Deconstructor of Matrix
	*/
	virtual ~Matrix();

	/** Set a single value at a position

		@param row		the row index of the position
		@param col		the column index of the position
		@param value	the value to be set
	*/
	virtual void setValueAt(int row, int col, double value);

	/** Retreive a value at a position
	
		@param row		the row index of the position
		@param col		the column index of the position
		@return			return the value at position row row, col column
	*/
	virtual double getValueAt(int row, int col) const;

	/** Conduct one row/column operation

		Each value of the target row/column will be updated by
		Line[target] = Line[left] * leftCoef + Line[target] * targetCoef
		where the `Line` is row or column controlled by parameter col.

		@param left			the row/column index of the line
		@param leftCoef		the coefficient applied to the so-called left line
		@param target		the row/column index of the target line
		@param targetCoef	the coefficient applied to the target line
		@param col			a boolean value indicating whther to use column as the line
		@param start		a integer value indicating the start position (either in row or column) of operation
	*/
	virtual void rocOp(int left, double leftCoef, int target, double targetCoef, bool col, int start);

	/** Swap two rows/columns

		Swap each value of the two rows/columns.

		@param index1		the row/column index of the one line
		@param index2		the row/column index of the other line
		@param col			a boolean value indicating whther to use column as the line
	*/
	virtual void swap(int index1, int index2, bool col);

	/** Matrix Addition
	
		@param rightMatrix	the matrix to add
		@param result		the result of the addition
	*/
	virtual void matAdd(const Matrix& rightMatrix, Matrix& result);

	/** Matrix Substraction

		@param rightMatrix	the matrix to minus
		@param result		the result of the substraction
	*/
	virtual void matMinus(const Matrix& rightMatrix, Matrix& result);

	/** Matrix Multiplication

		@param rightMatrix	the matrix to mupltiply
		@param result		the result of the multiplication
	*/
	virtual void matMul(const Matrix& rightMatrix, Matrix& result);

	/** Inverse the matrix

		@param result		the inversion of the matrix
	*/
	virtual bool inverse(Matrix& result);

	/** Transpose the matrix

		@param result_T		the transpose of the matrix
	*/
	void Transpose(Matrix& result_T);

	/** Apply row operations to produce the triangular upper matrix using Gaussian Elimination Method

		@param RHS			the RHS of the linear system
		@param result		the triangular upper matrix
	*/
	virtual bool GauEli(Matrix& RHS, Matrix& result);

	/** Output the matrix

		Print each value of the matrix by the order of index.
	*/
	virtual void printMatrix() const;

	/** Retreive the number of rows in the matrix
	*/
	int getRows() const;

	/** Retreive the number of columns in the matrix
	*/
	int getCols() const;

private:
	double* values_ptr; ///< values of the matrix

protected:
	bool preallocate; ///< preallocated or not
	int rows;	///< number of rows in the matrix
	int cols;	///< number of columns in the matrix
};
