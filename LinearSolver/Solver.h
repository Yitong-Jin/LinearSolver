#pragma once
#include "Matrix.h"
#include "CSRMatrix.h"

/**
    The class Solver provides various algorithms to solve a given linear system.
    The algorithms support two types of matrix class: Matrix and CSRMatrix.
    For the performance of different solvers and types, see on GitHub.
*/
class Solver {
public:

    /** Use Jacobi to solve the linear systems
        
        (i!=j) 
        The computation of xi^(k+1) requires each element in xj^(k) except xi^(k).
        We can't overwrite xi^(k) with xi^(k+1).
        The amount of storage is two vectors of size n.

        @param A       the reference of the given (n x n) symmetric LHS matrix A
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */   
    static bool Jacobi(Matrix& A, Matrix& b, Matrix& x);

    /** Use Jacobi to solve the linear systems (sparse)
        
        (i!=j)
        The computation of xi^(k+1) requires each element in xj^(k) except xi^(k).
        We can't overwrite xi^(k) with xi^(k+1).
        The amount of storage is two vectors of size n.

        @param A       the reference of the given (n x n) symmetric LHS matrix A
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool Jacobi(CSRMatrix& A, Matrix& b, Matrix& x);

    /** Use Gauss_Seidel to solve the linear systems
        
        (i!=j)
        The computation of xi^(k+1) uses the elements of xj^(k+1) that have already been computed.
        We just overwrite xi^(k) with xi^(k+1).
        The amount of storage is one vector of size n.

        @param A       the reference of the given (n x n) symmetric LHS matrix A
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool Gauss_Seidel(Matrix& A, Matrix& b, Matrix& x);

    /** Use Gauss_Seidel to solve the linear systems (sparse)

        (i!=j)
        The computation of xi^(k+1) uses the elements of xj^(k+1) that have already been computed.
        We just overwrite xi^(k) with xi^(k+1).
        The amount of storage is one vector of size n.

        @param A       the reference of the given (n x n) symmetric LHS matrix A
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool Gauss_Seidel(CSRMatrix& A, Matrix& b, Matrix& x);
    
    /** Using the matrix inversion to solve linear systems
        
        Compute the inversion of LHS matrix A using Gaussian Elimination Method.
        Then, multiply the A^-1 by b.

        @param A       the reference of the given (n x n) LHS matrix A
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool GaussElimination(Matrix& A, Matrix& b, Matrix& x);
    
    /** Use LU decomposition (without using partial pivoting) to solve linear systems
    
        @param A       the reference of the given (n x n) LHS matrix A
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool LU(Matrix& A, Matrix& b, Matrix& x);

    /** Use LU decomposition (without using partial pivoting) to solve linear systems (sparse)

        @param A       the reference of the given (n x n) LHS matrix A in Sparse Matrix format
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool LU(CSRMatrix& A, Matrix& b, Matrix& x);

    /** Use LU decomposition (using partial pivoting) to solve linear systems
    
        @param A       the reference of the given (m x n) LHS matrix A
        @param b       the reference of the given (m x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool LU_pp(Matrix& A, Matrix& b, Matrix& x);

    /** Use LU decomposition (using partial pivoting) to solve linear systems (sparse)

        @param A       the reference of the given (n x n) LHS matrix A in Sparse Matrix format
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool LU_pp(CSRMatrix& A, Matrix& b, Matrix& x);

    /** Use Cholesky Fractorization to solve the linear systems
    
        @param A       the reference of the given (n x n) symmetric LHS matrix A
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool Cholesky(Matrix& A, Matrix& b, Matrix& x);

    /** Use Cholesky Fractorization to solve the linear systems (sparse)

        @param A       the reference of the given (n x n) symmetric LHS matrix A in Sparse Matrix format
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
	static bool Cholesky(CSRMatrix& A, Matrix& b, Matrix& x);

private:

    /** Develop a method which decomposes or factorises the LHS matrix A without using partial pivoting

        This decomposition involves a lower-(L) and an upper-(U) triangular matrix
        and in such a way that it is cheap to compute a new solution vector x for any given RHS vector b.
        Notice that after call this function, the A matrix becomes to the upper-(U) triangular matrix
        which was mentioned above, while the lower-(L) triangular matrix will be constructed by constructing a new object.

        @param A       the reference of the given (n x n) LHS matrix A
        @param L       the reference of an empty matrix in which we want to populate the values of the lower-(L) triangular matrix
        @return        indicating whether the linear system is solvable
    */
    static bool LU_decomposition(Matrix& A, Matrix& L);

    /** Develop a method which decomposes or factorises the LHS matrix A without using partial pivoting (sparse)

        This decomposition involves a lower-(L) and an upper-(U) triangular matrix
        and in such a way that it is cheap to compute a new solution vector x for any given RHS vector b.
        Notice that after call this function, the A matrix becomes to the upper-(U) triangular matrix
        which was mentioned above, while the lower-(L) triangular matrix will be constructed by constructing a new object.

        @param A       the reference of the given (n x n) LHS matrix A in Sparse Matrix format
        @param L       the reference of an empty matrix in which we want to populate the values of the lower-(L) triangular matrix
        @return        indicating whether the linear system is solvable
    */
    static bool LU_decomposition(CSRMatrix& A, Matrix& L);

    /** Function to perform backward subsitution on the system Ax=b

        This function assumes that A is already an upper triangular matrix.

        @param A       the reference of the upper triangular LHS matrix A (m x n)
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool backward_substitution(Matrix& A, Matrix& b, Matrix& x);

    /** Function to perform backward subsitution on the system Ax=b (sparse)

        This function assumes that A is already an upper triangular matrix.

        @param A       the reference of the upper triangular LHS matrix A (m x n) in Sparse Matrix format
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool backward_substitution(CSRMatrix& A, Matrix& b, Matrix& x);

    /** Function to perform forward subsitution on the system Ax=b

        This function assumes that A is already an lower triangular matrix.

        @param A       the reference of the lower triangular LHS matrix A (m x n)
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool forward_substitution(Matrix& A, Matrix& b, Matrix& x);

    /** Function to perform forward subsitution on the system Ax=b (sparse)

        This function assumes that A is already an lower triangular matrix.

        @param A       the reference of the lower triangular LHS matrix A (m x n) in Sparse Matrix format
        @param b       the reference of the given (n x 1) RHS vector b
        @param x       the reference of the (n x 1) vector x that we want to find
        @return        indicating whether the linear system is solvable
    */
    static bool forward_substitution(CSRMatrix& A, Matrix& b, Matrix& x);

    /** Develop a method which decomposes or factorises the LHS matrix A, using partial pivoting

        Same method with the function "LU_decomposition" but additionally keep track of row swaps
        which can help us appropriately avoid the problems of dividing by zeros and round off errors.
        Record each pivoting operation as simple matrix multiplication, and at the end these additional
        matrices can be 'pulled through' the matrix multiplication to form a single permutation matrix P.
        Here we call P as P_ to emphasise that this is the permutation matrix on the LHS.

        @param A       the reference of the given (n x n) LHS matrix A
        @param L       the reference of an empty matrix in which we want to populate the values of the lower-(L) triangular matrix
        @param P_      the reference of the single permutation matrix which keeps track of every pivoting operations.
        @return        indicating whether the linear system is solvable
    */
    static bool LU_decomposition_pp(Matrix& A, Matrix& L, Matrix& P_);

    /** Develop a method which decomposes or factorises the LHS matrix A, using partial pivoting

        Same method with the function "LU_decomposition" but additionally keep track of row swaps
        which can help us appropriately avoid the problems of dividing by zeros and round off errors.
        Record each pivoting operation as simple matrix multiplication, and at the end these additional
        matrices can be 'pulled through' the matrix multiplication to form a single permutation matrix P.
        Here we call P as P_ to emphasise that this is the permutation matrix on the LHS.

        @param A       the reference of the given (n x n) LHS matrix A in Sparse Matrix format
        @param L       the reference of an empty matrix in which we want to populate the values of the lower-(L) triangular matrix
        @param P_      the reference of the single permutation matrix which keeps track of every pivoting operations.
        @return        indicating whether the linear system is solvable
    */
    static bool LU_decomposition_pp(CSRMatrix& A, Matrix& L, Matrix& P_);

    /** The decomposition process of Cholesky Fractorization, using The Cholesky�CBanachiewicz and Cholesky�CCrout algorithms
    
        Compute the lower-(L) triangular matix which satisfies A = LL^T
        
        @param A       the reference of the given (n x n) symmetric LHS matrix A
        @param L       the reference of an empty matrix in which we want to populate the values of the lower-(L) triangular matrix
        @return        indicating whether the linear system is solvable
    */
    static bool CholeskyDecomposition(Matrix& A, Matrix& L);

    /** The decomposition process of Cholesky Fractorization in sparse matrix format, using The Cholesky�CBanachiewicz and Cholesky�CCrout algorithms

        Compute the lower-(L) triangular matix which satisfies A = LL^T

        @param A       the reference of the given (n x n) symmetric LHS matrix A in Sparse Matrix format
        @param L       the reference of an empty matrix in which we want to populate the values of the lower-(L) triangular matrix
        @return        indicating whether the linear system is solvable
    */
	static bool CholeskyDecomposition(CSRMatrix& A, Matrix& L);

};