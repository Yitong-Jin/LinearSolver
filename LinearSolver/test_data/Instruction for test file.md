## Instruction for test file

For testing our linear solver, we wrote a test function which can check the consistency between the result from our solvers and that from the numpy.linalg.solve( ) function. 

From the txt file including several linear systems $A\boldsymbol{x}=\boldsymbol{b}$ , the test function read in the linear systems which have been solved by numpy solver. Every linear systems should be separated by a line break. And in each linear system, the data should be stored in the specific format. For a $n\times n$ LHS matrix $A$, the whole system is stored in ($n$+3) lines in the txt file. The first line should be the number of rows and number of cols, separated by a space. The data of LHS matrix $A$ is stored from the second line to the ($n$+1)th line, corresponding to each line of the $A$ from top to bottom, while the data of RHS vector $\boldsymbol{b}$ and vector $\boldsymbol{x}$ is stored in the ($n$+2)th and ($n$+3)th line separately and every elements should be separated by a space. 

For example, a $4\times 4$ linear system: $A\boldsymbol{x}=\boldsymbol{b}$

$$
\left(
  \begin{array}{rr}
    9929 & 924 & 948 & 537\\
    924 & 9543 & 94 & 270\\
    948 & 94 & 9975 & 120 \\
    537 & 270 & 120 & 9867
  \end{array}
\right)\left(
  \begin{array}{c}
    x1 \\
    x2 \\
    x3 \\
    x4 
  \end{array}
\right) = \left(
  \begin{array}{c}
    329 \\
    570 \\
    883 \\
    613   
  \end{array}
\right)
$$

Here, we need to get the exact value of vector $\boldsymbol{x}$ by numpy.linalg.solve( ) in advance, because for testing reasons, we want to check whether or not the result we obtain by our solvers are consistent with we get from numpy solver.

In this case, vector  $\boldsymbol{x}$ is:

$$
x =\left(
  \begin{array}{c}
    0.02877076 \\
    0.05368908 \\
    0.09264927 \\
    0.06115476 
  \end{array}
 \right)
 $$

 Above linear system in the txt file to be tested by our test function should in the format as follows:

 $$
\begin{array}{rr}
    4 & 4 \\
    9929 & 924 & 948 & 537\\
    924 & 9543 & 94 & 270\\
    948 & 94 & 9975 & 120 \\
    537 & 270 & 120 & 9867 \\
    329 & 570 & 883 & 613 \\
    0.02877076 & 0.05368908 & 0.09264927 & 0.06115476
  \end{array}$$

After reading in the data, we can call different solvers to solve all of the linear systems in the file in a row and compare the results with numpy solver. When the difference of two results is smaller than a certain accuracy which we can set manually, the test function will print "$Correct$", otherwise, it will print "$Wrong$".