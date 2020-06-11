#include <stdio.h>
#include <string.h>
#include <iostream>
#include "Matrix.h"
#include "Solver.h"
#include "CSRMatrix.h"
#include "Interface.h"

using namespace std;

int main(int argc, char** argv) {
	 Interface* stdInterface = new Interface();
	 stdInterface->start();
	 delete stdInterface;
	return 0;
}
