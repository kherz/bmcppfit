//!  BlochMcConnellSolver.h 
/*!
Solver function for BlochMcConnell equation

Kai Herz, 2019
kai.herz@tuebingen.mpg.de

**********************************
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
**********************************

*/

#pragma once
#include <Eigen/Eigen>
using namespace Eigen;


//! Solver function for BlochMcConnell equation
/*!
	M = (Mi + A^-1 * C) * exp(A*t) - A^-1 * C
	The matrix exponent is calculated with the Pade approximation
	see: Golub and Van Loan, Matrix Computations, Algorithm 11.3-1.
	Template function for arbitryry number of pools
	\param Mi Nx1 init magnetization vector
	\param A Nx7 matrix with relaxation and exchange rates + pulse parameters
	\param C Nx1 vector with the relaxation parameters
	\param t duration of the simulation [s]
	\param numApprox number of approximations for pade method
	\return the resulting new Nx1 magnetization vector
*/
template<int size> Matrix<double, size, 1> SolveBlochEquation(const Matrix<double, size, 1> &Mi, const Matrix<double, size, size> &A, const Matrix<double, size, 1> &C, double& t, int numApprox = 6)
{
	typedef Matrix<double, size, 1> VectorNd; // 
	typedef Matrix<double, size, size> MatrixNd; // 

	VectorNd AInvT = A.inverse()*C; // helper variable A^-1 * C
	MatrixNd At = A * t;			// helper variable A * t
	//solve exponential with pade method
	int infExp; //infinity exponent of the matrix
	int j;
	std::frexp(At.template lpNorm<Infinity>(), &infExp); // pade method is only stable if ||A||inf / 2^j <= 0.5
	j = std::max(0, infExp + 1);
	At = At * (1.0 / (pow(2, j)));
	//the algorithm usually starts with D = X = N = Identity and c = 1
	// since c is alway 0.5 after the first loop, we can start in the second round and init the matrices corresponding to that
	MatrixNd X(At); // X = A after first loop
	double c = 0.5; // c = 0.5 after first loop
	MatrixNd N(At);
	N.setIdentity();
	MatrixNd D = N - c * At;
	N += c * At;
	bool p = true; // D +- cX is dependent from (-1)^k, fastest way is with changing boolean in the loop
	double q = numApprox;
	MatrixNd cX; // helper variable for c * X
	// run the approximation
	for (int k = 2; k <= q; k++)
	{
		c = c * (q - k + 1) / (k*(2 * q - k + 1));
		X = At * X;
		cX = c * X;
		N = N + cX;
		if (p)
			D = D + cX;
		else
			D = D - cX;
		p = !p;
	}
	MatrixNd F = D.inverse()*N; // solve D*F = N for F
	for (int k = 1; k <= j; k++)
	{
		F *= F;
	}
	return F * (Mi + AInvT) - AInvT;
}