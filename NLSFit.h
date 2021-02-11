//!  NLSFit.h
/*!
Contains the gsl fitting stuff

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

#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include "SimulationParameters.h"

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
	std::cout << "iter: " << iter << std::endl;
}

class NLSFit
{
public:

	struct FitPoint
	{
		double current;
		double lower;
		double upper;
	};

	NLSFit(int dp);
	~NLSFit();
	void AddFreeParameter(double init, double lower, double upper);
	void Init(int(*f) (const gsl_vector * x, void *data, gsl_vector * f), SimulationParameters &sp);
	void SetTolerances(double x, double g, double f);
	//void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);


private:
	const gsl_multifit_nlinear_type* type;
	gsl_multifit_nlinear_workspace* workspace;
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params;
	gsl_vector_view x;
	gsl_vector_view wts;
	gsl_vector *f;
	gsl_matrix *J;
	gsl_matrix *covar;

	const int numDataPoints;
	int numFreeParams;
	std::vector<FitPoint> fitParams;
	bool isInit;

	double xtol;
	double gtol;
	double ftol;

	double* x_init;
	double* weights;


};

NLSFit::NLSFit(int dp) : type(gsl_multifit_nlinear_trust),
                                 numDataPoints(dp)
{
	numFreeParams = 0;
	xtol = 1e-8;
	gtol = 1e-8;
	ftol = 0.0;
	fdf_params = gsl_multifit_nlinear_default_parameters();
	fdf_params.h_df = 0.8;
	isInit = false;
}

NLSFit::~NLSFit()
{
	if (isInit) {
		delete x_init;
		delete weights;
	}
}

void NLSFit::AddFreeParameter(double init, double lower, double upper)
{
	if (!isInit) {
		fitParams.push_back({ init, lower, upper });
		numFreeParams++;
	}
}


void NLSFit::Init(int(*f) (const gsl_vector * x, void *data, gsl_vector * f), SimulationParameters &sp)
{
	if (!isInit) {

		x_init = new double[numFreeParams];
		for (int i = 0; i < numFreeParams; i++)	{
			x_init[i] = fitParams.at(i).current;
		}
		weights = new double[numDataPoints];
		for (int i = 0; i < numDataPoints; i++) {
			weights[i] = 1.0;
		}
		x = gsl_vector_view_array(x_init, numFreeParams);
		wts = gsl_vector_view_array(weights, numDataPoints);
		fdf.f = f;
		fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
		fdf.fvv = NULL;     /* not using geodesic acceleration */
		fdf.n = numDataPoints;
		fdf.p = numFreeParams;
		fdf.params = &sp;


		/* allocate workspace with default parameters */
		workspace = gsl_multifit_nlinear_alloc(type, &fdf_params, numDataPoints, numFreeParams);

		/* initialize solver with starting point and weights */
		gsl_multifit_nlinear_winit(&x.vector, &wts.vector, &fdf, workspace);

		/* compute initial cost function */
		double chisq, chisq0;
		gsl_vector *f = gsl_multifit_nlinear_residual(workspace);
		gsl_blas_ddot(f, f, &chisq0);

		/* solve the system with a maximum of 100 iterations */
		int status, info;
		status = gsl_multifit_nlinear_driver(10, xtol, gtol, ftol,
			callback, NULL, &info, workspace);

		for (int i = 0; i < numFreeParams; i++) {
			fitParams.at(i).current = gsl_vector_get(workspace->x, i);
		}

		gsl_blas_ddot(f, f, &chisq);
	}
}


void NLSFit::SetTolerances(double x, double g, double f)
{
	xtol = x;
	gtol = g;
	ftol = f;
}

