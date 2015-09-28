/*                               -*- C++ -*-
 * Copyright (C) 2014 Felix Salfelder
 * Author: Felix Salfelder
 *
 * This file is part of "Gnucap", the Gnu Circuit Analysis Package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * gsl supplementary stuff
 */
#include "u_gsl.h"

// int dgelss_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int
//     *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int
//     *info);
//
extern "C" {
	/*
		SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND,
		RANK, WORK, LWORK, INFO )

		INTEGER        INFO, LDA, LDB, LWORK, M, N, NRHS, RANK

		DOUBLE         PRECISION RCOND

		DOUBLE         PRECISION A( LDA, * ), B( LDB, * ), S( *), WORK( * )
		*/
	extern void zgelss_(integer *m, integer *n, integer *nrhs,
			double *a, integer *lda, double *b, integer *ldb,
			double *s, double *rcond,
			integer *rank, double *work,
			integer *lwork, double *rwork,
			integer *info);
}

int
gsl_blas_zgelss ( gsl_matrix_complex * A,
		gsl_matrix_complex * B,
		gsl_vector * S,
		double* RCOND, integer* RANK, integer* INFO )
{
	// M: rows in A
	// N: cols in A
	// NRHS: number of cols in B and X
	integer M = (integer)A->size2;
	integer N = (integer)A->size1;
	integer NRHS = (integer)B->size1;

	// if (!NRHS) {
	// 	GSL_ERROR ("no rhs?", GSL_EBADLEN);
	// }else if((unsigned)M>S->size && (unsigned)N>S->size){
	// 	GSL_ERROR ("S too small rhs?", GSL_EBADLEN);
	// }else

	{


		integer Atda = (integer)A->tda;
		integer Btda = (integer)B->tda;
		integer LWORK = -1;
		double WORK_[2];
		trace5("zgelss", M, N, NRHS, Atda, Btda);
		zgelss_(&M, &N, &NRHS, A->data, &Atda,
				B->data, &Btda, S->data, RCOND,
				RANK, WORK_, &LWORK, NULL, INFO );
		LWORK = (int)WORK_[0];
		double WORK[LWORK*2];
		double RWORK[std::max(5*std::min(M,N)-4,(integer)1)];

		trace1("zgelss input matrix\n", *A);
		trace1("zgelss input vectors\n", *B);
		zgelss_(&M, &N, &NRHS, A->data, &Atda,
				B->data, &Btda, S->data, RCOND,
				RANK, WORK, &LWORK, RWORK, INFO );
		trace1("zgelss done\n", *B);
		return GSL_SUCCESS;
	}
}
