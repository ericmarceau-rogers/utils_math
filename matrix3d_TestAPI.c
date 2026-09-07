
///////////////////////////////////////
//
//  $Id: matrix3d_TestAPI.c,v 1.1 2026/09/07 20:51:32 ericthered Exp $
//
//	API demo and test program
//
///////////////////////////////////////


//
// 3rd-Party Libraries
//
#include <stdio.h>
#include <math.h>


//
// Locally Authored Libraries
//
#include "matrix3d.h"
#include "matrix2d.h"

void reportRealValues( char label[40], double val ){
	printf("\n\t %25s %10.6f\n", label, val ) ;
} ;

void report3VectorValues( char label[4], double U[3] ){
	printf("\n\t %25s [ ", label ) ;
	for( int j = 0 ; j < 3 ; j++ ){
		printf(" %10.6f ", U[j] ) ;
	} ;
	printf(" ]\n") ;
}

void report3MatrixValues( char label[4], double U[3][3] ){
	printf("\n") ;
	for( int i = 0 ; i < 3 ; i++ ){
		printf("\t %25s [ ", label ) ;
		for( int j = 0 ; j < 3 ; j++ ){
			printf(" %10.6f ", U[i][j] ) ;
		} ;
		printf(" ]\n") ;
	} ;
}

void report2VectorValues( char label[4], double U[2] ){
	printf("\n\t %25s [ ", label ) ;
	for( int j = 0 ; j < 2 ; j++ ){
		printf(" %10.6f ", U[j] ) ;
	} ;
	printf(" ]\n") ;
}

void report2MatrixValues( char label[4], double U[2][2] ){
	printf("\n") ;
	for( int i = 0 ; i < 2 ; i++ ){
		printf("\t %25s [ ", label ) ;
		for( int j = 0 ; j < 2 ; j++ ){
			printf(" %10.6f ", U[i][j] ) ;
		} ;
		printf(" ]\n") ;
	} ;
}

int main(){

//
//	Required Variables
//
	double mult = 2.5 ;
	
	//
	//	Vectors
	//
	double D[3] = { 9, 8, 7 } ;
	double E[3] = { 1, 2, 3 } ;
	double F[3] = { 2, 3, 4 } ;
	double G[3] = { 0, 0, 0 } ;


	double L[3] = { 9, 8 } ;
	double M[3] = { 1, 2 } ;
	double N[3] = { 0, 0 } ;


	//
	//	Matrices
	//
	double P[3][3] = { { 1, 4, 8 } ,
			   { 2, 6, 9 } ,
			   { 3, 5, 7 } } ;

	double Q[3][3] = { { 4, 7, 1 } ,
			   { 5, 8, 2 } ,
			   { 6, 9, 3 } } ;

	double R[3][3] = { { 0, 0, 0 } ,
			   { 0, 0, 0 } ,
			   { 0, 0, 0 } } ;


	double H[2][2] = { { 1, 4 } ,
			   { 3, 5 } } ;

	double K[2][2] = { { 4, 7 } ,
			   { 6, 9 } } ;

	double S[2][2] = { { 0, 0 } ,
			   { 0, 0 } } ;


//
//	Testing of 3D Context Functions
//
	printf("\n\n Results for Testing of 3D Context Functions\n" ) ;

	report3VectorValues( "F : ", F ) ;

	reportRealValues( "m3d_Rdot : ", m3d_Rdot( E, F ) ) ;

	reportRealValues( "m3d_Rlngt : ", m3d_Rlngt( E ) ) ;

	m3d_Vscale( mult, E, G ) ;
	report3VectorValues( "m3d_Vscale : ", G ) ;

	m3d_Vadd( E, F, G ) ;
	report3VectorValues( "m3d_Vadd : ", G ) ;
	
	m3d_Vsubtr( E, F, G ) ;
	report3VectorValues( "m3d_Vsubtr : ", G ) ;
	
	m3d_Vcross( E, F, G ) ;
	report3VectorValues( "m3d_Vcross : ", G ) ;
	
	m3d_VgetCentroid( D, E, F, G ) ;
	report3VectorValues( "m3d_VgetCentroid : ", G ) ;

	m3d_VmultMatr( P, E, G ) ;
	report3VectorValues( "m3d_VmultMatr : ", G ) ;

	m3d_Mscale( mult, P, R ) ;
	report3MatrixValues( "m3d_Mscale : ", R ) ;

	m3d_Madd( P, Q, R ) ;
	report3MatrixValues( "m3d_Madd : ", R ) ;

	m3d_Msubtr( P, Q, R ) ;
	report3MatrixValues( "m3d_Msubtr : ", R ) ;

	m3d_Mtranspose( P, R ) ;
	report3MatrixValues( "m3d_Mtranspose : ", R ) ;

	m3d_Mmult( P, Q, R ) ;
	report3MatrixValues( "m3d_Mmult : ", R ) ;

	reportRealValues( "m3d_Rdeterminant : ", m3d_Rdeterminant( P ) ) ;

	m3d_Vunit( E, G ) ;
	report3VectorValues( "m3d_Vscale : ", G ) ;

	reportRealValues( "m3d_RdistSquare : ", m3d_RdistSquare( E, F ) ) ;

	reportRealValues( "m3d_Rdist : ", m3d_Rdist( E, F ) ) ;

	m3d_MrotTransf( E, 180.0, R ) ;
	report3MatrixValues( "m3d_MrotTransf : ", R ) ;

	m3d_VrotVect( F, E, 5.0, G ) ;
	report3VectorValues( "m3d_VrotVect : ", G ) ;

	reportRealValues( "m3d_Rangle_T2T : ", m3d_Rangle_T2T( E, F ) ) ;


	//
	//  Testing of 2D Context Functions
	//
	printf("\n\n Results for Testing of 2D Context Functions\n" ) ;


	report2VectorValues( "L : ", L ) ;

	reportRealValues( "m2d_Rdot : ", m2d_Rdot( L, M ) ) ;

	reportRealValues( "m2d_Rlngt : ", m2d_Rlngt( L ) ) ;

	m2d_Vscale( mult, L, N ) ;
	report2VectorValues( "m2d_Vscale : ", N ) ;

	m2d_Vadd( L, M, N ) ;
	report2VectorValues( "m2d_Vadd : ", N ) ;
	
	m2d_Vsubtr( L, M, N ) ;
	report2VectorValues( "m2d_Vsubtr : ", N ) ;
	
	//
	
	m2d_VgetCentroid( L, M, N ) ;
	report2VectorValues( "m2d_VgetCentroid : ", N ) ;

	m2d_VmultMatr( H, L, N ) ;
	report2VectorValues( "m2d_VmultMatr : ", N ) ;

	m2d_Mscale( mult, H, S ) ;
	report2MatrixValues( "m2d_Mscale : ", S ) ;

	m2d_Madd( H, K, S ) ;
	report2MatrixValues( "m2d_Madd : ", S ) ;

	m2d_Msubtr( H, K, S ) ;
	report2MatrixValues( "m2d_Msubtr : ", S ) ;

	m2d_Mtranspose( H, S ) ;
	report2MatrixValues( "m2d_Mtranspose : ", S ) ;

	m2d_Mmult( H, K, S ) ;
	report2MatrixValues( "m2d_Mmult : ", S ) ;

	reportRealValues( "m2d_Rdeterminant : ", m2d_Rdeterminant( H ) ) ;

	m2d_Vunit( L, N ) ;
	report2VectorValues( "m2d_Vscale : ", N ) ;

	reportRealValues( "m2d_RdistSquare : ", m2d_RdistSquare( L, M ) ) ;

	reportRealValues( "m2d_Rdist : ", m2d_Rdist( L, M ) ) ;

	m2d_MrotTransf( 180.0, S ) ;
	report2MatrixValues( "m2d_MrotTransf : ", S ) ;

	m2d_VrotVect( L, 5.0, N ) ;
	report2VectorValues( "m2d_VrotVect : ", N ) ;

	reportRealValues( "m2d_Rangle_T2T : ", m2d_Rangle_T2T( L, M ) ) ;




	return 0 ;
}
