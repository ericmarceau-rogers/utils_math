
//////////////////////////////////////////////////////////////////////////////
//
//  $Id: matrix2d.h,v 1.1 2026/09/07 20:49:28 ericthered Exp $
//
//  Header library for 2D math operations related to Catesian Matrices (x,y)
//
//  IMPORTANT:  The assumption is that the operations are applied
//	        to unit-length vectors with start point at Origin (0,0)
//
//////////////////////////////////////////////////////////////////////////////


//	Dot Product
double m2d_Rdot( double U[2], double V[2] ){
	return	U[0]*V[0] + U[1]*V[1] ;
} ;

//	Length of a vector
double m2d_Rlngt( double U[2] ){
	return	sqrt( m2d_Rdot( U, U ) ) ;
} ;;

//	Scalar multiple of input vector
void m2d_Vscale( double factor, double U[2], double R[2] ){
	R[0] = factor * U[0] ;
	R[1] = factor * U[1] ;
} ;

//	Vector Addition
void m2d_Vadd( double U[2], double V[2], double R[2] ){
	R[0] = U[0]+V[0] ;
	R[1] = U[1]+V[1] ;
} ;

//	Vector Subtraction
void m2d_Vsubtr( double U[2], double V[2], double R[2] ){
	R[0] = U[0]-V[0] ;
	R[1] = U[1]-V[1] ;
} ;

//	Vector Centroid (midpoint)
void m2d_VgetCentroid( double U[2], double V[2], double R[2] ){
	R[0] = ( U[0] + V[0] )/3 ;
	R[1] = ( U[1] + V[1] )/3 ;
} ;

//	Multiply a Matrix with a Vector
void m2d_VmultMatr( double W[2][2], double P[2], double R[2] ){
	// R[i] = SUM( W[ij] * P[j] )
	R[0] = W[0][0]*P[0] + W[0][1]*P[1] ;
	R[1] = W[1][0]*P[0] + W[1][1]*P[1] ;
} ;

//	Scalar multiple of input vector
void m2d_Mscale( double factor, double U[2][2], double R[2][2] ){
	R[0][0] = factor * U[0][0] ;
	R[0][1] = factor * U[0][1] ;

	R[1][0] = factor * U[1][0] ;
	R[1][1] = factor * U[1][1] ;
} ;


//	Matrix addition
void m2d_Madd( double U[2][2], double V[2][2], double R[2][2] ){
	R[0][0] = U[0][0] + V[0][0] ;
	R[0][1] = U[0][1] + V[0][1] ;

	R[1][0] = U[1][0] + V[1][0] ;
	R[1][1] = U[1][1] + V[1][1] ;
} ;

//	Matrix subtraction
void m2d_Msubtr( double U[2][2], double V[2][2], double R[2][2] ){
	R[0][0] = U[0][0] - V[0][0] ;
	R[0][1] = U[0][1] - V[0][1] ;

	R[1][0] = U[1][0] - V[1][0] ;
	R[1][1] = U[1][1] - V[1][1] ;
} ;

//	Transpose of matrix
void m2d_Mtranspose( double U[2][2], double R[2][2] ){
	R[0][0] = U[0][0] ;
	R[0][1] = U[1][0] ;

	R[1][0] = U[0][1] ;
	R[1][1] = U[1][1] ;
} ;

//	Matrix Multiplication
void m2d_Mmult( double U[2][2], double V[2][2], double R[2][2] ){
	// R[ij] = SUM( U[ij] * V[ji] )
	R[0][0] = U[0][0]*V[0][0] + U[0][1]*V[1][0] ;
	R[0][1] = U[0][0]*V[0][1] + U[0][1]*V[1][1] ;

	R[1][0] = U[1][0]*V[0][0] + U[1][1]*V[1][0] ;
	R[1][1] = U[1][0]*V[0][1] + U[1][1]*V[1][1] ;
} ;

//	Determinant of degree 3 Matrix
double m2d_Rdeterminant(double U[2][2] ){
	return U[0][0] * U[1][1] - U[0][1] * U[1][0] ;
} ;





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




//	Normalize to unit vector
void m2d_Vunit( double U[2], double R[2] ){
	// returns the unit vector as R[]
	double factor = 1/m2d_Rlngt(U) ;
	m2d_Vscale( factor, U, R ) ;
} ;

//	Square of distance between 2 points
double m2d_RdistSquare( double U[2], double V[2] ){
	double C[2] ;
	m2d_Vsubtr( U, V, C ) ;
	return	m2d_Rdot( C, C ) ;
} ;

//	Distance between 2 points
double m2d_Rdist( double U[2], double V[2] ){
	double C[2] ;
	m2d_Vsubtr( U, V, C ) ;
	return	sqrt( m2d_Rdot( C, C ) ) ;
} ;

//	Define matrix for transformation by rotation about Origin in XY plane
void m2d_MrotTransf( double rot, double R[2][2] ){
	// Origin is the point of rotation
	// rot is the angle of rotation
	// Matrix R[] returned is the transformation/rotation matrix for the applicable angle
	double c = cos( rot ) ;
	double s = sin( rot ) ;
	R[0][0] = c ;
	R[0][1] = -1.0 * s ;
	R[1][0] = s ;
	R[1][1] = c ;
} ;

//	Rotate a point in space around a unit vector radiating from the Origin
void m2d_VrotVect( double P[2], double rot, double R[2] ){
	// Matrix R[] returned is the point rotated about the provided
	// axis vector by the specified angle of rotation
	double W[2][2] ;
	m2d_MrotTransf( rot, W ) ;
	m2d_VmultMatr( W, P, R ) ;
} ;

//	Angle between 2 intersecting vectors
double m2d_Rangle_T2T( double U[2], double V[2] ){
	// Returns the value of the angle as measured between the heads of 2 vectors
	return  acos( m2d_Rdot( U, V ) / ( m2d_Rlngt(U) * m2d_Rlngt(V) ) ) ; 
} ;

///////////////////////////////////////

