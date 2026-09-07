
//////////////////////////////////////////////////////////////////////////////
//
//  $Id: matrix3d.h,v 1.1 2026/09/07 20:50:35 ericthered Exp $
//
//  Header library for 3D math operations related to Catesian Matrices (x,y,z)
//
//  IMPORTANT:  The assumption is that the operations are applied
//	        to unit-length vectors with start point at Origin (0,0,0)
//
//////////////////////////////////////////////////////////////////////////////


//	Dot Product
double m3d_Rdot( double U[3], double V[3] ){
	return U[0]*V[0] + U[1]*V[1] + U[2]*V[2] ;
} ;

//	Length of a vector
double m3d_Rlngt( double U[3] ){
	return	sqrt( m3d_Rdot( U, U ) ) ;
} ;;

//	Scalar multiple of input vector
void m3d_Vscale( double factor, double U[3], double R[3] ){
	R[0] = factor * U[0] ;
	R[1] = factor * U[1] ;
	R[2] = factor * U[2] ;
} ;

//	Vector Addition
void m3d_Vadd( double U[3], double V[3], double R[3] ){
	R[0] = U[0]+V[0] ;
	R[1] = U[1]+V[1] ;
	R[2] = U[2]+V[2] ;
} ;

//	Vector Subtraction
void m3d_Vsubtr( double U[3], double V[3], double R[3] ){
	R[0] = U[0]-V[0] ;
	R[1] = U[1]-V[1] ;
	R[2] = U[2]-V[2] ;
} ;

//	Cross-product of 2 vectors
void m3d_Vcross( double U[3], double V[3], double R[3] ){
	R[0] = U[1]*V[2] - U[2]*V[1] ;
	R[1] = U[2]*V[0] - U[0]*V[2] ;
	R[2] = U[0]*V[1] - U[1]*V[0] ;
} ;

//	Vector Centroid (midpoint)
void m3d_VgetCentroid( double U[3], double V[3], double W[3], double R[3] ){
	R[0] = ( U[0] + V[0] + W[0] )/3 ;
	R[1] = ( U[1] + V[1] + W[1] )/3 ;
	R[2] = ( U[2] + V[2] + W[2] )/3 ;
} ;

//	Multiply a Matrix with a Vector
void m3d_VmultMatr( double W[3][3], double P[3], double R[3] ){
	// R[i] = SUM( W[ij] * P[j] )
	R[0] = W[0][0]*P[0] + W[0][1]*P[1] + W[0][2]*P[2] ;
	R[1] = W[1][0]*P[0] + W[1][1]*P[1] + W[1][2]*P[2] ;
	R[2] = W[2][0]*P[0] + W[2][1]*P[1] + W[2][2]*P[2] ;
} ;

//	Scalar multiple of input vector
void m3d_Mscale( double factor, double U[3][3], double R[3][3] ){
	R[0][0] = factor * U[0][0] ;
	R[0][1] = factor * U[0][1] ;
	R[0][2] = factor * U[0][2] ;

	R[1][0] = factor * U[1][0] ;
	R[1][1] = factor * U[1][1] ;
	R[1][2] = factor * U[1][2] ;

	R[2][0] = factor * U[2][0] ;
	R[2][1] = factor * U[2][1] ;
	R[2][2] = factor * U[2][2] ;
} ;

//	Matrix addition
void m3d_Madd( double U[3][3], double V[3][3], double R[3][3] ){
	R[0][0] = U[0][0] + V[0][0] ;
	R[0][1] = U[0][1] + V[0][1] ;
	R[0][2] = U[0][2] + V[0][2] ;

	R[1][0] = U[1][0] + V[1][0] ;
	R[1][1] = U[1][1] + V[1][1] ;
	R[1][2] = U[1][2] + V[1][2] ;

	R[2][0] = U[2][0] + V[2][0] ;
	R[2][1] = U[2][1] + V[2][1] ;
	R[2][2] = U[2][2] + V[2][2] ;
} ;

//	Matrix subtraction
void m3d_Msubtr( double U[3][3], double V[3][3], double R[3][3] ){
	R[0][0] = U[0][0] - V[0][0] ;
	R[0][1] = U[0][1] - V[0][1] ;
	R[0][2] = U[0][2] - V[0][2] ;

	R[1][0] = U[1][0] - V[1][0] ;
	R[1][1] = U[1][1] - V[1][1] ;
	R[1][2] = U[1][2] - V[1][2] ;

	R[2][0] = U[2][0] - V[2][0] ;
	R[2][1] = U[2][1] - V[2][1] ;
	R[2][2] = U[2][2] - V[2][2] ;
} ;

//	Transpose of matrix
void m3d_Mtranspose( double U[3][3], double R[3][3] ){
	R[0][0] = U[0][0] ;
	R[0][1] = U[1][0] ;
	R[0][2] = U[2][0] ;

	R[1][0] = U[0][1] ;
	R[1][1] = U[1][1] ;
	R[1][2] = U[2][1] ;

	R[2][0] = U[0][2] ;
	R[2][1] = U[1][2] ;
	R[2][2] = U[2][2] ;
} ;

//	Matrix Multiplication
void m3d_Mmult( double U[3][3], double V[3][3], double R[3][3] ){
	// R[ij] = SUM( U[ij] * V[ji] )
	R[0][0] = U[0][0]*V[0][0] + U[0][1]*V[1][0] + U[0][2]*V[2][0] ;
	R[0][1] = U[0][0]*V[0][1] + U[0][1]*V[1][1] + U[0][2]*V[2][1] ;
	R[0][2] = U[0][0]*V[0][2] + U[0][1]*V[1][2] + U[0][2]*V[2][2] ;

	R[1][0] = U[1][0]*V[0][0] + U[1][1]*V[1][0] + U[1][2]*V[2][0] ;
	R[1][1] = U[1][0]*V[0][1] + U[1][1]*V[1][1] + U[1][2]*V[2][1] ;
	R[1][2] = U[1][0]*V[0][2] + U[1][1]*V[1][2] + U[1][2]*V[2][2] ;

	R[2][0] = U[2][0]*V[0][0] + U[2][1]*V[1][0] + U[2][2]*V[2][0] ;
	R[2][1] = U[2][0]*V[0][1] + U[2][1]*V[1][1] + U[2][2]*V[2][1] ;
	R[2][2] = U[2][0]*V[0][2] + U[2][1]*V[1][2] + U[2][2]*V[2][2] ;
} ;

//	Determinant of degree 3 Matrix
double m3d_Rdeterminant(double U[3][3] ){
	// Matrix components U[i][j]
	return   ( U[0][0] * ( U[1][1]*U[2][2] - U[1][2]*U[2][1] ) +
		   U[0][1] * ( U[1][2]*U[2][0] - U[1][0]*U[2][2] ) +
		   U[0][2] * ( U[1][0]*U[2][1] - U[1][1]*U[2][0] ) ) ;;
} ;




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



//	Normalize to unit vector
void m3d_Vunit( double U[3], double R[3] ){
	// returns the unit vector as R[]
	double factor = 1/m3d_Rlngt(U) ;
	m3d_Vscale( factor, U, R ) ;
} ;

//	Square of distance between 2 points
double m3d_RdistSquare( double U[3], double V[3] ){
	double W[3] ;
	m3d_Vsubtr( U, V, W ) ;
	return	m3d_Rdot( W, W ) ;
} ;

//	Distance between 2 points
double m3d_Rdist( double U[3], double V[3] ){
	double W[3] ;
	m3d_Vsubtr( U, V, W ) ;
	return	sqrt( m3d_Rdot( W, W ) ) ;
} ;

//	Define matrix for transformation by rotation about
//	a specified vector for a specified angle or rotation
void m3d_MrotTransf( double U[3], double rot, double R[3][3] ){
	// Applying the Rogrigues rotation formula for matrices
	// Vector U[] is axis of rotation
	// rot is the angle of rotation
	// Matrix R[] returned is the transformation/rotation matrix for the applicable angle
	double c = cos( rot ) ;
	double s = sin( rot ) ;
	double W[3] ;

	m3d_Vunit( U, W ) ;

	R[0][0] = W[0]*W[0]*(1-c) + c ;
	R[0][1] = W[0]*W[1]*(1-c) - W[2]*s ;
	R[0][2] = W[0]*W[2]*(1-c) + W[1]*s ;

	R[1][0] = W[1]*W[0]*(1-c) + W[2]*s ;
	R[1][1] = W[1]*W[1]*(1-c) + c ;
	R[1][2] = W[1]*W[2]*(1-c) - W[0]*s ;

	R[2][0] = W[2]*W[0]*(1-c) - W[1]*s ;
	R[2][1] = W[2]*W[1]*(1-c) + W[0]*s ;
	R[2][2] = W[2]*W[2]*(1-c) + c ;
} ;

//	Rotate a point in space around a unit vector radiating from the Origin
void m3d_VrotVect( double P[3], double U[3], double rot, double R[3] ){
	// Matrix R[] returned is the point rotated about the provided
	// axis vector by the specified angle of rotation
	double W[3][3] ;
	m3d_MrotTransf( U, rot, W ) ;
	m3d_VmultMatr( W, P, R ) ;
} ;

//	Angle between 2 intersecting vectors
double m3d_Rangle_T2T( double U[3], double V[3] ){
	// Returns the value of the angle as measured between the heads of 2 vectors
	return  acos( m3d_Rdot( U, V ) / ( m3d_Rlngt(U) * m3d_Rlngt(V) ) ) ; 
} ;

///////////////////////////////////////

