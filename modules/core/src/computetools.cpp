
#include "core.hpp"

#define SETSIZE 4			//FIXME: Set these as externs!!
#define COST_TRAJ_SIZE 4	//FIXME: Set these as externs!!
#define MERGE_COST_SIZE 6	//FIXME: Set these as externs!!

namespace lpt {

double calculateDistance(const lpt::Particle& P1, const lpt::Particle& P2) {
	double diffx = (P1.x - P2.x);
	double diffy = (P1.y - P2.y);
	double diffz = (P1.z - P2.z);
	return sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
}

double calculateDistance(const lpt::Particle::Ptr P1, const lpt::Particle::Ptr P2) {
	double diffx = (P1->x - P2->x);
	double diffy = (P1->y - P2->y);
	double diffz = (P1->z - P2->z);
	return sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
}

double calculateDistance(const lpt::ParticleImage::Ptr p1, const lpt::ParticleImage::Ptr p2) {
	double diffx = (p1->x - p2->x);
	double diffy = (p1->y - p2->y);
	return sqrt(diffx * diffx + diffy * diffy);
}

void convertCameraCoordinatesToWorld(const lpt::Camera& cam, const array<double, 3>& pc, array<double, 3>& pw) {
	pw[0] = cam.R[0][0] * (pc[0] - cam.T[0]) + cam.R[1][0] * (pc[1] - cam.T[1]) + cam.R[2][0] * (pc[2] - cam.T[2]);
	pw[1] = cam.R[0][1] * (pc[0] - cam.T[0]) + cam.R[1][1] * (pc[1] - cam.T[1]) + cam.R[2][1] * (pc[2] - cam.T[2]);
	pw[2] = cam.R[0][2] * (pc[0] - cam.T[0]) + cam.R[1][2] * (pc[1] - cam.T[1]) + cam.R[2][2] * (pc[2] - cam.T[2]);
}

void transposeMatrix(double A[3][3] , double A_T[3][3])  {
     for(int i = 0; i < 3; ++i)  {
        for(int j = 0; j < 3; ++j)  {
            A_T[i][j] = A[j][i];
        }
    }
}

void multiplyMatrix(double A[3][3], double B[3][3], double C[3][3])  {
    for (int i = 0; i < 3; ++i)  {
        for (int j = 0; j < 3; ++j)  {
            C[i][j] = 0;
            for (int k = 0; k < 3; ++k)  {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

double Regression::eval(double *coeff, int degree, double x) {
	double y = 0;
	for(int i = 0; i < degree; i++)
		y += coeff[i] * (pow(x,i));
	return y;
}

double Regression::mult(double *A, double *B, size_t m) {
	int i;
	double C=0;

	for (i=0; i<m; i++)
		C += A[i]*B[i]; 

	return C;
}

double Regression::Norm(double *A, size_t m) {
	double norm, sum = 0;
	for(int i = 0; i < m; i++)
		sum += A[i] * A[i];

	norm = sqrt(sum);

	return norm;
}
void Regression::Householder(double **A, double *b, size_t m, size_t n) {
	//Transforms the Vandermonde matrix into upper diagonal form
	double* e, *vk, *aj, *ak;
	double alphak, Bk, yj, norm, sign;

	e = new double[m];
	vk = new double[m];
	aj = new double[m];
	ak = new double[m];

	for(int i = 0; i < m; i++) {
		e[i] = 0;
		vk[i] = 0;
	}

	for (int k = 0; k < n; k++) {
		e[k] = 1.0;

		if(A[k][k] > 0)
			sign = 1;
		else
			sign = -1;

		for(int i = k; i < m; i++) {
			vk[i] = A[i][k];
			ak[i] = A[i][k];
		}		
		norm = Norm(ak, m);
		alphak = -1.0 * sign * norm;
		e[k] = alphak;
		vk[k] -= e[k];
		Bk = mult(vk, vk, m);
		if (Bk != 0) {
			for(int j = k; j < n; j++) {
				for(int i = 0; i < m; i++)
					aj[i] = A[i][j];

				yj = mult(vk, aj, m);

				for(int i = 0; i < m; i++)
					A[i][j] -= (2.0 * yj / Bk) * vk[i];
			}
			yj = mult(vk, b, m);
			for(int i = 0; i < m; i++)
				b[i] -= (2 * yj / Bk) * vk[i];
		}
		ak[k] = 0;
		e[k] = 0;
		vk[k] = 0;
	}
	delete[] e;
	delete[] vk;
	delete[] aj;
	delete[] ak;
}

void Regression::Backsub(double** A, double* b, size_t n, double* X) {
	//solves the upper triangular matrix A using back substitution
	for (int j = static_cast<int>(n)-1; j >= 0; j--) {
		if(A[j][j] == 0) {
			cout << "matrix is singular \n";
			break;
		}
		X[j] = b[j] / A[j][j];
		for (int i = 0; i < j; i++)
			b[i] = b[i] - A[i][j] * X[j];
	}

}
void Regression::Polyfit(double* t, double* y, int m, int degree, double* X) {
	//double* coeff=new double[degree];               //degree = order of the polynomial plus 1;
	double** A = new double*[m];
	double* b = new double[m];

	//----forming the Vandermonde matrix-----  (NOTE: This only works upto degree = 3  !!!)
	for (int i = 0; i < m; i++) {
		b[i] = y[i];
		A[i] = new double[degree];
		for (int j = 0; j < degree; j++) {
			if(j == 0)
				A[i][j] = 1;
			if(j == 1)
				A[i][j] = t[i];
			if(j == 2)                    // NOTE:  ADD more if statements if Higher order regression is needed
				A[i][j] = t[i] * t[i];
		}
	}	
	Householder(A, b, m, degree);                  //Transform A into upper triangular form
	Backsub(A, b, degree, X);                    //Solve for regression coefficients X using back substitution
	for(int i = 0; i < m; ++i)
		delete [] A[i];
	delete[] A;
	delete[] b;
	//X=vector of regression coef sorting from low to high polynomial order
}

////******************************************************************************
//double fastRegression::eval(double *coeff, double x)
//{
//	double y = 0.0f;
//	y += coeff[0];
//	y += coeff[1] * x;
//
//	return y;
//}
//
//double fastRegression::mult(double *A, double *B)
//{
//	double C = 0.0f;
//	
// 	for (unsigned int i = 0; i < COST_TRAJ_SIZE; i++)
//		C += A[i] * B[i]; 
//	
//	return C;
//}
//
//double fastRegression::Norm(double *A)
//{
//	double norm, sum = 0.0f;
//	for(unsigned int i = 0; i < COST_TRAJ_SIZE; i++)
//		sum += A[i] * A[i];
//
//	norm = sqrtf(sum);
//
//	return norm;
//}
//void fastRegression::Householder(double *A, double *b) 
//{	//Transforms the Vandermonde matrix into upper diagonal form
//	
//	double alphak, Bk, yj, norm, sign;
//	
//	double e[COST_TRAJ_SIZE];
//	double vk[COST_TRAJ_SIZE];
//	double aj[COST_TRAJ_SIZE];
//	double ak[COST_TRAJ_SIZE];
//	
//	for(unsigned int i = 0; i < COST_TRAJ_SIZE; i++){
//			e[i] = 0.0f;
//			vk[i] = 0.0f;
//	}
//    			
//	for (unsigned int k = 0; k < 2; k++){
//		e[k] = 1.0f;
//		
//		if(A[k * 2 + k] > 0)
//			sign = 1;
//		else
//			sign = -1;
//				
//		for(unsigned int i = k; i < COST_TRAJ_SIZE; i++){
//			vk[i] = A[i * 2 + k];
//			ak[i] = A[i * 2 + k];
//		}
//
//		//NORM START		
//		norm = Norm(ak);
//		//NORM END
//
//		alphak = -1.0f * sign * norm;
//		
//		e[k] = alphak;
//		vk[k] -= e[k];
//
//		//MULT START
//		Bk = mult(vk, vk);
//		//MULT END
// 
//		if (Bk != 0){
//			for(unsigned int j = k; j < 2; j++){
//				for(unsigned int i = 0; i < COST_TRAJ_SIZE; i++)
//					aj[i] = A[i * 2 + j]; 
//
//				//MULT START
//				yj = mult(vk, aj);
//				//MULT END
//
//				for(unsigned int i = 0; i < COST_TRAJ_SIZE; i++)
//					A[i * 2 + j] -= (2.0f * yj / Bk) * vk[i];
//			}
//		//MULT START
//		yj = mult(vk, b); 
//		//MULT END
//
//        	for(unsigned int i = 0; i < COST_TRAJ_SIZE; i++)
//			b[i] -= (2.0f * yj / Bk) * vk[i];
//		}
//
//		ak[k] = 0.0f;
//		e[k] = 0.0f;
//		vk[k] = 0.0f;
//	}
//
//}
//
//void fastRegression::Backsub(double *A, double *b, double *X)
//{	//solves the upper triangular matrix A using back substitution
//
//	/*if(A[j][j] == 0){
//	  printf("matrix is singular \n"); //need an error from CUDA..is there a flag we can set?
//	  break;
//	}*/
//
//	X[1] = b[1] / A[3];                    // A[1][1] = A[1 * 2 + 1] = A[3] 
//	b[0] = b[0] - A[1] * X[1];     // A[0][1] = A[0 * 2 + 1] = A[1]
//	X[0] = b[0] / A[0];            // A[0][0] = A[0 * 2 + 0] = A[0]	
//	   
//	
//
//}
//void fastRegression::Polyfit(double *t, double *y, double *X)
//{	                             //degree = order of the polynomial plus 1;
//	double A[COST_TRAJ_SIZE * 2];          //COST_TRAJ_SIZE rows and 2 columns
//	double b[COST_TRAJ_SIZE];
//	
//	//----forming the Vandermonde matrix-----  (NOTE: This only works up to degree = 2  !!!)
//	for (unsigned int row = 0; row < COST_TRAJ_SIZE; row++){
//		b[row] = y[row];
//		A[row * 2 + 0] = 1.0f;
//		A[row * 2 + 1] = t[row];   // NOTE:  ADD more statements if Higher order regression is needed
//	}
//
//
//	Householder(A,b);                  //Transform A into upper triangular form
//
//	Backsub(A,b,X);                    //Solve for regression coeficients X using back substitution
//
//}
//
////******************************************************************************
//double mergeRegression::eval(double *coeff, double x)
//{
//	double y = 0;
//	y += coeff[0];
//	y += coeff[1] * x;
//
//	return y;
//}
//
//double mergeRegression::mult(double *A, double *B)
//{
//	double C = 0;
//	
// 	for (unsigned int i = 0; i < (MERGE_COST_SIZE - 1); i++)
//		C += A[i] * B[i]; 
//	
//	return C;
//}
//
//double mergeRegression::Norm(double *A)
//{
//	double norm, sum = 0;
//	for(unsigned int i = 0; i < (MERGE_COST_SIZE - 1); i++)
//		sum += A[i] * A[i];
//
//	norm = sqrtf(sum);
//
//	return norm;
//}
//void mergeRegression::Householder(double *A, double *b) 
//{	//Transforms the Vandermonde matrix into upper diagonal form
//	
//	double alphak, Bk, yj, norm, sign;
//	
//	double e[(MERGE_COST_SIZE - 1)];
//	double vk[(MERGE_COST_SIZE - 1)];
//	double aj[(MERGE_COST_SIZE - 1)];
//	double ak[(MERGE_COST_SIZE - 1)];
//	
//	for(unsigned int i = 0; i < (MERGE_COST_SIZE - 1); i++){
//			e[i] = 0;
//			vk[i] = 0;
//	}
//    			
//	for (unsigned int k = 0; k < 2; k++){
//		e[k] = 1;
//		
//		if(A[k * 2 + k] > 0)
//			sign = 1;
//		else
//			sign = -1;
//				
//		for(unsigned int i = k; i < (MERGE_COST_SIZE - 1); i++){
//			vk[i] = A[i * 2 + k];
//			ak[i] = A[i * 2 + k];
//		}
//
//		//NORM START		
//		norm = Norm(ak);
//		//NORM END
//
//		alphak = -1 * sign * norm;
//		
//		e[k] = alphak;
//		vk[k] -= e[k];
//
//		//MULT START
//		Bk = mult(vk, vk);
//		//MULT END
// 
//		if (Bk != 0){
//			for(unsigned int j = k; j < 2; j++){
//				for(unsigned int i = 0; i < (MERGE_COST_SIZE - 1); i++)
//					aj[i] = A[i * 2 + j]; 
//
//				//MULT START
//				yj = mult(vk, aj);
//				//MULT END
//
//				for(unsigned int i = 0; i < (MERGE_COST_SIZE - 1); i++)
//					A[i * 2 + j] -= (2 * yj / Bk) * vk[i];
//			}
//		//MULT START
//		yj = mult(vk, b); 
//		//MULT END
//
//        	for(unsigned int i = 0; i < (MERGE_COST_SIZE - 1); i++)
//			b[i] -= (2 * yj / Bk) * vk[i];
//		}
//
//		ak[k] = 0;
//		e[k] = 0;
//		vk[k] = 0;
//	}
//
//}
//
//void mergeRegression::Backsub(double *A, double *b, double *X)
//{	//solves the upper triangular matrix A using back substitution
//
//	/*if(A[j][j] == 0){
//	  printf("matrix is singular \n"); //need an error from CUDA..is there a flag we can set?
//	  break;
//	}*/
//
//	X[1] = b[1] / A[3];                    // A[1][1] = A[1 * 2 + 1] = A[3] 
//	b[0] = b[0] - A[1] * X[1];     // A[0][1] = A[0 * 2 + 1] = A[1]
//	X[0] = b[0] / A[0];            // A[0][0] = A[0 * 2 + 0] = A[0]	
//	   
//	
//
//}
//void mergeRegression::Polyfit(double *t, double *y, double *X)
//{	                             //degree = order of the polynomial plus 1;
//	double A[(MERGE_COST_SIZE - 1) * 2];          //MERGE_COST_SIZE rows and 2 columns
//	double b[(MERGE_COST_SIZE - 1)];
//	
//	//----forming the Vandermonde matrix-----  (NOTE: This only works up to degree = 2  !!!)
//	for (unsigned int row = 0; row < (MERGE_COST_SIZE - 1); row++){
//		b[row] = y[row];
//		A[row * 2 + 0] = 1;
//		A[row * 2 + 1] = t[row];   // NOTE:  ADD more statements if Higher order regression is needed
//	}
//
//
//	Householder(A,b);                  //Transform A into upper triangular form
//
//	Backsub(A,b,X);                    //Solve for regression coefficients X using back substitution
//	
//}

} /* NAMESPACE_PT */



