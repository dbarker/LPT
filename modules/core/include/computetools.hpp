
#ifndef COMPUTETOOLS_H_
#define COMPUTETOOLS_H_

namespace lpt {

extern const int FRAMESPERSEC;

template <class Object_T> class RegressionModel;
template <class Object_T> class FiniteDiffModel;

typedef lpt::RegressionModel<lpt::Particle3d> Regress3d;
typedef lpt::FiniteDiffModel<lpt::Particle3d> FiniteDiff3d;

struct ComputeMode {
	enum {REGRESSION = 0, FINITEDIFF = 1};
};
struct Degree {
	enum {LINEAR = 1, QUADRADIC = 2};
};

double calculateDistance(const lpt::Particle& P1, const lpt::Particle& P2);
double calculateDistance(const lpt::Particle::Ptr P1, const lpt::Particle::Ptr P2);
double calculateDistance(const lpt::ParticleImage::Ptr p1, const lpt::ParticleImage::Ptr p2);

void convertCameraCoordinatesToWorld(const lpt::Camera& cam, const array<double, 3>& pc, array<double, 3>& pw);

class Regression {
public:
	double eval(double *coeff, int degree, double x);
    double mult(double *A, double *B, size_t m);
    double Norm(double *A, size_t m);
    void Householder(double **A, double *b, size_t m, size_t n);
    void Backsub(double **A, double *b, size_t n, double *X);
	void Polyfit(double *t, double *y, int m, int degree, double *X);
};
//
//class FastRegression {
//public:
//	double eval(double *coeff, double x);
//	double mult(double *A, double *B);
//	double Norm(double *A);
//	void Householder(double *A, double *b);
//	void Backsub(double *A, double *b, double *X);
//	void Polyfit(double *t, double *y, double *X);
//};
//
//class MergeRegression {
//public:
//	double eval(double *coeff, double x);
//	double mult(double *A, double *B);
//	double Norm(double *A);
//	void Householder(double *A, double *b);
//	void Backsub(double *A, double *b, double *X);
//	void Polyfit(double *t, double *y, double *X);
//};

//////////////////////////////////QR Factorization////////////////////////////////////////////////////
//template <class Object_T>
//class QRfactorization {
//public:
//	typedef typename Object_T::float_type float_type;
//	float_type eval(float_type *coeff, int degree, float_type x);
//	float_type mult(float_type *A, float_type *B,int m);
//	float_type Norm(float_type *A, int m);
//	void Householder(float_type **A, float_type *b, int m, int n);
//	void Backsub(float_type **A, float_type *b, int n, float_type *X);
//	void Polyfit(float_type *t, float_type *y, int m, int degree, float_type *X);
//};
//
//template<class Object_T>
//typename Object_T::float_type QRfactorization<Object_T>::eval(float_type *coeff, int degree, float_type x) {
//	float_type y = 0;
//	for(int i = 0; i < degree; i++)
//		y += coeff[i] * (pow(x,i));
//	return y;
//}
//
//template<class Object_T>
//typename Object_T::float_type QRfactorization<Object_T>::mult(float_type *A, float_type *B, int m) {
//	int i;
//	float_type C=0;
//
//	for (i=0; i<m; i++)
//		C += A[i]*B[i];
//
//	return C;
//}
//
//template<class Object_T>
//typename Object_T::float_type QRfactorization<Object_T>::Norm(float_type *A, int m) {
//	float_type norm, sum = 0;
//	for(int i = 0; i < m; i++)
//		sum += A[i] * A[i];
//
//	norm = sqrt(sum);
//
//	return norm;
//}
//
//template<class Object_T>
//void QRfactorization<Object_T>::Householder(float_type **A, float_type *b, int m, int n) {
//	//Transforms the Vandermonde matrix into upper diagonal form
//	float_type* e, *vk, *aj, *ak;
//	float_type alphak, Bk, yj, norm, sign;
//
//	e = new float_type[m];
//	vk = new float_type[m];
//	aj = new float_type[m];
//	ak = new float_type[m];
//
//	for(int i = 0; i < m; i++) {
//		e[i] = 0;
//		vk[i] = 0;
//	}
//
//	for (int k = 0; k < n; k++) {
//		e[k] = 1.0;
//
//		if(A[k][k] > 0)
//			sign = 1;
//		else
//			sign = -1;
//
//		for(int i = k; i < m; i++) {
//			vk[i] = A[i][k];
//			ak[i] = A[i][k];
//		}
//		norm = Norm(ak, m);
//		alphak = -1.0 * sign * norm;
//		e[k] = alphak;
//		vk[k] -= e[k];
//		Bk = mult(vk, vk, m);
//		if (Bk != 0) {
//			for(int j = k; j < n; j++) {
//				for(int i = 0; i < m; i++)
//					aj[i] = A[i][j];
//
//				yj = mult(vk, aj, m);
//
//				for(int i = 0; i < m; i++)
//					A[i][j] -= (2.0 * yj / Bk) * vk[i];
//			}
//			yj = mult(vk, b, m);
//			for(int i = 0; i < m; i++)
//				b[i] -= (2 * yj / Bk) * vk[i];
//		}
//		ak[k] = 0;
//		e[k] = 0;
//		vk[k] = 0;
//	}
//	delete[] e;
//	delete[] vk;
//	delete[] aj;
//	delete[] ak;
//}
//
//template<class Object_T>
//void QRfactorization<Object_T>::Backsub(float_type** A, float_type* b, int n, float_type* X) {
//	//solves the upper triangular matrix A using back substitution
//	for (int j = n-1; j >= 0; j--) {
//		if(A[j][j] == 0) {
//			cout << "QRfactorization::Backsub-- matrix is singular" << endl;
//			break;
//		}
//		X[j] = b[j] / A[j][j];
//		for (int i = 0; i < j; i++)
//			b[i] = b[i] - A[i][j] * X[j];
//	}
//
//}
//
//template<class Object_T>
//void QRfactorization<Object_T>::Polyfit(float_type* t, float_type* y, int m, int degree, float_type* X) {
//	//float_type* coeff=new float_type[degree];               //degree = order of the polynomial plus 1;
//	float_type** A = new float_type*[m];
//	float_type* b = new float_type[m];
//
//	//----forming the Vandermonde matrix-----  (NOTE: This only works upto degree = 3  !!!)
//	for (int i = 0; i < m; i++) {
//		b[i] = y[i];
//		A[i] = new float_type[degree];
//		for (int j = 0; j < degree; j++) {
//			if(j == 0)
//				A[i][j] = 1;
//			if(j == 1)
//				A[i][j] = t[i];
//			if(j == 2)                    // NOTE:  ADD more if statements if Higher order regression is needed
//				A[i][j] = t[i] * t[i];
//		}
//	}
//	Householder(A, b, m, degree);                  //Transform A into upper triangular form
//	Backsub(A, b, degree, X);                    //Solve for regression coefficients X using back substitution
//	for(int i = 0; i < m; ++i)
//		delete [] A[i];
//	delete[] A;
//	delete[] b;
//	//X=vector of regression coef sorting from low to high polynomial order
//}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//template<class Object_T>
//class Kinematics {
//public:
//	typedef typename Object_T::float_type float_type;
//	typedef Frame_<Object_T> Frame;
//	typedef Trajectory_<Object_T> Trajectory;
//
//	virtual ~Kinematics(){}
//	virtual void calculate(const Trajectory& trajectory, const int frame_loc, vector<vector<float_type> >& kinetics)=0;
//};
//
//template<class Object_T>
//class RegressionModel : public Kinematics<Object_T> {
//private:
//	typedef Object_T Object;
//	typedef typename Object_T::float_type float_type;
//	typedef Frame_<Object> Frame;
//	typedef Trajectory_<Object> Trajectory;
//	static const int DIM = Object_T::dim;
//
//	const int DEGREE;
//	const int TRAJSIZE;
//	QRfactorization<Object_T> regress;
//
//public:
//	RegressionModel(const int degree, const int size) : DEGREE(degree), TRAJSIZE(size) {}
//	inline void calculate(const Trajectory& trajectory, const int frame_loc, vector<vector<float_type> >& kinetics);
//};
//
//template <class Object_T>
//void RegressionModel<Object_T>::calculate(const Trajectory& trajectory, const int frame_loc, vector<vector<float_type> >& kinetics) {
//	kinetics.resize(DIM);						//TODO: Change kinetics from a vector<vector<>> to something faster float_type** or array
//	for(int dim = 0; dim < DIM; ++dim){
//		kinetics[dim].resize(DEGREE + 1, 0);
//	}
//
//	int trajsize = TRAJSIZE; 
//	int degree = DEGREE + 1;
//
//	if( trajsize > trajectory.objects.size() )
//		trajsize = trajectory.objects.size();
//	if(trajsize < degree)
//		degree = trajsize;
//
//	float_type* X = new float_type[trajsize];
//	float_type* t = new float_type[trajsize];
//	float_type* coeff = new float_type[degree];
//	//cout << "objects size = " << trajectory.objects.size() << ",\t trajsize = " << trajsize << endl;
//	for (int dim = 0; dim < DIM; ++dim) {
//		for (int p = trajectory.objects.size() - trajsize; p < trajectory.objects.size(); ++p) {
//			t[p] = (float_type)trajectory.objects[p]->frame_index;
//			X[p] = trajectory.objects[p]->X[dim];
//		}
//		regress.Polyfit(t, X, trajsize, degree, coeff);
//		for (int deg = 0; deg < degree; ++deg)
//			kinetics[dim][deg] = coeff[deg];
//	}
//
//	delete [] X;
//	delete [] t;
//	delete [] coeff;
//	//cout << "calculating velocity and acceleration using Regression: size = " << trajsize << "degree = " << degree << " Dim = " << DIM << endl;
//}
//
//template<class Object_T>
//class FiniteDiffModel : public Kinematics<Object_T> {
//private:
//	typedef Object_T Object;
//	typedef std::shared_ptr<Object_T> Object_Ptr; 
//	typedef typename Object_T::float_type float_type;
//
//	typedef Frame_<Object> Frame;
//	typedef std::shared_ptr<Frame> Frame_Ptr; 
//	
//	typedef Trajectory_<Object> Trajectory;
//	typedef std::shared_ptr<Trajectory> Trajectory_Ptr; 
//	
//	const int DEGREE;
//	const int VEL_ORDER;
//	const int ACCEL_ORDER;
//	static const int DIM = Object_T::dim;
//
//public:
//	FiniteDiffModel(
//			const int degree,
//			const int vel_order = 1,
//			const int accel_order = 2) :
//				DEGREE(degree),
//				VEL_ORDER(vel_order),
//				ACCEL_ORDER(accel_order) {}
//	inline void calculate(const Trajectory& trajectory, const int frame_loc, vector<vector< float_type> >& kinetics);
//};
//
//template <class Object_T>
//void FiniteDiffModel<Object_T>::calculate(const Trajectory& trajectory, const int frame_loc, vector<vector<float_type> >& coefficients) {
//
//	//double dt = 1/FRAMESPERSEC;
//
//	const vector<Object_Ptr>& P = trajectory.objects;
//	int i = P.size() - 1;						// Particle frame index along trajectory
//	coefficients.resize(DIM);						//TODO: Change kinetics from a vector<vector<>> to something faster float_type** or array
//	for(int dim = 0; dim < DIM; ++dim){
//		coefficients[dim].resize(DEGREE + 1, 0);
//		coefficients[dim][0] = P[i]->X[dim];
//	}
//
//	int vel_order = VEL_ORDER;
//	if ( i < VEL_ORDER )
//		vel_order = i;
//
//	switch(vel_order){
//	case 1:
//	{
//		for (int dim =0; dim < DIM; ++dim)
//			coefficients[dim][1] = (P[i]->X[dim] - P[i-1]->X[dim]);
//		break;
//	}
//	case 2:
//	{
//		for (int dim =0; dim < DIM; ++dim)
//			coefficients[dim][1] = 1.0/2.0 * (P[i-2]->X[dim] - 4.0 * P[i-1]->X[dim] + 3.0 * P[i]->X[dim]);
//		break;
//	}
//	case 3:
//	{
//		for (int dim =0; dim < DIM; ++dim)
//			coefficients[dim][1] = (-1.0/(6.0)) * (2.0 * P[i-3]->X[dim] - 9.0 * P[i-2]->X[dim] + 18.0 * P[i-1]->X[dim] - 11.0*P[i]->X[dim]);
//		break;
//	}
//	default:
//		for (int dim =0; dim < DIM; ++dim)
//			coefficients[dim][1] = (P[i]->X[dim] - P[i-1]->X[dim]);
//		break;
//	}
//	//cout << "calculating velocity using FiniteDiff, vel order = " << VEL_ORDER << endl;
//
//	if (DEGREE > 1 && trajectory.objects.size() > 2 ) {
//		switch(ACCEL_ORDER){
//		case 2:
//		{
//			for (int dim =0; dim < DIM; ++dim)
//				coefficients[dim][2] = -1 * (P[i]->X[dim] - 2*P[i-1]->X[dim] + P[i-2]->X[dim]);//(dt*dt);
//			break;
//		}
//		default:
//			break;
//		}
//		//cout << "calculating acceleration using FiniteDiff" << endl;
//	}
//}


} /* NAMESPACE_PT */
#endif /* COMPUTETOOLS_H_ */
