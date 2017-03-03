#ifndef CORE_H_
#define CORE_H_

#define EVALRESULT 0  //FIXME: MOVE THESE AND MAKE EXTERN!!!
#define CHARM 0       //FIXME: MOVE THESE AND MAKE EXTERN!!!

#include <cstdio>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include <cmath>
#include <map>
#include <vector>
#include <list>
#include <array>
#include <queue>
#include <set>
#include <tuple>
#include <new>
#include <utility>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <memory>
#include <limits>

#include <opencv2/opencv.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/bind.hpp>

#if(CHARM)
#include <charm++.h>
#endif /* CHARM */

//********namespace lpt: all LPT functionality**************
namespace lpt {

using namespace std;

struct Dimensions {	enum {XY = 2, XYZ = 3}; };
enum CameraSystemType { OPTITRACK, VIRTUAL };

class Particle;
class Frame;
class Trajectory;
class ParticleImage;
class Match;
class ImageFrame;
class Camera;
class CameraPair;
class KalmanFilter;

typedef vector<ImageFrame> ImageFrameGroup;

template <typename FLOATTYPE, int DIM> class Particle_;
template <class Object_T> class Frame_;
template <class Object_T> class Trajectory_;

typedef Particle_<double, Dimensions::XYZ> Particle3d;
typedef std::shared_ptr<Particle3d> Particle3d_Ptr;

typedef Frame_<Particle3d> Frame3d;
typedef std::shared_ptr<Frame3d> Frame3d_Ptr;

typedef Trajectory_<Particle3d> Trajectory3d;
typedef std::shared_ptr<Trajectory3d> Trajectory3d_Ptr;

typedef Particle_<double, Dimensions::XY> Particle2d;
typedef Frame_<Particle2d> Frame2d;
typedef Trajectory_<Particle2d> Trajectory2d;

typedef Particle_<float, Dimensions::XYZ> Particle3f;
typedef Frame_<Particle3f> Frame3f;
typedef Trajectory_<Particle3f> Trajectory3f;

typedef Particle_<float, Dimensions::XY> Particle2f;
typedef Frame_<Particle2f> Frame2f;
typedef Trajectory_<Particle2f> Trajectory2f;

template <class Object_T> class Kinematics;

template<typename Data>
class concurrent_queue
{
private:
	std::queue<Data> m_queue;
	boost::mutex m_mutex;
	boost::condition_variable m_condition_variable;
	int m_capacity;

public:
	concurrent_queue(int capacity = 1000) : m_capacity(capacity) {}

	bool push(Data const& data)
	{
		boost::mutex::scoped_lock lock(m_mutex);
		bool ispushed = false;

		if (m_queue.size() < m_capacity)
		{
			m_queue.push(data);
			ispushed = true;
		}

		lock.unlock();
        m_condition_variable.notify_one();
		return ispushed;
	}

	inline size_t size() const	{ return m_queue.size();	}

	bool empty()
	{
		boost::mutex::scoped_lock lock(m_mutex);
		return m_queue.empty();
	}

	inline void clear() {
		boost::mutex::scoped_lock lock(m_mutex);
		m_queue = std::queue<Data>();
	}
	
   inline bool full()
	{
		boost::mutex::scoped_lock lock(m_mutex);
		return (m_queue.size() < m_capacity ? false : true);
	}

	inline void setCapacity(int capacity) {
		m_capacity = capacity;
	}

	bool try_pop(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(m_mutex);
		if( m_queue.empty() )
			return false;

		popped_value = m_queue.front();
		m_queue.pop();
		return true;
	}

	void wait_and_pop(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(m_mutex);
		while( m_queue.empty() ) 
			m_condition_variable.wait(lock);

		popped_value = m_queue.front();
		m_queue.pop();
	}

};

class ParticleImage {
public:
	typedef std::shared_ptr<ParticleImage> Ptr;
	static inline ParticleImage::Ptr create(
		int idnumber = -1, 
		double x_pixel = 0.0, 
		double y_pixel = 0.0, 
		double radius = 0.0,
		double intensity = 0.0
		)
        { return std::make_shared<ParticleImage>(idnumber, x_pixel, y_pixel, radius, intensity); }
	
	double x, y, radius, intensity;
	int id, match_count; 
	bool is_4way_matched;
	
	ParticleImage():id(-1), x(0), y(0), radius(0), intensity(0), is_4way_matched(false), match_count(0){}
	ParticleImage(int idnumber, double x_pixel, double y_pixel, double radius = 0, double intensity = 0):
		id(idnumber),x(x_pixel), y(y_pixel), radius(radius), intensity(intensity), is_4way_matched(false), match_count(0){}
};

class Match {
public:
	typedef std::shared_ptr<Match> Ptr;
    static inline Match::Ptr create() { return std::make_shared<Match>(); }

    vector < pair < ParticleImage::Ptr, size_t > > particles;
	double residual;
	//array<ParticleImage*,4> p;
	//array<int, 4> i;
	//int id;
	Match():residual(0)/*, id(0)*/{}
	
    inline void addParticle(lpt::ParticleImage::Ptr new_particle, size_t camera_id) {
        new_particle->match_count++;
        particles.push_back( std::move(std::make_pair(new_particle, camera_id)) );
		//p[id] = added_particle;
		//i[id] = camera_id;
	}

	inline bool isUnique(){
		bool unique = true;
		for(int i = 0; i < this->particles.size(); ++i) {
			unique &= particles[i].first->match_count == 1;
			particles[i].first->match_count = 0;
		}
		/*for(int i = 0; i < p.size(); ++i) {
			unique &= p[i]->match_count == 1;
			p[i]->match_count = 0;
		}
*/
		return unique;
	}
};

class ImageFrame {
public:
	int frame_index;
	cv::Mat image;
	vector<ParticleImage::Ptr> particles;
	vector<vector<cv::Point>> contours;

	ImageFrame():frame_index(-1){}
	ImageFrame(int f):frame_index(f){}
	ImageFrame(int f, cv::Mat image):frame_index(f), image(image){}
	inline void addParticle(ParticleImage::Ptr new_particle) { particles.push_back(new_particle);}
};

inline vector< cv::Mat > getImageVector (ImageFrameGroup& framegroup) {
	vector< cv::Mat > images(framegroup.size()); 
	for (int i = 0; i < framegroup.size(); ++i)
		images[i] = framegroup[i].image.clone();
	return images;
}

enum ObjectType{ THREEPOINTLINE, SINGLEPOINT };

class Object {
public:
	virtual bool find(lpt::ImageFrame& frame)=0;
	virtual void draw(vector<lpt::ParticleImage::Ptr>& particles, cv::Mat& image)=0;
	lpt::ObjectType object_type;
};

class SinglePoint : public Object {
	bool find(lpt::ImageFrame& frame);
	void draw(vector<lpt::ParticleImage::Ptr>& particles, cv::Mat& image);
};

class ThreePointLine : public Object {
public:
	bool find(lpt::ImageFrame& frame);
	void draw(vector<lpt::ParticleImage::Ptr>& particles, cv::Mat& image);
};

class ObjectFrame {
public:
	typedef std::shared_ptr<ObjectFrame> Ptr;
    static inline ObjectFrame::Ptr create( lpt::ImageFrameGroup& frames, int index = -1 )
    {
        return std::make_shared<ObjectFrame>(frames, index);
    }

    lpt::ImageFrameGroup camera_frames;
	int frame_index;
	vector<lpt::Particle3d_Ptr> particles;

    ObjectFrame() {}
    ObjectFrame(lpt::ImageFrameGroup& frames, int index = -1) : camera_frames(frames), frame_index(index) {}
};

class Camera {
public:
	int id;
	string name;
	string path;

	double R [3][3];        //Rotation matrix
	double r_vec_u[3];      //Rotation vector uncertainty

	double T [3];           //Translation matrix
	double T_u [3];         //Translation matrix uncertainty

	double f [2];           //focal length fx and fy
	double f_u [2];         //focal length uncertainty 

	double c [2];           //center of image cx and cy
	double c_u [2];         //center of image uncertainty in cx and cy

	double X_u[2];			//x and y image coordinate uncertainty

	double dist_coeffs [4]; //distortion coeffs [k1, k2, p1, p2] FIXME: make number of distortion coeffs adjustable
	double dist_coeffs_u [4]; //distortion coeffs uncertainty [k1, k2, p1, p2] 

	double centriod_loc_uncertainty; //pixels
	double sensor_size[2];      //mm [width, height]
	double pixel_size[2]; 		//mm [x, y]
	vector<lpt::ImageFrame> frames;  
	vector<string> imagelist;
	double avg_reprojection_error;
	//vector<float> reprojection_errors;
	
	Camera(int id_number = -1): id(id_number), avg_reprojection_error(0.0){
		
		c[0] = 1.0;
		c[1] = 1.0;
		f[0] = 1.0;
		f[1] = 1.0;

		dist_coeffs[0] = 0.0; 
		dist_coeffs[1] = 0.0;
		dist_coeffs[2] = 0.0;
		dist_coeffs[3] = 0.0;
		
		dist_coeffs_u[0] = 0.0; 
		dist_coeffs_u[1] = 0.0;
		dist_coeffs_u[2] = 0.0;
		dist_coeffs_u[3] = 0.0;
		
		X_u[0] = 0;
		X_u[1] = 0;
		
		f_u[0] = 0;
		f_u[1] = 0;
		
		c_u[0] = 0;
		c_u[1] = 0;

		for (int i = 0; i < 3; ++i ) {
			r_vec_u[i] = 0;
			T_u[i] = 0;
		}
		
		centriod_loc_uncertainty = 0;
	}

    inline cv::Mat getCameraMatrix() const {
		cv::Mat camera_matrix = cv::Mat::eye(3,3, CV_64F);
		camera_matrix.at<double>(0,0) = this->f[0];
		camera_matrix.at<double>(1,1) = this->f[1];
		camera_matrix.at<double>(0,2) = this->c[0];
		camera_matrix.at<double>(1,2) = this->c[1];
		return camera_matrix;
	}

    inline cv::Mat getDistCoeffs() const {
		cv::Mat dist = cv::Mat(4, 1, CV_64F);
		dist.at<double>(0) = this->dist_coeffs[0]; 
		dist.at<double>(1) = this->dist_coeffs[1];
		dist.at<double>(2) = this->dist_coeffs[2];
		dist.at<double>(3) = this->dist_coeffs[3];
		return dist;
	}
	
	void readImageList();
	void readFrames(int startframe = 0); 
	void writeFrames();
};

class CameraPair {
public:
    lpt::Camera& cam_A;
    lpt::Camera& cam_B;
    double F [3][3];
    double epipolar_error;
    double reprojection_error;

    CameraPair(lpt::Camera& A, lpt::Camera& B ) : cam_A(A), cam_B(B),
        epipolar_error(0.0),reprojection_error(0.0){}
};

class SharedObjects { 
public:
    typedef std::shared_ptr<SharedObjects> Ptr;
    static inline SharedObjects::Ptr create() {return std::make_shared<SharedObjects>();}

    SharedObjects() : frame_rate(0), input_path("./"), output_path("./"), KF_isOn(false), isRotation_Correction(false) { }
		
    vector<lpt::Camera> cameras;
    vector<lpt::CameraPair> camera_pairs;
    lpt::CameraSystemType camera_type;
    cv::Mat image_type;
	
    int frame_rate;
    string input_path;
    string output_path;
    bool KF_isOn;
	array<array<double, 3>, 3> S;
	array<double, 3> P;
	bool isRotation_Correction;
};

class Particle {
public:
	typedef std::shared_ptr<Particle> Ptr;
    static inline Particle::Ptr create(
        int id = -1, double x = 0.0, double y = 0.0, double z = 0.0)
    {
        return std::make_shared<Particle>(id, x, y, z);
    }

	double x, y, z, e;
	int id, frame_index;

	Particle(){}
	Particle(double x, double y, double z) : x(x), y(y), z(z) {}
	Particle(int id, double x, double y, double z) : id(id), x(x), y(y), z(z) {}
	
#if(EVALRESULT)
	double vel, accel, cost, dr;
#endif /* EVALRESULT */

#if(CHARM)
	void pup(PUP::er &p) {
		p | x;
		p | y;
		p | z;
		p | frame_index;
		p | id;
	}
#endif/* CHARM */
};

class Frame {
public:
	typedef std::shared_ptr<Frame> Ptr;
    static inline Frame::Ptr create(int index = -1) { return std::make_shared<Frame>(index); }

	int frame_index;
	vector<lpt::Particle::Ptr> particles;
	lpt::ImageFrameGroup camera_frames;

	Frame(){}
	Frame(int index): frame_index(index){}

#if(EVALRESULT)
	double S;
#endif /* EVALRESULT */

#if(CHARM)
	CkVec<Particle> _particles;
	void pup(PUP::er &p) {
		if (p.isSizing()) {
			_particles.clear();

			for (int i = 0; i < particles.size(); i++)
				_particles.push_back(*particles[i]);

		}

		p | frame_index;
		p | _particles;

		if (p.isUnpacking()) {
			particles.clear();

			for (int i = 0; i < _particles.size(); i++)
				particles.push_back(&_particles[i]);

		}
	}
#endif /* CHARM */
};

class Trajectory {
public:
	typedef std::shared_ptr<Trajectory> Ptr;
    static inline Trajectory::Ptr create() { return std::make_shared<Trajectory>(); }

	int id, gap, startframe;
	vector<Particle::Ptr> particles;
	vector<bool> matches;

#if(EVALRESULT)
	double R,Rh,Rl,C,Ch,Cl;
#endif /* EVALRESULT */

#if(CHARM)
	CkVec<Particle> _particles;

	Trajectory *trim(bool end) {
		Trajectory *t = new Trajectory();
		t->id = this->id;
		t->gap = this->gap;
		t->startframe = this->startframe;

		// assuming that particles it at least MERGE_COST_SIZE/2 of size
		if (end) {
			for (int i = particles.size() - MERGE_COST_SIZE/2; i < particles.size(); i++)
				t->particles.push_back(this->particles[i]);
		} else {
			for (int i = 0; i < MERGE_COST_SIZE/2; i++)
				t->particles.push_back(this->particles[i]);
		}

		return t;
	}

	// TODO: JL push this into the copy constructor
	// TODO:  is this fixParts() needed? or is this taken care of in the pup routine below?

	void fixParts() {
		particles.clear();
		for (int i = 0; i < _particles.size(); i++) {
			particles.push_back(&_particles[i]);
		}
	}

	void pup(PUP::er &p) {
		if (p.isSizing()) {namespace
			_particles.clear();

			for (int i = 0; i < particles.size(); i++)
				_particles.push_back(*particles[i]);
		}

		p | id;
		p | gap;
		p | _particles;
		p | startframe;

		if (p.isUnpacking()) {
			particles.clear();

			for (int i = 0; i < _particles.size(); i++)
				particles.push_back(&_particles[i]);
		}
	}
#endif /* CHARM */
};


//**************************************TEMPLATES************************

template <typename FLOATTYPE, int DIM>
class Particle_ {
public:
	typedef FLOATTYPE float_type;
	array<float_type, DIM> X;
    array<size_t, 4> camera_ids;
	int id, frame_index;
	static const int dim = DIM;
	
	Particle_() : id(-1), frame_index(-1) {
		for (int i = 0; i < DIM; ++i)
			X[i] = 0.0;
	}

    Particle_(array<float_type, DIM> coords, int id) : X(coords), id(id), frame_index(-1) {}

	static inline std::shared_ptr<Particle_<FLOATTYPE, DIM> > create() {  
        return std::make_shared<Particle_<FLOATTYPE, DIM>>();
	}
	static inline std::shared_ptr<Particle_<FLOATTYPE, DIM> > create(array<float_type, DIM> coords, int id = -1) { 
        return std::make_shared<Particle_<FLOATTYPE, DIM>>(coords, id);
	}
};

template <class Object_T>
class Frame_{
public:
	typedef Object_T Object;
	vector<std::shared_ptr<Object_T>> objects;	//Fixme Can we eliminate pointers here?
	lpt::ImageFrameGroup camera_frames;
	int frame_index;

	Frame_(){}
    Frame_(lpt::ImageFrameGroup& frames, int index) : camera_frames(frames), frame_index(index) {}

	static inline std::shared_ptr<Frame_<Object_T>> create() { 
        return std::make_shared<Frame_<Object_T>>();
	}
	static inline std::shared_ptr<Frame_<Object_T>> create( lpt::ImageFrameGroup& frames, int index = -1 ) { 
        return std::make_shared<Frame_<Object_T>>(frames, index);
	}
};

class TrajectoryVTKBase {
public:
	typedef std::shared_ptr<TrajectoryVTKBase> Ptr;
	TrajectoryVTKBase () {}
	virtual ~TrajectoryVTKBase() {}
};

class KalmanFilter {
public:
    typedef std::shared_ptr<KalmanFilter> Ptr;
    static inline KalmanFilter::Ptr create(cv::Mat stateIni)
    { return std::make_shared<KalmanFilter>(stateIni); }

    KalmanFilter(cv::Mat stateIni)
      : sigma_a(1E-5), sigma_z(0.1), dt(1.0/120)
    {
		KF = cv::KalmanFilter(9, 9, 0, CV_64F);
		KF.statePre = stateIni;
		cv::setIdentity(KF.measurementMatrix);
		cv::setIdentity(KF.errorCovPost);
		cv::setIdentity(KF.processNoiseCov, cv::Scalar::all(sigma_a));
		cv::setIdentity(KF.measurementNoiseCov, cv::Scalar::all(sigma_z));
    }

    ~KalmanFilter() {}

	cv::KalmanFilter KF;
private:

    double sigma_a;                     // accleration standard deviation
    double sigma_z;                     // measurement uncertanty
    double dt;                          // time step size (inverse of frame rate)

    //std::shared_ptr<SharedObjects> shared_objects;
};

template <class Object_T>
class Trajectory_{
public:
	typedef Object_T Object;
	typedef typename Object_T::float_type float_type;
	
	static inline std::shared_ptr<Trajectory_<Object_T>> create() { 
        return std::make_shared<Trajectory_<Object_T>>();
    }

    static int count;

    void initializeKF(cv::Mat stateIni, double sigma_a, double sigma_z) {
		kalman_filter = cv::KalmanFilter(9, 9, 0, CV_64F);
		kalman_filter.statePost = stateIni;
		cv::setIdentity(kalman_filter.transitionMatrix);
		for (int i=0; i<9; i++) {
			if (i < 6)
				kalman_filter.transitionMatrix.at<float_type>(i,i+3) = 1.0/120;
			if (i < 3)
				kalman_filter.transitionMatrix.at<float_type>(i,i+6) = 0.5/120/120;
		}
		cv::setIdentity(kalman_filter.measurementMatrix);
		cv::setIdentity(kalman_filter.errorCovPost);
		cv::setIdentity(kalman_filter.processNoiseCov, cv::Scalar::all(sigma_a));
		cv::setIdentity(kalman_filter.measurementNoiseCov, cv::Scalar::all(sigma_z));
	}

    vector<std::pair<std::shared_ptr<Object_T>, array<double, 9> > > objects;
    //vector<array<double, 9>> obj_state;
	
	array<float_type, Object_T::dim> V;
	array<float_type, Object_T::dim> A;

    vector<pair<int, float_type> > candidates_costs;
	int id, gap, startframe, last_cell_id;
	float_type search_radius;
	Object_T extrap_object;
	lpt::TrajectoryVTKBase::Ptr trajvtk_ptr;
    cv::KalmanFilter kalman_filter;

    Trajectory_() : search_radius(0.0), id(-1), gap(0), startframe(0), last_cell_id(-1) {}
    ~Trajectory_() {}
};

template<class Object_T>
typename Object_T::float_type calculateDistance(const Object_T& P1, const Object_T& P2) {
	typename Object_T::float_type distance = 0;
	for (int i = 0; i < Object_T::dim; ++i)
		distance += (P1.X[i] - P2.X[i]) * (P1.X[i] - P2.X[i]);
	return sqrt(distance);
}

template<class Object_T>
inline void extrapolatePosition(lpt::Trajectory_<Object_T>& traj)
{
	typedef typename Object_T::float_type float_type;

	Object_T& t2 = *(traj.objects.rbegin()->first);
	Object_T& t1 = *((++traj.objects.rbegin())->first);
	for (int d = 0; d < Object_T::dim; ++d) {
		traj.V[d] = t2.X[d] - t1.X[d];
		traj.extrap_object.X[d] = t2.X[d] + ( traj.V[d] );
		traj.V[d] *= 120; //FIXME: This may not be in use ....This should be varialble based on frame rate and number of frames between points -- (t2.frame_index -  t1.frame_index) );
		traj.A[d] = 0;
	}

	traj.search_radius = lpt::calculateDistance(t2, t1);  
}

} /* NAMESPACE_PT */

#include "utilities.hpp"
#include "computetools.hpp"

#endif /* CORE_H_ */

