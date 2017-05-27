#ifndef _CORRESPONDCUDA_H_
#define _CORRESPONDCUDA_H_

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/experimental/cuda/pinned_allocator.h> 
#include <core.hpp>
#include <correspond.hpp>

using namespace std;

namespace lpt {

//void findMatchesHost(lpt::ImageFrameGroup& cameragroup, vector<lpt::Match::Ptr>& matches); 

class MatchIDs {
public:
	int ids[NUM_MATCHES];	
};

class CameraPairCUDA {
public:
	float F[3][3];
	int cam_a_id;
	int cam_b_id;
};

template < typename T >
class KernelArray {
public:    
	KernelArray<T>( thrust::device_vector< T >& d_vec ) {
		data = thrust::raw_pointer_cast( &d_vec[0] );
		size = ( int ) d_vec.size();
	}
	
	T*  data;
    int size;
};

class PointMatcherCUDA : public lpt::Correspondence {
public:
	typedef std::shared_ptr<lpt::PointMatcherCUDA> Ptr;
	static inline lpt::PointMatcherCUDA::Ptr create() { return std::make_shared<lpt::PointMatcherCUDA>(); }

	class Parameters {
	public:
		Parameters() : match_thresh_level(10), match_threshold(0.5){}
		float match_threshold;
		int match_thresh_level;
	} params;

	PointMatcherCUDA();

	virtual void initialize();
	virtual void initializeEpipolarMatchThread(int thread_id);
	virtual void addControls();
    virtual void findUniqueMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches);
    virtual void findEpipolarMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap);
	virtual void find3WayMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches);

	friend void callbackMatchThreshcuda(int state, void* data) {
		PointMatcherCUDA* matcher = static_cast<PointMatcherCUDA*>(data);
		matcher->params.match_threshold = static_cast<float>(matcher->params.match_thresh_level) / 10.0;
	}

private:
	void findEpipolarMatchesStreams(lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap);
	void findEpipolarMatchesManyThreads(lpt::ImageFrameGroup& frame_group); //TODO: try to make this work!!
	int getNextComputeDeviceID();
	thrust::host_vector<CameraPairCUDA> camera_pairs_h;
	thrust::host_vector<int> num_matches_h;
    thrust::device_vector<int> num_matches_d;
	thrust::host_vector< MatchIDs > matches2way_h;
	thrust::device_vector< MatchIDs > matches2way_d;
	
	MatchIDs match_initializer;

	vector<cudaStream_t> streams;
		
	thrust::host_vector<float> particles_x_h;
	thrust::host_vector<float> particles_y_h;

	thrust::device_vector<float> particles_x_d;
	thrust::device_vector<float> particles_y_d;

	thrust::host_vector<int> num_particles_h;
    thrust::device_vector<int> num_particles_d;

	boost::mutex mutex;
	queue<int> compute_devices_available;
};

} // namespace lpt

#endif // _CORRESPONDCUDA_H_
