
#ifndef CORRESPOND_H_
#define CORRESPOND_H_

#include <core.hpp>
#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>

#define NUM_MATCHES 10

namespace lpt {

using namespace std;

typedef vector< vector < array <int, NUM_MATCHES> > > MatchMap;

class Correspondence {
public:	
	typedef	std::shared_ptr<Correspondence> Ptr;
	
	Correspondence() : 
		initial_max_particles_per_image(2048), 
		map_storage_size(500) { }
	
	void findMatches(lpt::ImageFrameGroup& frame_group, vector<lpt::Match::Ptr>& matches) {
		resetMatchMap(current_matchmap);
		findEpipolarMatches(frame_group, current_matchmap);
		findUniqueMatches(frame_group, current_matchmap, matches);
	}

	void  run(
		lpt::concurrent_queue<lpt::ImageFrameGroup>* in_queue, 
		lpt::concurrent_queue<std::pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > >* out_queue,
		int number_epipolarmatching_threads = 1,
		int number_uniquematching_threads = 1
		); 

	void stop();
	void runEpipolarMatching(int thread_id, lpt::concurrent_queue<lpt::ImageFrameGroup>* in_queue);
	void runUniqueMatching(int thread_id, lpt::concurrent_queue<std::pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > >* out_queue);	

	virtual void findEpipolarMatches(lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap)=0;
	virtual void findUniqueMatches(lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches)=0;

	virtual void initialize()=0;
	virtual void initializeEpipolarMatchThread(int thread_id)=0;
	virtual void addControls()=0;
		
	inline double calculateEpipolarResidual(double lineB[3], lpt::ParticleImage& pointB);
	inline void calculateEpiline( lpt::ParticleImage& pointA, const double F[3][3], double line[3] );
	
	inline void setSharedObjects( std::shared_ptr<SharedObjects> new_shared_objects) { shared_objects = new_shared_objects; }
	inline std::shared_ptr<SharedObjects> getSharedObjects() { return shared_objects; }
	void testMatches(lpt::ImageFrameGroup& cameragroup, vector<lpt::Match::Ptr>& matches);
	void printMatchMap(lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, string output_file_name);

protected:
	void resetMatchMap(lpt::MatchMap& matchmap);
	void initializeMatchMap();
	lpt::MatchMap init_matchmap;
	lpt::MatchMap current_matchmap;

	std::shared_ptr < lpt::SharedObjects > shared_objects; 
	
	vector<pair<lpt::ImageFrameGroup, lpt::MatchMap>> matchmap_storage;

	lpt::concurrent_queue< vector<pair<lpt::ImageFrameGroup, lpt::MatchMap>>::iterator > empty_maps_queue;
	lpt::concurrent_queue< vector<pair<lpt::ImageFrameGroup, lpt::MatchMap>>::iterator > full_maps_queue;
	lpt::concurrent_queue< pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > > final_match_queue;

	boost::thread_group epipolar_match_threads;
	boost::thread_group unique_match_threads;

	int initial_max_particles_per_image;
	int map_storage_size;
};

class PointMatcher : public Correspondence {
public:
	typedef	std::shared_ptr<lpt::PointMatcher> Ptr;
	static inline lpt::PointMatcher::Ptr create() { return std::make_shared<lpt::PointMatcher>(); }

	class Parameters {
	public:
		Parameters() : match_thresh_level(5), match_threshold(0.5){}
		double match_threshold;
		int match_thresh_level;
	} params;

	PointMatcher() { cout << "Epipolor Point matcher created " << endl; }
		
	virtual void initialize(){
		this->initializeMatchMap();
	}
	virtual void initializeEpipolarMatchThread(int thread_id){}

	virtual void addControls() {
		void* matcher_void_ptr = static_cast<void*> ( this );
		cv::createTrackbar("Match Thresh", string() , &params.match_thresh_level, 50, callbackMatchThresh, matcher_void_ptr);
	}
	virtual void findEpipolarMatches(lpt::ImageFrameGroup& cameragroup, lpt::MatchMap& matchmap);
	virtual void findUniqueMatches(lpt::ImageFrameGroup& frame_group, lpt::MatchMap& MatchMap, vector<lpt::Match::Ptr>& matches);

	friend void callbackMatchThresh(int state, void* data) {
		PointMatcher* matcher = static_cast<PointMatcher*>(data);
		matcher->params.match_threshold = matcher->params.match_thresh_level / 10.0;
	}

private:
	inline void removeNonUniqueMatches(vector<Match::Ptr>& matches);
	//void find3wayMatches(lpt::ParticleImage::Ptr Pa, int camB_id, int camC_id );
	//virtual void find3wayMatches(tuple<int,int,int>& triple, lpt::ImageFrameGroup& cameragroup, vector<lpt::Match::Ptr>& matches);
	//virtual void find4WayMatches(tuple<int,int,int,int>& quad, lpt::ImageFrameGroup& cameragroup, vector<lpt::Match::Ptr>& matches);
};

class Reconstruct3D {
public:
	typedef std::shared_ptr<Reconstruct3D> Ptr;
	static inline Reconstruct3D::Ptr create(vector<lpt::Camera>& cams) { return std::make_shared<Reconstruct3D>(); }

	Reconstruct3D() { }
	virtual void reconstruct3DFrame(vector<lpt::Match::Ptr>& matches, lpt::Frame3d& frame);
	virtual void draw(lpt::Frame3d& frame);

	inline void setSharedObjects( std::shared_ptr<SharedObjects> new_shared_objects) { shared_objects = new_shared_objects; }
	inline std::shared_ptr<SharedObjects> getSharedObjects() { return shared_objects; }
protected:
	std::shared_ptr < lpt::SharedObjects > shared_objects;
private:
	lpt::Regression solver;	
};

class Recontstuct3DwithSVD : public Reconstruct3D {
public:
	typedef std::shared_ptr<Recontstuct3DwithSVD> Ptr;
	static inline Recontstuct3DwithSVD::Ptr create() { return std::make_shared<lpt::Recontstuct3DwithSVD>(); }

	Recontstuct3DwithSVD() : Reconstruct3D() {}
	virtual void reconstruct3DFrame(vector<lpt::Match::Ptr>& matches, lpt::Frame3d& frame);
private:
	cv::SVD svd;
};

} /* NAMESPACE_PT */

#endif /* CORRESPOND_H_ */
