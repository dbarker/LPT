/**
 * @file correspond.hpp
 * @brief Correspondence module declaration
 * The module contains correspondence and 3D reconstruction classes
 */

#ifndef CORRESPOND_H_
#define CORRESPOND_H_

#include <core.hpp>
#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>

#define NUM_MATCHES 10

namespace lpt
{

using namespace std;


//typedef accumulator_set< double, features<tag::mean, tag::variance, tag::count, tag::max, tag::min> > boost_accumulator;

/**
 * @brief MatchMap
 * Bookkeeping tool to identify 4-way unique matches from epipolar matches
 * Size: # total particles by # cameras by # max_matches
 */
typedef vector< vector < array <int, NUM_MATCHES> > > MatchMap;

/**
 * @brief The Correspondence class
 * Base class for correspondence functionality
 * Completes correspondence in two steps: 1) Epipolar matching 2) 4-way unique matching
 * Also defines the streaming pipeline functions
 */
class Correspondence
{
public:	
	typedef	std::shared_ptr<Correspondence> Ptr;
	
    Correspondence();

    virtual ~Correspondence();
	
    /**
     * @brief Completes correspondence in an image frame group
     * Invokes resetMatchMap, findEpipolarMathes, and findUniqueMatches
     *
     * @param frame_group Input image frame group
     * @param matches Output vector of Match pointers
     */
    void findMatches(const ImageFrameGroup &frame_group, vector<lpt::Match::Ptr>& matches);

    /**
     * @brief Run correspondence
     *
     * @param in_queue Input queue
     * @param out_queue Output queue
     * @param number_epipolarmatching_threads
     * @param number_uniquematching_threads
     */
	void  run(
		lpt::concurrent_queue<lpt::ImageFrameGroup>* in_queue, 
		lpt::concurrent_queue<std::pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > >* out_queue,
		int number_epipolarmatching_threads = 1,
		int number_uniquematching_threads = 1
		); 

    /**
     * @brief Stop correspondence
     */
	void stop();

    /**
     * @brief Run Epipolar matching
     * @param thread_id
     * @param in_queue Input queue
     */
	void runEpipolarMatching(int thread_id, lpt::concurrent_queue<lpt::ImageFrameGroup>* in_queue);

    /**
     * @brief Run 4-way unique matching
     *
     * @param thread_id
     * @param out_queue Output queue
     */
	void runUniqueMatching(int thread_id, lpt::concurrent_queue<std::pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > >* out_queue);	

    /**
     * @brief Find Epipolar matches
     * Epipolar matches are created by calculating Epiploar lines and Epipolar residuals
     * MatchMap is generated for use in finding unique matches
     *
     * @param frame_group
     * @param matchmap
     */
	virtual void findEpipolarMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap)=0;

    /**
     * @brief Find 4-way unique matches
     * A match must be found in at least 4 cameras to be considered a unique match
     * Unique matches are added to the final match list
     *
     * @param frame_group
     * @param matchmap
     * @param matches The final match list
     */
	virtual void findUniqueMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches)=0;

	/**
	 * @brief Find 3-way unique matches
	 *
	 * @param frame_group
     * @param matchmap
     * @param matches The final match list
     */
	virtual void find3WayMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches)=0;

    /**
     * @brief Initialize MatchMap
     */
	virtual void initialize()=0;
	virtual void initializeEpipolarMatchThread(int thread_id)=0;

    /**
     * @brief addControls
     */
	virtual void addControls()=0;

    /**
     * @brief Calculate the residual between an Epipolar line and a point
     * @param line The Epipolar line
     * @param point The Given point
     * @return Epiolar residual
     */
    inline double calculateEpipolarResidual(double line[3], const lpt::ParticleImage& point) const;

    /**
     * @brief Calculate the Epipolar line
     * @param point Starting point
     * @param F Fundamental matrix
     * @param line The resulting Epipolar line
     */
    inline void calculateEpiline( const lpt::ParticleImage& point, const double F[3][3], double line[3] ) const;
	
    /**
     * @brief Set SharedObjects
     * @param new_shared_objects
     */
	inline void setSharedObjects( std::shared_ptr<SharedObjects> new_shared_objects) { shared_objects = new_shared_objects; }
    /**
     * @brief Get SharedObjects
     * @return SharedObjects
     */
    inline std::shared_ptr<SharedObjects> getSharedObjects() const { return shared_objects; }
    /**
     * @brief Test matches
     * @param cameragroup
     * @param matches
     */
    void testMatches(const lpt::ImageFrameGroup& cameragroup, const vector<lpt::Match::Ptr>& matches) const;
    /**
     * @brief Print MatchMap
     * @param frame_group
     * @param matchmap
     * @param output_file_name
     */
    void printMatchMap(const lpt::ImageFrameGroup& frame_group, string output_file_name) const;

	void printMatchMap(const lpt::ImageFrameGroup& frame_group, const lpt::MatchMap& match_map, string output_file_name) const;

protected:
    /**
     * @brief Reset MatchMap
     * @param matchmap
     */
	void resetMatchMap(lpt::MatchMap& matchmap);
    /**
     * @brief Initialize MatchMap
     */
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

/**
 * @brief The PointMatcher class
 * Derived class from the Correspondence class
 */
class PointMatcher : public Correspondence
{
public:
	typedef	std::shared_ptr<lpt::PointMatcher> Ptr;
	static inline lpt::PointMatcher::Ptr create() { return std::make_shared<lpt::PointMatcher>(); }

	class Parameters {
	public:
        Parameters() : match_threshold(0.5), match_thresh_level(5) {}
		float match_threshold;
		int match_thresh_level;
	} params;

    PointMatcher();
    virtual ~PointMatcher();
		
    virtual void initialize();
	virtual void initializeEpipolarMatchThread(int thread_id){}
    virtual void addControls();

	virtual void findEpipolarMatches(const lpt::ImageFrameGroup& cameragroup, lpt::MatchMap& matchmap);
	virtual void findUniqueMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& MatchMap, vector<lpt::Match::Ptr>& matches);
	virtual void find3WayMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& MatchMap, vector<lpt::Match::Ptr>& matches);

	friend void callbackMatchThresh(int state, void* data) {
		PointMatcher* matcher = static_cast<PointMatcher*>(data);
		matcher->params.match_threshold = static_cast<float>(matcher->params.match_thresh_level) / 10.0;
	}

private:
    /**
     * @brief Remove NonUnique Matches
     * @param matches
     */
	inline void removeNonUniqueMatches(vector<Match::Ptr>& matches);
	//void find3wayMatches(lpt::ParticleImage::Ptr Pa, int camB_id, int camC_id );
	//virtual void find3wayMatches(tuple<int,int,int>& triple, lpt::ImageFrameGroup& cameragroup, vector<lpt::Match::Ptr>& matches);
	//virtual void find4WayMatches(tuple<int,int,int,int>& quad, lpt::ImageFrameGroup& cameragroup, vector<lpt::Match::Ptr>& matches);
};

class Reconstruct3D
{
public:
	typedef std::shared_ptr<Reconstruct3D> Ptr;
    static inline Reconstruct3D::Ptr create() { return std::make_shared<Reconstruct3D>(); }

    Reconstruct3D();
    virtual ~Reconstruct3D();

    /**
     * @brief Reconstruct 3D frame for matched particles
     * @param matches Match list
     * @param frame Output 3D frame
     */
	virtual void reconstruct3DFrame(vector<lpt::Match::Ptr>& matches, lpt::Frame3d& frame);

    /**
     * @brief Draw reconstruced objects
     * @param frame
     */
	virtual void draw(lpt::Frame3d& frame);

	inline void setSharedObjects( std::shared_ptr<SharedObjects> new_shared_objects) { shared_objects = new_shared_objects; axis.open(shared_objects->output_path + "axis.txt"); }
    inline std::shared_ptr<SharedObjects> getSharedObjects() const { return shared_objects; }
protected:
	std::shared_ptr < lpt::SharedObjects > shared_objects;
private:
	lpt::Regression solver;	
	ofstream axis;
	vector<array<double, 3>> positions;
	size_t frame_count;
};

class Reconstruct3DwithSVD : public Reconstruct3D
{
public:
    typedef std::shared_ptr<Reconstruct3DwithSVD> Ptr;
    static inline Reconstruct3DwithSVD::Ptr create() { return std::make_shared<lpt::Reconstruct3DwithSVD>(); }

    Reconstruct3DwithSVD();
    virtual ~Reconstruct3DwithSVD();

	virtual void reconstruct3DFrame(vector<lpt::Match::Ptr>& matches, lpt::Frame3d& frame);
private:
	cv::SVD svd;
};

} /* NAMESPACE_PT */

#endif /* CORRESPOND_H_ */
