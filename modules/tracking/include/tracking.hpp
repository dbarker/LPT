/**
 * @file tracking.hpp
 * The tracking module declaration
 * Includes CostCalculator class and Tracker class
 */
#ifndef TRACKING_H
#define TRACKING_H

#include "core.hpp"

namespace lpt {

using namespace std;

/**
 * @brief The CostCalculator class
 * Base class for various cost calculators
 */
class CostCalculator {
public:
	typedef shared_ptr<CostCalculator> Ptr;
    virtual double calculate (const lpt::Trajectory3d& traj, const lpt::Particle3d& particle) const = 0;
    virtual ~CostCalculator();
};

/**
 * @brief The CostNearestNeighbor class
 * Derived class from the CostCalculator class
 * Calculates candidate cost based on a 2-frame nearest neighbor criteria
 * Minimizes the displacement within the 2-frame interval
 */
class CostNearestNeighbor : public CostCalculator {
public:
	typedef shared_ptr<CostNearestNeighbor> Ptr;
	static inline CostNearestNeighbor::Ptr create() { return std::make_shared<lpt::CostNearestNeighbor>(); }
	
    CostNearestNeighbor();
    virtual ~CostNearestNeighbor();

    /**
     * @brief Calculates the candidate cost
     * @param traj The trajectory
     * @param particle The candidate particle
     * @return The cost of adding current candidate to this trajectory
     */
    double calculate(const lpt::Trajectory3d& traj, const lpt::Particle3d& particle) const;
};

/**
 * @brief The CostMinimumAcceleration class
 * Derived class from the CostCalculator class
 * Calculates candidate cost based on a 3-frame minimum accelearation algorithm
 * Minimizes the acceleration within the 3-frame interval
 */
class CostMinimumAcceleration : public CostCalculator {
public:
    typedef shared_ptr<CostMinimumAcceleration> Ptr;
    static inline CostMinimumAcceleration::Ptr create () { return std::make_shared<lpt::CostMinimumAcceleration>(); }

    CostMinimumAcceleration();
    ~CostMinimumAcceleration();

    /**
     * @brief Calculates the candidate cost
     * @param traj The trajectory
     * @param particle The candidate particle
     * @return The cost of adding current candidate to this trajectory
     */
    double calculate(const Trajectory3d &traj, const Particle3d &particle) const;
};

/**
 * @brief The Tracker class
 * Manages the tracking parameters and algorithms
 */
class Tracker {
public:
	typedef	shared_ptr<Tracker> Ptr;
	static Tracker::Ptr create() { return std::make_shared<lpt::Tracker>(); }

	friend void callbackSetAlpha(int state, void* data);
	friend void callbackSetRMax(int state, void* data);
	friend void callbackSetRMin(int state, void* data);
    friend inline bool compareCost(pair<int,double>& one, pair<int,double>& two) {return one.second > two.second;}
	
	class Parameters {
	public:
        Parameters(double max_radius = 20, double min_radius = 2, double alpha = 1.0, double sigma_a = 1E-5, double sigma_z = 1E-1)
            : max_radius(max_radius), min_radius(min_radius), alpha(alpha) , KF_sigma_a(sigma_a), KF_sigma_z(sigma_z) {
			max_radius_level = static_cast<int>(max_radius);
			min_radius_level = static_cast<int>(min_radius);
			alpha_level = static_cast<int>(alpha) * 100;
		}
		
		double max_radius;
		int max_radius_level;

		double min_radius;
		int min_radius_level;
		
		double alpha;
		int alpha_level;

        double KF_sigma_a;
        double KF_sigma_z;
	} params;

    /**
     * @brief Tracker constructor
     */
    Tracker();

    /**
     * @brief Tracker destructor
     */
    ~Tracker();

    /**
     * @brief Track particles across consecutive frames
     * @param frame1 The current frame [n]
     * @param frame2 The future frame [n+1]
     */
	void trackFrames(lpt::Frame3d& frame1, lpt::Frame3d& frame2 );
	
    /**
     * @brief addControls
     */
	void addControls();
	
    /**
     * @brief setTrajectoryViews
     * @param cameras
     * @param image_type
     */
	void setTrajectoryViews(vector<lpt::Camera>& cameras, cv::Mat image_type);
	
    /**
     * @brief Set SharedObjects
     * @param new_shared_objects
     */
	inline void setSharedObjects( std::shared_ptr<SharedObjects> new_shared_objects) { shared_objects = new_shared_objects; }

    /**
     * @brief Set CostCalculator
     * @param new_cost_calculator
     */
    inline void setCostCalculator( lpt::CostCalculator::Ptr new_cost_calculator ) { cost_calculator = new_cost_calculator; }
	
    /**
     * @brief Get SharedObjects
     * @return SharedObjects
     */
    inline std::shared_ptr<SharedObjects> getSharedObjects() const { return shared_objects; }
	
    /**
     * @brief Get ActiveTrajectories
     * @return The list of active trajectories
     */
    inline list<lpt::Trajectory3d_Ptr>& getActiveTrajectories() { return active_trajs; }
	
    /**
     * @brief drawTrajectories
     * @param cameras
     */
	void drawTrajectories(vector<lpt::Camera>& cameras);

    void saveTrajectory(ofstream& output, lpt::Trajectory3d_Ptr traj, int id);

	bool clear_drawings;

private:
    lpt::CostCalculator::Ptr					cost_calculator;
	list<lpt::Trajectory3d_Ptr>					trajectories;
	list<lpt::Trajectory3d_Ptr>					active_trajs;
	list<lpt::Particle3d_Ptr>					unmatched_particles;
	vector<cv::Mat>								trajectory_views;
    std::shared_ptr<lpt::SharedObjects>         shared_objects;
    bool                                        isSave;
	size_t										particle_count;
}; 

/**
 * @brief The TestTracker class
 * Tests tracking algorithms and prints results
 */
class TestTracker{
public:
	vector<Trajectory3d_Ptr>& trajectories;
	vector<Trajectory3d_Ptr>& gold_trajectories;

	double correct_ratio;
	double cover_ratio;
	size_t links_created;
	size_t links_total;
	size_t links_correct;

    TestTracker(vector<Trajectory3d_Ptr>& trajs, vector<Trajectory3d_Ptr>& gold_trajs);

    void testTrajectories(size_t maxframes);
    void printTestResults();
};

} /* NAMESPACE_PT */

#endif /* TRACKING_H */
