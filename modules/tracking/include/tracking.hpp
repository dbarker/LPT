#ifndef TRACKING_H
#define TRACKING_H

#include "core.hpp"

namespace lpt {

using namespace std;

class CostCalculator {
public:
	typedef shared_ptr<CostCalculator> Ptr;
	virtual double calculate (lpt::Trajectory3d& traj, lpt::Particle3d& particle)=0;
};

class CostNearestNeighbor : public CostCalculator {
public:
	typedef shared_ptr<CostNearestNeighbor> Ptr;
	static inline CostNearestNeighbor::Ptr create() { return std::make_shared<lpt::CostNearestNeighbor>(); }
	
	CostNearestNeighbor() {}
	double calculate(lpt::Trajectory3d& traj, lpt::Particle3d& particle) {
        return lpt::calculateDistance(*(traj.objects.back()), particle);
	}
};

class CostMinimumAcceleration : public CostCalculator {
public:
    typedef shared_ptr<CostMinimumAcceleration> Ptr;
    static inline CostMinimumAcceleration::Ptr create () { return std::make_shared<lpt::CostMinimumAcceleration>(); }

    CostMinimumAcceleration() {}
    double calculate(lpt::Trajectory3d& traj, lpt::Particle3d& particle) {
        return lpt::calculateDistance(traj.extrap_object, particle);
    }
};

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
		Parameters(double max_radius = 20, double min_radius = 2, double alpha = 1.0) : max_radius(max_radius), min_radius(min_radius), alpha(alpha) {
			max_radius_level = max_radius;
			min_radius_level = min_radius;
			alpha_level = alpha * 100;
		}
		
		double max_radius;
		int max_radius_level;

		double min_radius;
		int min_radius_level;
		
		double alpha;
		int alpha_level;
	} params;

	Tracker() : clear_drawings(false) { cout << "tracker created " << endl; }
	void trackFrames(lpt::Frame3d& frame1, lpt::Frame3d& frame2 );
	
	void addControls();
	
	void setTrajectoryViews(vector<lpt::Camera>& cameras, cv::Mat image_type);
	
	inline void setSharedObjects( std::shared_ptr<SharedObjects> new_shared_objects) { shared_objects = new_shared_objects; }

    inline void setCostCalculator( lpt::CostCalculator::Ptr new_cost_calculator ) { cost_calculator = new_cost_calculator; }
	
	inline std::shared_ptr<SharedObjects> getSharedObjects() { return shared_objects; }
	
	inline list<lpt::Trajectory3d_Ptr>& getActiveTrajectories() { return active_trajs; }
	
	void drawTrajectories(vector<lpt::Camera>& cameras);

	bool clear_drawings;

private:
    lpt::CostCalculator::Ptr					cost_calculator;
	list<lpt::Trajectory3d_Ptr>					trajectories;
	list<lpt::Trajectory3d_Ptr>					active_trajs;
	list<lpt::Particle3d_Ptr>					unmatched_particles;
	vector<cv::Mat>								trajectory_views;
    std::shared_ptr<lpt::SharedObjects>         shared_objects;
}; 


class TestTracker{
public:
	vector<Trajectory3d_Ptr>& trajectories;
	vector<Trajectory3d_Ptr>& gold_trajectories;

	double correct_ratio;
	double cover_ratio;
	int links_created;
	int links_total;
	int links_correct;
	TestTracker(
			vector<Trajectory3d_Ptr>& trajs,
			vector<Trajectory3d_Ptr>& gold_trajs
			):
				trajectories(trajs),
				gold_trajectories(gold_trajs),
				correct_ratio(0.0),
				cover_ratio(0.0),
				links_created(0),
				links_total(0),
				links_correct(0) {}
	void testTrajectories(int maxframes);
	void printTestResults();
};

} /* NAMESPACE_PT */

#endif /* TRACKING_H */
