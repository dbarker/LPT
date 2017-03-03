#include <core.hpp>
#include <tracking.hpp>
#include <boost/chrono.hpp>

using namespace std;

int main(int argc, char** argv) {

	string input_file;
	string output_file;

	if (argc > 1) {
		input_file = argv[1];
	}
	else{
		input_file = "../../../data/tests/synthetic.txt";
	}

	if (argc > 2) {
		output_file = argv[2];
	}
	else{
		output_file = "../../../data/output/";
	}

	cout << input_file << endl;

	boost::chrono::duration<double> tracking_time;
	boost::chrono::system_clock::time_point start, stop;

	//----Input frames from known trajectory data----
	lpt::Input in;

	vector<lpt::Trajectory3d_Ptr> gold_trajectories; 
	vector<lpt::Frame3d_Ptr> frames;
	
	in.readTrajectoryFile(input_file, gold_trajectories);
	
	cout << "trajectories size " << gold_trajectories.size() << endl;

	in.convertTrajectoriesToFrames(gold_trajectories, frames);

	cout << "frames size " << frames[0]->objects.size() << endl;
		
	cout << "--Tracking Started:" << endl;
	lpt::Tracker tracker;
	auto shared_objects = lpt::SharedObjects::create();
	shared_objects->frame_rate = 60;
	shared_objects->KF_isOn = true;
	shared_objects->output_path = output_file;
	tracker.setSharedObjects(shared_objects);
	tracker.setCostCalculator(lpt::CostMinimumAcceleration::create());
	tracker.params.KF_sigma_a = 1E-5;
	tracker.params.KF_sigma_z = 1E-2;

	tracker.params.alpha = 1.0;  	  //0.8 for pivSTD352       //0.8 for CFD 1000p4096f100fps         //1.0 for scaling data sets
	tracker.params.max_radius = 50.0; 
	tracker.params.min_radius = 1.0;	 

	//--Initialize Timing function
	start = boost::chrono::system_clock::now();		
	
	for (int f = 0; f < frames.size() - 1; ++f) {
		tracker.trackFrames(*frames[f], *frames[f+1]);
	}
	
	stop = boost::chrono::system_clock::now();		
	
	tracking_time = stop - start;

	auto trajectory_list = tracker.getActiveTrajectories();
	vector<lpt::Trajectory3d_Ptr> trajectories( trajectory_list.begin(), trajectory_list.end() );
	
	cout << "--Tracking Complete: " << endl;
	cout << "\tTime = " << tracking_time.count() << " seconds, " << frames.size() / tracking_time.count() << " fps" << endl;
	cout << "\tTotal number of trajectories tracked = " << trajectories.size() << endl;

	//----Run Final Test----
	lpt::TestTracker evaluate(gold_trajectories, trajectories);
	evaluate.testTrajectories( frames.size() );
	evaluate.printTestResults();
	return 0;
}
