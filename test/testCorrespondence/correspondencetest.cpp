#include <core.hpp>
#include <correspond.hpp>
#include <boost/chrono.hpp>

#ifdef USE_CUDA
	#include <correspondcuda.h>
#endif

using namespace std;

int main(int argc, char** argv){

	string project_path;
	
	if (argc > 1) {
        project_path = argv[1];
    }
    else{
		project_path = "../../../data/tests/pivstd/6cams/";
    }

	cout << "using project path: " << project_path << endl;
    
    size_t number_of_frames = 0;
    auto shared_objects = lpt::SharedObjects::create();
	auto& cameras = shared_objects->cameras;
	auto& camera_pairs = shared_objects->camera_pairs;

	vector<lpt::ImageFrameGroup> framegroups;

    string full_cameras_filename = project_path + "cameras.yaml";
    lpt::readCamerasFile(full_cameras_filename, cameras);
    
    for (size_t c = 0; c < cameras.size(); ++c){
    	cout << cameras[c];
    	stringstream file_name;
    	file_name << cameras[c].id << "_pts.yaml";
    	string full_frame_filename = project_path + file_name.str();
    	lpt::readImageFramesFile(full_frame_filename, cameras[c].frames);	
        cout << "Number of frames = " << cameras[c].frames.size() << endl;
    }

	number_of_frames = cameras[0].frames.size();
    
	string pairs_filename = project_path + "camera_pairs.yaml";
	lpt::readCameraPairsFile(pairs_filename, cameras, camera_pairs);
	framegroups.resize(number_of_frames);
    for (size_t f = 0; f < framegroups.size(); ++f) {
        for (size_t c = 0; c < cameras.size(); ++c) {
			std::random_shuffle(cameras[c].frames[f].particles.begin(), cameras[c].frames[f].particles.end()); 
			framegroups[f].push_back(std::move(cameras[c].frames[f]));
		}
	}

    for (size_t i = 0; i < camera_pairs.size(); ++i)
    	cout << camera_pairs[i];
	
	double match_threshold = 0.05;
	
	lpt::PointMatcher::Ptr host_matcher = lpt::PointMatcher::create();
	host_matcher->setSharedObjects(shared_objects);
	host_matcher->params.match_threshold = match_threshold;
	
	cout << "Match threshold = " << host_matcher->params.match_threshold << endl;

	lpt::Correspondence::Ptr host_matcher_ptr = host_matcher;

	host_matcher_ptr->initialize();
	cout << "Host correspondence solver initialization complete" << endl;
	
	boost::chrono::duration<double> hosttime;
	boost::chrono::system_clock::time_point start, stop;
	
	vector<lpt::Match::Ptr> hostmatches;
	
	start = boost::chrono::system_clock::now();	
    for (size_t f = 0; f < framegroups.size(); ++f) {
		host_matcher_ptr->findMatches(framegroups[f], hostmatches);
        cout << "Frame " << f << ":";
		host_matcher_ptr->testMatches(framegroups[f], hostmatches);
		//cout << hostmatches.size() << endl;
	}
	stop = boost::chrono::system_clock::now();
	hosttime = stop - start;
	
	//host_matcher_ptr->testMatches(framegroups[144], hostmatches);
    host_matcher_ptr->printMatchMap(framegroups[144], "matchmap.txt");
	cout << "Host correspondence solution complete: " << endl;	
	cout << "Host time = " << hosttime.count() << " seconds, " << framegroups.size() / hosttime.count() << " fps" << endl;

    for (size_t f = 0; f < framegroups.size(); ++f) {
        for (size_t c = 0; c < framegroups[f].size(); ++c) {
            for (size_t p = 0; p < framegroups[f][c].particles.size(); ++p) {
				framegroups[f][c].particles[p]->is_4way_matched = false;
				framegroups[f][c].particles[p]->match_count = 0;
			}
		}
	}

#ifdef USE_CUDA
	lpt::PointMatcherCUDA::Ptr cuda_matcher = lpt::PointMatcherCUDA::create();
	cuda_matcher->setSharedObjects(shared_objects);
	cuda_matcher->params.match_threshold = static_cast<float>(match_threshold); 

	lpt::Correspondence::Ptr cuda_matcher_ptr = cuda_matcher;

	cuda_matcher_ptr->initialize();
	cuda_matcher_ptr->initializeEpipolarMatchThread(0);

	vector<lpt::Match::Ptr> cudamatches;

	boost::chrono::duration<double> cudatime;

	start = boost::chrono::system_clock::now();		
    for (size_t f = 0; f < framegroups.size(); ++f) {
		cuda_matcher_ptr->findMatches(framegroups[f], cudamatches);
	}

	stop = boost::chrono::system_clock::now();
	cudatime = stop - start;
	
	cuda_matcher_ptr->testMatches(framegroups[144], cudamatches);
	
	cout << endl << "GPU (CUDA) correspondence solution complete" << endl;
	cout << "CUDA time = " << cudatime.count() << " seconds, " << framegroups.size() / cudatime.count() << " fps" << endl;
	cout << "speed up = " << hosttime.count() / cudatime.count() << "x " << endl;
#endif
	cout << endl << "Finished" << endl;
}
