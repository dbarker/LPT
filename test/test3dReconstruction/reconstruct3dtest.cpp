
#include "core.hpp"
#include "correspond.hpp"
#include "datagen.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace std;
using namespace boost::accumulators;

typedef accumulator_set< double, features<tag::mean, tag::variance, tag::count, tag::max, tag::min > > boost_accumulator;

int main(int argc, char** argv){
	string trajectories_file;
	string data_path;

	if (argc > 1) {
		trajectories_file = argv[1];
	}
	else{
		trajectories_file = "../../../data/tests/pivstd/300p145fpiv352";
	}

	if (argc > 2) {
		data_path = argv[2];
	}
	else{
		data_path = "../../../data/tests/pivstd/6cams/";
	}

	auto shared_objects = std::make_shared<lpt::SharedObjects>();
	auto& cameras = shared_objects->cameras;
	auto& camera_pairs = shared_objects->camera_pairs;
    
	string cameras_file = data_path + "cameras.yaml";
   lpt::readCamerasFile(cameras_file, cameras);
	lpt::generateCameraPairs(cameras, camera_pairs);
	
	auto image_creator = std::make_shared<lpt::ImageCreator>();
	image_creator->radius = 1;
	image_creator->intensity = 255;
	//image_creator->image_type = cv::Mat::zeros( cv::Size(1280, 1024), CV_8UC1 );

	lpt::DataSetGenerator generate;
	generate.setSharedObjects(shared_objects);
	cout << "camera pairs size: " << camera_pairs.size() << endl;
	generate.setDataPath(data_path);
	generate.setImageCreator(image_creator);

	generate.read3DTrajectoryFile(trajectories_file, lpt::PLAINTEXT);
	vector< lpt::Frame3d > gold_frames( generate.frames.size() );
	vector< vector<lpt::Match::Ptr> > globalmatches( gold_frames.size() );
	
	for (int f = 0; f < generate.frames.size(); ++f) {
		auto& points = generate.frames[f].second;
		for (int p = 0; p < points.size(); ++p) {
			auto newparticle = std::make_shared<lpt::Particle3d>();
			newparticle->X[0] = points[p].x;
			newparticle->X[1] = points[p].y;
			newparticle->X[2] = points[p].z;

			gold_frames[f].objects.push_back(newparticle);		
	 	}
    }

	cout << "Projecting 3D points to 2D image planes" << endl;
	generate.project3DFramesTo2D();
	generate.showImages();
    
    for (int f = 0; f < cameras[0].frames.size(); ++f) {
		for (int p = 0; p < cameras[0].frames[f].particles.size(); ++p) {
			auto newmatch = std::make_shared<lpt::Match>();
    		for (int c = 0; c < cameras.size(); ++c) {
    			newmatch->addParticle(cameras[c].frames[f].particles[p].get(), cameras[c].id);
    		}
    		globalmatches[f].push_back(newmatch);
    	}
    }

   // double t1, t2;
   // t1 = lpt::get_clock();

   vector<lpt::Frame3d> frames(gold_frames.size());
	lpt::Reconstruct3D recon;
	recon.setSharedObjects(shared_objects);

	for (int f = 0; f < frames.size(); ++f) {
		recon.reconstruct3DFrame(globalmatches[f], frames[f]);
	}

	boost_accumulator error_stats;
	cout << "Calculating reprojection error on " << frames.size() << " frames " << endl;
	for (int f = 0; f < frames.size(); ++f) {
		auto& gold_points3D = gold_frames[f].objects;
		auto& points3D = frames[f].objects;
		for (int p = 0; p < gold_points3D.size(); ++p) {

			//3D reprojection error
			double error_3d = 
				sqrt(
				( gold_points3D[p]->X[0] - points3D[p]->X[0]) * (gold_points3D[p]->X[0] - points3D[p]->X[0]) +
				( gold_points3D[p]->X[1] - points3D[p]->X[1]) * (gold_points3D[p]->X[1] - points3D[p]->X[1]) +
				( gold_points3D[p]->X[2] - points3D[p]->X[2]) * (gold_points3D[p]->X[2] - points3D[p]->X[2])
				);
			error_stats(error_3d);
		}
	}

	double mean = extract_result<tag::mean>(error_stats);
	double stdev = extract_result<tag::variance>(error_stats);
	stdev = (stdev > 0 ? sqrt(stdev) : 0 ); 
	double min = extract_result<tag::min>(error_stats);
	double max = extract_result<tag::max>(error_stats);
	cout << "3D pos error = " << mean << " +- " << stdev << ", Max = " << max << ", Min = " << min << endl;
	cout << "Finished" << endl;	
    return 0;
}

