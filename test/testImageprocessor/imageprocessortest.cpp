#include <core.hpp>
#include <imageproc.hpp>
#include <datagen.hpp>

using namespace std;

int main(int argc, char** argv){
   string cameras_path;
   string trajectories_file;

   if (argc > 1)
      cameras_path = argv[1];
   else
      cameras_path = "../../../data/tests/pivstd/6cams/";
   if (argc > 2) {
		trajectories_file = argv[1];
	}
	else{
		trajectories_file = "../../../data/tests/pivstd/300p145fpiv352";
	}

   cout << cameras_path << endl;
   cout << trajectories_file << endl;

   auto shared_objects = std::make_shared<lpt::SharedObjects>(); 

	auto& cameras = shared_objects->cameras;
	auto& camera_pairs = shared_objects->camera_pairs;

   lpt::ImageProcessor processor;
   processor.addProcess( lpt::GaussianBlur(5) );
   processor.addProcess( lpt::Threshold(20) );

   auto fc_detector = std::make_shared<lpt::FindContoursDetector>();
   fc_detector->params.max_contour_area = 150;

   //auto gftt_detector = std::make_shared<lpt::GoodFeaturesToTrackDetector>();

   lpt::Detector::Ptr detector = fc_detector;  //or gftt_detector;

   string cameras_file = cameras_path + "cameras.yaml";
   cout << cameras_file << endl;
   lpt::readCamerasFile(cameras_file, cameras);

   auto image_creator = std::make_shared<lpt::ImageCreator>();
	image_creator->radius = 0;
	image_creator->intensity = 0;
	image_creator->object_intensity = 5E8;
	image_creator->object_size = 3;
	image_creator->blur_ksize = 3;
	//creator.image_type = cv::Mat::zeros( cv::Size(1280, 1024), CV_8UC1 );

	lpt::DataSetGenerator generate;

	generate.setSharedObjects(shared_objects);
	cout << "camera pairs size: " << shared_objects->camera_pairs.size() << endl;
	//generate.setDataPath(output_path);
	generate.setImageCreator(image_creator);
	generate.read3DTrajectoryFile(trajectories_file, lpt::PLAINTEXT);
	cout << "Projecting" << endl;
	generate.project3DFramesTo2D();
	generate.showImages();

   cout << endl << "Running Test" << endl;
   for (int c = 0; c < cameras.size(); ++c){
      stringstream file_name;
      file_name << cameras[c].id << "_pts.yaml";
      vector<lpt::ImageFrame> true_frames;
      string imageframe_points_file = cameras_path + file_name.str();
      lpt::readImageFramesFile(imageframe_points_file, true_frames);
      cout << true_frames[0].particles.size() << endl;
      
      for (int f = 0; f < true_frames.size(); ++f) {
         cout << "Frame " << f << ":";
         lpt::testDetectedParticles(true_frames[f].particles, cameras[c].frames[f].particles);
      }
   }
   cout <<"Finished" <<endl;
   return 0;

}
