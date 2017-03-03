#include <core.hpp>
#include <imageproc.hpp>
#include <datagen.hpp>

using namespace std;

int main(int argc, char** argv)
{
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

    auto shared_objects = lpt::SharedObjects::create();

	auto& cameras = shared_objects->cameras;
    auto& camera_pairs = shared_objects->camera_pairs;

    auto processor = lpt::ImageProcessor::create();
    processor->addProcess( lpt::GaussianBlur::create(5) );
    processor->addProcess( lpt::Threshold::create(20) );

    auto fc_detector = lpt::FindContoursDetector::create();
    fc_detector->params.max_contour_area = 150;

    //auto gftt_detector = std::make_shared<lpt::GoodFeaturesToTrackDetector>();

    lpt::Detector::Ptr detector = fc_detector;  //or gftt_detector;

    string cameras_file = cameras_path + "cameras.yaml";
    cout << cameras_file << endl;
    lpt::readCamerasFile(cameras_file, cameras);

    string camera_pairs_file = cameras_path + "camera_pairs.yaml";
    cout << camera_pairs_file << endl;
    lpt::readCameraPairsFile(camera_pairs_file, cameras, camera_pairs);

    auto image_creator = std::make_shared<lpt::ImageCreator>();
	image_creator->radius = 0;
	image_creator->intensity = 0;
	image_creator->object_intensity = static_cast<int>(5E8);
	image_creator->object_size = 3;
	image_creator->blur_ksize = 3;
	//creator.image_type = cv::Mat::zeros( cv::Size(1280, 1024), CV_8UC1 );

	lpt::DataSetGenerator generate;

	generate.setSharedObjects(shared_objects);
    cout << "camera size: " << cameras.size() << endl;
    cout << "camera pairs size: " << camera_pairs.size() << endl;
	//generate.setDataPath(output_path);
	generate.setImageCreator(image_creator);
	generate.read3DTrajectoryFile(trajectories_file, lpt::PLAINTEXT);
	cout << "Projecting" << endl;
	generate.project3DFramesTo2D();
	generate.showImages();

    cout << endl << "Running Test" << endl;
    for (size_t c = 0; c < cameras.size(); ++c){
        stringstream file_name;
        file_name << cameras[c].id << "_pts.yaml";
        vector<lpt::ImageFrame> true_frames;
        string imageframe_points_file = cameras_path + file_name.str();
        lpt::readImageFramesFile(imageframe_points_file, true_frames);
        cout << true_frames[0].particles.size() << endl;
      
        for (size_t f = 0; f < true_frames.size(); ++f) {
            cout << "Frame " << f << ":";
            processor->processImage(cameras[c].frames[f].image);
            detector->detectFeatures(cameras[c].frames[f].image, cameras[c].frames[f].particles, cameras[c].frames[f].contours);
            lpt::testDetectedParticles(true_frames[f].particles, cameras[c].frames[f].particles);
        }
   }

   cout <<"Finished" <<endl;
   return 0;
}
