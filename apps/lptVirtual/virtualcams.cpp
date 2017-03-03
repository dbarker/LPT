#include <core.hpp>
#include <imageproc.hpp>
#include <datagen.hpp>
#include <dataaquisition.hpp>
#include <tracking.hpp>

using namespace std;

int main(int argc, char** argv) {
	
	string input = (argc > 1 ? argv[1] : "../../../data/input/");
	string output = (argc > 2 ? argv[2] : "../../../data/output/");
	string camerasfile = "../../../data/tests/pivstd/6cams/cameras.yaml"; //"./output/cameras.yaml"; //"../../../data/pivstd/4cams/cameras.yaml"; 
	string campairsfile = "../../../data/tests/pivstd/6cams/camera_pairs.yaml"; //"./output/camera_pairs.yaml"; //
	string traj_file = "../../../data/tests/pivstd/300p145fpiv352"; //"../../../data/grid.txt"; 
	
	lpt::StreamingPipeline pipeline;
	pipeline.setQueueCapacity(10);

	lpt::ImageProcessor::Ptr processor = lpt::ImageProcessor::create();
    lpt::ImageProcess::Ptr blur = lpt::GaussianBlur::create(3);
    lpt::ImageProcess::Ptr thresh = lpt::Threshold::create(20);
    processor->addProcess( blur );
	processor->addProcess( thresh );
	
	lpt::FindContoursDetector::Ptr detector = lpt::FindContoursDetector::create();
	
	auto image_creator = std::make_shared<lpt::ImageCreator>();
	image_creator->radius = 0;
	image_creator->intensity = 0;
	image_creator->object_intensity = 5E8;
	image_creator->object_size = 3;
	image_creator->blur_ksize = 3;
	auto camera_system = lpt::VirtualCameras::create(camerasfile, traj_file);
	camera_system->getGenerator()->setImageCreator(image_creator);
	
	lpt::PointMatcher::Ptr matcher = lpt::PointMatcher::create();
	matcher->params.match_threshold = 2.0;
	matcher->params.match_thresh_level = 20;
	
	lpt::PointMatcherCUDA::Ptr matcher_cuda =  lpt::PointMatcherCUDA::create();
	matcher_cuda->params.match_threshold = 2.0; //pixels
	matcher_cuda->params.match_thresh_level = 20;

	lpt::Tracker::Ptr tracker = lpt::Tracker::create();
    tracker->setCostCalculator(lpt::CostMinimumAcceleration::create());
	tracker->params.min_radius = 4.0; //mm
	tracker->params.min_radius_level = 4;
	tracker->params.max_radius = 25.0; //mm
	tracker->params.max_radius_level = 25;
	tracker->params.KF_sigma_a = 1E-5;
	tracker->params.KF_sigma_z = 1E-1;

	bool KalmanFilter = false;

	lpt::Visualizer::Ptr visualizer = lpt::Visualizer::create();
	visualizer->getVolumeGrid()->setGridOrigin(0,0,0);
	visualizer->getVolumeGrid()->setGridDimensions(100,100,100); // mm
	visualizer->getVolumeGrid()->setGridCellCounts(10,10,10);
	visualizer->params.queue_capacity = 10;
	pipeline.setInputDataPath(input);
	pipeline.setOutputDataPath(output);
	pipeline.setKalmanFilter(KalmanFilter);
	pipeline.attachCameraSystem(camera_system);
	pipeline.attachImageProcessor(processor);
	pipeline.attachDetector(detector);
	pipeline.attachMatcher(matcher_cuda);
	pipeline.attachTracker(tracker);
	pipeline.attachVisualizer(visualizer);

	bool on = pipeline.initialize();

	if (on){ 
		camera_system->loadCameraParams(camerasfile);
		camera_system->loadCameraPairParams(campairsfile);
		pipeline.run();
	} else 
		cout << "System could not initialize: Shutting down" << endl;
	
	cout << "Finished Stream Data" << endl;
	return 0;
}

