/*
Real-time particle tracking system:
This is the main code for the real-time particle tracking system. 
This creates the streaming pipeline of tasks including:
	1) data aquisition from camera
	2) image processing
	3) object detection in each image
	4) camera correspondence solver
	5) 3D reconstruction
	6) Temporal tracking
	7) Visualization

 copyright: Douglas Barker 2011
*/

#include "core.hpp"
#include "imageproc.hpp"
#include "datagen.hpp"
#include "dataaquisition.hpp"
#include "correspond.hpp"
#include "tracking.hpp"
#include "visualization.hpp"
//--All particle tracking functionality is loaded in namespace lpt::

using namespace std;

int main(int argc, char** argv) {
	
	string input = (argc > 1 ? argv[1] : "../../../data/input/");
	string output = (argc > 2 ? argv[2] : "../../../data/output/");
	string cameras_file = input + "cameras.yaml"; //"../../../data/pivstd/4cams/cameras.yaml"; //
	string camera_pairs_file = input + "camera_pairs.yaml"; //"../../../data/pivstd/4cams/camera_pairs.yaml"; //
	
	lpt::StreamingPipeline pipeline;
	pipeline.setQueueCapacity(1000);

	auto processor = std::make_shared<lpt::ImageProcessor>();
	lpt::GaussianBlur blur(3);
	lpt::Threshold thresh(20);
	processor->addProcess( blur );
	//processor->addProcess( thresh );
	
	auto detector = std::make_shared<lpt::FindContoursDetector>();
	
	auto camera_system = std::make_shared<lpt::Optitrack>();

	auto matcher = std::make_shared<lpt::PointMatcher>();

	auto matcher_cuda =  std::make_shared<lpt::PointMatcherCUDA>();
	matcher_cuda->params.match_threshold = 2.0; //pixels
	matcher_cuda->params.match_thresh_level = 20;

	auto tracker = std::make_shared<lpt::Tracker>();
	tracker->params.min_radius = 4.0; //mm
	tracker->params.min_radius_level = 4;
	tracker->params.max_radius = 25.0; //mm
	tracker->params.max_radius_level = 25;

	auto visualizer = std::make_shared<lpt::Visualizer>();
	double grid_side_length = 300; // mm
	int cell_counts[3] = {21, 31, 31};
	visualizer->getVolumeGrid()->setGridOrigin(137,-1.0 * grid_side_length/2.0 + 83,  -1.0*grid_side_length/2.0);
	visualizer->getVolumeGrid()->setGridCellCounts(cell_counts[0], cell_counts[1], cell_counts[2]);
	visualizer->getVolumeGrid()->setGridDimensions( cell_counts[0]*grid_side_length/cell_counts[2], grid_side_length, grid_side_length); // mm
	
	pipeline.setInputDataPath(input);
	pipeline.setOutputDataPath(output);
	pipeline.attachCameraSystem(camera_system);
	pipeline.attachImageProcessor(processor);
	pipeline.attachDetector(detector);
	pipeline.attachMatcher(matcher_cuda);
	pipeline.attachTracker(tracker);
	pipeline.attachVisualizer(visualizer);
	
	bool on = pipeline.initialize();

	if (on){ 
		camera_system->loadCameraParams(cameras_file);
		camera_system->loadCameraPairParams(camera_pairs_file);
		pipeline.run();
	} else 
		cout << "System could not initialize: Shutting down" << endl;
	
	cout << "Finished" << endl;
	return 0;
}

