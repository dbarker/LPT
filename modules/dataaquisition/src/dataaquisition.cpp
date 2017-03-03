/*
Data aquisition module implementation:
These classes provide data aquisition funcionality to the particle tracking system and the overall multi-threaded streaming framework.

copyright: Douglas Barker 2011
*/

#include "dataaquisition.hpp"

namespace lpt {

using namespace std;
using namespace CameraLibrary;

/*****lpt::Video class implementation*****/

void Video::writeAVI(int playback_fps, int codec) {
	cv::Size frame_size = video_frames[0].size();
	bool color = ( video_frames[0].channels() == 3 ) ? true : false;
	cv::VideoWriter writer(file_path, codec, playback_fps, frame_size, color);
	for (int f = 0; f < video_frames.size(); ++f) 
		writer.write(video_frames[f]);
	cout << "Wrote video to file: " << file_path << endl;
}

void Video::writeImages(string basename, string filetype) {
	for (int f = 0; f < video_frames.size(); ++f) {
		stringstream file_name;
		file_name << file_path << basename << "_" << f << filetype;
		cv::imwrite(file_name.str(), video_frames[f]);
	}
}

void Video::readAVI() {
	if ( ! file_path.empty() ) {
		cout << "reading " << file_path << endl;
		cv::VideoCapture capture(file_path);
		cv::Mat frame;
		while ( capture.read(frame) ) {
			cv::cvtColor(frame, frame, CV_RGB2GRAY);
			this->addImageToVideo( frame.clone() );
		}
	}
}

void Video::readImages(string basename, int numframes, string filetype ) {
	video_frames.resize(numframes);
	for (int f = 0; f < video_frames.size(); ++f) {
		stringstream file_name;
		file_name << file_path << basename << "_" << f << filetype;
		video_frames[f] = cv::imread( file_name.str() );
	}
}


/*****lpt::Recorder class implementation*****/

void Recorder::takeSnapShot( vector<cv::Mat>& frames ) {
	snapshot_id++;
	for (int camera_id = 0; camera_id < frames.size(); camera_id++) { 
		stringstream filename;
		filename << path << camera_id << "_" << snapshot_id << ".jpg";
		cv::imwrite(filename.str(), frames[camera_id]);	
		imagelists[camera_id].push_back( std::move( filename.str() ) );
	}
	snapshot_requested = false;
	cout << "Snapshot taken: ID = " << snapshot_id << endl;
}

void Recorder::writeSnapShotImageLists() {
	stringstream mainlistname;
	mainlistname << path << "mainlist.txt";
	ofstream mainlistfile( mainlistname.str().c_str() );
	for (int camera_id = 0; camera_id < imagelists.size(); camera_id++) {
		stringstream listname; 
		listname << path << "imagelist" << camera_id << ".txt";
		mainlistfile << listname.str() << endl;
		ofstream imagelistfile( listname.str().c_str() );
		for (int image_id = 0; image_id < imagelists[camera_id].size(); image_id++) {
			imagelistfile << imagelists[camera_id][image_id] << endl;
		}
		imagelistfile.close();
	}
	mainlistfile.close();
}

void Recorder::createVideos() {
	video_id++;
	for (int cam_id = 0; cam_id < number_of_cameras; cam_id++) {
		stringstream filename;
		filename << path << cam_id << "_" << video_id << ".avi";
		Video new_video( filename.str() );
		videos[cam_id].push_back( std::move(new_video) );
	}
	record_video = true;
	cout << "Recording video clip (" << clip_length << " frames): ID = " << video_id << endl;
}

void Recorder::addFramesToVideos(vector<cv::Mat>& frames) {
	if (videos[0][video_id].getVideoLength() < clip_length) {
		for (int cam_id = 0; cam_id < frames.size(); ++cam_id) 
			videos[cam_id][video_id].addImageToVideo(frames[cam_id]);
	} 
	else {
		record_video = false;
		cout << "Finished recording video clip: ID = " << video_id << endl;
	}
}

void Recorder::writeVideosToFile(int fps, int codec) {
	for (int cam_id = 0; cam_id < number_of_cameras; cam_id++)
		for (int count = 0; count < videos[cam_id].size(); count++)
			videos[cam_id][count].writeAVI(fps, codec);
}

void callbackRecordVideo(int state, void* data) {
	Recorder* recorder = static_cast<Recorder*> (data);
	if ( !recorder->isVideoRecording() )
		recorder->createVideos();
}

void callbackTakeSnapshot(int state, void* data) {
	Recorder* recorder = static_cast<Recorder*> (data);
	recorder->requestSnapShot();
}


/*****lpt::CameraSystem class implementation*****/

void CameraSystem::loadCameraParams(const string filepath) {
	auto& cameras = shared_objects->cameras;
	if ( cameras.empty() ) {
		cout << "pt_cameras is empty: cannot load camera data" << endl;
		return;
	}
	vector<lpt::Camera> cameras_data;
	lpt::readCamerasFile(filepath, cameras_data);
	 
	if ( cameras_data.size() == cameras.size() ) {
		for (int a = 0; a < cameras.size(); ++a) {
			for (int b = 0; b < cameras_data.size(); ++b) {  
				if ( cameras[a].name == cameras_data[b].name ) {
					//The cameras match - Load camera object with matched data from file
					auto frames = cameras[a].frames;
					cameras[a] = cameras_data[b];
					cameras[a].frames = frames;
					break;
				}
			}
		}
		cout << "Loaded camera parameters for " << cameras_data.size() << " cameras: " << endl;
	}
}

void CameraSystem::loadCameraPairParams(const string filepath) {
	auto& cameras = shared_objects->cameras;
	auto& camera_pairs = shared_objects->camera_pairs;
	if ( cameras.empty() ) {
		cout << "cameras is empty: cannot load camera pair data" << endl;
		return;
	}

	lpt::readCameraPairsFile(filepath, cameras, camera_pairs);  //TODO: This will fail if the file contains pairs for a greater number of cameras than pt_cameras has
	cout << "Loaded camera parameters for " << camera_pairs.size() << " camera pairs: " << endl;
}


/*****lpt::VirtualCameras class implementation*****/

bool VirtualCameras::initializeCameras(){
	cout << "Virtual Camera System initializing" << endl;

	auto& cameras = shared_objects->cameras;
	auto& camera_pairs = shared_objects->camera_pairs;
	
	shared_objects->camera_type = lpt::VIRTUAL;
	lpt::readCamerasFile(cameras_file, cameras); 

	this->generator->setSharedObjects(this->shared_objects);		
	
	generator->setDataPath(shared_objects->output_path);
	generator->read3DTrajectoryFile(trajectory_file);
	generator->project3DFramesTo2D();
	shared_objects->frame_rate = 60;
	shared_objects->image_type = this->generator->getImageCreator()->image_type; 
	cout << shared_objects->image_type.size().height << endl << endl;
	if ( ! cameras[0].frames.empty() ) {
		this->setCamerasStatus(true);
		return true; 
	} else
		return false;
}

void VirtualCameras::addControls() {
	string null = "";
	cv::createTrackbar("FrameRate", null , &frame_rate_level, 50000, 0, 0);
	cv::createTrackbar("Light intensity", null , &this->generator->getImageCreator()->object_intensity, 10E10, 0, 0);
//	cv::createTrackbar("Object size (mm)", null , &this->generator->getImageCreator()->object_size, 20, 0, 0);   //TODO: allow object size to be adjusted during simulation
	cout << "Virtual Camera Controls added" << endl; 
}

bool VirtualCameras::grabFrameGroup(lpt::ImageFrameGroup& frame_group){
	auto& cameras = shared_objects->cameras;

	for ( int c = 0; c < cameras.size(); ++c) {
		lpt::ImageFrame& frame = cameras[c].frames[frame_index];
		frame_group[c].frame_index = frame.frame_index;
		frame_group[c].image = frame.image.clone();			 
	}

	if ( frame_index < cameras[0].frames.size() - 1 )
		++frame_index;
	else 
		frame_index = 0;

	boost::posix_time::microseconds sleeptime(this->frame_rate_level);
	boost::this_thread::sleep(sleeptime);

	return true;
}

void VirtualCameras::shutdown(){
	cout << "\t-------VirtualCameras Shutdown Complete" << endl;
}


#ifdef USE_NP_CAMERASDK 
//---Using Natural Point camera SDK---

/*****lpt::Optitrack class implementation*****/

bool Optitrack::initializeCameras() {
	this->shared_objects->camera_type = lpt::OPTITRACK; 
	auto& cameras = shared_objects->cameras;
	cout << "Searching for attached cameras" << endl;
	CameraLibrary::CameraManager::X().WaitForInitialization();
	cout << ( CameraLibrary::CameraManager::X().AreCamerasInitialized() ? 
			"\t --Complete" : "\t --Failed" ) << endl;
	CameraLibrary::CameraList list;
	
	if(list.Count() <= 0) {
		cout << "No Cameras detected:  Exiting" << endl;
		return false;
	}
	else {
		optitrack_cameras.resize( list.Count() );
		cameras.resize( list.Count() );
		cout << "Number of cameras detected: " << optitrack_cameras.size() << endl; 
	}

	for( int i = 0; i < optitrack_cameras.size(); i++ ) {
		CameraLibrary::Camera* newcamera = CameraLibrary::CameraManager::X().GetCamera( list[i].UID() );
		int camera_id = newcamera->CameraID() - 1;
		optitrack_cameras[camera_id] = newcamera;
		cameras[camera_id].id = camera_id;
		stringstream cameraname;
		cameraname <<  list[i].Name();
		cameras[camera_id].name = cameraname.str();
		cameras[camera_id].sensor_size[0] = optitrack_cameras[camera_id]->ImagerWidth();
		cameras[camera_id].sensor_size[1] = optitrack_cameras[camera_id]->ImagerHeight();
		cameras[camera_id].pixel_size[0] = optitrack_cameras[camera_id]->ImagerWidth() / optitrack_cameras[camera_id]->PhysicalPixelWidth();
		cameras[camera_id].pixel_size[1] = optitrack_cameras[camera_id]->ImagerHeight() / optitrack_cameras[camera_id]->PhysicalPixelHeight();
		
		cout << "Camera " <<  camera_id << ": " << list[i].Name() << endl;
		if(!optitrack_cameras[camera_id]) {
			cout << "Fail! camera does not exist " << camera_id << endl;
			return false;
		}
		//Print camera capabilities
		cout << "\t --Filter Switcher: " << (optitrack_cameras[camera_id]->IsFilterSwitchAvailable() ? "Avaliable" : "Not Avaliable") << endl;
		cout << "\t --MJPEG Mode: " << (optitrack_cameras[camera_id]->IsMJPEGAvailable() ? "Avaliable" : "Not Avaliable") << endl;
		cout << "\t --Default FrameRate = " << optitrack_cameras[camera_id]->ActualFrameRate() << endl;

		//Set some initial camera operating parameters
		optitrack_cameras[camera_id]->SetFrameRate(100);   // Percentage for V120-SLIM Cameras
		optitrack_cameras[camera_id]->SetVideoType(Core::ObjectMode);
		optitrack_cameras[camera_id]->SetLateMJPEGDecompression(false); 
		optitrack_cameras[camera_id]->SetAEC(false);
		optitrack_cameras[camera_id]->SetAGC(true);
		optitrack_cameras[camera_id]->SetExposure(40);
		optitrack_cameras[camera_id]->SetThreshold(40);
		optitrack_cameras[camera_id]->SetIntensity(0);   //IR LED intensity if available
		optitrack_cameras[camera_id]->SetTextOverlay(false);
		optitrack_cameras[camera_id]->SetObjectColor(255);
		optitrack_cameras[camera_id]->SetIRFilter(true);
	}

	shared_objects->frame_rate = optitrack_cameras[0]->ActualFrameRate();
	cout << "\t Actual Frame Rate = " << this->getFrameRate() << endl;
	
	sync = CameraLibrary::cModuleSync::Create();
	for ( int i = 0; i < optitrack_cameras.size(); ++i) 
		sync->AddCamera(optitrack_cameras[i]);

	sensor_dim = cv::Size(optitrack_cameras[0]->Width(), optitrack_cameras[0]->Height() );
	cout << "Camera sensor size [" << sensor_dim.width << " x " << sensor_dim.height << "]" << endl;

	shared_objects->image_type = cv::Mat(sensor_dim, CV_8UC1);
	
	optitrack_frames.resize( sync->CameraCount() );
	this->setCamerasStatus(true);

	for ( int i = 0; i < optitrack_cameras.size(); ++i) 
		optitrack_cameras[i]->Start();

	return true;
}

void Optitrack::shutdown() {
	
	cout << "------------Shutting Down Cameras------------" << endl;
	sync->RemoveAllCameras();
	cModuleSync::Destroy( sync );

	for ( int i = 0; i < optitrack_cameras.size(); ++i) {
		optitrack_cameras[i]->Release();
		//cout << i << " ";
	}
	CameraLibrary::CameraManager::X().Shutdown();
	cout << "\t---Optitrack Shutdown Complete" << endl;
}

bool Optitrack::grabFrameGroup(lpt::ImageFrameGroup& frame_group) {
			
	CameraLibrary::FrameGroup* native_frame_group = sync->GetFrameGroup();
	if(native_frame_group && native_frame_group->Count() == shared_objects->cameras.size()) {
		auto& image_type = shared_objects->image_type;
		for ( int camera_id = 0; camera_id < native_frame_group->Count(); ++camera_id ) { 
			optitrack_frames[camera_id] = native_frame_group->GetFrame(camera_id);

			frame_group[camera_id].image = image_type.clone();

			frame_group[camera_id].frame_index = optitrack_frames[camera_id]->FrameID();

			if( this->collect_particledata_from_camera ) {
				frame_group[camera_id].particles.clear();
				for (int object_id = 0; object_id < optitrack_frames[camera_id]->ObjectCount(); object_id++) {
					auto object = optitrack_frames[camera_id]->Object(object_id);
					frame_group[camera_id].particles.push_back
						(
							std::move( 
								lpt::ParticleImage::create( object_id, object->X(), object->Y(), object->Radius() ) 
							) 
						);
				}
			}
			else {
				cv::Mat temp = image_type.clone();
				optitrack_frames[camera_id]->Rasterize(image_type.size().width, image_type.size().height, 
					static_cast<unsigned int>(temp.step), 8, temp.data);
				cv::Rect ROI(0, 0, image_type.size().width/2, image_type.size().height/2);
				cv::resize(temp(ROI), frame_group[camera_id].image, image_type.size());
			}
			optitrack_frames[camera_id]->Release();
		}
		native_frame_group->Release();
		if(sync->LastFrameGroupMode() != CameraLibrary::FrameGroup::Hardware) {
			for ( int camera_id = 0; camera_id < native_frame_group->Count(); ++camera_id ) 
				frame_group[camera_id].image = image_type.clone();
			cout << "\t Cameras NOT Synchronized: Frame # = " << frame_count << endl;  
		}
		++frame_count;
		return true;
	} else 
		return false;
	
}

void Optitrack::addControls() {
	void* optitrack_void_ptr = static_cast<void*> ( this );
	void* cameras_void_ptr = static_cast<void*> ( &optitrack_cameras );

	string null;
	auto camera = optitrack_cameras[0];
	init_video_mode = 2;		// {0 = Precision, 1 = Segment, 2 = Object, 3 = MJPEG Mode}
	init_threshold = 40;
	int max_threshold = camera->MaximumThreshold();
	int min_threshold = camera->MinimumThreshold();
	init_exposure = 40;
	int max_exposure = camera->MaximumExposureValue();
	int min_exposure = camera->MinimumExposureValue();
	init_intensity = 5;
	init_framerate_mode = 2;    //Mode number: {2 = 100%, 1 = 50%, 0 = 25%} for V120:SLIM cameras 
	
	cv::createButton("IR Filter", callbackSetIRFilter, cameras_void_ptr ,CV_CHECKBOX, 1 );
	cv::createButton("Automatic Gain Control", callbackSetAGC, cameras_void_ptr, CV_CHECKBOX, 1 );
	cv::createButton("Autoamtic Exposure Control", callbackSetAEC, cameras_void_ptr, CV_CHECKBOX, 0 );
	cv::createButton("Text Overlay", callbackSetTextOverlay, cameras_void_ptr, CV_CHECKBOX, 0 );
	cv::createTrackbar("VideoMode", null , &init_video_mode, 3, callbackSetVideoType, optitrack_void_ptr);
	cv::createTrackbar("Threshold", null , &init_threshold, max_threshold-min_threshold, callbackSetThreshold, cameras_void_ptr);
	cv::createTrackbar("Exposure", null , &init_exposure, max_exposure-min_exposure, callbackSetExposure, cameras_void_ptr);
	cv::createTrackbar("FrameRateMode", null , &init_framerate_mode, 2, callbackSetFrameRate, optitrack_void_ptr);
}

//-------------OPTITRACK CALL BACK FUNCTIONS ------------------------------
void callbackSetVideoType( int mode, void* data )
{
	Optitrack* system = static_cast< Optitrack*> (data);
	vector<CameraLibrary::Camera*>& optitrack_cameras = system->optitrack_cameras;
	Core::eVideoMode new_mode;

	cout << "Setting Camera Mode: ";
	switch (mode) {
	case 0:
		cout << "Precision" << endl;
		new_mode = Core::PrecisionMode;
		system->setParticleCollectionFromCamera(true);
		break;
	case 1:
		cout << "Segment" << endl;
		new_mode = Core::SegmentMode;
		system->setParticleCollectionFromCamera(true);
		break;
	case 2:
		cout << "Object" << endl;
		new_mode = Core::ObjectMode;
		system->setParticleCollectionFromCamera(true);
		break;
	case 3:
		cout << "MJPEG" << endl;
		new_mode = Core::MJPEGMode;
		system->setParticleCollectionFromCamera(false);
		break;
	default:
		new_mode = Core::ObjectMode;
		break;
	}

	for (int id = 0; id < optitrack_cameras.size(); id++) { 
		if (new_mode == Core::MJPEGMode && optitrack_cameras[id]->IsMJPEGAvailable() == false)
			cout << "Camera " << id << " does not support MJPEG Mode" << endl;
		else
			optitrack_cameras[id]->SetVideoType(new_mode);
	}
}

void callbackSetThreshold( int value, void* data) {
	vector<CameraLibrary::Camera*>& optitrack_cameras = *static_cast< vector<CameraLibrary::Camera*>* > (data);
	
	value += optitrack_cameras[0]->MinimumThreshold();
	cout << "Setting Threshold value: " << value << endl;

	for (int id = 0; id < optitrack_cameras.size(); id++)
		optitrack_cameras[id]->SetThreshold( value );
}

void callbackSetIRFilter(int state, void* data) {
	vector<CameraLibrary::Camera*>& optitrack_cameras = *static_cast< vector<CameraLibrary::Camera*>* > (data);
	cout << "Setting IR Filter: " << (state ? "on":"off") << endl;

	for (int id = 0; id < optitrack_cameras.size(); id++) {
		if (optitrack_cameras[id]->IsFilterSwitchAvailable() ) 
			optitrack_cameras[id]->SetIRFilter( state ? true : false );
	}
}

void callbackSetAEC(int state, void* data) {
	// Enable or Disable autmatic exposure control AEC
	vector<CameraLibrary::Camera*>& optitrack_cameras = *static_cast< vector<CameraLibrary::Camera*>* > (data);
	cout << "Setting Auto Exposure Control (AEC):" << (state ? "on":"off") << endl;

	for (int id = 0; id < optitrack_cameras.size(); id++) 
		optitrack_cameras[id]->SetAEC(state ? true : false);
}

void callbackSetAGC(int state, void* data) {
	// Enable or Disable autmatic gain control AGC
	vector<CameraLibrary::Camera*>& optitrack_cameras = *static_cast< vector<CameraLibrary::Camera*>* > (data);
	cout << "Setting Auto Gain Control (AGC): " << (state ? "on":"off") << endl;

	for (int id = 0; id < optitrack_cameras.size(); id++) 
		optitrack_cameras[id]->SetAGC(state ? true : false);		
}

void callbackSetTextOverlay(int state, void* data) {
	// Enable or Disable text overlay
	vector<CameraLibrary::Camera*>& optitrack_cameras = *static_cast< vector<CameraLibrary::Camera*>* > (data);
	cout << "Setting Text Overlay: " << (state ? "on":"off") << endl;

	for (int id = 0; id < optitrack_cameras.size(); id++) 
		optitrack_cameras[id]->SetTextOverlay(state ? true : false);		
}

void callbackSetExposure(int value, void* data) {
	// Sets Exposure manually when AEC is not activated
	vector<CameraLibrary::Camera*>& optitrack_cameras = *static_cast< vector<CameraLibrary::Camera*>* > (data);

	value += optitrack_cameras[0]->MinimumExposureValue();
	cout << "Setting Exposure: Desired = " << value;

	for (int id = 0; id < optitrack_cameras.size(); id++) 
		optitrack_cameras[id]->SetExposure(value);
	Sleep(2);
	//int actual = optitrack_cameras[0]->Exposure();
	//cout << ", Actual = "<<  actual << "("<< 17.13 * actual << " micro seconds)" << endl; //Only valid for V120:SLIM
}

void callbackSetIntensity(int value, void* data) {
	// Sets IR illumination intensity for optitrack_cameras that support this feature
	vector<CameraLibrary::Camera*>& optitrack_cameras = *static_cast< vector<CameraLibrary::Camera*>* > (data);
	cout << "Setting Intensity: " << "Desired = "<< value;

	for (int id = 0; id < optitrack_cameras.size(); id++) 
		optitrack_cameras[id]->SetIntensity(value);
	Sleep(2);
	cout << ", Actual = "<< optitrack_cameras[0]->Intensity()<< endl;
}

void callbackSetFrameRate(int value, void* data) {
	lpt::Optitrack* camera_system = static_cast< lpt::Optitrack* > ( data );
	if (camera_system) {
		int fps;
		switch(value) {
		case 2:
			fps = 120;
			break;
		case 1:
			fps = 60;
			break;
		case 0:
			fps = 30;
			break;
		default:
			fps = 120;
			break;
		}
		cout << "Setting Frame Rate: Desired  = " << fps;
		vector<CameraLibrary::Camera*>& optitrack_cameras = camera_system->getOptitrackCameras();
		for (int id = 0; id < optitrack_cameras.size(); id++) 
			optitrack_cameras[id]->SetFrameRate(fps);
		Sleep(2);
		cout << ", Actual FPS = " << optitrack_cameras[0]->ActualFrameRate() << endl;
		camera_system->shared_objects->frame_rate = optitrack_cameras[0]->ActualFrameRate();
	}
}

#endif /*USE_NP_CAMERASDK*/

void callbackSetImageViewStatus( int state, void* data ) {
	StreamingPipeline* system = static_cast< StreamingPipeline*> (data);
	if (system)
		system->setImageViewStatus( (state != 0) );
	else
		cout << "---INVALID pointer to CameraSystem:  Cannot Set ImageViewStatus" << endl;

}

void callbackSetCompositeView(int state, void* data) {
	StreamingPipeline* system = static_cast< StreamingPipeline*> (data);
	if (system)
		system->setCompositeView( (state != 0) );
	else
		cout << "---INVALID pointer to CameraSystem:  Cannot Set CompositeView" << endl;
}

void callbackSetDetectionView(int state, void* data) {
	StreamingPipeline* system = static_cast< StreamingPipeline*> (data);
	if (system)
		system->setDetectionView( (state != 0) );
	else
		cout << "---INVALID pointer to CameraSystem:  Cannot Set DetectionView" << endl;
}

void callbackSetReprojectionView( int state, void* data ) {
	StreamingPipeline* system = static_cast< StreamingPipeline*> (data);
	if (system)
		system->setReprojectionView( (state != 0) );
	else
		cout << "---INVALID pointer to CameraSystem:  Cannot Set ReprojectionView" << endl;
}

void callbackSetTrajectoryView( int state, void* data ) {
	StreamingPipeline* system = static_cast< StreamingPipeline*> (data);
	if (system) {
		system->tracker->clear_drawings = true;
		system->setTrajectoryView( (state != 0) );
	}
	else
		cout << "---INVALID pointer to CameraSystem:  Cannot Set TrajectoryView" << endl;
}

void callbackFlushFrameQueue(int state, void* data) {
	StreamingPipeline* system = static_cast< StreamingPipeline*> (data);
	if (system) {
		system->frame_queue.clear();
	}
	else
		cout << "---INVALID pointer to CameraSystem:  Cannot flush queue!!" << endl;
}

void callbackFlushProcessedQueue(int state, void* data) {
	StreamingPipeline* system = static_cast< StreamingPipeline*> (data);
	if (system) {
		system->processed_queue.clear();
	}
	else
		cout << "---INVALID pointer to CameraSystem:  Cannot flush queue!!" << endl;
}

void callbackStopCameras(int state, void* data) {
	StreamingPipeline* system = static_cast< StreamingPipeline*> (data);
	if (system)
		system->getCameraSystem()->setCamerasStatus(false);
	else
		cout << "---INVALID pointer to CameraSystem:  Cannot Stop!!" << endl;
}


/*****lpt::StreamingPipeline class implementation*****/

bool StreamingPipeline::initialize() {	
	cout << "Initializing streaming pipeline " << endl;
	frame_queue.setCapacity(queue_capacity); 
	processed_queue.setCapacity(queue_capacity);
	match_queue.setCapacity(queue_capacity);
	frame3D_queue.setCapacity(queue_capacity);
	monitor_queue.setCapacity(queue_capacity);
	
	cout << "Concurrent Queue capacities set to " << queue_capacity << endl;
	
	bool cameras_ok = false;
	
	if (camera_system)
		cameras_ok = camera_system->initializeCameras();
	else
		cout << "No camera system found. Make sure to call StreamingPipeline::attachCameraSystem()" << endl;

	if (cameras_ok) {
		auto& cameras = shared_objects->cameras;
		auto& camera_pairs = shared_objects->camera_pairs;
		auto& output_path = shared_objects->output_path;
		if (! visualizer)
			this->visualizer = std::make_shared < lpt::Visualizer > ();
		this->visualizer->setSharedObjects(this->shared_objects);
		this->visualizer->initialize();

		this->calibrator = std::make_shared < lpt::Calibrator > (cameras, camera_pairs, camera_displayed);
		this->calibrator->setCalibViews(shared_objects->image_type);
		this->matcher->initialize();
		this->tracker->setTrajectoryViews(cameras, shared_objects->image_type);
		this->recorder = std::make_shared < lpt::Recorder > (cameras.size(), output_path );
		this->reconstructor = std::make_shared < lpt::Reconstruct3D >(); //std::make_shared < lpt::Reconstruct3DwithSVD > ();
		this->reconstructor->setSharedObjects(this->shared_objects);
				
		//if ( camera_pairs.empty() )       
		//	lpt::generateCameraPairs(this->pt_cameras, this->camera_pairs);            
		
		initializeControlWindow();
	}

	return cameras_ok;
}

void StreamingPipeline::initializeControlWindow() {
	auto& cameras = shared_objects->cameras;
	
	if (cameras.empty() ) {
		cout << "Could not initialize control window: No cameras found" << endl;
		return;
	}
	cout << "Initializing Control Window with " << cameras.size() << " Cameras" << endl;
	// Set up display window using opencv
	string null;  
	camera_displayed = 0;       // index of initial camera to be displayed in opencv window
	
	image_view_status = true;
	composite_view_status = false;
	detection_view_status = false;
	reprojection_view_status = false;
	trajectory_view_status = false;

	void* system_void_ptr = static_cast<void*> ( this );
	void* recorder_void_ptr = static_cast<void*> ( this->recorder.get() );
	cv::namedWindow( camera_system->getWindowName() );
	cv::createTrackbar("Camera", camera_system->getWindowName() , &camera_displayed, static_cast<int>(cameras.size() - 1), 0);
	
	cv::createButton("Record Video Clip", callbackRecordVideo, recorder_void_ptr, CV_PUSH_BUTTON, 0 );
	cv::createButton("Take Snapshot", callbackTakeSnapshot, recorder_void_ptr, CV_PUSH_BUTTON, 0 );
	cv::createButton("Stop Cameras", callbackStopCameras, system_void_ptr , CV_PUSH_BUTTON, 0 );
	cv::createButton("Show ImageView", callbackSetImageViewStatus, system_void_ptr , CV_CHECKBOX, 1 );
	cv::createButton("Show Composite", callbackSetCompositeView, system_void_ptr , CV_CHECKBOX, 0 );
	cv::createButton("Show Detected", callbackSetDetectionView, system_void_ptr , CV_CHECKBOX, 0 );
	cv::createButton("Reproject 3D", callbackSetReprojectionView, system_void_ptr , CV_CHECKBOX, 0 );

	cv::createButton("Flush raw queue", callbackFlushFrameQueue, system_void_ptr, CV_PUSH_BUTTON, 0 );
	cv::createButton("Flush proc queue", callbackFlushProcessedQueue, system_void_ptr, CV_PUSH_BUTTON, 0 );
	this->camera_system->addControls();
	this->calibrator->addControls();
	this->processor->addControls();
	this->detector->addControls();
	this->matcher->addControls();
	this->tracker->addControls();
	this->visualizer->addControls();
			
	cv::waitKey(50);
}

void StreamingPipeline::runControlWindow() {
	cout << "Thread: Displaying images" << endl;
	auto& cameras = shared_objects->cameras;
	string window_name = camera_system->getWindowName();
	if (cameras.empty())
		return;

	// Displays images and controls camera parameters
	if (calibrator) {
		cout<< "calibrator is good" << endl;	
	}
	if (recorder)
		cout <<"recorder is good" << endl;

	int last_frame_index = 0;
	long int count = 0, count2 = 0;

	while( camera_system->areCamerasRunning() ) {
		lpt::Frame3d_Ptr frame3d;

		if( monitor_queue.try_pop( frame3d ) ) {

			lpt::ImageFrameGroup& camera_frames = frame3d->camera_frames;

			if( recorder->isSnapShotRequested() )
				recorder->takeSnapShot( lpt::getImageVector( camera_frames ) );
			if( recorder->isVideoRecording() )
				recorder->addFramesToVideos( lpt::getImageVector( camera_frames ) );

			if ( calibrator->isStereoDataCollecting() ) 
				if ( count % calibrator->getFrameInterval() == 0) 
					calibrator->findCalibObject(camera_frames );

			if ( calibrator->isIntParamDataCollecting() ) 
				if ( count % calibrator->getFrameInterval() == 0) 
					calibrator->findCalibBoard( camera_frames[camera_displayed] );

			if ( calibrator->isSettingGlobalReference() ) {
				bool ref_found = calibrator->findGlobalReference( camera_frames );
				if (ref_found)
					calibrator->setGlobalReference(false);
			}

			if( this->showReprojectionView() && ! frame3d->objects.empty() ) {
				reconstructor->draw( *frame3d );
			}

			if ( count2 >= shared_objects->frame_rate / 20 ) {
				count2 = 0;
				stringstream capturedetails;
				if ( this->showCompositeView() ) {
					for (int c = 0; c < cameras.size(); ++c) {
						capturedetails.str("");
						capturedetails << "Objects = " << camera_frames[c].particles.size() << "\t3D Objects = " << frame3d->objects.size() << "\t monitor_queue size = " << monitor_queue.size();
						cv::imshow(cameras[c].name, camera_frames[c].image);
						cv::displayStatusBar(cameras[c].name, capturedetails.str(), 1000);
					}
				}
				else {
					capturedetails << "Camera #" << cameras[camera_displayed].id << ": " << cameras[camera_displayed].name <<
						"\n\t2D# = " << camera_frames[camera_displayed].particles.size() <<
						"\t3D# = " << frame3d->objects.size() << "\tvis_queue = " << visualizer->getQueueSize() << "\tren_queue = " << visualizer->getRenderQueueSize();

					cv::displayStatusBar(window_name, capturedetails.str(), 1000);
					if ( this->getImageViewStatus() )
						cv::imshow(window_name, camera_frames[camera_displayed].image);
				}
				capturedetails.str("");
				capturedetails << "\t frame_queue size = " << frame_queue.size() 
					<< "\t processed_queue size = " << processed_queue.size()
					<< "\t match_queue size = " << match_queue.size()
					<< "\t frame3D_queue size = " << frame3D_queue.size()
					<< "\t monitor_queue size = " << monitor_queue.size();
				cv::displayOverlay(window_name, capturedetails.str(), 100);			

				if (count % 60 == 0)
					cv::waitKey(1);
			}
			++count2;
			++count;
		} else
			cv::waitKey(1);
			// Check for frame queue overload:  this depends on the type of camera since frameID is defined differently
			//if( monitor_queue.front()[camera_displayed].frame_index - last_frame_index != 1) {
			//	stringstream overloaddisplay;
			//	overloaddisplay << "FrameQueue Overload: Queue size = " << monitor_queue.size();
			//	cv::displayOverlay(window_name, overloaddisplay.str(), 100);
			//}
			//last_frame_index = monitor_queue.front()[camera_displayed].frame_index;
	}
	cout << "Control Window thread done" << endl;
	cv::destroyAllWindows();
}

void StreamingPipeline::aquireImageData() {
	auto& cameras = shared_objects->cameras;
	cout << "------------Starting " << cameras.size() << " cameras-------------" << endl;
	this->setFrameRate( camera_system->getFrameRate() );
	boost::posix_time::microseconds sleeptime(100);

	while( camera_system->areCamerasRunning() ) {
		lpt::ImageFrameGroup frame_group( cameras.size() );
		bool good_frames = camera_system->grabFrameGroup( frame_group );

		if ( good_frames && ! frame_queue.full() ) {
			frame_queue.push( std::move(frame_group) );
		}
		else
			boost::this_thread::sleep(sleeptime);
	}
}

void StreamingPipeline::processImages(int index) {
	auto& cameras = shared_objects->cameras;
	auto& image_type = shared_objects->image_type;

	cout << "Thread #" << index << " processing images" << endl;
//	boost::chrono::system_clock::time_point start, stop;
	cv::Mat map1, map2;
	cv::Mat temp_image;
	cv::Mat cam_mat = cameras[index].getCameraMatrix();
	cv::Mat dist_coeffs = cameras[index].getDistCoeffs();
	cv::Size image_size = image_type.size();
	//cv::initUndistortRectifyMap(cam_mat, dist_coeffs, cv::Mat(), cv::Mat(), image_size, CV_32FC1, map1, map2);   
	bool first = true;
	while( camera_system->areCamerasRunning() ) {
		
		boost::unique_lock<boost::mutex> lock(imageproc_mutex, boost::try_to_lock);
		imageproc_barrier->wait();
		if ( lock.owns_lock() ) {
			//stop = boost::chrono::system_clock::now();
			/*if ( ! first) {
				this->imageproc_time += stop - start;
			} else 
				first = false;*/

			if ( !processed_queue.full() ) 
				processed_queue.push( std::move(imageproc_frames) );	
			frame_queue.wait_and_pop(imageproc_frames);
			//start = boost::chrono::system_clock::now();
	
		}
		imageproc_barrier->wait();	
		//start = boost::chrono::system_clock::now();		
		if ( this->camera_system->getParticleCollectionFromCamera() == false ) {
			temp_image = imageproc_frames[index].image.clone();
			//cout << "empty? " << temp_image.empty() << endl;
			processor->processImage( temp_image );
			detector->detectFeatures( temp_image, imageproc_frames[index].particles, imageproc_frames[index].contours );
		}
		if ( this->showDetectionView() ) {
			detector->drawResult( imageproc_frames[index] );
			//detector->drawContours( imageproc_frames[index].image, imageproc_frames[index].contours );
		}
		
		lpt::undistortPoints(cameras[index], imageproc_frames[index] );
	
		boost::this_thread::interruption_point();	
	}
	
	cout << index << " Image Processor thread done:" <<endl; //<< this->imageproc_time << endl;
}

void StreamingPipeline::solveCorrespondence() {
	cout << "Correspondence Thread started" << endl;
//	boost::chrono::system_clock::time_point start, stop;
//	start = boost::chrono::system_clock::now();
	matcher->run(&processed_queue, &match_queue, 1,1); 
	
//	stop = boost::chrono::system_clock::now();
//	this->correspond_time = stop - start;
	cout << "Matching thread done: " << endl; //this->correspond_time << endl;
}

void StreamingPipeline::reconstuct3DObjects() {
	cout << "Reconstruct 3D Objects Thread started" << endl;
//	boost::chrono::system_clock::time_point start, stop;
	
	while( camera_system->areCamerasRunning() ) {	
		std::pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > newpair;
		match_queue.wait_and_pop(newpair);
		//start = boost::chrono::system_clock::now();
		lpt::ImageFrameGroup& frame_group = newpair.first;
		vector<lpt::Match::Ptr>& matches =  newpair.second;

		if (this->visualizer->getTakeImageMeasurement() ) {
			this->visualizer->accumulateCentroidDetectionUncertainty( matches );
		}

		lpt::Frame3d_Ptr newframe3d = lpt::Frame3d::create(frame_group, frame_group[0].frame_index);
		reconstructor->reconstruct3DFrame( matches, *newframe3d ); 

		if( !frame3D_queue.full() )
			frame3D_queue.push( std::move(newframe3d) );
		//stop = boost::chrono::system_clock::now();
		//this->recon_time += stop - start;
		boost::this_thread::interruption_point();	
	}
	//stop = boost::chrono::system_clock::now();
	cout << "3D reconstruction thread done: " << endl;//this->recon_time <<  endl;
}

void StreamingPipeline::trackObjects() {
	cout << "Track Objects Thread started" << endl;
	boost::posix_time::microseconds sleeptime(500);
//	boost::chrono::system_clock::time_point start, stop;
	int count = 0;
	int last_frame_index = 0;
	lpt::Frame3d_Ptr frame1, frame2;

	while( camera_system->areCamerasRunning() ) {
		
        if( !frame1 ) {
			frame3D_queue.wait_and_pop(frame1);
        }
		frame3D_queue.wait_and_pop(frame2);

		//start = boost::chrono::system_clock::now();		
		if (last_frame_index != frame1->frame_index ) {
            cout << "!!!!!!!!!!Frame skipped!!!!!!!!!!!!! " << frame1->frame_index << endl;
			last_frame_index = frame2->frame_index;
			frame1 = frame2;
			continue;
		}

		last_frame_index = frame2->frame_index;

        tracker->trackFrames(*frame1, *frame2);
		
		if ( visualizer->getVisualizationStatus() ) 
			visualizer->addTrajectoriesToQueue( tracker->getActiveTrajectories() );
		
		if ( !monitor_queue.full() ) 	
			monitor_queue.push( std::move( frame1 ) );	

		frame1 = frame2;
		++count;
		//stop = boost::chrono::system_clock::now();	
		//this->tracking_time += stop - start;
	
		boost::this_thread::interruption_point();	
	}
	//tracker->finalizeTrajectoryList();
	cout << "tracking thread done: " << endl;// this->tracking_time << endl;
}

void StreamingPipeline::run() {
	auto& cameras = shared_objects->cameras;

	imagegrabber_thread = boost::thread(&StreamingPipeline::aquireImageData, this); 
	
	imageproc_barrier = std::make_shared<boost::barrier>( cameras.size() );
	
	lpt::ImageFrame init_frame;
	init_frame.image = shared_objects->image_type;

	imageproc_frames.resize(cameras.size(), init_frame);
	for (int c = 0; c < cameras.size(); ++c ) 
		imageproc_workers_thread_group.create_thread( boost::bind(&StreamingPipeline::processImages, this, c) );
		
	matcher_thread = boost::thread(&StreamingPipeline::solveCorrespondence, this);
	reconstructor_thread = boost::thread(&StreamingPipeline::reconstuct3DObjects, this);
	tracker_thread = boost::thread(&StreamingPipeline::trackObjects, this);
	visualizer->setCameras( cameras );
	visualizer_thread = boost::thread(&StreamingPipeline::runVisualizer, this);
	
#ifdef WIN32
	// Priority settings available: THREAD_PRIORITY_TIME_CRITICAL, THREAD_PRIORITY_HIGHEST, THREAD_PRIORITY_ABOVE_NORMAL,  
	//	  THREAD_PRIORITY_NORMAL, THREAD_PRIORITY_BELOW_NORMAL, THREAD_PRIORITY_LOWEST

	SetThreadPriority(imagegrabber_thread.native_handle(), THREAD_PRIORITY_TIME_CRITICAL);

	HANDLE mainthread = GetCurrentThread();
	SetThreadPriority(mainthread, THREAD_PRIORITY_TIME_CRITICAL);
	SetThreadPriority(matcher_thread.native_handle(), THREAD_PRIORITY_HIGHEST);
	SetThreadPriority(reconstructor_thread.native_handle(), THREAD_PRIORITY_HIGHEST);
	SetThreadPriority(tracker_thread.native_handle(), THREAD_PRIORITY_HIGHEST);
	
#else 
		cout << "need to set thread priorities for Linux systems" << endl;
#endif

	boost::this_thread::sleep( boost::posix_time::seconds(3) );
	this->runControlWindow();
	// All threads running...stop called when cameras are stopped
	this->stop();	
}

void StreamingPipeline::stop() {
	
	imagegrabber_thread.join();

	imageproc_workers_thread_group.interrupt_all();
	imageproc_workers_thread_group.join_all();

	matcher_thread.interrupt();
	matcher_thread.join();
	matcher->stop();

	reconstructor_thread.interrupt();
	reconstructor_thread.join();

	tracker_thread.interrupt();
	tracker_thread.join();
	this->camera_system->shutdown();

	visualizer_thread.join();
}

void StreamingPipeline::load_Rotation_Matrix() 
{
	ifstream S_in, P_in;
	string S = this->shared_objects->input_path + "S.txt";
	string P = this->shared_objects->input_path + "P.txt";
	S_in.open(S);
	P_in.open(P);

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			S_in >> this->shared_objects->S[i][j];
		}
		P_in >> this->shared_objects->P[i];
	}
	this->shared_objects->isRotation_Correction = true;
}

} /*NAMESPACE PT*/
