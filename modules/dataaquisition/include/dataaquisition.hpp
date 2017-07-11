/*
Data aquisition module headers:
These classes provide data aquisition funcionality to the particle tracking system and the overall multi-threaded streaming framework.

copyright: Douglas Barker 2011
*/

#ifndef DATAAQUSITION_H_
#define DATAAQUSITION_H_

#include <core.hpp>
#include <calib.hpp>
#include <datagen.hpp>
#include <imageproc.hpp>
#include <correspond.hpp>

#ifdef USE_CUDA
#include <correspondcuda.h>
#endif

#include <tracking.hpp>
#include <visualization.hpp>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/chrono.hpp>

#ifdef USE_NP_CAMERASDK
#include <windows.h>
#include "cameralibrary.h"
#endif


namespace lpt {

	using namespace std;
	using namespace CameraLibrary;
		
	class Video;
	class Recorder;
	class CameraSystem;
	class VirtualCameras;
	class Optitrack;
	class StreamingPipeline;
	
    /*
     * This class provides the video object that the Recorder class will use to store images and convert into avi files
     */
	class Video {
	public: 
		typedef std::shared_ptr<Video> Ptr;
		static inline Video::Ptr create(string path = "./") { return Video::Ptr( new Video(path) ); }

		vector<cv::Mat> video_frames;
		Video(string path = "./") : file_path(path) {}
		inline void setFilePath(string path) { file_path = path; }
		inline void addImageToVideo(cv::Mat image) { video_frames.push_back(image); }
		inline size_t getVideoLength() const { return video_frames.size(); }
		void writeAVI(int playback_fps, int codec);
		void writeImages(string basename, string filetype = ".jpg");
		void readAVI();
		void readImages(string basename, int numframes, string filetype = ".jpg");
		
		virtual ~Video() {}

	private:
		string file_path;
	};

    /*
     * This class provides the functionality to record avi files from the camera system images
     */
	class Recorder {
	public:
		typedef std::shared_ptr<lpt::Recorder> Ptr;
		static inline lpt::Recorder::Ptr create(int num_cams, string path = "./", int length = 900) { 
			return lpt::Recorder::Ptr( new lpt::Recorder(num_cams, path, length) ); 
		}

		Recorder(int num_cams, string path = "./", int length = 900) 
			: number_of_cameras(num_cams), path(path), clip_length(length)
		{ 
			videos.resize( number_of_cameras ); 
			imagelists.resize( number_of_cameras );
			snapshot_requested = false;
			record_video = false;
			snapshot_id = -1;
			video_id = -1;
		}
		
		virtual ~Recorder() {
			if ( getVideoCount() > 0 ) {
				cout << "Writing " << getVideoCount() << " video sequence to file for each camera" << endl;
				int playback_rate;
				cout << "Enter desired playback frame rate" << endl;
				cin >> playback_rate;
				writeVideosToFile(playback_rate);
			}
			if ( getSnapShotCount() > 0 )
				writeSnapShotImageLists(); 
		}

		inline void setClipLength(int length) { clip_length = length;}
		inline void setNumCameras(int num) {
			number_of_cameras = num; 
			videos.resize( number_of_cameras ); 
			imagelists.resize( number_of_cameras ); 
		} 
		inline void requestSnapShot() { snapshot_requested = true; }
		inline bool isSnapShotRequested() const { return snapshot_requested; }
		inline bool isVideoRecording() const { return record_video; }
		inline size_t getVideoCount() const { return videos[0].size(); }
		inline size_t getSnapShotCount() const { return imagelists[0].size(); }
		void takeSnapShot( vector<cv::Mat>& frames );
		void writeSnapShotImageLists();
		void createVideos();
		void addFramesToVideos(vector<cv::Mat>& frames);
		void writeVideosToFile(int playback_fps = 10, int codec = -1);
	
	private:
		bool snapshot_requested;
		bool record_video;
		int video_id;
		int snapshot_id;
		int clip_length;
		int number_of_cameras;
		vector<vector<string> > imagelists;
		vector<vector<Video> > videos;
		string path;
	
	};

    /*
     * This is a generic base class for the camera system used in particle tracking
     */
	class CameraSystem {
	public:	
		typedef std::shared_ptr < lpt::CameraSystem > Ptr;

		CameraSystem() : collect_particledata_from_camera(true) { }

		virtual void addControls()=0;
		virtual bool initializeCameras()=0;
		virtual bool grabFrameGroup(lpt::ImageFrameGroup& frame_group)=0;
		virtual void shutdown()=0;	
		virtual ~CameraSystem() {}
		bool initialize();	
		
		void loadCameraParams(const string filepath);
		void loadCameraPairParams(const string filepath);

		inline shared_ptr < lpt::SharedObjects > getSharedObjects() { return shared_objects; }
		inline void setSharedObjects(shared_ptr < lpt::SharedObjects > new_shared_objects) { shared_objects = new_shared_objects; }
		
		inline void setInputDataPath(string& path) { shared_objects->input_path = path; }
		inline string getInputDataPath() { return shared_objects->input_path; }
		inline void setOutputDataPath(string& path) { shared_objects->output_path = path; }
		inline string getOutputDataPath() { return shared_objects->output_path; }		

		inline void setImageType(cv::Mat image_type) {shared_objects->image_type = image_type;}
		inline cv::Mat getImageType() {return shared_objects->image_type;}
	
		inline void setParticleCollectionFromCamera(bool state) { collect_particledata_from_camera = state; }
		inline bool getParticleCollectionFromCamera() { return collect_particledata_from_camera; }
		inline string getWindowName() { return window_name; }
		inline bool areCamerasRunning() const { return cameras_status;}
		inline void setCamerasStatus(bool state) { cameras_status = state; cout << "Cameras " << ( state ? "on" : "off" ) << endl; }
		inline int getFrameRate() {return shared_objects->frame_rate;}
		
	protected:
		shared_ptr < lpt::SharedObjects > shared_objects;
		
		bool cameras_status;
		bool collect_particledata_from_camera;
		string window_name;
	};

    /*
     * This class provides virtual cameras to allow simulation of particle images based on synthetic trajectories given in a file
     */
	class VirtualCameras : public CameraSystem {
	public:
		typedef std::shared_ptr<VirtualCameras> Ptr;
		static VirtualCameras::Ptr create( string cams_file, string traj_file ) { return std::make_shared<VirtualCameras> (cams_file, traj_file); }

		VirtualCameras( string cams_file, string traj_file ) : cameras_file(cams_file), trajectory_file(traj_file), frame_rate_level(10000), frame_index(0) {
			this->generator = std::make_shared<lpt::DataSetGenerator>();
			window_name = "Virtual Camera Control Window";
			cout << "Created Virtual Camera System" << endl;
		} 
		
		virtual bool initializeCameras();
		virtual void addControls();
		virtual bool grabFrameGroup(lpt::ImageFrameGroup& frame_group);
		void shutdown();
		inline lpt::DataSetGenerator::Ptr getGenerator() { return generator; }
		
		int frame_rate_level;
	private:
		string cameras_file;
		string trajectory_file;
		lpt::DataSetGenerator::Ptr generator;
		int frame_index;
	};

#ifdef USE_NP_CAMERASDK
	
    /*
     * This class identifies, controls, and grabs images/data from Optitrack motion capture cameras using the Natural Point C++ camera SDK
     */
	class Optitrack : public CameraSystem {
	public:		
		typedef std::shared_ptr<Optitrack> Ptr;
		static Optitrack::Ptr create() { return make_shared< lpt::Optitrack >(); }

		Optitrack() : frame_count (0) 
		{ 
			window_name = "Optitrack Control Window"; 
		}

		bool initializeCameras();
		void addControls();
		bool grabFrameGroup(lpt::ImageFrameGroup& frame_group);	
		void shutdown();
		
		vector<CameraLibrary::Camera*>& getOptitrackCameras() { return optitrack_cameras; }

		//----CALL BACK FUNCTIONS-----
		friend void callbackSetVideoType( int mode, void* data );
		friend void callbackSetThreshold( int value, void* data );
		friend void callbackSetIRFilter(int state, void* data );
		friend void callbackSetAEC( int state, void* data );
		friend void callbackSetAGC( int state, void* data );
		friend void callbackSetTextOverlay( int state, void* data );
		friend void callbackSetExposure( int value, void* data );
		friend void callbackSetIntensity( int value, void* data );
		friend void callbackSetFrameRate( int value, void* data );

	private:
		cModuleSync* sync;
		cv::Size sensor_dim;
		vector<CameraLibrary::Camera*> optitrack_cameras;
		vector<CameraLibrary::Frame*> optitrack_frames;
		
		int init_video_mode;		// {0 = Precision, 1 = Segment, 2 = Object, 3 = MJPEG Mode}
		int init_threshold; 
		int init_exposure;
		int init_intensity; 
		int init_framerate_mode;    //Mode number: {2 = 100%, 1 = 50%, 0 = 25%} for V120:SLIM cameras
		int frame_count;
	};

#endif /*USE_NP_CAMERASDK*/
	
    /*
     * This class provides the multi threaded streaming pipeline structure of the particle tracking system
     */
	class StreamingPipeline {
	public:
		typedef std::shared_ptr<lpt::StreamingPipeline> Ptr;
		static StreamingPipeline::Ptr create() { return make_shared< lpt::StreamingPipeline >(); }

		StreamingPipeline() : queue_capacity(100) { shared_objects = std::make_shared < lpt::SharedObjects > (); }
		virtual ~StreamingPipeline() {}

		virtual bool initialize();
		virtual void initializeControlWindow();		
		virtual void run();
		virtual void stop();
		
		virtual void aquireImageData();
		virtual void processImages(int index);
		virtual void solveCorrespondence();
		virtual void reconstuct3DObjects();
		virtual void trackObjects();
		virtual void runControlWindow();
		virtual void runVisualizer() { 
			this->visualizer->start();
			cout << "Visualizer thread done"<<endl;
			this->visualizer->stop();
		}
		
		inline void attachCameraSystem(lpt::CameraSystem::Ptr system) { 
			this->camera_system.reset();
			this->camera_system = system; 
			camera_system->setSharedObjects(shared_objects); 
			cout << "Camera System attached to pipeline " << endl;
		}
		inline void attachImageProcessor(lpt::ImageProcessor::Ptr processor) { 
			this->processor.reset();
			this->processor = processor;
			cout << "Image Processor attached to pipeline " << endl;
		}
		
		inline void attachDetector(lpt::Detector::Ptr detector) { 
			this->detector.reset();
			this->detector = detector;
			cout << "Object Detector attached to pipeline " << endl;
		}
		
		inline void attachMatcher(lpt::Correspondence::Ptr matcher) { 
			this->matcher.reset();
			this->matcher = matcher;
			matcher->setSharedObjects(shared_objects);
			cout << "Correspondence solver attached to pipeline " << endl;
		}
		
		inline void attachTracker(lpt::Tracker::Ptr tracker ) {
			this->tracker.reset();
			this->tracker = tracker; 
			tracker->setSharedObjects(shared_objects);
			cout << "Tracker attached to pipeline " << endl;
		}
		
		inline void attachVisualizer(lpt::Visualizer::Ptr visualizer ) {
			this->visualizer.reset();
			this->visualizer = visualizer; 
			visualizer->setSharedObjects(shared_objects);
			cout << "Visualizer attached to pipeline " << endl;
		}
		
		inline void setSharedObjects( std::shared_ptr<SharedObjects> new_shared_objects) { shared_objects = new_shared_objects; }
		inline std::shared_ptr<SharedObjects> getSharedObjects() { return shared_objects; }
		
		inline std::shared_ptr<lpt::CameraSystem> getCameraSystem() { return camera_system; }
		inline std::shared_ptr<lpt::Calibrator> getCalibrator() { return calibrator; }
					
		inline void setInputDataPath(string& path) { shared_objects->input_path = path; }
		inline string getInputDataPath() { return shared_objects->input_path; }

		inline void setOutputDataPath(string& path) { shared_objects->output_path = path; }
		inline string getOutputDataPath() { return shared_objects->output_path; }
		
		inline void setImageType(cv::Mat image_type) {image_type = image_type;}
		inline cv::Mat getImageType() {return shared_objects->image_type;}

        inline void setFrameRate(int fps) { shared_objects->frame_rate = fps; }
        inline int getFrameRate() { return shared_objects->frame_rate; }

        inline void setKalmanFilter(bool state) { shared_objects->KF_isOn = state; }
        inline bool getKalmanFilter() { return shared_objects->KF_isOn; }
		
		inline void setQueueCapacity(int capacity) {queue_capacity = capacity; }
		
		inline void setCamerasStatus(bool state) { cameras_status = state; }
		inline bool getImageViewStatus() const { return image_view_status; }
		inline void setImageViewStatus(bool state) { image_view_status = state; }
		inline bool showCompositeView() const { return composite_view_status; }
		inline void setCompositeView(bool state) { composite_view_status = state; }
		inline bool showDetectionView() const { return detection_view_status; }
		inline void setDetectionView(bool state) { detection_view_status = state; }
		inline bool showReprojectionView() const { return reprojection_view_status; }
		inline void setReprojectionView(bool state) { reprojection_view_status = state; }
		inline bool showTrajectoryView() const { return trajectory_view_status; }
		inline void setTrajectoryView(bool state) { trajectory_view_status = state; }
	 
		friend void callbackRecordVideo( int state, void* data );
		friend void callbackTakeSnapshot( int state, void* data );
		friend void callbackSetImageViewStatus( int state, void* data ); 
		friend void callbackSetCompositeView( int state, void* data ); 
		friend void callbackStopCameras( int state, void* data );
		friend void callbackFlushFrameQueue(int state, void* data);
		friend void callbackFlushProcessedQueue(int state, void* data);
		friend void callbackSetDetectionView( int state, void* data );
		friend void callbackSetReprojectionView( int state, void* data );
		friend void callbackSetTrajectoryView( int state, void* data );

		void load_Rotation_Matrix();

	protected:

		std::shared_ptr < lpt::SharedObjects >	shared_objects;
		std::shared_ptr < lpt::CameraSystem >	camera_system;
		std::shared_ptr < lpt::Visualizer >		visualizer;
		std::shared_ptr < lpt::Calibrator >		calibrator;
		std::shared_ptr < lpt::Recorder >		recorder;
		std::shared_ptr < lpt::ImageProcessor >	processor;
		std::shared_ptr < lpt::Detector >		detector;
		std::shared_ptr < lpt::Correspondence >	matcher;
		std::shared_ptr < lpt::Reconstruct3D >	reconstructor;
		std::shared_ptr < lpt::Tracker >			tracker;

		// Concurrent Queues
		lpt::concurrent_queue < lpt::ImageFrameGroup >	frame_queue;
		lpt::concurrent_queue < lpt::ImageFrameGroup >	processed_queue;
		lpt::concurrent_queue < std::pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > > 	match_queue;
		lpt::concurrent_queue < std::shared_ptr<lpt::Frame3d > >		frame3D_queue;
		lpt::concurrent_queue < std::shared_ptr<lpt::Frame3d > >		monitor_queue;
	
		boost::thread	imagegrabber_thread;
		boost::thread	matcher_thread;
		boost::thread	reconstructor_thread;
		boost::thread	tracker_thread;
		boost::thread	visualizer_thread;
		boost::thread_group		imageproc_workers_thread_group;
		
		shared_ptr<boost::barrier> imageproc_barrier;
		boost::mutex imageproc_mutex;

		lpt::ImageFrameGroup imageproc_frames;
		
		int camera_displayed;       // index of initial camera to be displayed in opencv window
		int view_point;
		int queue_capacity;
		
		bool cameras_status;
		bool composite_view_status;
		bool image_view_status;
		bool detection_view_status;
		bool reprojection_view_status;
		bool trajectory_view_status;
		bool run_calibration;
	};

	class Process {
	public:
		Process() {}
		virtual void initialize(){}
		virtual void addControls(){}
		virtual void run(){}
		virtual void shutdown(){}
		virtual void setInputQueue(){}
		virtual void setOuputQueue(){}
		virtual ~Process(){}
	protected:
		boost::thread process_thread;
	};

	class Process1 : public Process {
	public:
		Process1() { }
		virtual void run() { }
		virtual ~Process1() { }
	private:

	};

	class Process2 : public Process {
	public:
		Process2() { }
		virtual void run() { }
		virtual ~Process2() { }
	private:

	};

	class ProcessPipeline {
	public:
		ProcessPipeline(){}
		virtual ~ProcessPipeline(){}

		virtual void addProcess(std::shared_ptr<lpt::Process>& process) {
			vector<std::shared_ptr<lpt::Process>> new_process;
			new_process.push_back(process);
			processes.push_back(new_process);
		}
		virtual void addParallelProcess(vector<std::shared_ptr<lpt::Process> >& process_vector) {
			processes.push_back(process_vector);
		}
	private:
		vector<vector<std::shared_ptr<lpt::Process>>> processes;
	};

} /* NAMESPACE_PT */
#endif /*DATAAQUSITION_H_*/
