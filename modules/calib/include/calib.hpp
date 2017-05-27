
#ifndef CALIB_H_
#define CALIB_H_

#include "core.hpp"  //LPT core module
#include "opencv2/opencv.hpp"

#include <boost/thread.hpp>
#include <boost/random.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/variance.hpp>

namespace lpt {

using namespace std;
using namespace boost::accumulators;

typedef accumulator_set< double, features<tag::mean, tag::variance, tag::count, tag::max, tag::min > > boost_accumulator;

enum CalibBoardType{ CHESSBOARD, CIRCLESGRID, ASYMMETRIC_CIRCLESGRID };

class CalibrationBoard;
class Calibrator;
class Calibration;

void compute3DPositionUncertainties(vector<lpt::Camera>& cameras, vector<lpt::Particle3d_Ptr>& particles, vector<double>& uncertainties, array<double,12>& stats);

void convertImagePoints(const vector<ImageFrame>& frame, vector<vector<cv::Point2f> >& image_points);
void convertImagePoints(const vector<vector<cv::Point2f> >& image_points, vector<ImageFrame>& frame);
void convertFrame(const vector<cv::Point2f>& image_points, ImageFrame& frame, int frameindex);
void convertFrame(const ImageFrame& frame, vector<cv::Point2f>& image_points);
void findCommonFramesAndConvert(
		const vector<ImageFrame>& framesA,
		const vector<ImageFrame>& framesB,
		vector<vector<cv::Point2f> >& common_pointsA,
		vector<vector<cv::Point2f> >& common_pointsB);
void convertToMat(Camera& cam, cv::Mat& camera_matrix, cv::Mat& dist_coeffs);
void operator>>(cv::Mat& M, double N[3][3]);

class CalibrationBoard {
    public:
		lpt::CalibBoardType object_type;
		cv::Size board_size;
        double square_size;
        vector<cv::Point3f> object_points;
		vector<cv::Point2f> image_points;

		virtual bool find(lpt::ImageFrame& frame)=0;
		void draw(cv::Mat& image);
};

class Chessboard : public CalibrationBoard {
public:
	Chessboard(cv::Size board_size, double square_size);
	bool find(lpt::ImageFrame& frame);
	int find_chessboard_flags;
};

class CirclesGrid : public CalibrationBoard {
public:
	CirclesGrid(cv::Size board_size = cv::Size(4, 11), double square_size = 20.0);
	bool find(lpt::ImageFrame& frame);
	int find_circlesgrid_flags;
};

class Calibrator {
public:
	typedef std::shared_ptr<lpt::Calibrator> Ptr;
	static inline lpt::Calibrator::Ptr create(vector<lpt::Camera>& pt_cameras, vector<lpt::CameraPair>& pt_camera_pairs, int& cam_id) { 
		return lpt::Calibrator::Ptr(new lpt::Calibrator(pt_cameras, pt_camera_pairs, cam_id)); 
	} 

	friend void callbackSetStereoDataCollection(int state, void* data);
	friend void callbackSetIntParamDataCollection(int state, void* data);
	friend void callbackClearStereoData(int state, void* data);
	friend void callbackClearIntParamData(int state, void* data);
	friend void callbackRunStereoCalibration(int state, void* data);
	friend void callbackRunIntParamCalibration(int state, void* data);
	friend void callbackFindGlobalReference(int state, void* data);

	Calibrator(vector<lpt::Camera>& pt_cameras, vector<lpt::CameraPair>& pt_camera_pairs, int& cam_id)
		: cameras(pt_cameras), camera_pairs(pt_camera_pairs), window_name("calibration"), frame_interval(20), stereo_frame_max(200), current_camera(cam_id) 
		{ 
			stereo_data_collection = false;
			intparam_data_collection = false;
			get_global_reference = false;
			updated = false;
			object = new lpt::ThreePointLine(); 
			//board = new lpt::Chessboard(cv::Size(9, 6), 25.6);
			board = new lpt::CirclesGrid();
			calib_flags = CV_CALIB_FIX_ASPECT_RATIO | CV_CALIB_FIX_K3 | CV_CALIB_FIX_K4 | CV_CALIB_FIX_K5;
		}

	virtual ~Calibrator() {
		if (updated) {
			bool savefiles;
			cout << "Camera Parameters have been updated" << endl;
			cout << "save files? 0 or 1" << endl;
			cin >> savefiles;
			if (savefiles) {
				string camerasfile, campairsfile;
				cout << "Enter new camera file name" << endl;
				cin >> camerasfile;
				lpt::writeCamerasFile(camerasfile, cameras);
				cout << "Wrote " << cameras.size() << " cameras to file" << endl;
				cout << "Enter new camera pairs file name" << endl;
				cin >> campairsfile;
				lpt::writeCameraPairsFile(campairsfile, camera_pairs);
				cout << "Wrote " << camera_pairs.size() << " camera pairs to file" << endl;
			}
		}
	}
	
	inline void setBoard(CalibrationBoard& calib_board) {
		if(board)
			delete board;
		board = &calib_board; 
	}

	inline void setObject(lpt::Object& calib_object) {
		if(object)
			delete object;
		object = &calib_object; 
	}

	inline void setCalibViews(cv::Mat image_type) {
		calib_views.clear();
		stereo_views.clear();
		ref_views.clear();
		calib_views.resize( cameras.size() );
		stereo_views.resize( cameras.size() );
		ref_views.resize( cameras.size() );
		for( int i = 0; i < cameras.size(); ++i ) {
			calib_views[i] = cv::Mat::zeros(image_type.size(), CV_8UC3);
			stereo_views[i] = cv::Mat::zeros(image_type.size(), CV_8UC3);
			ref_views[i] = cv::Mat::zeros(image_type.size(), CV_8UC3); 
		}
	}
	
	inline void setStereoDataCollection(bool state) { stereo_data_collection = state; }
	inline bool isStereoDataCollecting() { return stereo_data_collection; }
	inline void setIntParamDataCollection(bool state) { intparam_data_collection = state; }
	inline bool isIntParamDataCollecting() { return intparam_data_collection; }
	inline void setGlobalReference(bool state) { get_global_reference = state; }
	inline bool isSettingGlobalReference(){ return get_global_reference;}
	inline int getFrameInterval() { return frame_interval; }
	void clearStereoData() { 
		stereo_data_frames.clear();
		for( int i = 0; i < cameras.size(); ++i ) 
			stereo_views[i] = cv::Scalar(0,0,0);
	}
	void clearIntParamData() { 
		cameras[current_camera].frames.clear();
		calib_views[current_camera] = cv::Scalar(0,0,0);
	}
	void addControls();
	bool findCalibObject(lpt::ImageFrameGroup& group);
	bool findCalibBoard(lpt::ImageFrame& frame);
	void findFundamentalMatrices();
	void calibrateCamera();
	bool findGlobalReference(lpt::ImageFrameGroup& group);
	array<double,2> computeReprojectionErrors(
        const vector<vector<cv::Point3f> >& object_points,
        const vector<vector<cv::Point2f> >& image_points,
        const vector<cv::Mat>& rvecs, const vector<cv::Mat>& tvecs,
        const cv::Mat& camera_matrix, const cv::Mat& dist_coeffs,
        vector<double>& per_view_errors );
	void Calibrator::storeCameraParameters(lpt::Camera& cam,
		cv::Mat& camera_matrix, cv::Mat& dist_coeffs,
		double avg_reprojection_error);
	double checkStereoCalibration( vector<cv::Point2f>& pointsA, vector<cv::Point2f>& pointsB, cv::Mat& F);
	void calcFundamentalMatrices( vector<lpt::CameraPair>& camera_pairs);
	vector<cv::Mat> calib_views;
	vector<cv::Mat> stereo_views;
	vector<cv::Mat> ref_views;
private:
	int& current_camera;
	int stereo_frame_max;
	int frame_interval;
	int calib_flags;
	bool stereo_data_collection;
	bool intparam_data_collection;
	bool get_global_reference;
	bool updated;
	string window_name;
	vector<lpt::Camera>& cameras;
	vector<lpt::CameraPair>& camera_pairs;
	lpt::Object* object;
	CalibrationBoard* board;
	vector<lpt::ImageFrameGroup> stereo_data_frames;
};



//**********************************************************************************
class Calibration {
    public:
		int calib_flags;
		int calib_flags_stereo;
		int find_chessboard_flags;
		cv::Size image_size;
        CalibrationBoard &calib_board;
        vector < vector < cv::Point3f > > object_points;

        Calibration(CalibrationBoard& board) : calib_board(board) {}
        bool findCalibrationPoints(vector<Camera>& cameras);
        bool runCalibration(vector<Camera>& camera);
        array<double, 2> computeReprojectionErrors(
                const vector<vector<cv::Point3f> >& objectPoints,
                const vector<vector<cv::Point2f> >& imagePoints,
                const vector<cv::Mat>& rvecs, const vector<cv::Mat>& tvecs,
                const cv::Mat& cameraMatrix, const cv::Mat& distCoeffs,
                vector<double>& per_view_errors);
        void storeCameraParameters(Camera& cam,
        		vector<vector<cv::Point2f> >& image_points,
        		vector<cv::Mat>& rotation_vecs, vector<cv::Mat>& translation_vecs,
        		cv::Mat& camera_matrix, cv::Mat& dist_coeffs, array<double, 2>& error_stats);
        void calcFundamentalMatrices(vector<Camera>& cameras, vector<CameraPair>& camera_pairs);
        double checkStereoCalibration(
        		vector<vector<cv::Point2f> >& image_pointsA,
        		vector<vector<cv::Point2f> >& image_pointsB,
        		cv::Mat& camera_matrixA, cv::Mat& camera_matrixB,
        		cv::Mat& dist_coeffsA, cv::Mat& dist_coeffsB, cv::Mat& F);
        void writeObjectPoints(string file_path);
        virtual ~Calibration();
};

} /* NAMESPACE_PT */

#endif /* CALIB_H_ */
