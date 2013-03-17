
#ifndef DATAGEN_H_
#define DATAGEN_H_
#include "core.hpp"
#include "opencv2/opencv.hpp"

namespace lpt {

using namespace std;

enum InputFormat{ PLAINTEXT, YAMLTEXT, BINARY }; //TODO: Read trajectories from yaml and binary files

class ImageCreator;
class DataSetGenerator;

void convertFrame(const vector<cv::Point2f> &image_points,
		lpt::ImageFrame &frame, vector<int>& particleIDs);
void convertCamParameters2CV(const lpt::Camera& cam, cv::Mat& camera_matrix,
		cv::Mat& dist_coeffs, cv::Mat& rotation_matrix, cv::Mat& translation_vec);
void setCameraRotation(lpt::Camera& cam, double angle[3]);
void setCameraTranslation(lpt::Camera& cam, double trans[3]);
void setCameraIntrinsics(lpt::Camera& cam,
		double focal_length, double pixel_width,
		double aspect_ratio, int image_width,
		int image_height, double dist_coeffs[4]);

void calcFundamentalMatrices(vector<lpt::CameraPair>& camera_pairs);

class ImageCreator {
public:
	ImageCreator(
		cv::Mat type = cv::Mat::zeros( cv::Size(640, 480), CV_8UC1 ), 
		double radius = 1.0, int intensity = 255,double object_size = 2.0, int object_intensity = 255, int blur_ksize = 3 ) 
		: image_type(type), radius(radius), intensity(intensity), object_size(object_size), object_intensity(object_intensity), blur_ksize(blur_ksize) {}
	void createImage(lpt::ImageFrame& frame);

	cv::Mat image_type;
	double radius;
	int intensity;
	int object_intensity;
	double object_size;
	int blur_ksize;
};
	
class DataSetGenerator {
public:
	typedef std::shared_ptr<DataSetGenerator> Ptr;
	map< int, pair< vector<int>, vector<cv::Point3f> > > frames;
	
	DataSetGenerator()
	{
		image_creator = std::make_shared<lpt::ImageCreator>();
		cout << "Dataset generator created" << endl;
	}

	inline void setDataPath(string path ) {	data_path = path; }
	inline void setImageCreator(std::shared_ptr<lpt::ImageCreator> new_image_creator ) { 
		image_creator.reset(); 
		image_creator = new_image_creator; 
	}  
	inline std::shared_ptr<lpt::ImageCreator> getImageCreator() { return image_creator; }
	inline void setSharedObjects( std::shared_ptr<SharedObjects> new_shared_objects) { shared_objects = new_shared_objects; }
	inline std::shared_ptr<SharedObjects> getSharedObjects() { return shared_objects; }

	void read3DTrajectoryFile(string filename, lpt::InputFormat format = lpt::PLAINTEXT);
	void project3DFramesTo2D();
	void writeImageFramePoints(string data_path, string basename);
	void writeCameraPairs(string filename);
	void showImages();
	/*TODO*/
	void createSpiralTrajectories(
			int number_of_frames, int number_of_particles,
			int d, double theta);
	void createOpenFOAMTrajectories(int number_of_frames, int number_of_particles);
	virtual ~DataSetGenerator();

protected:
	std::shared_ptr<lpt::SharedObjects> shared_objects;
	std::shared_ptr<lpt::ImageCreator> image_creator;
	string data_path;
};

} /* NAMESPACE_PT */
#endif /* DATAGEN_H_ */
