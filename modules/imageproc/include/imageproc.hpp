
#ifndef IMAGEPROC_H_
#define IMAGEPROC_H_

#include "core.hpp"  //LPT core module
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
namespace lpt {

using namespace std;

class ImageProcess;
class ImageProcessor;
class SubtractBackground;
class Threshold;
class Erode;
class GaussianBlur;

class Detector;
class FindContoursDetector;
class GoodFeaturesToTrackDetector;

void undistortPoints(lpt::Camera& camera, lpt::ImageFrame& frame);

class ImageProcess {
public:
	virtual void addControls()=0;
	virtual void process(cv::Mat& image)=0;
};

class ImageProcessor {
public:
	typedef std::shared_ptr<ImageProcessor> Ptr;
	static inline ImageProcessor::Ptr create() { return std::make_shared<lpt::ImageProcessor>(); }

	ImageProcessor() { cout << "Image processor created " << endl; }

	void processImage( cv::Mat& image ) {
		for(int i = 0; i < processes.size(); ++i)
			processes[i]->process(image);
	}
	void addControls() {
		cout << "\t--Adding Image Processor Controls to window: " << endl;
		for(int i = 0; i < processes.size(); ++i) {
			if (processes[i])
				processes[i]->addControls();
			else
				cout << "process " << i << " invalid pointer" << endl;
		}
	}
	void addProcess( ImageProcess& process) { processes.push_back(&process); cout << "Image process added " << endl; }
	void showResults( const string& window_name );
private:
	vector<ImageProcess*> processes; 
};

class SubtractBackground : public ImageProcess {
public:
	void addControls() {
		cout << "Add controls for Background Subtraction" << endl;
	}
	inline void process(cv::Mat& image) {
		image = image - background;
	}

private:
	cv::Mat background;
};

class Threshold : public ImageProcess {
public:
	Threshold(int threshold = 100, int max = 255) 
		: threshold_value(threshold), max_threshold(max) {}
	void addControls() {
		cv::createTrackbar("ThresholdImage", string(), &threshold_value, max_threshold, 0, 0 );
		cout << "\t\tAdded Threshold process to Control Window" << endl;
	}
	inline void process(cv::Mat& image) {
		cv::threshold(image, image, threshold_value, max_threshold, cv::THRESH_BINARY);
	}
private:
	int threshold_value;
	int max_threshold;
};

class Erode : public ImageProcess {
public:
	Erode(int iter = 1, int max_iter = 5) : iterations(iter),  max_iterations(max_iter) {}
	void addControls() {
		cv::createTrackbar("Erode iterations", string(), &iterations, max_iterations, 0, 0 );
	}
	inline void process(cv::Mat& image) {
		cv::erode(image, image, cv::Mat(), cv::Point(-1,-1), iterations);
	}
private:
	int iterations;
	int max_iterations;
};

class EqualizeHistogram : public ImageProcess {
public:
	EqualizeHistogram () {}
	void addControls(){}
	inline void process(cv::Mat& image){
		cv::equalizeHist(image, image);
	}
};

class Dilate : public ImageProcess {
public:
	Dilate(int iter = 1, int max_iter = 5) : iterations(iter),  max_iterations(max_iter) {}
	void addControls() {
		cv::createTrackbar("Dilate iterations", string(), &iterations, max_iterations, 0, 0 );
	}
	inline void process(cv::Mat& image) {
		cv::dilate(image, image, cv::Mat(), cv::Point(-1,-1), iterations);
	}
private:
	int iterations;
	int max_iterations;
};

class GaussianBlur : public ImageProcess {
public:
	GaussianBlur(int ksize = 5, double sigma1 = 0.0, double sigma2 = 0.0, int boarder = 4 ) 
		: kernel_size(ksize), sigma1(sigma1), sigma2(sigma2), boarder_type(boarder) {}
	void addControls() {
		//cv::createTrackbar("GaussBlur ksize", window_name, &kernel_size, 9, 0, 0 );
		cout << "\t\tAdded GaussianBlur to Control Window" << endl;
	}
	inline void process(cv::Mat& image) {
		cv::GaussianBlur( image, image, cv::Size(kernel_size, kernel_size), sigma1, sigma2, boarder_type);
	}

private:
	int kernel_size;
	double sigma1;
	double sigma2;
	int boarder_type;
};

class Detector {
public:
	typedef std::shared_ptr<Detector> Ptr;

	virtual void detectFeatures( cv::Mat& image, vector<lpt::ParticleImage::Ptr>& features )=0;
	virtual void addControls()=0;

	void drawResult(ImageFrame& frame) {
		cv::cvtColor(frame.image, frame.image, CV_GRAY2BGR);
		for (int i = 0; i < frame.particles.size(); i++) {
			cv::circle( frame.image, cv::Point( (int) (frame.particles[i]->x + 0.5), (int) (frame.particles[i]->y) + 0.5), 5, cv::Scalar(0, 255, 0), 1, 8, 2);
		}
	}
};

class FindContoursDetector : public Detector {
	
public:
	typedef std::shared_ptr<FindContoursDetector> Ptr;
	static inline FindContoursDetector::Ptr create() { return std::make_shared<lpt::FindContoursDetector>(); }

	class Parameters {
	public:
		Parameters(int mode = cv::RETR_EXTERNAL, int method = cv::CHAIN_APPROX_NONE, int min_area = 2, int max_area = 200 )
			: mode(mode), method(method), min_contour_area(min_area), max_contour_area(max_area) {}
		int mode;
		int method;
		int min_contour_area;
		int max_contour_area;
	} params;
	
	FindContoursDetector() { cout << " Find Contours detector created " << endl; }

	void detectFeatures( cv::Mat& image, vector<lpt::ParticleImage::Ptr>& features ) {
		vector<vector< cv::Point > > contours;
		cv::findContours(image, contours, params.mode, params.method);
		double area;
		for(int c = 0; c < contours.size(); ++c) {
			//area = cv::contourArea( contours[c] );
			if( contours[c].size() > (double)params.min_contour_area  && contours[c].size() < (double)params.max_contour_area) {
			//if( area > (double)params.min_contour_area  && area < (double)params.max_contour_area) {
				cv::Moments mom = cv::moments( cv::Mat(contours[c]) );
				if( mom.m00 > 0 ) {
					features.push_back( ParticleImage::create(features.size(), mom.m10 / mom.m00, mom.m01 / mom.m00) );
				}
			}
		}
	}

	void addControls() {
		cout << "\t--Adding Find Contours Detector Controls to window: " << endl;
		cv::createTrackbar("Min Area", string(), &params.min_contour_area, 500, 0, 0 );
		cv::createTrackbar("Max Area", string(), &params.max_contour_area, 1000, 0, 0 );
	}
	void drawContours(cv::Mat& result_image, vector<vector< cv::Point >> contours) {
		cv::drawContours( result_image, contours, -1, cv::Scalar(0), 1 );
	}
};

class GoodFeaturesToTrackDetector : public Detector {
public:
	typedef std::shared_ptr<GoodFeaturesToTrackDetector> Ptr;
	static inline GoodFeaturesToTrackDetector::Ptr create() { return std::make_shared<lpt::GoodFeaturesToTrackDetector>(); }

	class Parameters {
    public:
        Parameters ( int max_corners = 10000, double quality = 0.01,
                double min_distance = 5.0, cv::Mat mask = cv::Mat(), int blocksize = 6,
                bool use_harris = false, double k = 0.04 )
				: max_corners(max_corners), quality_level(quality), min_distance(min_distance),
				mask(mask), neighborhood_size(blocksize), use_harris(use_harris), k(k) {}
        int max_corners;
        double quality_level;
        double min_distance;
        cv::Mat mask;
		int neighborhood_size;
        bool use_harris;
        double k;
    } params;

	GoodFeaturesToTrackDetector() { cout << " Good Features To Track (GFTT) detector created " << endl; }
	
	void detectFeatures( cv::Mat& image, vector<lpt::ParticleImage::Ptr>& features ) {
		
		//CV_Assert(image.depth() != sizeof(uchar));
		if (! image.isContinuous())
			cout << "Mat is not continuous "<< image.size().width << " " << image.size().height << endl;
		
		vector<cv::Point2f> corners;
		cv::goodFeaturesToTrack(image, corners, params.max_corners, params.quality_level, params.min_distance, 
					 params.mask, params.neighborhood_size, params.use_harris, params.k);
		int index = 0;
		for( int particle_count = 0; particle_count < corners.size(); ++particle_count ) {
			lpt::ParticleImage::Ptr newparticle = lpt::ParticleImage::create(particle_count, corners[particle_count].x, corners[particle_count].y);
			features.push_back(newparticle);

			////image boundary check
			//int row_start = corners[particle_count].y - params.neighborhood_size / 2 > 0 ?
			//		corners[particle_count].y - params.neighborhood_size / 2 : 0;
			//int row_end = corners[particle_count].y + params.neighborhood_size / 2 < image.size().height ?
			//		corners[particle_count].y + params.neighborhood_size / 2 : image.size().height - 1;
			//int col_start = corners[particle_count].x - params.neighborhood_size / 2 > 0 ?
			//		corners[particle_count].x - params.neighborhood_size / 2 : 0;
			//int col_end = corners[particle_count].x + params.neighborhood_size / 2 < image.size().width ?
			//		corners[particle_count].x + params.neighborhood_size / 2 : image.size().width - 1;

			////find max intensity pixel location of detected corner
			//int intensity_max = 0;
			//cv::Point2i max(0, 0);
			//for (int row = row_start; row <= row_end; ++row) {
			//	for (int col = col_start; col <= col_end; ++col) {
			//		if (image.at<uchar>(row, col) > intensity_max) {
			//			max.x = col;
			//			max.y = row;
			//			intensity_max = image.at<uchar>(row, col);
			//		}
			//	}
			//}
			//if (intensity_max > 0) {
			//	double temp_x = (double)max.x + (log((double)image.at<uchar>(max.y, max.x - 1))
			//			- log((double)image.at<uchar>(max.y, max.x + 1)) ) /
			//			( 2.0 * log((double) image.at<uchar>(max.y , max.x - 1)) -
			//					4.0 * log((double) image.at<uchar>(max.y , max.x)) +
			//					2.0 * log((double) image.at<uchar>(max.y, max.x + 1)) );
			//	double temp_y = (double)max.y + (log((double)image.at<uchar>(max.y - 1, max.x))
			//					- log((double)image.at<uchar>(max.y + 1, max.x)) ) /
			//			(2.0 * log((double)image.at<uchar>(max.y - 1, max.x)) -
			//					4.0 * log((double) image.at<uchar>(max.y, max.x)) +
			//					2.0 * log((double) image.at<uchar>(max.y + 1, max.x)) );
			//	features.push_back( std::move(ParticleImage::create(index, temp_x, temp_y)) );
			//	++index;
			//}
		}
	}

	void addControls() { cout << "GFTT Detector Controls Added" << endl;}
};


void processImages(
		Camera& camera,
		ImageProcessor& processor,
		Detector& detector);
		
void testDetectedParticles(
		vector<ParticleImage::Ptr>& true_particles,
		vector<ParticleImage::Ptr>& detected_particles);

}/* NAMESPACE_PT */

#endif /*IMAGEPROC_H_*/

