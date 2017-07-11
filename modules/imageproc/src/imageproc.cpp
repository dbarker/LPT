/**
 * @file imageproc.cpp
 * Image processing module definition
 */

#include "imageproc.hpp"

namespace lpt {

using namespace std;

void undistortPoints(const lpt::Camera& camera, lpt::ImageFrame& frame) {
	if (!frame.particles.empty()) {
		vector<cv::Point2d> image_points(frame.particles.size());
		
		for (int j = 0; j < frame.particles.size(); ++j) {
			image_points[j].x = frame.particles[j]->x;
			image_points[j].y = frame.particles[j]->y;
		}

		cv::undistortPoints(image_points, image_points, camera.getCameraMatrix(), camera.getDistCoeffs(), cv::Mat(), camera.getCameraMatrix() );
		for (int j = 0; j < frame.particles.size(); ++j) { 
			frame.particles[j]->x = image_points[j].x;
			frame.particles[j]->y = image_points[j].y;	
		}
	}
}

ImageProcess::~ImageProcess()
{
    cout << "ImageProcess destructed" << endl;
}

ImageProcessor::ImageProcessor()
{
    cout << "ImageProcessor constructed" << endl;
}

void ImageProcessor::processImage(cv::Mat &image)
{
    for(int i = 0; i < m_processes.size(); ++i) {
        m_processes[i]->process(image);
    }
}

void ImageProcessor::addControls()
{
    cout << "\t--Adding Image Processor Controls to window: " << endl;
    for(int i = 0; i < m_processes.size(); ++i) {
        if (m_processes[i])
            m_processes[i]->addControls();
        else
            cout << "process " << i << " invalid pointer" << endl;
    }
}

void ImageProcessor::addProcess(ImageProcess::Ptr process)
{
    m_processes.push_back(process);
    cout << "Image process added " << endl;
}

SubtractBackground::SubtractBackground()
{
    cout << "SubtractBackground created" << endl;
}

void SubtractBackground::addControls()
{
    cout << "Add controls for background subtraction" << endl;
}

SubtractBackground::~SubtractBackground()
{
    cout << "SubtractBackgournd destructed" << endl;
}

Threshold::Threshold(int threshold, int max_threshold)
  : m_threshold(threshold), m_max_threshold(max_threshold)
{
    cout << "Threshold constructed" << endl;
}

void Threshold::addControls()
{
    cv::createTrackbar("ThresholdImage", string(), &m_threshold, m_max_threshold, 0, 0 );
    cout << "\t\tAdded Threshold process to Control Window" << endl;
}

Threshold::~Threshold()
{
    cout << "Threshold destructed" << endl;
}

Erode::Erode(int iterations, int max_iterations)
  : m_iterations(iterations), m_max_iterations(max_iterations)
{
    cout << "Erode constructed" << endl;
}

void Erode::addControls()
{
    cv::createTrackbar("Erode iterations", string(), &m_iterations, m_max_iterations, 0, 0 );
}

Erode::~Erode()
{
    cout << "Erode desturcted" << endl;
}

EqualizeHistogram::EqualizeHistogram()
{
    cout << "EqualizeHistogram constructed" << endl;
}

void EqualizeHistogram::addControls()
{
    cout << "Add controls for EqualizeHistogram" << endl;
}

EqualizeHistogram::~EqualizeHistogram()
{
    cout << "EqualizeHistogram destructed" << endl;
}

Dilate::Dilate(int iterations, int max_iterations)
  : m_iterations(iterations), m_max_iterations(max_iterations)
{
    cout << "Dilate constructed" << endl;
}

void Dilate::addControls()
{
    cv::createTrackbar("Dilate iterations", string(), &m_iterations, m_max_iterations, 0, 0 );
}

Dilate::~Dilate()
{
    cout << "Dilate destructed" << endl;
}

GaussianBlur::GaussianBlur(int kernel_size, double sigma1, double sigma2, int boarder)
  : m_kernel_size(kernel_size), m_sigma1(sigma1), m_sigma2(sigma2), m_boarder_type(boarder)
{
    cout << "Gaussian blur constructed" << endl;
}

void GaussianBlur::addControls()
{
    //cv::createTrackbar("GaussBlur ksize", window_name, &kernel_size, 9, 0, 0 );
    cout << "\t\tAdded GaussianBlur to Control Window" << endl;
}

GaussianBlur::~GaussianBlur()
{
    cout << "Gaussian blur destructed" << endl;
}

void Detector::drawResult(ImageFrame &frame)
{
    //cv::cvtColor(frame.image, frame.image, CV_GRAY2BGR);
    for (int i = 0; i < frame.particles.size(); i++) {
		auto particle = frame.particles[i];
        cv::circle( frame.image, cv::Point( static_cast<int>(particle->x), static_cast<int>(particle->y) ), 
			static_cast<int>(particle->radius), 200, 1);
    }
}

Detector::~Detector()
{
    cout << "Detector destructed" << endl;
}

FindContoursDetector::FindContoursDetector()
{
    cout << "Find Contours detector constructed" << endl;
}

void FindContoursDetector::detectFeatures(const cv::Mat &image, vector<ParticleImage::Ptr> &features, vector<vector<cv::Point>>& contours)
{
    //vector<vector< cv::Point > > contours;
    cv::findContours(image, contours, params.mode, params.method);
    //double area;
    for(int c = 0; c < contours.size(); ++c) {
        //area = cv::contourArea( contours[c] );
        if( contours[c].size() > (double)params.min_contour_area  && contours[c].size() < (double)params.max_contour_area) {
        //if( area > (double)params.min_contour_area  && area < (double)params.max_contour_area) {
            cv::Moments mom = cv::moments( cv::Mat(contours[c]) );
            if( mom.m00 > 0 ) {
                features.push_back( ParticleImage::create(static_cast<int>(features.size()), mom.m10 / mom.m00, mom.m01 / mom.m00) );
            }
        }
    }
}

//void FindContoursDetector::detectFeatures(const cv::Mat &image, vector<ParticleImage::Ptr> &features, vector<vector<cv::Point>>& contours)
//{
//    cv::findContours(image, contours, params.mode, params.method);
//    //double area;
//    for(int c = 0; c < contours.size(); ++c) {
//        //area = cv::contourArea( contours[c] );
//        if( contours[c].size() > (double)params.min_contour_area  && contours[c].size() < (double)params.max_contour_area) {
//        //if( area > (double)params.min_contour_area  && area < (double)params.max_contour_area) {
//            cv::Moments mom = cv::moments( cv::Mat(contours[c]) );
//            if( mom.m00 > 0 ) {
//                features.push_back( ParticleImage::create(static_cast<int>(features.size()), mom.m10 / mom.m00, mom.m01 / mom.m00) );
//            }
//        }
//    }
//}

void FindContoursDetector::addControls()
{
    cout << "\t--Adding Find Contours Detector Controls to window: " << endl;
    cv::createTrackbar("Min Area", string(), &params.min_contour_area, 500, 0, 0 );
    cv::createTrackbar("Max Area", string(), &params.max_contour_area, 1000, 0, 0 );
}

void Detector::drawContours(cv::Mat &image, vector<vector<cv::Point> > contours)
{
    cv::drawContours( image, contours, -1, cv::Scalar(0, 255, 0), 1 );
}

FindContoursDetector::~FindContoursDetector()
{
    cout << "Find Contours detector destructed" << endl;
}

GoodFeaturesToTrackDetector::GoodFeaturesToTrackDetector()
{
    cout << "Good Features To Track (GFTT) detector constructed" << endl;
}

void GoodFeaturesToTrackDetector::detectFeatures(const cv::Mat &image, vector<ParticleImage::Ptr> &features, vector<vector<cv::Point>>& contours)
{
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

void GoodFeaturesToTrackDetector::addControls()
{
    cout << "GoodFeaturesToTrack detector controls added" << endl;
}

GoodFeaturesToTrackDetector::~GoodFeaturesToTrackDetector()
{
    cout << "GoodFeaturesToTrack detector destructed" << endl;
}

void processImages( Camera& camera, ImageProcessor& processor, Detector& detector )
{
	cout << "Camera "<< camera.id << ": --Processing Images" << endl;
    	
	stringstream result_window;
	result_window << camera.name << ": detected particles";
	
	processor.addControls();
	detector.addControls();
	cv::waitKey(10);

	size_t number_of_frames = camera.imagelist.size();
	camera.frames.resize( number_of_frames );
	for (int i = 0; i < number_of_frames; ++i) {
		cv::Mat temp_image = camera.frames[i].image.clone();
		processor.processImage( temp_image );
		cv::imshow("processed image", temp_image );
		detector.detectFeatures( temp_image, camera.frames[i].particles, camera.frames[i].contours);
		detector.drawResult( camera.frames[i] );
		cv::imshow(result_window.str(), camera.frames[i].image);
	}
	cv::destroyWindow( result_window.str() );
}

void testDetectedParticles(
		vector<ParticleImage::Ptr>& true_particles,
		vector<ParticleImage::Ptr>& detected_particles)
{
	double total_residual = 0;
	int total_matches = 0;
	for (int m = 0; m < true_particles.size(); m++) {
		int number_of_matches = 0;
		vector<double> residuals(true_particles.size(), 0);
		for(int n = 0; n < detected_particles.size(); n++) {
			double diff_x = true_particles[m]->x - detected_particles[n]->x;
			double diff_y = true_particles[m]->y - detected_particles[n]->y;
			double r = sqrt( diff_x * diff_x + diff_y * diff_y );
            if(r <= true_particles[m]->radius) {   //FIXME: make this a function of actual particle radius
				residuals[m] += r;
				number_of_matches++;
				break;                                    //FIXME: consider if more than one particle matches
			}
		}
		total_matches += number_of_matches;
		total_residual += residuals[m];
	}
	//TODO: print the results to a file and plot residuals (histogram)
	cout << "Correct Ratio: " << (double)total_matches / detected_particles.size();
	cout << "\tCover Ratio: " << (double)total_matches / true_particles.size();
	cout << "\tAvg Frame Residual: "
			<< (total_matches > 0 ? total_residual / total_matches : 0) << endl;
}

} /*NAMESPACE_PT_*/

