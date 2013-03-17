
#include "imageproc.hpp"

namespace lpt {

using namespace std;

void undistortPoints(lpt::Camera& camera, lpt::ImageFrame& frame) {
	if (frame.particles.size() > 0) {
		vector<cv::Point2f> image_points(frame.particles.size());
		
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


void processImages( Camera& camera, ImageProcessor& processor, Detector& detector )
{
	cout << "Camera "<< camera.id << ": --Processing Images" << endl;
    	
	stringstream result_window;
	result_window << camera.name << ": detected particles";
	
	processor.addControls();
	detector.addControls();
	cv::waitKey(10);

	int number_of_frames = camera.imagelist.size();
	camera.frames.resize( number_of_frames );
	for (int i = 0; i < number_of_frames; ++i) {
		cv::Mat temp_image = camera.frames[i].image.clone();
		processor.processImage( temp_image );
		cv::imshow("processed image", temp_image );
		detector.detectFeatures( temp_image, camera.frames[i].particles);
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
			if(r <= 1.2/*true_particles[m]->radius*/) {   //FIXME: make this a function of actual particle radius
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

