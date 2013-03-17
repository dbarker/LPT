
#include "core.hpp"

namespace lpt {

using namespace std;

bool ThreePointLine::find(lpt::ImageFrame& frame) {
	bool found = false;
	if (frame.particles.size() == 3) {
		lpt::ParticleImage::Ptr points[3];
		double dist[3];
		int count = 0;
		double maxdist = 0.0;
		for (int i = 0; i < 3; ++i) {
			points[i] = frame.particles[i];
			for (int j = i + 1; j < 3; ++j) {
				dist[count] = lpt::calculateDistance(frame.particles[i], frame.particles[j]);
				count++;
			}
		}

		if (dist[0] < dist[1] && dist[0] < dist[2] && dist[2] < dist[1] ) {
			frame.particles[0] = points[0];
			frame.particles[1] = points[1];
			frame.particles[2] = points[2];
		} 
		else if (dist[1] < dist[0] && dist[1] < dist[2] && dist[2] < dist[0] ) {
			frame.particles[0] = points[0];
			frame.particles[1] = points[2];
			frame.particles[2] = points[1];
		}
		else if (dist[2] < dist[1] && dist[2] < dist[0] && dist[1] < dist[0] ) {
			frame.particles[0] = points[1];
			frame.particles[1] = points[2];
			frame.particles[2] = points[0];
		}
		else if (dist[0] < dist[1] && dist[0] < dist[2] && dist[1] < dist[2] ) {
			frame.particles[0] = points[1];
			frame.particles[1] = points[0];
			frame.particles[2] = points[2];
		}
		else if ( dist[0] < dist[1] && dist[2] < dist[0] && dist[2] < dist[1] ) {
			frame.particles[0] = points[2];
			frame.particles[1] = points[1];
			frame.particles[2] = points[0];
		}
		else if (dist[1] < dist[0] && dist[0] < dist[2] && dist[1] < dist[2] ) {
			frame.particles[0] = points[2];
			frame.particles[1] = points[0];
			frame.particles[2] = points[1];
		}
		
		double mid_pt_dist = abs( (points[2]->x - points[0]->x)*(points[0]->y - points[1]->y) - 
			(points[0]->x - points[1]->x)*(points[2]->y - points[0]->y) ) / 
			sqrt( (points[2]->x - points[0]->x)*(points[2]->x - points[0]->x) + 
			(points[2]->y - points[0]->y)*(points[2]->y - points[0]->y) );

		if ( mid_pt_dist < 1.5 )
			found = true;
		else
			frame.particles.clear();
	}
	return found;
}

void ThreePointLine::draw(vector<lpt::ParticleImage::Ptr>& particles, cv::Mat& image) {
	vector<cv::Point2f> points( particles.size() );
	for (int i = 0; i < points.size(); ++i)
		points[i] = cv::Point2f(particles[i]->x, particles[i]->y);

	cv::line( image, points[0], points[2], cv::Scalar(0,255,0), 1, 8 );
	cv::circle( image, points[0], 3, cv::Scalar(0,0,255), -1 );
	cv::circle( image, points[1], 3, cv::Scalar(255,0,0), -1 );
	cv::circle( image, points[2], 3, cv::Scalar(0,255,0), -1 );			
}

/***********************************************************
 * Camera
 ***********************************************************/

void Camera::readImageList() {
	if (path.empty() )
		path = "./";
	stringstream list;
	list << path << "imagelist" << id << ".txt";
	std::ifstream fin( list.str().c_str() );
	string line;
	while ( getline(fin, line) )
		imagelist.push_back(line);
}

void Camera::readFrames(int startframe) {
	if ( imagelist.empty() )
		readImageList();
	
	if (! imagelist.empty() ) {
		frames.resize( imagelist.size() );
		for (int f = 0; f < imagelist.size(); ++f){
			frames[f].frame_index = f + startframe;
			frames[f].image = cv::imread(imagelist[f], 0); //NOTE: Reads only greyscale images
		}
	} else
		cout << "Imagelist could not be read" << endl;
}

void Camera::writeFrames(){
//TODO: Write particle x and y pixel locations during calibration
}


} /* NAMESPACE_PT */
