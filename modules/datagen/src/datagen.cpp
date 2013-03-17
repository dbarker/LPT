
#include "datagen.hpp"

namespace lpt {

using namespace std;

void convertFrame(const vector<cv::Point2f>& image_points, lpt::ImageFrame& frame, vector<int>& particleIDs) {
	for (int j = 0; j < image_points.size(); ++j) {
		lpt::ParticleImage::Ptr newparticle = lpt::ParticleImage::create(particleIDs[j], image_points[j].x, image_points[j].y);
		frame.particles.push_back(newparticle);
	}
}

void convertCamParameters2CV(const lpt::Camera& cam, cv::Mat& camera_matrix, cv::Mat& dist_coeffs,
		cv::Mat& rotation_matrix, cv::Mat& translation_vec)
{
	camera_matrix =  cv::Mat::eye(3, 3, CV_64F);
	rotation_matrix = cv::Mat::zeros(3, 3, CV_64F);
	translation_vec = cv::Mat::zeros(3, 1, CV_64F);
	dist_coeffs = cv::Mat::zeros(4, 1, CV_64F);

	camera_matrix.at<double>(0,0) = cam.f[0];
	camera_matrix.at<double>(1,1) = cam.f[1];
	camera_matrix.at<double>(0,2) = cam.c[0];
	camera_matrix.at<double>(1,2) = cam.c[1];

	dist_coeffs.at<double>(0) = cam.dist_coeffs[0];
	dist_coeffs.at<double>(1) = cam.dist_coeffs[1];
	dist_coeffs.at<double>(2) = cam.dist_coeffs[2];
	dist_coeffs.at<double>(3) = cam.dist_coeffs[3];

	for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				rotation_matrix.at<double>(i,j) = cam.R[i][j];

	translation_vec.at<double>(0) = cam.T[0];
	translation_vec.at<double>(1) =	cam.T[1];
	translation_vec.at<double>(2) =	cam.T[2];
}

void setCameraRotation(lpt::Camera& cam, double angle[3]) {
	double alpha = angle[0];
	double beta = angle[1];
	double gamma = angle[2];

	double rx[] =
	        {
	            1.0,	0.0,	0.0,
	            0.0,	cos(alpha), -sin(alpha),
	            0.0,	sin(alpha), cos(alpha)
	        };
	cv::Mat Rx(3, 3, CV_64F, rx);

	double ry[] =
	        {
	            cos(beta),	0.0,	sin(beta),
	            0.0,	1.0, 	0.0,
	            -sin(beta),	0.0, cos(beta)
	        };
	cv::Mat Ry(3, 3, CV_64F, ry);

	double rz[] =
	        {
	            cos(gamma),	-sin(gamma),	0.0,
	            sin(gamma),	cos(gamma),		0.0,
	            0.0,	0.0, 	1.0
	        };
	cv::Mat Rz(3, 3, CV_64F, rz);
	cv::Mat Rmat = Rx * Ry * Rz;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			cam.R[i][j] = Rmat.at<double>(i,j);
}

void setCameraTranslation(lpt::Camera& cam, double trans[3]) {
	for (int i = 0; i < 3; ++i)
		cam.T[i] = trans[i];
}

void setCameraIntrinsics(lpt::Camera& cam,
		double focal_length, double pixel_width,
		double aspect_ratio, int image_width,
		int image_height, double dist_coeffs[4])
{
	cam.f[0] = focal_length / pixel_width;
	cam.f[1] = focal_length / (pixel_width * aspect_ratio );
	cam.c[0] = (double)image_width / 2.0 - 0.5;
	cam.c[1] = (double)image_height / 2.0 - 0.5;
	cam.dist_coeffs[0] = dist_coeffs[0];
	cam.dist_coeffs[1] = dist_coeffs[1];
	cam.dist_coeffs[2] = dist_coeffs[2];
	cam.dist_coeffs[3] = dist_coeffs[3];
}

void calcFundamentalMatrices( vector<lpt::CameraPair>& camera_pairs) {
	
	for (vector<lpt::CameraPair>::iterator pair_it = camera_pairs.begin();
			pair_it != camera_pairs.end(); ++pair_it)
	{
		int cam_a_id = pair_it->cam_A.id;
		int cam_b_id = pair_it->cam_B.id;
		cout << "Calculating F for " << cam_a_id << " and " << cam_b_id << endl;
		cv::Mat Ma, Mb;
		cv::Mat Ra, Rb, Ta, Tb;
		cv::Mat notused;
		lpt::convertCamParameters2CV(pair_it->cam_A, Ma, notused, Ra, Ta);
		lpt::convertCamParameters2CV(pair_it->cam_B, Mb, notused, Rb, Tb);

		cv::Mat R = cv::Mat::zeros(3,3,CV_64F);
		cv::Mat T = cv::Mat::zeros(3,1,CV_64F);
		cv::Mat E = cv::Mat::zeros(3,3,CV_64F);
		cv::Mat F = cv::Mat::zeros(3,3,CV_64F);
		cv::Mat Fs;

		R = Rb * Ra.t();
		T = Ta - (R.t() * Tb);

		double s[] = {
					 0.0, 	-T.at<double>(2),	T.at<double>(1),
					 T.at<double>(2), 	0.0,   -T.at<double>(0),
					-T.at<double>(1), 	T.at<double>(0), 	0.0
					};
		cv::Mat S(3, 3, CV_64F, s);

		E = R * S;
		F = Mb.inv().t() * E * Ma.inv();

		double scale = fabs( F.at<double>(2,2) ) > 0 ? 1.0 / F.at<double>(2,2) : 1.0;
		F.convertTo(Fs, CV_64F, scale);

		for ( int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				pair_it->F[i][j] = Fs.at<double>(i,j);

//		for (int c = 0; c < cameras[cam_a].frames[0].particles.size(); c++) {
//			cv::Mat pa (3, 1, CV_64F);
//			pa.at<double>(0) = cameras[cam_a].frames[0].particles[c]->x;
//			pa.at<double>(1) = cameras[cam_a].frames[0].particles[c]->y;
//			pa.at<double>(2) = 1.0;
//
//			cv::Mat pb (3, 1, CV_64F);
//			pb.at<double>(0) = cameras[cam_b].frames[0].particles[c]->x;
//			pb.at<double>(1) = cameras[cam_b].frames[0].particles[c]->y;
//			pb.at<double>(2) = 1.0;
//
//			cout << "resid ab = " << c << "  " << pa.t() * F * pb << endl;
//			cout << "resid ba = " << pb.t() * F * pa << endl;
//		}
	}
}

void ImageCreator::createImage(lpt::ImageFrame& frame) {
	frame.image = image_type.clone();
	for ( int p = 0; p < frame.particles.size(); ++p) {
		double x = frame.particles[p]->x;
		double y = frame.particles[p]->y;
		double r, I;
		
		if (this->radius > 0)
			r = this->radius;
		else
			r = frame.particles[p]->radius;

		if (this->intensity > 0)
			I = this->intensity;
		else
			I = frame.particles[p]->intensity;		
		if ( r >= 1 ) {
			cv::Point center(x, y);
			cv::circle(frame.image, center, r, cv::Scalar( I, I, I ), -1, 8 );
		} else {
			frame.image.at<uchar>(y,x) = I;
		}
	}
	cv::GaussianBlur( frame.image, frame.image, cv::Size( blur_ksize, blur_ksize ), 0, 0 );
}

void DataSetGenerator::read3DTrajectoryFile(string filename, lpt::InputFormat format) {

	switch (format) {
	case lpt::BINARY: //TODO: ADD BINARY FILE SUPPORT FOR TRAJ INPUT
		cout << "BINARY format support needs to be coded for trajectory input" <<endl; //TODO: add YAML support
		break;
	case lpt::PLAINTEXT:
		{
			lpt::Input in;
			vector<lpt::Trajectory::Ptr> trajs = in.trajinput( filename.c_str() );
			cout << "(datagen.cpp) Number of Trajectories = " << trajs.size() << endl;
			for (int i = 0; i < trajs.size(); ++i){
				for (int j = 0; j < trajs[i]->particles.size(); ++j){
					double x = trajs[i]->particles[j]->x;
					double y = trajs[i]->particles[j]->y;
					double z = trajs[i]->particles[j]->z;
					int id = trajs[i]->particles[j]->id;
					int frame_index = trajs[i]->particles[j]->frame_index;
					cv::Point3f newparticle(x,y,z);
					this->frames[frame_index].first.push_back(id);             //vector of particle ids
					this->frames[frame_index].second.push_back(newparticle);   //vector of particles (cv::Point3f)
				}
			}
		}
		break;
	case lpt::YAMLTEXT: //TODO: ADD YAML FILE SUPPORT FOR TRAJ INPUT
		cout << "YAML format support needs to be coded for trajectory input" <<endl; //TODO: add YAML support
		break;
	default:
		cout << "3D trajectory format selected is not supported" << endl;
		break;
	}
}

void DataSetGenerator::project3DFramesTo2D() {
	auto& cameras = this->shared_objects->cameras;
	for (int camera_index = 0; camera_index < cameras.size(); ++camera_index ) {
		lpt::Camera& cam = cameras[camera_index];
		cv::Mat camera_matrix;
		cv::Mat rotation_matrix;
		cv::Mat translation_vec;
		cv::Mat dist_coeffs;
		lpt::convertCamParameters2CV(cam, camera_matrix, dist_coeffs,
				rotation_matrix, translation_vec);

		cv::Mat rotation_vec;
		cv::Rodrigues(rotation_matrix, rotation_vec);
		array<double, 3> p_camera, p_world;
		
		p_camera[0] = cam.sensor_size[0] / 2;
		p_camera[0] = cam.sensor_size[1] / 2;
		p_camera[0] = 0.0;

		lpt::convertCameraCoordinatesToWorld(cam, p_camera, p_world);
		double focal_length = cam.f[0] * cam.pixel_size[0];

		map<int, pair<vector<int>, vector<cv::Point3f> > >::iterator frame_it;
		for( frame_it = frames.begin(); frame_it != frames.end(); ++frame_it ) {
			int frame_index = frame_it->first;
			
			vector<int>& point_ids = frame_it->second.first;     // this is a vector<int> particleIDs
			vector<cv::Point3f>& points3D = frame_it->second.second;  // this is a vector<cv::Point3f> frame of 3D particles
			vector<cv::Point2f> image_points;
			cv::projectPoints(cv::Mat(points3D), rotation_vec, translation_vec,
					camera_matrix, dist_coeffs, image_points);
			lpt::ImageFrame newframe(frame_index);
		
			//lpt::convertFrame(image_points, newframe, point_ids);

			for (int p = 0; p < image_points.size(); ++p) {
				double distance = 
					sqrt(
						  (points3D[p].x - p_world[0]) * (points3D[p].x - p_world[0]) 
						+ (points3D[p].y - p_world[1]) * (points3D[p].y - p_world[1]) 
						+ (points3D[p].z - p_world[2]) * (points3D[p].z - p_world[2]) 
					);
	
				double radius = focal_length * image_creator->object_size / ( abs( distance - focal_length ) ) / 2.0 / cam.pixel_size[0];
				if (radius < 1.0 && radius >= 0.5) 
					radius = 0.5;
				double intensity =  image_creator->object_intensity / ( distance * distance );
				//cout << "radius = " << radius << " intensity = " << intensity << " distance " << distance << endl;
				lpt::ParticleImage::Ptr newparticle = lpt::ParticleImage::create(point_ids[p], image_points[p].x, image_points[p].y, radius, intensity);
				newframe.particles.push_back(newparticle);
			}
			
			//!!!!!!!!!!!!!!!
			//std::random_shuffle(newframe.particles.begin(), newframe.particles.end()); 
			//!!!!!!!!!!!!!!!
			image_creator->createImage(newframe);
			cam.frames.push_back(newframe);
		}
    }
	frames.clear(); 
	frames.swap(map<int, pair<vector<int>, vector<cv::Point3f> > >());
}

void DataSetGenerator::writeImageFramePoints(string data_path, string basename) {
	ofstream fout;
	auto& cameras = this->shared_objects->cameras;
	for (int c = 0; c < cameras.size(); ++c){
		YAML::Emitter out;
		out << cameras[c].frames;
		stringstream framesfile;
		framesfile << data_path << cameras[c].id << basename;
		fout.open( framesfile.str().c_str() );
		fout << out.c_str();
		fout.close();
	}
}

void DataSetGenerator::writeCameraPairs(string filename) {
	ofstream fout;
	auto& camera_pairs = this->shared_objects->camera_pairs;

	YAML::Emitter pairs_out;
	pairs_out << camera_pairs;
	fout.open( filename.c_str() );
	fout << pairs_out.c_str();
	fout.close();
}

void DataSetGenerator::createSpiralTrajectories(
        int number_of_frames, int number_of_particles,
        int d, double theta)
{  /* TODO: port matlab code here
    vector<Frame> frames(number_of_frames);
    for (int f; f < frames.size(); ++f){
        double xo, yo, zo;

        Particle newParticle;
        newParticle.x =
        newParticle.y =
        newParticle.z =
        frames.particles

    }
   */
}

void DataSetGenerator::showImages() {
	auto& cameras = this->shared_objects->cameras;
	stringstream capturedetails;
	for (int i = 0; i < cameras.size(); i++) {
		for (int f = 0; f < cameras[i].frames.size(); f++) {
			capturedetails.str("");
			capturedetails << cameras[i].name << "  --  Frame " << f;
			cv::imshow("show", cameras[i].frames[f].image);
			cv::displayStatusBar("show", capturedetails.str(), 1);
			cv::waitKey(2);
		}
	}
}

void DataSetGenerator::createOpenFOAMTrajectories(int number_of_frames, int number_of_particles){
	//TODO: Incorporate OpenFOAM simulation to create trajectories
}


DataSetGenerator::~DataSetGenerator() {
    // TODO Auto-generated destructor stub
}

} /* NAMESPACE_PT */
