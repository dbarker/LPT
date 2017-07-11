#include "calib.hpp"

namespace lpt {

using namespace std;

void compute3DPositionUncertainties(vector<lpt::Camera>& cameras, vector<lpt::Particle3d_Ptr>& particles, vector<double>& uncertainties, array<double,12>& stats){ 
	
	lpt::boost_accumulator uncertainty_stats, error_stats, condA_stats;
	vector<lpt::boost_accumulator> uncertainty_components;
	boost::mt19937 rng;
	boost::normal_distribution<double> normal_dist(0, 2.5);
	boost::uniform_real<double> uniform_dist(-2.0, 2.0);
	
	stringstream filename;
	filename << "output/";
	int num_dist_coeffs = 4;
	uncertainties.resize(particles.size(), 0);

	vector<vector<cv::Point2f>> image_points_distorted( cameras.size() );
	vector<vector<cv::Point2f>> image_points( cameras.size() );
	vector<vector< vector < cv::Point2d > > > distortion_sensitivities( cameras.size() ); 

	vector<cv::Point3d> points3D( particles.size() );
	for( int i = 0; i < points3D.size(); ++i ) {
		points3D[i].x = particles[i]->X[0];
		points3D[i].y = particles[i]->X[1];
		points3D[i].z = particles[i]->X[2];
	}

	for (int camera_index = 0; camera_index < cameras.size(); ++camera_index ) {
		lpt::Camera& cam = cameras[camera_index];
		filename << cam.id << "_";
		cv::Mat camera_matrix = cam.getCameraMatrix();
		cv::Mat rotation_matrix(3,3, CV_64F, cam.R); 		
		cv::Mat translation_vec(3,1, CV_64F, cam.T);
		cv::Mat dist_coeffs = cam.getDistCoeffs();

		cv::Mat rotation_vec;
		cv::Rodrigues(rotation_matrix, rotation_vec);

		cv::projectPoints(cv::Mat(points3D), rotation_vec, translation_vec, camera_matrix, dist_coeffs, image_points_distorted[camera_index]);
		cv::undistortPoints( image_points_distorted[camera_index], image_points[camera_index], camera_matrix, dist_coeffs, cv::Mat(), camera_matrix);

		distortion_sensitivities[camera_index].resize( num_dist_coeffs ); 

		for(int dist_id = 0; dist_id < 4; ++dist_id) {
			
			distortion_sensitivities[camera_index][dist_id].resize(image_points[camera_index].size() ) ;

			cv::Mat dist_coeffs_perturbed = dist_coeffs.clone();
			dist_coeffs_perturbed.at<double>(dist_id) += cam.dist_coeffs_u[dist_id];
			
			vector<cv::Point2f> image_points_perturbed;
			
			cv::undistortPoints( image_points_distorted[camera_index], image_points_perturbed, camera_matrix, dist_coeffs_perturbed, cv::Mat(), camera_matrix);

			for ( int p = 0; p < image_points[camera_index].size(); ++p ) {
				double x_sensitivity = 0;
				double y_sensitivity = 0;

				if ( cam.dist_coeffs_u[dist_id] != 0 ) { 
					x_sensitivity = ( image_points[camera_index][p].x - image_points_perturbed[p].x ) / cam.dist_coeffs_u[dist_id]; 
					y_sensitivity = ( image_points[camera_index][p].y - image_points_perturbed[p].y ) / cam.dist_coeffs_u[dist_id]; 
				}
				
				distortion_sensitivities[camera_index][dist_id][p].x = x_sensitivity;
				distortion_sensitivities[camera_index][dist_id][p].y = y_sensitivity;
			}
			
		}
		
	}
	
	// uncomment to add noise to the particle positions
	/*for (int camera_index = 0; camera_index < cameras.size(); ++camera_index ) {
		for ( int p = 0; p < image_points[camera_index].size(); ++p ) {
			image_points[camera_index][p].x += uniform_dist(rng);
			image_points[camera_index][p].y += uniform_dist(rng);
		}
	}*/

	cout << "Cameras " << filename.str() << endl;
	
	filename << ".txt";
	ofstream fout( filename.str().c_str() );
	fout << "3D Error\t";
	fout << "cond(A)\t";
	fout << "R L2 norm\t";
	fout << "Theta\t";
	fout << "u(x')^2\t";
	fout << "u(y')^2\t";
	fout << "u(T[0])^2\t";
	fout << "u(T[1])^2\t";
	fout << "u(T[2])^2\t";
	fout << "u(rvec[0])^2\t";
	fout << "u(rvec[1])^2\t";
	fout << "u(rvec[2])^2\t";
	fout << "uc\tx\ty\tz" << endl;
	
	for( int i = 0; i < points3D.size(); ++i ) {
		cv::Mat A = cv::Mat::zeros( 2 * static_cast<int>(cameras.size()), 3, CV_64F );
		cv::Mat B = cv::Mat::zeros( 2 * static_cast<int>(cameras.size()), 1, CV_64F );
		cv::Mat X = cv::Mat::zeros( 3, 1, CV_64F );
		//cout << "\n" << i << "\t";
		for(int c = 0; c < cameras.size(); ++c) {
			cv::Point2f& P = image_points[c][i];

			double x = ( P.x - cameras[c].c[0] ) / cameras[c].f[0];
			double y = ( P.y - cameras[c].c[1] ) / cameras[c].f[1];
		
			int s = c * 2;
			int e = s + 1;
			
			A.at<double>(s, 0) = x * cameras[c].R[2][0] - cameras[c].R[0][0];   
			A.at<double>(s, 1) = x * cameras[c].R[2][1] - cameras[c].R[0][1];
			A.at<double>(s, 2) = x * cameras[c].R[2][2] - cameras[c].R[0][2];

			A.at<double>(e, 0) = y * cameras[c].R[2][0] - cameras[c].R[1][0];
			A.at<double>(e, 1) = y * cameras[c].R[2][1] - cameras[c].R[1][1];
			A.at<double>(e, 2) = y * cameras[c].R[2][2] - cameras[c].R[1][2];

			B.at<double>(s) = cameras[c].T[0] - x * cameras[c].T[2];
			B.at<double>(e) = cameras[c].T[1] - y * cameras[c].T[2];
		}

		cv::SVD svd;
		svd(A);
		svd.backSubst(B,X);
		array<double,3> X_true;

		X_true[0] = X.at<double>(0);
		X_true[1] = X.at<double>(1);
		X_true[2] = X.at<double>(2);
		
		double error_3d = 
			sqrt(
			( X_true[0] - points3D[i].x) * (X_true[0] - points3D[i].x) +
			( X_true[1] - points3D[i].y) * (X_true[1] - points3D[i].y) +
			( X_true[2] - points3D[i].z) * (X_true[2] - points3D[i].z) 
			);

		//C1: 3D reprojection error
		error_stats(error_3d);
		fout << error_3d << "\t";

		double max_svalue = 0, min_svalue = 0; 
		cv::minMaxIdx(svd.w, &min_svalue, &max_svalue);
		double cond_A = max_svalue / min_svalue;
		double normL2_A = max_svalue;
		double normL2_B = cv::norm(B);
		double normL2_X = cv::norm(X);
		
		cv::Mat AX = A * X;
		cv::Mat R = B - AX;
		double normL2_R = cv::norm(R); 
		double theta =  atan(normL2_R / cv::norm(AX) );
		//C2: cond(A)
		fout << cond_A << "\t";
		condA_stats(cond_A);
		//C3: residual norm
		fout << normL2_R << "\t";
		//C4: theta
		fout << theta * 180/3.14159 << "\t";

		for(int cam_id = 0; cam_id < cameras.size(); ++cam_id) {
			//cout << cam_id << "\t";
			double x_u2 = 0; 
			double y_u2 = 0;

			for(int dist_id = 0; dist_id < 4; ++dist_id) {
				double sensitivity = distortion_sensitivities[cam_id][dist_id][i].x;
				double uncertainty = cameras[cam_id].dist_coeffs_u[dist_id];
				x_u2 += sensitivity*sensitivity * uncertainty*uncertainty;

				sensitivity = distortion_sensitivities[cam_id][dist_id][i].y;
				uncertainty = cameras[cam_id].dist_coeffs_u[dist_id];
				y_u2 += sensitivity*sensitivity * uncertainty*uncertainty;
			
			}
		
			x_u2 += 1.0 * cameras[cam_id].centriod_loc_uncertainty * cameras[cam_id].centriod_loc_uncertainty; // 0.25 = 0.5 * 0.5 = squared uncertainty of centriod pixel coordinate
			y_u2 += 1.0 * cameras[cam_id].centriod_loc_uncertainty * cameras[cam_id].centriod_loc_uncertainty; // 0.25 = 0.5 * 0.5 = squared uncertainty of centriod pixel coordinate
			
			//cout << "Image point uncertainty from distortion parameters and centroid localization " << sqrt(x_u2) << " " << sqrt(y_u2) << endl;

			cv::Point2f& P = image_points[cam_id][i];
			double x = ( P.x - cameras[cam_id].c[0] ) / cameras[cam_id].f[0];
			double y = ( P.y - cameras[cam_id].c[1] ) / cameras[cam_id].f[1];

			double dximdx = 1.0 / cameras[cam_id].f[0];
			double dximdc = -1.0 / cameras[cam_id].f[0];
			double dximdf = ( cameras[cam_id].c[0] - P.x ) / ( cameras[cam_id].f[0] * cameras[cam_id].f[0] );
			
			cameras[cam_id].X_u[0] = 
				sqrt( 
						  dximdx * dximdx * x_u2
						+ dximdc * dximdc * cameras[cam_id].c_u[0] * cameras[cam_id].c_u[0]
						+ dximdf * dximdf * cameras[cam_id].f_u[0] * cameras[cam_id].f_u[0]
					);

			double dyimdy =  1.0 / cameras[cam_id].f[1];
			double dyimdc = -1.0 / cameras[cam_id].f[1];
			double dyimdf = ( cameras[cam_id].c[1] - P.y ) / ( cameras[cam_id].f[1] * cameras[cam_id].f[1] );
			
			cameras[cam_id].X_u[1] = 
				sqrt( 
						  dyimdy * dyimdy * y_u2 
						+ dyimdc * dyimdc * cameras[cam_id].c_u[1] * cameras[cam_id].c_u[1]
						+ dyimdf * dyimdf * cameras[cam_id].f_u[1] * cameras[cam_id].f_u[1]
					);
			
			// 3D position uncertainty due to x' and y' image coordinates
			if ( cameras[cam_id].X_u[0] != 0 && cameras[cam_id].X_u[1] != 0) {
				for (int dim = 0; dim < 2; ++dim) {

					cv::Mat delta_B = cv::Mat::zeros( 2 * static_cast<int>(cameras.size()), 1, CV_64F );
					cv::Mat E = cv::Mat::zeros( 2 * static_cast<int>(cameras.size()), 3, CV_64F );

					int s = (dim == 0 ? cam_id * 2 : cam_id * 2 + 1);

					E.at<double>(s, 0) = cameras[cam_id].X_u[dim] * cameras[cam_id].R[2][0];   
					E.at<double>(s, 1) = cameras[cam_id].X_u[dim] * cameras[cam_id].R[2][1];
					E.at<double>(s, 2) = cameras[cam_id].X_u[dim] * cameras[cam_id].R[2][2];
					delta_B.at<double>(s) = -1.0 * cameras[cam_id].X_u[dim] * cameras[cam_id].T[2];

					svd(E);
					double max_s = 0, min_s = 0; 
					cv::minMaxIdx(svd.w, &min_s, &max_s);
					double normL2_E = max_s;
					double normL2_deltaB = cv::norm(delta_B);
					double E_sen = (cond_A * cond_A * tan(theta) + cond_A ) * normL2_E / normL2_A * normL2_X / cameras[cam_id].X_u[dim];
					double deltaB_sen = cond_A / cos(theta) * normL2_deltaB / normL2_B * normL2_X / cameras[cam_id].X_u[dim];
					double squared_uncertainty_homogenious_coordinate =	
						( E_sen * E_sen + deltaB_sen * deltaB_sen) * cameras[cam_id].X_u[dim] * cameras[cam_id].X_u[dim];
					fout << squared_uncertainty_homogenious_coordinate << "\t";
					uncertainties[i] += squared_uncertainty_homogenious_coordinate;
				}
			}
			// 3D position uncertainty to translation vector T
			for (int t = 0; t < 3; ++t) {
				double T_u[3] = {0};
				T_u[t] = cameras[cam_id].T_u[t];

				cv::Mat delta_B = cv::Mat::zeros( 2 * static_cast<int>(cameras.size()), 1, CV_64F );

				int s = cam_id * 2;
				int e = s + 1;
			
				delta_B.at<double>(s) = T_u[0] - x * T_u[2];
				delta_B.at<double>(e) = T_u[1] - y * T_u[2];

				double normL2_deltaB = cv::norm(delta_B);
				double deltaB_sen = 0;
				if ( T_u[t] != 0)
					deltaB_sen = cond_A / cos(theta) * normL2_deltaB / normL2_B * normL2_X / cameras[cam_id].T_u[t];
				double squared_uncertainty_tvec =	deltaB_sen * deltaB_sen * T_u[t] * T_u[t];
				
				fout << squared_uncertainty_tvec << "\t";
				uncertainties[i] += squared_uncertainty_tvec;		
			}

			// 3D position uncertainty due to Rotation vector r_vec
		
			
			cv::Mat A_u = cv::Mat::zeros( 2 * static_cast<int>(cameras.size()), 3, CV_64F );
			cv::Mat B_u = cv::Mat::zeros( 2 * static_cast<int>(cameras.size()), 1, CV_64F );
			cv::Mat X_u = cv::Mat::zeros( 3, 1, CV_64F );
			//cout << "\n" << i << "\t";
			for(int c = 0; c < cameras.size(); ++c) {
				cv::Point2f& Pc = image_points[c][i];

				double xc = ( Pc.x - cameras[c].c[0] ) / cameras[c].f[0];
				double yc = ( Pc.y - cameras[c].c[1] ) / cameras[c].f[1];

				cv::Mat R_u = cv::Mat(3, 3, CV_64F, cameras[c].R).clone(); 		
				cv::Mat r_vec(3, 1, CV_64F);
				cv::Rodrigues(R_u, r_vec);

				r_vec.at<double>(0) += cameras[c].r_vec_u[0];
				r_vec.at<double>(1) += cameras[c].r_vec_u[1];
				r_vec.at<double>(2) += cameras[c].r_vec_u[2];

				cv::Rodrigues(r_vec, R_u);

				int s = c * 2;
				int e = s + 1;

				A_u.at<double>(s, 0) = xc * R_u.at<double>(2,0) - R_u.at<double>(0,0);   
				A_u.at<double>(s, 1) = xc * R_u.at<double>(2,1) - R_u.at<double>(0,1);
				A_u.at<double>(s, 2) = xc * R_u.at<double>(2,2) - R_u.at<double>(0,2);

				A_u.at<double>(e, 0) = yc * R_u.at<double>(2,0) - R_u.at<double>(1,0);
				A_u.at<double>(e, 1) = yc * R_u.at<double>(2,1) - R_u.at<double>(1,1);
				A_u.at<double>(e, 2) = yc * R_u.at<double>(2,2) - R_u.at<double>(1,2);

				B_u.at<double>(s) = cameras[c].T[0] - xc * cameras[c].T[2];
				B_u.at<double>(e) = cameras[c].T[1] - yc * cameras[c].T[2];
			}

			svd(A_u);
			svd.backSubst(B_u, X_u);
			double dX = 
				sqrt(
				(X_true[0] - X_u.at<double>(0)) * (X_true[0] - X_u.at<double>(0)) +
				(X_true[1] - X_u.at<double>(1)) * (X_true[1] - X_u.at<double>(1)) +
				(X_true[2] - X_u.at<double>(2)) * (X_true[2] - X_u.at<double>(2)) 
				);

			double squared_uncertainty_rvec = dX*dX;

			fout << squared_uncertainty_rvec << "\t";
			uncertainties[i] += squared_uncertainty_rvec;					

		/*	for (int r = 0; r < 3; ++r) {
				cv::Mat R_u = cv::Mat::zeros( 3, 3, CV_64F);

				cv::Mat r_u_vec = cv::Mat::zeros( 1, 3, CV_64F);
				r_u_vec.at<double>(r) = cameras[cam_id].r_vec_u[r];

				cv::Rodrigues(r_u_vec, R_u);

				cv::Mat E = cv::Mat::zeros( 2 * cameras.size(), 3, CV_64F );

				int s = cam_id * 2;
				int e = s + 1;
				
				E.at<double>(s, 0) = x * R_u.at<double>(2, 0) -  R_u.at<double>(0, 0);   
				E.at<double>(s, 1) = x * R_u.at<double>(2, 1) -  R_u.at<double>(0, 1);
				E.at<double>(s, 2) = x * R_u.at<double>(2, 2) -  R_u.at<double>(0, 2);

				E.at<double>(e, 0) = y * R_u.at<double>(2, 0) -  R_u.at<double>(1, 0);   
				E.at<double>(e, 1) = y * R_u.at<double>(2, 1) -  R_u.at<double>(1, 1);
				E.at<double>(e, 2) = y * R_u.at<double>(2, 2) -  R_u.at<double>(1, 2);

				svd(E);
				double max_s = 0, min_s = 0; 
				cv::minMaxIdx(svd.w, &min_s, &max_s);
				double normL2_E = max_s;
				double E_sen =  (cond_A * cond_A * tan(theta) + cond_A ) * normL2_E / normL2_A * normL2_X / cameras[cam_id].r_vec_u[r];
				double squared_uncertainty_rvec = E_sen * E_sen * cameras[cam_id].r_vec_u[r] * cameras[cam_id].r_vec_u[r];

				fout << squared_uncertainty_rvec << "\t" << cameras[cam_id].r_vec_u[r] << "\t";
				uncertainties[i] += squared_uncertainty_rvec;					
			}*/
		}
		uncertainties[i] = (uncertainties[i] > 0 ? sqrt(uncertainties[i]) : 0 );
		uncertainty_stats( uncertainties[i] );
		fout << uncertainties[i] << "\t" << X_true[0] << "\t" << X_true[1] << "\t" << X_true[2] << "\t" << particles[i]->id << endl;
	}
	double mean = extract_result<tag::mean>(uncertainty_stats);
	double stdev = extract_result<tag::variance>(uncertainty_stats);
	stdev = (stdev > 0 ? sqrt(stdev) : 0 ); 

	double min = extract_result<tag::min>(uncertainty_stats);
	double max = extract_result<tag::max>(uncertainty_stats);
	
	cout << "\t3D pos std uncertainty = " << mean << " +- " << stdev << endl;
	stats[0] = mean;
	stats[1] = stdev;
	stats[2] = min;
	stats[3] = max;

	mean = extract_result<tag::mean>(error_stats);
	stdev = extract_result<tag::variance>(error_stats);
	stdev = (stdev > 0 ? sqrt(stdev) : 0 ); 
	min = extract_result<tag::min>(error_stats);
	max = extract_result<tag::max>(error_stats);
	cout << "\t3D pos error = " << mean << " +- " << stdev << endl;
	stats[4] = mean;
	stats[5] = stdev;
	stats[6] = min;
	stats[7] = max;

	mean = extract_result<tag::mean>(condA_stats);
	stdev = extract_result<tag::variance>(condA_stats);
	stdev = (stdev > 0 ? sqrt(stdev) : 0 ); 
	min = extract_result<tag::min>(condA_stats);
	max = extract_result<tag::max>(condA_stats);
	cout << "\tcond(A) = " << mean << " +- " << stdev << endl;
	stats[8] = mean;
	stats[9] = stdev;
	stats[10] = min;
	stats[11] = max;

	fout.close();
}

Chessboard::Chessboard(cv::Size board_size, double square_size) {
	this->board_size = board_size;
	this->square_size = square_size;
	this->object_type = lpt::CHESSBOARD;
	for( int i = 0; i < board_size.height; ++i ) {
		for( int j = 0; j < board_size.width; ++j ) {
			object_points.push_back(cv::Point3f(float(j * square_size),
				float(i * square_size), 0));
		}
	}
	find_chessboard_flags = CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK | CV_CALIB_CB_NORMALIZE_IMAGE;
}

bool Chessboard::find(lpt::ImageFrame& frame) {
	frame.particles.clear();
	//vector<cv::Point2f> image_points;
	image_points.clear();
	bool found = cv::findChessboardCorners( frame.image, board_size, image_points,
							find_chessboard_flags);
	if (found) {
		cv::cornerSubPix(frame.image, image_points, cv::Size(11,11), cv::Size(-1,-1),
		cv::TermCriteria( CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1 ));
		for (int j = 0; j < image_points.size(); ++j) {
			ParticleImage::Ptr newparticle = ParticleImage::create(j, image_points[j].x, image_points[j].y);
			frame.particles.push_back(newparticle);
		}
	}
	return found;
}

CirclesGrid::CirclesGrid(cv::Size board_size, double square_size) {
	this->board_size = board_size;
	this->square_size = square_size;
	this->object_type = lpt::CIRCLESGRID;
	find_circlesgrid_flags = cv::CALIB_CB_ASYMMETRIC_GRID; 
	cout << "Calib Board: height = " << board_size.height << " , width = " << board_size.width << endl; 
	for( int i = 0; i < board_size.height; i++ ) {
		for( int j = 0; j < board_size.width; j++ ) {
			object_points.push_back(cv::Point3f(double( (2 * j + i % 2) * square_size),
				double(i * square_size), 0));
		}
	}
	find_circlesgrid_flags = cv::CALIB_CB_ASYMMETRIC_GRID | cv::CALIB_CB_CLUSTERING;
}

bool CirclesGrid::find(lpt::ImageFrame& frame) {
	frame.particles.clear();
	image_points.clear();
	bool found = cv::findCirclesGrid(frame.image, board_size, image_points, find_circlesgrid_flags );
	if (found) {
		for (int j = 0; j < image_points.size(); ++j) {
			ParticleImage::Ptr newparticle = ParticleImage::create(j, image_points[j].x, image_points[j].y);
			frame.particles.push_back(newparticle);
		}
	}
	return found;
}

void CalibrationBoard::draw(cv::Mat& image) {
	cv::drawChessboardCorners( image, board_size, cv::Mat(image_points), true );
}

void Calibrator::addControls() {
	void* calibrator_void_ptr = static_cast<void*> ( this );
	cv::createButton("Get Int Data", callbackSetIntParamDataCollection, calibrator_void_ptr, CV_CHECKBOX, 0 );
	cv::createButton("Clear Int Data", callbackClearIntParamData, calibrator_void_ptr, CV_PUSH_BUTTON, 0 );
	cv::createButton("Calc Int Params", callbackRunIntParamCalibration, calibrator_void_ptr, CV_PUSH_BUTTON, 0 );

	cv::createButton("Get Stereo Data", callbackSetStereoDataCollection, calibrator_void_ptr, CV_CHECKBOX, 0 );
	cv::createButton("Clear Stereo Data", callbackClearStereoData, calibrator_void_ptr, CV_PUSH_BUTTON, 0 );
	cv::createButton("Calc F Mats", callbackRunStereoCalibration, calibrator_void_ptr, CV_PUSH_BUTTON, 0 );

	cv::createButton("Find Glob_Ref", callbackFindGlobalReference, calibrator_void_ptr, CV_PUSH_BUTTON, 0 );
	
}

bool Calibrator::findCalibBoard(lpt::ImageFrame& frame) {
	bool found = this->board->find(frame);
	if (found) {
		board->draw( calib_views[current_camera] );
		cameras[current_camera].frames.push_back(frame);
		cout << cameras[current_camera].name << " -- " << cameras[current_camera].frames.size() << endl;
		cv::imshow( cameras[current_camera].name, calib_views[current_camera] );
	}
	
	return found;
}
	
bool Calibrator::findCalibObject(lpt::ImageFrameGroup& group) {
	bool allfound = true;
	for (int i = 0; i < group.size(); ++i) {
		allfound &= object->find( group[i] );
		if ( ! allfound )
			break;
	}
	if ( allfound ) {
		for (int i = 0; i < group.size(); ++i) {
			object->draw( group[i].particles, stereo_views[i] );
			object->draw( group[i].particles, group[i].image );
			cv::imshow( cameras[i].name, stereo_views[i] );
		}

		stereo_data_frames.push_back(group);
		cout << "Stereo data frames = " << stereo_data_frames.size() << endl;
	}
	return allfound;
}

double Calibrator::checkStereoCalibration(
		vector<cv::Point2f>& pointsA,
		vector<cv::Point2f>& pointsB,
		cv::Mat& F)
{
    double err = 0;
    size_t npoints = pointsA.size();
    vector<vector<cv::Vec3f> > lines(2);
    
    cv::Mat imgpt[2];
    
	imgpt[0] = cv::Mat(pointsA).clone();
    cv::computeCorrespondEpilines(imgpt[0], 1, F, lines[0]);
    
	imgpt[1] = cv::Mat(pointsB).clone();
    cv::computeCorrespondEpilines(imgpt[1], 2, F, lines[1]);
	
	for( int j = 0; j < npoints; j++ ) {
       	err +=
     			fabs( pointsA[j].x * lines[1][j][0] +
      				  pointsA[j].y * lines[1][j][1] + lines[1][j][2]) +
       			fabs( pointsB[j].x * lines[0][j][0] +
       				  pointsB[j].y * lines[0][j][1] + lines[0][j][2] );
    }   
    return err/npoints;
}

void Calibrator::calibrateCamera() {
	vector< vector<cv::Point2f> > image_points;
	vector< vector<cv::Point3f> > object_points;
	lpt::convertImagePoints(cameras[current_camera].frames, image_points);
	object_points.resize(image_points.size(), board->object_points);

	cv::Size image_size = cameras[current_camera].frames[0].image.size();
	vector<double> reprojection_errors;
	cv::Mat camera_matrix, dist_coeffs;
	vector<cv::Mat> rotation_vecs( image_points.size() );
	vector<cv::Mat> translation_vecs( image_points.size() );
	double total_avg_err = 0;
	double rms = cv::calibrateCamera(object_points, image_points,
		image_size, camera_matrix, dist_coeffs, rotation_vecs,
		translation_vecs, calib_flags);

	cout << "RMS error reported by cv::calibrateCamera(): " << rms << endl;

	bool ok = cv::checkRange(camera_matrix) && cv::checkRange(dist_coeffs);

	array<double,2> error_stats = computeReprojectionErrors(object_points, image_points,
		rotation_vecs, translation_vecs, camera_matrix, dist_coeffs, reprojection_errors);

	cout << ( ok ? "Calibration succeeded" : "Calibration failed" ) <<  " avg reprojection error = " << error_stats[0] << " +- " << error_stats[1] << endl;

	if (ok){
		storeCameraParameters(cameras[current_camera], camera_matrix, dist_coeffs, total_avg_err );
		//double aperature_width, aperature_height;
		double fovx, fovy, focal_length, aspect_ratio;
		cv::Point2d principal_point;
		cout << "camera_matrix" << endl;
		cout << cameras[current_camera].getCameraMatrix() << endl;
		cout << "dist_coeffs" << endl;
		cout << cameras[current_camera].getDistCoeffs() << endl;
		cv::calibrationMatrixValues(camera_matrix, image_size,
			cameras[current_camera].sensor_size[0], cameras[current_camera].sensor_size[1], fovx, fovy,
			focal_length, principal_point, aspect_ratio);
		cout << "Focal Length = " << focal_length << " mm" << endl;
		cout << "Field of View = " << fovx << ", " << fovy << " degrees"<< endl;
		cout << "Principal Point (pixel x, y) = " << principal_point.x << ", " << principal_point.y << endl;
		cout << "Aspect Ratio = " << aspect_ratio << endl << endl;
		updated = true;

		vector<int> sample_ids;
		for (int i = 0; i < image_points.size(); ++i)
			sample_ids.push_back(i);
		
		int number_of_samples = 25;
		array<lpt::boost_accumulator, 2> f_samples;
		array<lpt::boost_accumulator, 2> c_samples;
		array<lpt::boost_accumulator, 4> dist_samples; //FIXME allow variable number of dist coeffs!!
		lpt::boost_accumulator rms_samples;
		cv::Mat optimal_camera_matrix = camera_matrix.clone();
		for (int i = 0; i < number_of_samples; ++i) {
			cout << "\tcalculating parameters for sample " << i << endl;
			std::random_shuffle( sample_ids.begin(), sample_ids.end() );
			vector<vector<cv::Point2f>> sample_2d_points;
			vector<vector<cv::Point3f>> sample_3d_points;
			
			for (int j = 0; j < image_points.size() * 3.0 / 4.0; ++j) {
				sample_2d_points.push_back(image_points[ sample_ids[j] ] );
				sample_3d_points.push_back(object_points[ sample_ids[j] ] );
			}
			
			camera_matrix = optimal_camera_matrix.clone();

			double rms = cv::calibrateCamera(sample_3d_points, sample_2d_points,
							image_size, camera_matrix, dist_coeffs, rotation_vecs,
							translation_vecs, calib_flags | CV_CALIB_USE_INTRINSIC_GUESS | CV_CALIB_FIX_PRINCIPAL_POINT);

			bool okay = cv::checkRange(camera_matrix) && cv::checkRange(dist_coeffs);
			if (okay) {
				f_samples[0]( camera_matrix.at<double>(0,0) );
				f_samples[1]( camera_matrix.at<double>(1,1) );
				c_samples[0]( camera_matrix.at<double>(0,2) );
				c_samples[1]( camera_matrix.at<double>(1,2) );

				dist_samples[0]( dist_coeffs.at<double>(0) );
				dist_samples[1]( dist_coeffs.at<double>(1) );
				dist_samples[2]( dist_coeffs.at<double>(2) );
				dist_samples[3]( dist_coeffs.at<double>(3) );

				rms_samples(rms);
			}
		}

		double f_means[2] = {0};
		f_means[0] = extract_result<tag::mean>( f_samples[0] );
		f_means[1] = extract_result<tag::mean>( f_samples[1] );

		cameras[current_camera].f_u[0] = sqrt( extract_result<tag::variance>( f_samples[0] ) );
		cameras[current_camera].f_u[1] = sqrt( extract_result<tag::variance>( f_samples[1] ) );
		
		double c_means[2] = {0};
		c_means[0] = extract_result<tag::mean>( c_samples[0] );
		c_means[1] = extract_result<tag::mean>( c_samples[1] );
		cameras[current_camera].c_u[0] = sqrt( extract_result<tag::variance>( c_samples[0] ) );
		cameras[current_camera].c_u[1] = sqrt( extract_result<tag::variance>( c_samples[1] ) );
		 
		double dist_means[4] = {0};
		dist_means[0] = extract_result<tag::mean>( dist_samples[0] );
		dist_means[1] = extract_result<tag::mean>( dist_samples[1] );
		dist_means[2] = extract_result<tag::mean>( dist_samples[2] );
		dist_means[3] = extract_result<tag::mean>( dist_samples[3] );

		cameras[current_camera].dist_coeffs_u[0] = sqrt( extract_result<tag::variance>( dist_samples[0] ) );
		cameras[current_camera].dist_coeffs_u[1] = sqrt( extract_result<tag::variance>( dist_samples[1] ) );
		cameras[current_camera].dist_coeffs_u[2] = sqrt( extract_result<tag::variance>( dist_samples[2] ) );
		cameras[current_camera].dist_coeffs_u[3] = sqrt( extract_result<tag::variance>( dist_samples[3] ) );

		cout << "f optimal \t mean \t stdev " << endl;
		cout << cameras[current_camera].f[0] << "\t" << f_means[0] << "\t+- " << cameras[current_camera].f_u[0] << endl;
		cout << cameras[current_camera].f[1] << "\t" << f_means[1] << "\t+- " << cameras[current_camera].f_u[1] << endl;
		
		cout << "c optimal \t mean \t stdev " << endl;
		cout << cameras[current_camera].c[0] << "\t" << c_means[0] << "\t+- " << cameras[current_camera].c_u[0] << endl;
		cout << cameras[current_camera].c[1] << "\t" << c_means[1] << "\t+- " << cameras[current_camera].c_u[1] << endl;

		cout << "dist coeff optimal \t mean \t stdev " << endl;
		cout << cameras[current_camera].dist_coeffs[0] << "\t" << dist_means[0] << "\t" << cameras[current_camera].dist_coeffs_u[0] << endl;
		cout << cameras[current_camera].dist_coeffs[1] << "\t" << dist_means[1] << "\t" << cameras[current_camera].dist_coeffs_u[1] << endl;
		cout << cameras[current_camera].dist_coeffs[2] << "\t" << dist_means[2] << "\t" << cameras[current_camera].dist_coeffs_u[2] << endl;
		cout << cameras[current_camera].dist_coeffs[3] << "\t" << dist_means[3] << "\t" << cameras[current_camera].dist_coeffs_u[3] << endl;

		cout << "RMS error of samples " << endl;
		cout << extract_result<tag::mean>(rms_samples) << " +- " << sqrt( extract_result<tag::variance>(rms_samples) ) << endl;
	}
}

void Calibrator::calcFundamentalMatrices( vector<lpt::CameraPair>& camera_pairs) {
	
	for (vector<lpt::CameraPair>::iterator pair_it = camera_pairs.begin();
			pair_it != camera_pairs.end(); ++pair_it)
	{
		int cam_a_id = pair_it->cam_A.id;
		int cam_b_id = pair_it->cam_B.id;
		auto& camA = pair_it->cam_A;
		auto& camB = pair_it->cam_B;

		cout << "Calculating F for " << cam_a_id << " and " << cam_b_id << endl;
		cv::Mat Ma, Mb;
		cv::Mat Ra, Rb, Ta, Tb;
		cv::Mat notused;

		Ma =  camA.getCameraMatrix();
		Ra = cv::Mat(3, 3, CV_64F, camA.R);
		Ta = cv::Mat(3, 1, CV_64F, camA.T);
		
		Mb =  camB.getCameraMatrix();
		Rb = cv::Mat(3, 3, CV_64F, camB.R);
		Tb = cv::Mat(3, 1, CV_64F, camB.T);
		
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
	}
}

bool Calibrator::findGlobalReference(lpt::ImageFrameGroup& group) {
	bool allfound = true;
	vector< vector<cv::Point2f > > image_points( group.size() );
	for (int c = 0; c < group.size(); ++c) {
		allfound &= board->find(group[c]);
		if ( allfound )
			image_points[c] = board->image_points;
		else
			break;
	}
	
	if (allfound) {
		//ofstream fout("matrix.txt");
		for (int c = 0; c < group.size(); ++c) {
			ref_views[c] = group[c].image.clone();
			board->image_points = image_points[c];
			board->draw( ref_views[c] );
			cv::imshow(cameras[c].name, ref_views[c] ); 
			cv::Mat r_vec, t_vec;
			cv::Mat camera_matrix = cameras[c].getCameraMatrix();
			cv::Mat dist_coeffs = cameras[c].getDistCoeffs();

			cv::solvePnP(board->object_points, board->image_points, camera_matrix, dist_coeffs, r_vec, t_vec, false);
			cv::Mat R = cv::Mat::zeros(3, 3, CV_64F); 
			cv::Rodrigues(r_vec, R);

			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					cameras[c].R[i][j] = R.at<double>(i,j);
				}
				cameras[c].T[i] = t_vec.at<double>(i);
			}

			vector<cv::Point2f> image_points2;
			
			cv::Mat jacobian;
			cv::projectPoints(cv::Mat(board->object_points), r_vec, t_vec, camera_matrix, dist_coeffs, image_points2, jacobian);
			
			vector<double> point_residuals( board->object_points.size() * 2 );
			auto resid_iter = point_residuals.begin();
			for (int i = 0; i < board->object_points.size(); ++i) {
				*resid_iter = abs(board->image_points[i].x - image_points2[i].x);
				++resid_iter;
				*resid_iter = abs(board->image_points[i].y - image_points2[i].y);
				++resid_iter;
			}

			vector<double> point_residual_mean;
			vector<double> point_residual_stdev;
			cv::meanStdDev(point_residuals, point_residual_mean, point_residual_stdev);
							
			cout << cameras[c].name << endl;
			cout << "Average pixel reprojection residual = " << point_residual_mean[0] << " +- " << point_residual_stdev[0] << endl;
			
			array<lpt::boost_accumulator, 3> rvec_samples;
			array<lpt::boost_accumulator, 3> tvec_samples;

			vector<int> point_ids;
			for (int ii = 0; ii < board->image_points.size(); ++ii) 
				point_ids.push_back(ii);
			
			int number_of_points = 8;
			for(int u = 0; u < 100; ++u) {
				vector<cv::Point2f> points(number_of_points);
				vector<cv::Point3f> object_points(number_of_points);

				std::random_shuffle( point_ids.begin(), point_ids.end() );

				for (int n = 0; n < number_of_points; ++n) {
					object_points[n] = board->object_points[point_ids[n]]; 
					points[n] = board->image_points[point_ids[n]];
				}
				
				cv::solvePnP(object_points, points, camera_matrix, dist_coeffs, r_vec, t_vec, false);

				rvec_samples[0]( r_vec.at<double>(0) );
				rvec_samples[1]( r_vec.at<double>(1) );
				rvec_samples[2]( r_vec.at<double>(2) );

				tvec_samples[0]( t_vec.at<double>(0) );
				tvec_samples[1]( t_vec.at<double>(1) );
				tvec_samples[2]( t_vec.at<double>(2) );
			}

			double rvec_means[3] = {0};
			rvec_means[0] = extract_result<tag::mean>( rvec_samples[0] );
			rvec_means[1] = extract_result<tag::mean>( rvec_samples[1] );
			rvec_means[2] = extract_result<tag::mean>( rvec_samples[2] );

			double rvec_stdevs[3] = {0};
			rvec_stdevs[0] = sqrt( extract_result<tag::variance>( rvec_samples[0] ) );
			rvec_stdevs[1] = sqrt( extract_result<tag::variance>( rvec_samples[1] ) );
			rvec_stdevs[2] = sqrt( extract_result<tag::variance>( rvec_samples[2] ) );
			
			cameras[c].r_vec_u[0] = rvec_stdevs[0];
			cameras[c].r_vec_u[1] = rvec_stdevs[1];
			cameras[c].r_vec_u[2] = rvec_stdevs[2];
		
			double tvec_means[3] = {0};
			tvec_means[0] = extract_result<tag::mean>( tvec_samples[0] );
			tvec_means[1] = extract_result<tag::mean>( tvec_samples[1] );
			tvec_means[2] = extract_result<tag::mean>( tvec_samples[2] );

			double tvec_stdevs[3] = {0};
			tvec_stdevs[0] = sqrt( extract_result<tag::variance>( tvec_samples[0] ) );
			tvec_stdevs[1] = sqrt( extract_result<tag::variance>( tvec_samples[1] ) );
			tvec_stdevs[2] = sqrt( extract_result<tag::variance>( tvec_samples[2] ) );
		
			cameras[c].T_u[0] = tvec_stdevs[0];
			cameras[c].T_u[1] = tvec_stdevs[1];
			cameras[c].T_u[2] = tvec_stdevs[2];
			
			cout << "R vec optimal \t mean \t stdev " << endl;
			cout << r_vec.at<double>(0) << "\t" << rvec_means[0] << "\t+- " << rvec_stdevs[0] << endl;
			cout << r_vec.at<double>(1) << "\t" << rvec_means[1] << "\t+- " << rvec_stdevs[1] << endl;
			cout << r_vec.at<double>(2) << "\t" << rvec_means[2] << "\t+- " << rvec_stdevs[2] << endl;
	
			cout << "\n\nT vec optimal \t mean \t stdev " << endl;
			cout << cameras[c].T[0] << "\t" << tvec_means[0] << "+- " << tvec_stdevs[0] << endl;
			cout << cameras[c].T[1] << "\t" << tvec_means[1] << "+- " << tvec_stdevs[1] << endl;
			cout << cameras[c].T[2] << "\t" << tvec_means[2] << "+- " << tvec_stdevs[2] << endl;

			//cv::Mat dxdrot = jacobian.colRange( 0, 3 );
			//cv::Mat dxdt = jacobian.colRange( 3, 6 );
			//cv::Mat dxdf = jacobian.colRange( 6, 8 );
			//cv::Mat dxdc = jacobian.colRange( 8, 10 );
			//cv::Mat dxddist = jacobian.colRange( 10, 14 ); //FIXME make dist coef number variable

			//cv::Mat jacobian_matrix( jacobian.size().width, jacobian.size().width, CV_64F ); 
			//
			//cv::Mat A = cv::Mat::zeros(jacobian.size().height, 8, CV_64F );
			//cv::Mat B = cv::Mat::zeros(jacobian.size().height, 6, CV_64F );
			//cv::Mat roi;

			//roi = cv::Mat(A, cv::Rect( cv::Point(0,0), dxdf.size() ) ); 
			//dxdf.copyTo(roi);

			//roi = cv::Mat(A, cv::Rect( cv::Point(2,0), dxdc.size() ) );
			//dxdc.copyTo(roi);

			//roi = cv::Mat(A, cv::Rect( cv::Point(4,0), dxddist.size() ) );
			//dxddist.copyTo(roi);

			//roi = cv::Mat(B, cv::Rect( cv::Point(0,0), dxdrot.size() ) );
			//dxdrot.copyTo(roi);

			//roi = cv::Mat(B, cv::Rect( cv::Point(3,0), dxdt.size() ) );
			//dxdt.copyTo(roi);

			//A = A.t();
			//B = B.t();

			//cv::Mat AA_t = A * A.t();
			//cv::Mat BB_t = B * B.t();

			//roi = cv::Mat( jacobian_matrix, cv::Rect( cv::Point(0 , 0), AA_t.size() ) );
			//AA_t.copyTo(roi);

			//roi = cv::Mat( jacobian_matrix, cv::Rect( cv::Point(8 , 8), BB_t.size() ) );
			//BB_t.copyTo(roi);

			//cv::Mat AB = A * B.t();

			//roi = cv::Mat( jacobian_matrix, cv::Rect( cv::Point(8 , 0), AB.size() ) );
			//AB.copyTo(roi);

			//roi = cv::Mat( jacobian_matrix, cv::Rect( cv::Point(0 , 8), AB.t().size() ) );
			//AB = AB.t();
			//AB.copyTo(roi);

			//cv::Mat jacobian_matrix_inverse;
			//cv::invert( jacobian_matrix, jacobian_matrix_inverse);
			//
			//cv::Mat sigmas;
			//cv::sqrt(jacobian_matrix_inverse.diag(0), sigmas);
			//
			//sigmas = 3 * sigmas * point_residual_stdev[0];
			//
			//cout << "Uncertainties = \n" << sigmas << endl << endl;

			//
			//fout <<"Jacobian \n" << jacobian << endl << endl;
			//fout <<"A \n " << A << endl << endl;
			//fout <<"B \n " << B << endl << endl;
			//fout <<"Jacobian matrix \n" << jacobian_matrix << endl << endl;
			//fout <<"Jacobian matrix inv \n" << jacobian_matrix_inverse << endl;
			//fout <<"Pixel residual mean = " << point_residual_mean[0] << " +- " << point_residual_stdev[0] << endl;
			//fout <<"Uncertainties \n" << sigmas << endl;

		
			////cout << inverse << endl;
			//
			////cout << sigmas << endl << endl;

			//vector<double> sum_squared_resid( jacobian.size().width, 0 );
			//vector<double> avg_jacobian( jacobian.size().width, 0);
			//for (int i = 0; i < image_points.size(); ++i) {
			//	cv::Point2f point_residual;
			//	point_residual.x = board->image_points[i].x - image_points2[i].x;
			//	point_residual.y = board->image_points[i].y - image_points2[i].y;
			//	
			//	for(int j = 0; j < jacobian.size().width; ++j) {
			//		if (jacobian.type() == CV_64F) {
			//			double sensitivity_x = jacobian.at<double>(i*2, j);
			//			double param_resid_x = ( sensitivity_x != 0 ? point_residual.x / sensitivity_x : 0 );
			//			sum_squared_resid[j] +=  param_resid_x * param_resid_x;
			//			
			//			double sensitivity_y = jacobian.at<double>(i*2 + 1, j);
			//			double param_resid_y = ( sensitivity_y != 0 ? point_residual.y / sensitivity_y : 0 );
			//			sum_squared_resid[j] +=  param_resid_y * param_resid_y;

			//			avg_jacobian[j] += sensitivity_x + sensitivity_y;
			//		}
			//		/*else if (jacobian.type() == CV_32F) {
			//			double sensitivity_x = jacobian.at<float>(i*2, j);
			//			double param_resid_x = ( sensitivity_x != 0 ? point_residual.x / sensitivity_x : 0 );
			//			sum_squared_resid[j] +=  param_resid_x * param_resid_x;
			//			
			//			double sensitivity_y = jacobian.at<float>(i*2 + 1, j);
			//			double param_resid_y = ( sensitivity_y != 0 ? point_residual.y / sensitivity_y : 0 );
			//			sum_squared_resid[j] +=  param_resid_y * param_resid_y;
			//		}*/
			//	}	
			//}

			//for(int j = 0; j < jacobian.size().width; ++j) {
			//	avg_jacobian[j] /= 2 * image_points.size();
			//	cout << avg_jacobian[j] << " ";
			//}
			//cout << endl;
			//cv::Mat r_vec_uncertainty = cv::Mat::zeros(3,1, CV_64F);
			//r_vec_uncertainty.at<double>(0,0) = sqrt( sum_squared_resid[0] / ( 2 * image_points2.size() ) );
			//r_vec_uncertainty.at<double>(1,0) = sqrt( sum_squared_resid[1] / ( 2 * image_points2.size() ) );
			//r_vec_uncertainty.at<double>(2,0) = sqrt( sum_squared_resid[2] / ( 2 * image_points2.size() ) );

			//cv::Mat R_u = cv::Mat::zeros(3, 3, CV_64F); 
			//cv::Rodrigues(r_vec_uncertainty, R_u);
			//			
			//cout << "R relative uncertainty = " << endl;
			//for (int i = 0; i < 3; ++i) {
			//	for (int j = 0; j < 3; ++j) {
			//		cameras[c].R_u[i][j] = R_u.at<double>(i,j);
			//		cout << abs( cameras[c].R_u [i][j] / cameras[c].R[i][j] ) << " ";
			//	}
			//	cout << endl;
			//}

			//cameras[c].T_u[0] = sqrt( sum_squared_resid[3] / ( 2 * image_points2.size() ) );
			//cameras[c].T_u[1] = sqrt( sum_squared_resid[4] / ( 2 * image_points2.size() ) );
			//cameras[c].T_u[2] = sqrt( sum_squared_resid[5] / ( 2 * image_points2.size() ) );

			//cout << "T uncertainty = " << abs( cameras[c].T_u[0] / cameras[c].T[0] ) << " " << 
			//	abs( cameras[c].T_u[1] / cameras[c].T[1] ) << " " <<  abs( cameras[c].T_u[2] / cameras[c].T[2] ) << endl;
		
		}
		this->calcFundamentalMatrices(this->camera_pairs);
		for (auto pair_it = camera_pairs.begin(); pair_it != camera_pairs.end(); ++pair_it) {
			pair_it->epipolar_error = checkStereoCalibration( image_points[pair_it->cam_A.id], image_points[pair_it->cam_B.id], cv::Mat(3,3, CV_64F, pair_it->F) ); 
			cout << "Fund Mat: " << pair_it->cam_A.id << " " << pair_it->cam_B.id << " epi error = " << pair_it->epipolar_error << endl;
		}
		//fout.close();
	}
	updated = true;
	return allfound;
}

array<double,2> Calibrator::computeReprojectionErrors(
        const vector<vector<cv::Point3f> >& object_points,
        const vector<vector<cv::Point2f> >& image_points,
        const vector<cv::Mat>& rvecs, const vector<cv::Mat>& tvecs,
        const cv::Mat& camera_matrix, const cv::Mat& dist_coeffs,
        vector<double>& per_view_errors )
{
    vector<cv::Point2f> image_points2;
    int i, total_points = 0;
    
    per_view_errors.resize(object_points.size());
	lpt::boost_accumulator error_stats;
    for( i = 0; i < (int)object_points.size(); ++i )
    {
    	cv::projectPoints(cv::Mat(object_points[i]), rvecs[i], tvecs[i],
                camera_matrix, dist_coeffs, image_points2);
       
		for (int j = 0; j < object_points[i].size(); ++j) {
			double pixel_distance =
				sqrt(
				  (image_points[i][j].x - image_points2[j].x) * (image_points[i][j].x - image_points2[j].x) 
				+ (image_points[i][j].y - image_points2[j].y) * (image_points[i][j].y - image_points2[j].y) 
			    );
			error_stats(pixel_distance);
		}
    }

	array<double, 2> error;
	error[0] = extract_result<tag::mean>(error_stats);
	error[1] = extract_result<tag::variance>(error_stats);
	error[1] = error[1] > 0 ? sqrt(error[1]) : 0;
    return error;
}

void Calibrator::storeCameraParameters(lpt::Camera &cam,
		cv::Mat &camera_matrix, cv::Mat &dist_coeffs,
		double avg_reprojection_error)
{
	cam.f[0] = camera_matrix.at<double>(0,0);
	cam.f[1] = camera_matrix.at<double>(1,1);
	cam.c[0] = camera_matrix.at<double>(0,2);
	cam.c[1] = camera_matrix.at<double>(1,2);

	cam.dist_coeffs[0] = dist_coeffs.at<double>(0);
	cam.dist_coeffs[1] = dist_coeffs.at<double>(1);
	cam.dist_coeffs[2] = dist_coeffs.at<double>(2);
	cam.dist_coeffs[3] = dist_coeffs.at<double>(3);

	cam.avg_reprojection_error = avg_reprojection_error;
}

void Calibrator::findFundamentalMatrices() {
	cout << "calculating F mats " << camera_pairs.size() << endl;
	vector<lpt::CameraPair>::iterator pair_it;
	for (pair_it = camera_pairs.begin(); pair_it != camera_pairs.end(); ++pair_it) {
		vector<cv::Point2f> pointsA, pointsB;
		int a_id = pair_it->cam_A.id;
		int b_id = pair_it->cam_B.id;
		for (int f = 0; f < stereo_data_frames.size(); ++f) {
			lpt::ImageFrameGroup& group = stereo_data_frames[f]; 
			if (group[a_id].particles.size() == 3 && group[b_id].particles.size() == 3) {
				for (int p = 0; p < 3; ++p) {
					pointsA.push_back(cv::Point2d(group[a_id].particles[p]->x, group[a_id].particles[p]->y) );
					pointsB.push_back(cv::Point2d(group[b_id].particles[p]->x, group[b_id].particles[p]->y) );
				}
			}
		}

		if (pointsA.size() > 8 && pointsB.size() > 8 ) {
			cv::Mat F = cv::findFundamentalMat(pointsA, pointsB, CV_FM_LMEDS, 1.0, 1.0-1E-6);
			F >> pair_it->F;
			pair_it->epipolar_error = checkStereoCalibration( pointsA, pointsB, F ); 
		}
		cout << "Fund Mat: " << pair_it->cam_A.id << " " << pair_it->cam_B.id << " epi error = " << pair_it->epipolar_error << endl;
	}
	updated = true;
}

void callbackSetStereoDataCollection(int state, void* data) {
	lpt::Calibrator* calibrator = static_cast<lpt::Calibrator*>(data);
	calibrator->setStereoDataCollection( (state != 0) );
	cout << "started stereo data collection" << endl;
}

void callbackSetIntParamDataCollection(int state, void* data) {
	lpt::Calibrator* calibrator = static_cast<lpt::Calibrator*>(data);
	calibrator->setIntParamDataCollection( (state != 0) );
	cout << "started internal param data collection" << endl;
}

void callbackClearStereoData(int state, void* data) {
	lpt::Calibrator* calibrator = static_cast<lpt::Calibrator*>(data);
	calibrator->clearStereoData();
}

void callbackClearIntParamData(int state, void* data) {
	lpt::Calibrator* calibrator = static_cast<lpt::Calibrator*>(data);
	calibrator->clearIntParamData();
}

void callbackRunStereoCalibration(int state, void* data) {
	lpt::Calibrator* calibrator = static_cast<lpt::Calibrator*>(data);
	calibrator->findFundamentalMatrices();
}

void callbackRunIntParamCalibration(int state, void* data) {
	lpt::Calibrator* calibrator = static_cast<lpt::Calibrator*>(data);
	if ( !calibrator->cameras[calibrator->current_camera].frames.empty() )
		calibrator->calibrateCamera();
}

void callbackFindGlobalReference(int state, void* data) {
	lpt::Calibrator* calibrator = static_cast<lpt::Calibrator*>(data);
	calibrator->setGlobalReference(true);
}

void operator>>(cv::Mat& M, double N[3][3]){
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			N[i][j] = M.at<double>(i,j);
}

void convertImagePoints(const vector<ImageFrame>& frames, vector<vector<cv::Point2f> >& image_points) {
	for (int i = 0, f = 0; i < frames.size(); ++i) {
		if (frames[i].particles.size() > 0) {
			vector<cv::Point2f> empty;
			image_points.push_back( empty );
			convertFrame(frames[i], image_points[f]);
			f++;
		}
	}
}

void convertImagePoints(const vector<vector<cv::Point2f> >& image_points, vector<ImageFrame>& frames) {
	frames.resize( image_points.size() );
	for (int i = 0; i < image_points.size(); ++i) {
		convertFrame(image_points[i], frames[i], i);
	}
}

void convertFrame(const vector<cv::Point2f>& image_points, ImageFrame& frame, int frameindex) {
	frame.frame_index = frameindex;
	for (int j = 0; j < image_points.size(); ++j) {
		ParticleImage::Ptr newparticle = ParticleImage::create(j, image_points[j].x, image_points[j].y);
		frame.particles.push_back(newparticle);
	}
}

void convertFrame(const ImageFrame &frame, vector<cv::Point2f> &image_points) {
	for (int j = 0; j < frame.particles.size(); ++j) {
		cv::Point2d newparticle(frame.particles[j]->x, frame.particles[j]->y );
		image_points.push_back(newparticle);
	}
}

void findCommonFramesAndConvert(const vector<ImageFrame> &framesA, const vector<ImageFrame> &framesB,
		vector<vector<cv::Point2f> > &common_pointsA, vector<vector<cv::Point2f> > &common_pointsB)
{
	for (int i = 0; i < framesA.size(); ++i) {
		for (int j = 0; j < framesB.size(); ++j) {
			if ( framesA[i].frame_index == framesB[j].frame_index ) {
				if (framesA[i].particles.size() > 0 && framesB[j].particles.size() > 0 ){
					vector<cv::Point2f> commonA, commonB;

					convertFrame(framesA[i], commonA);
					common_pointsA.push_back(commonA);

					convertFrame(framesB[j], commonB);
					common_pointsB.push_back(commonB);
				}
			}
		}
	}
}

void convertToMat(Camera &cam, cv::Mat &camera_matrix, cv::Mat &dist_coeffs) {
	camera_matrix.at<double>(0,0) = cam.f[0];
	camera_matrix.at<double>(1,1) = cam.f[1];
	camera_matrix.at<double>(0,2) = cam.c[0];
	camera_matrix.at<double>(1,2) = cam.c[1];

	dist_coeffs.at<double>(0) = cam.dist_coeffs[0];
	dist_coeffs.at<double>(1) = cam.dist_coeffs[1];
	dist_coeffs.at<double>(2) = cam.dist_coeffs[2];
	dist_coeffs.at<double>(3) = cam.dist_coeffs[3];
}

bool Calibration::findCalibrationPoints(vector<Camera> &cameras) {
	bool ok = true;
	cv::namedWindow( "Calibration View" );
	for (int c = 0; c < cameras.size(); ++c){
		vector<vector<cv::Point2f> > image_points(cameras[c].imagelist.size());
		for (int n = 0; n < cameras[c].imagelist.size(); ++n){
			string input_filename = cameras[c].imagelist[n];
			cout << input_filename << endl;
			cv::Mat image, gray_image;
			image = cv::imread(input_filename);
			cout << image.size().height << " " << image.size().width << endl;
			if (image.size().height > 0 ){
				if (image.channels() > 1)
					cv::cvtColor(image, gray_image, CV_BGR2GRAY);
				else 
					gray_image = image.clone();

				image_size = image.size();
				
				cv::imshow("Calibration View", gray_image);
				cv::waitKey(50);
				bool found = false;

				switch (calib_board.object_type){
				case CHESSBOARD:
					found = cv::findChessboardCorners( gray_image, calib_board.board_size, image_points[n],
							find_chessboard_flags);
					if (found)
						cv::cornerSubPix(gray_image, image_points[n], cv::Size(11,11), cv::Size(-1,-1),
								cv::TermCriteria( CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1 ));
					break;
				case CIRCLESGRID:
					found = cv::findCirclesGrid(gray_image, calib_board.board_size, image_points[n] /*TODO add flags*/);
					break;
				case ASYMMETRIC_CIRCLESGRID:
					found = cv::findCirclesGrid(gray_image, calib_board.board_size, image_points[n], cv::CALIB_CB_ASYMMETRIC_GRID  );
					break;
				default:
					cerr << "calibration object type not recognized" << endl;
					exit(1);
					break;
				}
				if (found) {
					cv::drawChessboardCorners( image, calib_board.board_size, cv::Mat(image_points[n]), found );
					cv::imshow("Calibration View",image);
					cv::waitKey(100);
					cout << ":  " << image_points[n].size() << " calibration points found" << endl;
				} else {
					cout << ":  calibration point identification unsuccessful" << endl;
					image_points[n].clear();
				}
			}
		}
		convertImagePoints(image_points, cameras[c].frames);
    }
    return ok;
}

bool Calibration::runCalibration( vector<Camera> &cameras) {
	bool ok = false;
	for (int c = 0; c < cameras.size(); ++c) {
		vector< vector<cv::Point2f> > image_points;
		convertImagePoints(cameras[c].frames, image_points);
		object_points.resize(image_points.size(), calib_board.object_points);

		vector<double> reprojection_errors;
		cv::Mat camera_matrix, dist_coeffs;
		vector<cv::Mat> rotation_vecs( image_points.size() );
		vector<cv::Mat> translation_vecs( image_points.size() );
	
		double rms = cv::calibrateCamera(object_points, image_points,
				image_size, camera_matrix, dist_coeffs, rotation_vecs,
				translation_vecs, calib_flags);

		printf("RMS error reported by cv::calibrateCamera(): %g\n", rms);

		bool ok = cv::checkRange(camera_matrix) && cv::checkRange(dist_coeffs);

		array<double,2> error_stats = computeReprojectionErrors(object_points, image_points,
				rotation_vecs, translation_vecs, camera_matrix, dist_coeffs, reprojection_errors);

		cout << (ok ? "Calibration succeeded" : "Calibration failed") << " avg reprojection error = " << error_stats[0] << " +- " << error_stats[1];
		if (ok){
			storeCameraParameters(cameras[c], image_points,
					rotation_vecs, translation_vecs, camera_matrix,
					dist_coeffs, error_stats );
			//double aperature_width, aperature_height;
			double fovx, fovy, focal_length, aspect_ratio;
			cv::Point2d principal_point;
			cv::calibrationMatrixValues(camera_matrix, image_size,
					cameras[c].sensor_size[0], cameras[c].sensor_size[1], fovx, fovy,
					focal_length, principal_point, aspect_ratio);
			cout << "Focal Length = " << focal_length << " mm" << endl;
			cout << "Field of View = " << fovx << ", " << fovy << " degrees"<< endl;
			cout << "Principal Point (pixel x, y) = " << principal_point.x << ", " << principal_point.y << endl;
			cout << "Aspect Ratio = " << aspect_ratio << endl << endl;
		}
	}
    return ok;
}


/*void Calibration::TransformCoordinates(
        const cv::Mat& rvec1, const cv::Mat& tvec1,
        const cv::Mat& rvec2, const cv::Mat& tvec2,
        cv::Mat& rvec3, cv::Mat& tvec3)
{
//FIXME:  Complete code to transform all chess board corner coordinates to
    //a single coordinate system.
    composeRT(rvec1, tvec1, rvec2, tvec2, rvec3, tvec3);
}*/

array<double,2> Calibration::computeReprojectionErrors(
        const vector<vector<cv::Point3f> >& object_points,
        const vector<vector<cv::Point2f> >& image_points,
        const vector<cv::Mat>& rvecs, const vector<cv::Mat>& tvecs,
        const cv::Mat& camera_matrix, const cv::Mat& dist_coeffs,
        vector<double>& per_view_errors )
{
    vector<cv::Point2f> image_points2;
    int i, total_points = 0;
    
    per_view_errors.resize(object_points.size());
	lpt::boost_accumulator error_stats;
    for( i = 0; i < (int)object_points.size(); ++i )
    {
    	cv::projectPoints(cv::Mat(object_points[i]), rvecs[i], tvecs[i],
                camera_matrix, dist_coeffs, image_points2);
       
		for (int j = 0; j < object_points[i].size(); ++j) {
			double pixel_distance =
				sqrt(
				  (image_points[i][j].x - image_points2[j].x) * (image_points[i][j].x - image_points2[j].x) 
				+ (image_points[i][j].y - image_points2[j].y) * (image_points[i][j].y - image_points2[j].y) 
			    );
			error_stats(pixel_distance);
		}
    }

	array<double, 2> error;
	error[0] = extract_result<tag::mean>(error_stats);
	error[1] = extract_result<tag::variance>(error_stats);
	error[1] = error[1] > 0 ? sqrt(error[1]) : 0;
    return error;
}

void Calibration::calcFundamentalMatrices(vector<Camera> &cameras, vector<CameraPair> &camera_pairs){
    for (vector<CameraPair>::iterator pair_it = camera_pairs.begin();
    		pair_it != camera_pairs.end(); ++pair_it)
    {
    	
    	vector<vector<cv::Point2f> > common_pointsA, common_pointsB;
		findCommonFramesAndConvert(pair_it->cam_A.frames, pair_it->cam_B.frames,
    			common_pointsA, common_pointsB);

    	object_points.resize(common_pointsA.size(), calib_board.object_points);

    	cv::Mat camera_matrixA = cv::Mat::eye(3,3, CV_64F);
    	cv::Mat camera_matrixB = cv::Mat::eye(3,3, CV_64F);
    	cv::Mat dist_coeffsA = cv::Mat::zeros(4,1, CV_64F);  //FIXME: change this to a variable to include more than 4 dist coeffs
    	cv::Mat dist_coeffsB = cv::Mat::zeros(4,1, CV_64F);

		convertToMat(pair_it->cam_A, camera_matrixA, dist_coeffsA );
    	convertToMat(pair_it->cam_B, camera_matrixB, dist_coeffsB );

    	cv::Mat R, T, E, F;
    	printf("Analyzing Cameras %d and %d\n",pair_it->cam_A.id, pair_it->cam_B.id);

    	double rms = cv::stereoCalibrate(object_points, common_pointsA, common_pointsB,
    			camera_matrixA, dist_coeffsA, camera_matrixB, dist_coeffsB,
    			image_size, R, T, E, F, cv::TermCriteria(CV_TERMCRIT_ITER + CV_TERMCRIT_EPS, 100, 1e-5),
    			calib_flags_stereo);
    	printf("Stereo system analyzed: Calibration RMS re-projection error = %f\n", rms);

    	double epipolar_error = checkStereoCalibration(common_pointsA, common_pointsB,
    			camera_matrixA, camera_matrixB, dist_coeffsA, dist_coeffsB, F);

    	printf("Fundamental Matrix Analyzed: average deviation from epipolar lines = %f\n", epipolar_error);
    	printf("Storing fundamental Matrix for %d, %d\n\n",pair_it->cam_A.id, pair_it->cam_B.id);
    	F >> pair_it->F;
    	pair_it->reprojection_error = rms;
    	pair_it->epipolar_error = epipolar_error;

    }
}

double Calibration::checkStereoCalibration(
		vector<vector<cv::Point2f> > &image_pointsA,
		vector<vector<cv::Point2f> > &image_pointsB,
		cv::Mat &camera_matrixA, cv::Mat &camera_matrixB,
		cv::Mat &dist_coeffsA, cv::Mat &dist_coeffsB, cv::Mat &F)
{

    double err = 0;
    int npoints = 0;
    vector<vector<cv::Vec3f> > lines(2);
    for( int i = 0; i < image_pointsA.size(); i++ )
    {
        int npt = (int)image_pointsA[i].size();
        cv::Mat imgpt[2];
        imgpt[0] = cv::Mat(image_pointsA[i]).clone();
        cv::undistortPoints(imgpt[0], imgpt[0], camera_matrixA, dist_coeffsA, cv::Mat(), camera_matrixA);
        cv::computeCorrespondEpilines(imgpt[0], 1, F, lines[0]);
        imgpt[1] = cv::Mat(image_pointsB[i]).clone();
        cv::undistortPoints(imgpt[1], imgpt[1], camera_matrixB, dist_coeffsB, cv::Mat(), camera_matrixB);
        cv::computeCorrespondEpilines(imgpt[1], 2, F, lines[1]);

        for( int j = 0; j < npt; j++ )
        {
        	double errij =
        			fabs( image_pointsA[i][j].x * lines[1][j][0] +
        				  image_pointsA[i][j].y * lines[1][j][1] + lines[1][j][2]) +
        			fabs( image_pointsB[i][j].x * lines[0][j][0] +
        				  image_pointsB[i][j].y*lines[0][j][1] + lines[0][j][2] );
        	err += errij;
        }
        npoints += npt;
    }
    return err/npoints;
}

void Calibration::storeCameraParameters(Camera &cam,
		vector<vector<cv::Point2f> > &image_points,
		vector<cv::Mat> &rotation_vecs, vector<cv::Mat> &translation_vecs,
		cv::Mat &camera_matrix, cv::Mat &dist_coeffs,
		array<double,2>& error_stats)
{
	cam.f[0] = camera_matrix.at<double>(0,0);
	cam.f[1] = camera_matrix.at<double>(1,1);
	cam.c[0] = camera_matrix.at<double>(0,2);
	cam.c[1] = camera_matrix.at<double>(1,2);
	cv::Mat R;
	cv::Rodrigues(rotation_vecs[0], R);               //FIXME: This assumes that the desired coord sys is index 0
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			cam.R[i][j] = R.at<double>(i,j);

	cam.T[0] = translation_vecs[0].at<double>(0);    //FIXME: This assumes that the desired coord sys is index 0
	cam.T[1] = translation_vecs[0].at<double>(1);    //FIXME: This assumes that the desired coord sys is index 0
	cam.T[2] = translation_vecs[0].at<double>(2);    //FIXME: This assumes that the desired coord sys is index 0

	cam.dist_coeffs[0] = dist_coeffs.at<double>(0);
	cam.dist_coeffs[1] = dist_coeffs.at<double>(1);
	cam.dist_coeffs[2] = dist_coeffs.at<double>(2);
	cam.dist_coeffs[3] = dist_coeffs.at<double>(3);

	cam.avg_reprojection_error = error_stats[0];
	cam.centriod_loc_uncertainty = error_stats[1];
}

void Calibration::writeObjectPoints(string file_path) {
    ofstream fout;
    stringstream pointfilename;
    pointfilename << file_path << "obj.pts";
    fout.open( pointfilename.str().c_str() );
    fout.flags(ios::scientific);
    if (fout.is_open()){
        for(int j = 0; j < object_points[0].size(); j++) {
            fout << object_points[0][j].x << "\t" << object_points[0][j].y << "\t" << object_points[0][j].z << endl;
        }
        fout.close();
    }
}

Calibration::~Calibration() {
    // TODO Auto-generated destructor stub
}

}
