
#include "core.hpp"
#include "datagen.hpp"
#include "calib.hpp"
#include "boost/random.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace std;
using namespace boost::accumulators;

typedef accumulator_set< double, features<tag::mean, tag::variance, tag::count, tag::max, tag::min > > boost_accumulator;

int main(int argc, char** argv) {

	string trajectories_file;
	string output_path;

	if (argc > 1) {
		trajectories_file = argv[1];
	}
	else{
		trajectories_file = "../../../data/tests/pivstd/300p145fpiv352";
	}

	if (argc > 2) {
		output_path = argv[2];
	}
	else{
		output_path = "../../../data/tests/pivstd/6cams/";
	}
	cout << output_path << endl;
	cout << trajectories_file << endl;

	int number_of_cameras = 6;
	auto shared_objects = std::make_shared<lpt::SharedObjects>(); 

	auto& cameras = shared_objects->cameras;
	auto& camera_pairs = shared_objects->camera_pairs;
	
	int image_width = 640;//1280;//
	int image_height = 480;//1024;//
	double pixel_width = 0.006;//0.0048;//mm
	double aspect_ratio = 1;
	double focal_length = 4.5;//5; //mm
	double dist_coeffs[4] = {0};       //FIXME:: Assumes only 4 distortion coeffs
	double sensor_width = pixel_width * image_width; 	//mm
	double sensor_height = pixel_width * aspect_ratio * image_height;  //mm

	double camera_spacing = 500;  //mm
	double camera_plain_dist = 500; //mm
	double pi = 4.0 * atan(1.0);
	double alpha = 45 * pi / 180 / 2.0;//asin(( camera_spacing / 2.0 )/ camera_plain_dist );
	
	for (int i = 0; i < number_of_cameras; i++) {
		cameras.push_back( lpt::Camera(i) );
		cameras[i].path = output_path;
		stringstream name;
		name << "Synthetic Camera-" << i;
		cameras[i].name = name.str();
		cameras[i].pixel_size[0] = pixel_width;
		cameras[i].pixel_size[1] = pixel_width * aspect_ratio;
		cameras[i].sensor_size[0] = sensor_width;
		cameras[i].sensor_size[1] = sensor_height;
      std::cout << name.str() << std::endl;
	}

	//////
	//double angle = 0;
	//vector<double> angles (cameras.size(), 0); 
	//for (int i = 0; i < cameras.size(); ++i) {
	//	double a[] = {0, angle, 0};
	//	double T[] = {0, 0, camera_plain_dist};
	//	lpt::setCameraRotation(cameras[i], a);
	//	lpt::setCameraTranslation(cameras[i], T);
	//	lpt::setCameraIntrinsics(cameras[i], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);
	//	angles[i] = angle;
	//	angle += pi / 4.0; //pi / static_cast<double>( cameras.size() );// 
	//}
	////////

	double angles0[] = {-alpha, -alpha, 0};
	double trans0[] = {0.0 , 0.0, camera_plain_dist};

	double angles1[] = {-alpha, alpha, 0};
	double trans1[] =  {0.0 , 0.0, camera_plain_dist};

	double angles2[] = {alpha, alpha, 0};
	double trans2[] =  {0.0, 0.0, camera_plain_dist};

	double angles3[] = {alpha, -alpha, 0};
	double trans3[] =  {0.0, 0.0, camera_plain_dist};
	
	double angles4[] = {2*alpha, 0, 0};
	double trans4[] =  {0.0, 0.0, camera_plain_dist};
	
	
	double angles5[] = {-2*alpha, 0, 0};
	double trans5[] =  {0.0, 0.0, camera_plain_dist};
	/*
	
	double angles6[] = {0, 0, 0};
	double trans6[] =  {0.0 , 0.0, camera_plain_dist};

	double angles7[] = {0, 2*alpha, 0};
	double trans7[] =  {0.0, 0.0, camera_plain_dist};

	double angles8[] = {alpha, 3 * alpha, 0};
	double trans8[] =  {0.0, 0.0, camera_plain_dist};

	double angles9[] = {-alpha, 3 * alpha, 0};
	double trans9[] =  {0.0, 0.0, camera_plain_dist};

	double angles10[] = {-2 * alpha, 2 * alpha, 0};
	double trans10[] =  {0.0, 0.0, camera_plain_dist};

	double angles11[] = {2 * alpha, 2 * alpha, 0};
	double trans11[] =  {0.0, 0.0, camera_plain_dist};
	*/
	
	lpt::setCameraRotation(cameras[0], angles0);
	lpt::setCameraTranslation(cameras[0], trans0);
	lpt::setCameraIntrinsics(cameras[0], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);

	lpt::setCameraRotation(cameras[1], angles1);
	lpt::setCameraTranslation(cameras[1], trans1);
	lpt::setCameraIntrinsics(cameras[1], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);

	lpt::setCameraRotation(cameras[2], angles2);
	lpt::setCameraTranslation(cameras[2], trans2);
	lpt::setCameraIntrinsics(cameras[2], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);

	lpt::setCameraRotation(cameras[3], angles3);
	lpt::setCameraTranslation(cameras[3], trans3);
	lpt::setCameraIntrinsics(cameras[3], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);

	lpt::setCameraRotation(cameras[4], angles4);
	lpt::setCameraTranslation(cameras[4], trans4);
	lpt::setCameraIntrinsics(cameras[4], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);

	lpt::setCameraRotation(cameras[5], angles5);
	lpt::setCameraTranslation(cameras[5], trans5);
	lpt::setCameraIntrinsics(cameras[5], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);
/*		
	lpt::setCameraRotation(cameras[6], angles6);
	lpt::setCameraTranslation(cameras[6], trans6);
	lpt::setCameraIntrinsics(cameras[6], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);
	
	lpt::setCameraRotation(cameras[7], angles7);
	lpt::setCameraTranslation(cameras[7], trans7);
	lpt::setCameraIntrinsics(cameras[7], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);

	lpt::setCameraRotation(cameras[8], angles8);
	lpt::setCameraTranslation(cameras[8], trans8);
	lpt::setCameraIntrinsics(cameras[8], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);
	
	lpt::setCameraRotation(cameras[9], angles9);
	lpt::setCameraTranslation(cameras[9], trans9);
	lpt::setCameraIntrinsics(cameras[9], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);

	lpt::setCameraRotation(cameras[10], angles10);
	lpt::setCameraTranslation(cameras[10], trans10);
	lpt::setCameraIntrinsics(cameras[10], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);

	lpt::setCameraRotation(cameras[11], angles11);
	lpt::setCameraTranslation(cameras[11], trans11);
	lpt::setCameraIntrinsics(cameras[11], focal_length, pixel_width, aspect_ratio, image_width, image_height, dist_coeffs);*/

	string cameras_file = output_path + "cameras.yaml";
	string pairs_file = output_path + "camera_pairs.yaml";
	string points_basename = "_pts.yaml";
	//lpt::Input in;
	//lpt::Output out;
	//boost::mt19937 rng;
	//boost::normal_distribution<double> normal_dist(0, 2.5);
	//boost::uniform_real<double> uniform_dist(-2.5, 2.5);
	//boost_accumulator accum, accum2;

	//for (int i = 0; i < 100000;++i ) {
	//	accum(normal_dist(rng));
	//	accum2(uniform_dist(rng));
	//}
//	cout << "mean = " << extract_result<tag::mean>(accum) << " +- " << sqrt(extract_result<tag::variance>(accum)) << endl;
//	cout << "mean2 = " << extract_result<tag::mean>(accum2) << " +- " << sqrt(extract_result<tag::variance>(accum2)) << endl;

	//auto trajs = in.trajinput(trajectories_file);
	//out.fprintTrajs("300p145fpiv352", trajs);

	double x_factor = 1.0/1000.0;
	double xx_factor = 1.0;

	lpt::Camera cam1; 
	double tt = 90 * pi / 180 * x_factor;
	double dd [] = {tt, tt, tt};
	lpt::setCameraRotation(cam1, dd);
	
	cv::Mat rvec = cv::Mat::zeros(3,1,CV_64F);
	
	cv::Rodrigues(cv::Mat(3,3,CV_64F, cam1.R),rvec);
	double rvec_u = rvec.at<double>(1);
	
	for (int i = 0; i < cameras.size(); i++) {
		cameras[i].f_u[0] = x_factor * focal_length / cameras[i].pixel_size[0];// 0.01 * x_factor; //cameras[i].f[0] / 100.0;
		cameras[i].f_u[1] = x_factor * focal_length / cameras[i].pixel_size[1];//0.01 * x_factor; //cameras[i].f[1] / 100.0;
		
		cameras[i].c_u[0] = x_factor * image_width / 2.0; //0.01 * x_factor; // cameras[i].c[0] / 1000.0;
		cameras[i].c_u[1] = x_factor * image_height / 2.0;//0.01 * x_factor; //cameras[i].c[1] / 1000.0;

		cameras[i].T_u[0] = x_factor * camera_plain_dist / 10;//0.01 * x_factor;//camera_plain_dist / 100.0;
		cameras[i].T_u[1] = x_factor * camera_plain_dist / 10;//0.01 * x_factor;//camera_plain_dist / 100.0;
		cameras[i].T_u[2] = x_factor * camera_plain_dist;//0.1 * x_factor;
		
		cameras[i].centriod_loc_uncertainty = x_factor * image_width;//0.01 * x_factor;
		for (int d = 0; d < 3; ++d ) {
			cameras[i].r_vec_u[d] = rvec_u;//0.0001 * xx_factor;
		}
	}
	
	for (int cb = 1; cb < cameras.size(); ++cb) {
			lpt::CameraPair newpair(cameras[0], cameras[cb]);
			camera_pairs.push_back( newpair );
	}

   std::cout << "generating camera pairs" << std::endl;
	lpt::generateCameraPairs(cameras,camera_pairs);
	lpt::calcFundamentalMatrices( camera_pairs);
	
	lpt::writeCamerasFile(cameras_file, cameras);
	lpt::writeCameraPairsFile(pairs_file, camera_pairs);

	/*double length = 300;
	array<double, 3> origin = {{-length/2,-length/2,-length/2}};
	int num_points = 15;
	double spacing = length / num_points;
	vector<lpt::Particle3d_Ptr> objects;
	int id = 0;
	for (int i = 0; i <= num_points; i++ ) {
		for (int j = 0; j <= num_points; j++ ) {		
			for (int k = 0; k <= num_points; k++ ) {
				array<double,3> coords =  {{i * spacing + origin[0], j * spacing  + origin[1], k * spacing  + origin[2]}};
				if (coords[0] != 0 && coords[1] != 0 && coords[2] != 0) {
					objects.push_back( lpt::Particle3d::create(coords, id) );
					++id;
				}
			}
		}
	}
	ofstream fout;
	string gridfile = "../../../data/grid.txt";*/
	/*fout.open(gridfile);

	for (int p = 0; p < objects.size(); ++p) {
		for (int f = 0; f < 100; f++) {
			fout << objects[p]->id << "\t" << f << "\t" <<objects[p]->X[0] << "\t" <<objects[p]->X[1] << "\t" <<objects[p]->X[2] << endl; 
		}
		fout << endl;
	}
	fout.close();
	*/
	//cout << "Particles size = " << objects.size() << endl;
	//fout.open("stats_3mm_2_0n_45deg.txt");
	//int i = 1;
	//for (auto iter = camera_pairs.begin(); iter != camera_pairs.end(); ++iter, ++i) {
	//	vector<lpt::Camera> cams;
	//	cams.push_back(iter->cam_A);
	//	cams.push_back(iter->cam_B);
	//	vector<double> uncertainties;
	//	array<double, 12> stats = {{0,0,0,0,0,0,0,0,0,0,0,0}};
	//	compute3DPositionUncertainties(cams, objects, uncertainties, stats);
	//	fout << iter->cam_A.id << "-" << iter->cam_B.id <<"\t"; //angles[i]*180/pi <<"\t"; //
	//	for (int s = 0; s < 12; ++s) {
	//		fout << stats[s] << "\t";
	//	}
	//	fout << endl;
	//}

	
	/*for (int ca = 0; ca < cameras.size() - 1; ++ca) {
		for (int cb = ca + 1; cb < cameras.size() - 1; ++cb) {
			for (int cc = cb + 1; cc < cameras.size(); ++cc) {
				vector<lpt::Camera> cams3combos;
				cams3combos.push_back(cameras[ca]);
				cams3combos.push_back(cameras[cb]);
				cams3combos.push_back(cameras[cc]);
				fout << ca << "-"<< cb << "-" << cc <<"\t";
				vector<double> uncertainties;
				array<double, 12> stats = {{0,0,0,0,0,0,0,0,0,0,0,0}};
				compute3DPositionUncertainties(cams3combos, objects, uncertainties, stats);
				for (int s = 0; s < 12; ++s) {
					fout << stats[s] << "\t";
				}
				fout << endl;

			}
		}
	}*/

	/*vector<lpt::Camera> cams4combos;
	cams4combos.push_back(cameras[0]);
	cams4combos.push_back(cameras[1]);
	cams4combos.push_back(cameras[2]);
	cams4combos.push_back(cameras[3]);
	fout << 0 << "-"<< 1 << "-" << 2 << "-" << 3 << "\t";
	vector<double> uncertainties;
	array<double, 12> stats = {{0,0,0,0,0,0,0,0,0,0,0,0}};
	compute3DPositionUncertainties(cams4combos, objects, uncertainties, stats);
	for (int s = 0; s < 12; ++s) {
		fout << stats[s] << "\t";
	}
	fout << endl;*/

	//fout << "\n\n" << 
	//	cameras[0].f_u[0] << "\t" <<
	//	cameras[0].f_u[1] << "\t" <<
	//	
	//	cameras[0].c_u[0] << "\t" <<
	//	cameras[0].c_u[1] << "\t" <<

	//	cameras[0].T_u[0] << "\t" <<
	//	cameras[0].T_u[1] << "\t" <<
	//	cameras[0].T_u[2] << "\t" <<
	//	
	//	cameras[0].centriod_loc_uncertainty << "\t" <<
	//	
	//	cameras[0].r_vec_u[0] << "\t" << endl;

	//fout.close();
	auto image_creator = std::make_shared<lpt::ImageCreator>();
	image_creator->radius = 0;
	image_creator->intensity = 0;
	image_creator->object_intensity = static_cast<int>(5E8);
	image_creator->object_size = 3;
	image_creator->blur_ksize = 3;
	//creator.image_type = cv::Mat::zeros( cv::Size(1280, 1024), CV_8UC1 );

	lpt::DataSetGenerator generate;

	generate.setSharedObjects(shared_objects);
	cout << "camera pairs size: " << shared_objects->camera_pairs.size() << endl;
	generate.setDataPath(output_path);
	generate.setImageCreator(image_creator);

	generate.read3DTrajectoryFile(trajectories_file, lpt::PLAINTEXT);
	cout << "Projecting" << endl;
	generate.project3DFramesTo2D();
	generate.showImages();

	generate.writeImageFramePoints(output_path, points_basename);
	cout << "camera pairs size: " << shared_objects->camera_pairs.size() << endl;
		
	cout << endl << "Finished" << endl;
	return 0;
}
