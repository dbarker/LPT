
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

	string input_path;
	string output_path;

	if (argc > 1) {
		input_path = argv[1];
	}
	else{
		input_path = "../../../data/input/";
	}

	if (argc > 2) {
		output_path = argv[2];
	}
	else{
		output_path = "../../../data/output/";
	}
	cout << input_path << endl;
	cout << output_path << endl;
			
	vector<lpt::Camera> cameras;
	vector<lpt::CameraPair> camera_pairs;
	string cameras_file = input_path + "cameras.yaml";
	string pairs_file = input_path + "camera_pairs.yaml";

	lpt::readCamerasFile(cameras_file, cameras);
	lpt::readCameraPairsFile(pairs_file, cameras, camera_pairs);
	int sample_size = 25;
	double denominator = std::sqrt(static_cast<double>(sample_size));
	for (int i = 0; i < cameras.size(); ++i) {
		cameras[i].centriod_loc_uncertainty = 0.1 / 10;
		cameras[i].c_u[0] /= denominator;
		cameras[i].c_u[1] /= denominator;
		cameras[i].f_u[0] /= denominator;
		cameras[i].f_u[1] /= denominator;
		cameras[i].dist_coeffs_u[0] /= denominator;
		cameras[i].dist_coeffs_u[1] /= denominator;
		cameras[i].dist_coeffs_u[2] /= denominator;
		cameras[i].dist_coeffs_u[3] /= denominator;

		cameras[i].r_vec_u[0] /= denominator;
		cameras[i].r_vec_u[1] /= denominator;
		cameras[i].r_vec_u[2] /= denominator;

		cameras[i].T_u[0] /= denominator;
		cameras[i].T_u[1] /= denominator;
		cameras[i].T_u[2] /= denominator;
	}

	string file = input_path + "points.txt";
	ifstream fin(file);
	string line;
	if(!fin){
		cerr << "Cannot find/open the input file \n" << file << "\n Program exiting" << endl;
		exit(1);
	}
	vector<lpt::Particle3d_Ptr> objects;
	int id = 0;
	while (getline(fin,line)){ 
		stringstream ss(line);
		auto particle = lpt::Particle3d::create();
		ss >> particle->X[0];
		ss >> particle->X[1];
		ss >> particle->X[2];
		particle->id = id++;
		objects.push_back(particle);
	}

	string outfile = output_path + "uncertainty_stats.txt";
	ofstream fout(outfile);
	cout << " Particles size " << objects.size() << endl;
	for (int ca = 0; ca < cameras.size() - 1; ++ca) {
		for (int cb = ca + 1; cb < cameras.size() - 1; ++cb) {
			for (int cc = cb + 1; cc < cameras.size(); ++cc) {
				for (int cd = cc + 1; cd < cameras.size(); ++cd) {
					vector<lpt::Camera> cams4combos;
					cams4combos.push_back(cameras[ca]);
					cams4combos.push_back(cameras[cb]);
					cams4combos.push_back(cameras[cc]);
					cams4combos.push_back(cameras[cd]);
					
					fout << ca << "-"<< cb << "-" << cc <<"-" << cd << "\t";
					vector<double> uncertainties;
					array<double, 12> stats = {{0,0,0,0,0,0,0,0,0,0,0,0}};
					compute3DPositionUncertainties(cams4combos, objects, uncertainties, stats);
					for (int s = 0; s < 12; ++s) {
						fout << stats[s] << "\t";
					}
					fout << endl;
					cout << ca << " " << cb << " " << cc << " " << cd << endl; 
				}
			}
		}
	}
	
	cout << endl << "Finished" << endl;
	return 0;
}
