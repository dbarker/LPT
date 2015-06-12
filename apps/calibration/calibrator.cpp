
#include "core.hpp"
#include "calib.hpp"

using namespace std;

int main(int argc, char** argv) {
    
  string input_path;
  string output_path;
     
     if (argc > 1) {
    	 input_path = argv[1];
     } 
     else{
          input_path = "../../../data/calib/jai1402/mainlist.txt";
     }
	
     if (argc > 2) {
      	  output_path = argv[2];
     }
     else{
    	 output_path = "../../../data/output/";
     }

	cout << input_path << endl;
    int corners_height = 9;//11;//
    int corners_width = 6;//4;//
    float squaresize = static_cast<float>(30.3);//10.5; //mm
	lpt::CalibrationBoard* board = new lpt::Chessboard(cv::Size(corners_height, corners_width), squaresize);
    //board.readParameters(filepath); TODO: Make this work (ALSO write a YAML function to read and write a CalibrationBoard object!!

    ifstream fin;
    ofstream fout;
    fin.open( input_path.c_str() );

    double sensor_width = 6E-3 * 640; //mm
    double sensor_height = 6E-3 * 480; //mm
	
    vector< lpt::Camera > cameras;
    string listfile;
    int id = 0;
    while( getline(fin, listfile) ){
    	lpt::Camera newcam(id);
		newcam.path = "../../../data/calib/jai1402/";
		cout << listfile << endl;
    	newcam.readImageList();
    	//newcam.sensor_size[0] = sensor_width;
    	//newcam.sensor_size[1] = sensor_height;
    	cameras.push_back(newcam);
    	listfile.clear();
    	id++;
    }
    fin.close();
	
	for (int c = 0; c < cameras.size(); c++) {
		cout << c << " list size " << cameras[c].imagelist.size() << endl; 
		for (int i = 0; i < cameras[c].imagelist.size(); i++) {
			cout << c <<" images = " << cameras[c].imagelist[i] << endl;
		}
	}
    cout << "size of cameras = " << cameras.size() << endl;
	
    lpt::Calibration calib(*board);
	calib.calib_flags = CV_CALIB_FIX_ASPECT_RATIO| CV_CALIB_FIX_K3 |
		CV_CALIB_FIX_K4 | CV_CALIB_FIX_K5;
	
	calib.calib_flags_stereo = calib.calib_flags | CV_CALIB_FIX_INTRINSIC; /*| CV_CALIB_USE_INTRINSIC_GUESS*/

    calib.find_chessboard_flags = CV_CALIB_CB_ADAPTIVE_THRESH | /*CV_CALIB_CB_FAST_CHECK |*/
		CV_CALIB_CB_FILTER_QUADS | CV_CALIB_CB_NORMALIZE_IMAGE;

    bool ok = calib.findCalibrationPoints(cameras);
    ok |= calib.runCalibration(cameras);
	cout << "calibration done" << endl;
    for (int c = 0; c < cameras.size(); ++c){
    	YAML::Emitter cameras_output_buffer;
    	cameras_output_buffer << cameras[c];
    	stringstream full_cameras_file_name;
    	full_cameras_file_name << output_path << "camera" << cameras[c].id << ".yaml";
        fout.open( full_cameras_file_name.str().c_str() );
    	fout << cameras_output_buffer.c_str();
    	fout.close();

    	YAML::Emitter frames_output_buffer;
    	frames_output_buffer << cameras[c].frames;
    	stringstream full_frames_file_name;
    	full_frames_file_name << output_path << cameras[c].id << "_pts.yaml";
    	fout.open( full_frames_file_name.str().c_str() );
    	fout << frames_output_buffer.c_str();
    	fout.close();
    }

    if (cameras.size() > 1) {
    	YAML::Emitter allcameras_output_buffer;
    	allcameras_output_buffer << cameras;
    	stringstream full_file_name;
    	full_file_name << output_path << "cameras_out.yaml";
    	fout.open( full_file_name.str().c_str() );
    	fout << allcameras_output_buffer.c_str();
    	fout.close();

    	vector<lpt::CameraPair> camera_pairs;
    	lpt::generateCameraPairs(cameras, camera_pairs);
    	calib.calcFundamentalMatrices(cameras, camera_pairs);

    	stringstream cam_pairs_file;
    	cam_pairs_file << output_path << "stereo_pairs_test.yaml";
		lpt::writeCameraPairsFile(cam_pairs_file.str(), camera_pairs);
    }

    calib.writeObjectPoints(output_path);

	cout << "Finished" << endl;

    return 0;
}
