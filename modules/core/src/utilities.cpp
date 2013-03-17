
#include "core.hpp"

namespace lpt {

using namespace std;

ostream& operator<<(ostream& os, const Particle* particle)
{
	os.flags (ios::scientific);
	os << particle->id <<'\t' << particle->frame_index <<'\t';
    os << particle->x << '\t' << particle->y << '\t' << particle->z << endl;
    return os;
}

ostream& operator<<(ostream& os, const vector<Particle*>& particles)
{
	for (int i = 0; i < particles.size(); i++){
		os.flags (ios::scientific);
        os << particles[i];
	}
    return os;
}

ostream& operator<< (ostream& os, const Match* match) {
	os << "Match : [\t";
	for (int i = 0; i < match->particles.size(); ++i) {
		os << match->particles[i].second << "-" << match->particles[i].first->id << "\t";
	}
	os << "]"<<endl;
	return os;
}

ostream& operator<< (ostream& os, const ImageFrame& frame) {
	os << "Frame Index = "<< frame.frame_index << endl;
	os << "Particles: size = " << frame.particles.size() << endl;
	for (int i = 0; i < frame.particles.size(); ++i) {
		os << frame.particles[i]->id << ": [";
		os << frame.particles[i]->x << "\t";
		os << frame.particles[i]->y << "]" << endl;
	}
	return os;
}

ostream& operator<<(ostream& os, const CameraPair& pair)
{
	os << "CameraA: "<<  pair.cam_A.id<< endl;
	os << "CameraB: "<<  pair.cam_B.id<< endl;
	os.flags (ios::scientific);
	os << "Fundamental Matrix:  " << endl;
	for (int i = 0; i < 3; ++i) {
		os << "[" << "\t";
		for (int j = 0; j < 3; ++j) {
			os << pair.F[i][j] <<'\t';
		}
		os << "]" << endl;
	}
	os << endl;
    return os;
}

ostream& operator<<(ostream& os, const Camera& cam)
{
	os << "CameraID: "<<  cam.id << endl;
	os << "Camera Name: " << cam.name << endl;
	os.flags (ios::scientific);
	os << "Rotation Matrix:  " << endl;
	for (int i = 0; i < 3; ++i) {
		os << "[" << "\t";
		for (int j = 0; j < 3; ++j) {
			os << cam.R[i][j] <<'\t';
		}
		os << "]" << endl;
	}
	os << "Translation Vector = "<< endl;
	os << "[" << cam.T[0] << "\t" << cam.T[1] << "\t" << cam.T[2] << "]" << endl;
	os << "f = " <<"["<< cam.f[0] << "\t" << cam.f[1] <<"]" << endl;
	os << "c = " <<"["<< cam.c[0] << "\t" << cam.c[1] <<"]" << endl;
	os << "DistCoeffs = " << endl;
	os << "[" << cam.dist_coeffs[0] << "\t"<< cam.dist_coeffs[1] << "\t"
			<< cam.dist_coeffs[2] << "\t"<< cam.dist_coeffs[3] << "]" <<endl;

	os << "sensor size (mm) = [" << cam.sensor_size[0] << "\t" << cam.sensor_size[1] << "]" << endl;
	os << "pixel size (mm) = [" << cam.pixel_size[0] << "\t" << cam.pixel_size[1] << "]" <<endl;
	os << "centroid location uncertainty = \t" << cam.centriod_loc_uncertainty << endl;
	os << "avg_reprojection_error = \t" << cam.avg_reprojection_error << endl;
	os << endl;
    return os;
}

YAML::Emitter& operator<<(YAML::Emitter& out, const Particle* particle) {
	out << YAML::Flow << YAML::BeginSeq << particle->id << particle->frame_index;
	out << particle->x << particle->y << particle->z << YAML::EndSeq;
	return out;
}

YAML::Emitter& operator<<(YAML::Emitter& out, const vector<Particle*>& particles) {
	out << YAML::BeginMap;
	out << YAML::Key << "Trajectory" << YAML::Value << YAML::BeginSeq;
	for (int i = 0; i < particles.size(); i++){
		out << particles[i];
	}
	out << YAML::EndSeq << YAML::EndMap;
	return out;
}

YAML::Node& operator>> (YAML::Node& node, vector<Particle*>& particles) {
	const YAML::Node* particles_node = node.FindValue("Trajectory");
	for( YAML::Iterator it = (*particles_node).begin(); it != (*particles_node).end(); ++it){
		Particle* new_particle = new Particle;
		(*it)[0] >> new_particle->id;
		(*it)[1] >> new_particle->frame_index;
		(*it)[2] >> new_particle->x;
		(*it)[3] >> new_particle->y;
		(*it)[4] >> new_particle->z;
		particles.push_back(new_particle);
	}
	return node;
}

YAML::Emitter& operator<< (YAML::Emitter& out, const ImageFrame& frame) {
	out << YAML::BeginMap << YAML::Key << "FrameIndex";
	out << YAML::Value << frame.frame_index;
	out << YAML::Key << "Particles" << YAML::Value << YAML::BeginSeq;
	for (int i = 0; i < frame.particles.size(); ++i){
		out << YAML::Flow << YAML::BeginSeq << frame.particles[i]->id;
		out << frame.particles[i]->x << frame.particles[i]->y << YAML::EndSeq;
	}
	out << YAML::EndSeq << YAML::EndMap;
	return out;
}
YAML::Emitter& operator<< (YAML::Emitter& out, const vector<ImageFrame>& frames) {
    out << YAML::Comment("--x and y image pixel locations") << YAML::Newline;
    for (int f = 0; f < frames.size(); ++f)
    	out << YAML::BeginDoc << frames[f] << YAML::EndDoc;
    return out;
}

YAML::Node& operator>> (YAML::Node& node, ImageFrame& frame) {
	node["FrameIndex"] >> frame.frame_index;
    const YAML::Node* p_node = node.FindValue("Particles");
	for( YAML::Iterator part_it = (*p_node).begin(); part_it != (*p_node).end(); ++part_it){
		ParticleImage::Ptr new_particle = ParticleImage::create();
		(*part_it)[0] >> new_particle->id;
		(*part_it)[1] >> new_particle->x;
		(*part_it)[2] >> new_particle->y;
		frame.particles.push_back(new_particle);
	}
	return node;
}

YAML::Emitter& operator<< (YAML::Emitter& out, const Frame* frame) {
	out << YAML::BeginMap << YAML::Key << "FrameIndex";
	out << YAML::Value << frame->frame_index;
	out << YAML::Key << "Particles" << YAML::Value << YAML::BeginSeq;
	for (int i = 0; i < frame->particles.size(); ++i){
		out << YAML::Flow << YAML::BeginSeq << frame->particles[i]->id;
		out << frame->particles[i]->x << frame->particles[i]->y << frame->particles[i]->z << YAML::EndSeq;
	}
	out << YAML::EndSeq << YAML::EndMap;
	return out;
}
YAML::Emitter& operator<< (YAML::Emitter& out, const vector<Frame*>& frames) {
    out << YAML::Comment("--3D particle locations") << YAML::Newline;
    for (int f = 0; f < frames.size(); ++f)
    	out << YAML::BeginDoc << frames[f] << YAML::EndDoc;
    return out;
}

YAML::Node& operator>> (YAML::Node& node, Frame& frame) {
	node["FrameIndex"] >> frame.frame_index;
    const YAML::Node* p_node = node.FindValue("Particles");
	for( YAML::Iterator part_it = (*p_node).begin(); part_it != (*p_node).end(); ++part_it){
		Particle::Ptr new_particle = Particle::create();
		(*part_it)[0] >> new_particle->id;
		(*part_it)[1] >> new_particle->x;
		(*part_it)[2] >> new_particle->y;
		(*part_it)[2] >> new_particle->z;
		frame.particles.push_back(new_particle);
	}
	return node;
}

YAML::Emitter& operator<< (YAML::Emitter& out, const Camera& cam) {

	out << YAML::BeginMap;
	out << YAML::Key << "CameraID" << YAML::Value << cam.id;
	out << YAML::Key << "Name" << YAML::Value << cam.name;
	out << YAML::Key << "path"<< YAML::Value << cam.path.c_str();

	out << YAML::Key << "R" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j)
			out << cam.R[i][j];
	}
	out << YAML::EndSeq;

	out << YAML::Key << "r_vec_u" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	for (int i = 0; i < 3; ++i)
		out << cam.r_vec_u[i];
	out << YAML::EndSeq;

	out << YAML::Key << "T" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	for (int i = 0; i < 3; ++i)
			out << cam.T[i];
	out << YAML::EndSeq;

	out << YAML::Key << "T_u" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	for (int i = 0; i < 3; ++i)
			out << cam.T_u[i];
	out << YAML::EndSeq;

	out << YAML::Key << "f" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	out << cam.f[0] << cam.f[1] << YAML::EndSeq;

	out << YAML::Key << "f_u" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	out << cam.f_u[0] << cam.f_u[1] << YAML::EndSeq;

	out << YAML::Key << "c" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	out << cam.c[0] << cam.c[1] << YAML::EndSeq;

	out << YAML::Key << "c_u" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	out << cam.c_u[0] << cam.c_u[1] << YAML::EndSeq;

	out << YAML::Key << "distortion" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	for (int i = 0; i < 4; ++i)      //FIXME: assumes only four distortion coeffs!! Make this a parameter
		out << cam.dist_coeffs[i];
	out << YAML::EndSeq;

	out << YAML::Key << "distortion_uncertainty" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	for (int i = 0; i < 4; ++i)      //FIXME: assumes only four distortion coeffs!! Make this a parameter
		out << cam.dist_coeffs_u[i];
	out << YAML::EndSeq;

	out << YAML::Key << "centriod_loc_uncertainty" << YAML::Value << cam.centriod_loc_uncertainty;
	
	out << YAML::Key << "sensor_size" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	out << cam.sensor_size[0] << cam.sensor_size[1] << YAML::EndSeq;
	
	out << YAML::Key << "pixel_size" << YAML::Value << YAML::Flow << YAML::BeginSeq;
	out << cam.pixel_size[0] << cam.pixel_size[1] << YAML::EndSeq;

	out << YAML::Key << "avg_reprojection_error" << YAML::Value << cam.avg_reprojection_error;
	out << YAML::Key << "imagelist"<< YAML::Value << YAML::BeginSeq;
	for (int i = 0; i < cam.imagelist.size(); ++i)
		out << cam.imagelist[i];
	out << YAML::EndSeq;
	out << YAML::EndMap;

	return out;
}

YAML::Emitter& operator<< (YAML::Emitter& out, const vector<Camera>& cameras) {

	for (int i = 0; i < cameras.size(); ++i)
		out << YAML::BeginDoc << cameras[i] << YAML::EndDoc;

	return out;
}

YAML::Node& operator>> (YAML::Node& camera_node, Camera& cam) {

	const YAML::Node* id_node = camera_node.FindValue("CameraID");
	(*id_node) >> cam.id;

	const YAML::Node* name_node = camera_node.FindValue("Name");
	(*name_node) >> cam.name;

	const YAML::Node* path_node = camera_node.FindValue("path");
	(*path_node) >> cam.path;

	const YAML::Node* r_node = camera_node.FindValue("R");
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j)
			(*r_node)[j + 3*i] >> cam.R[i][j];
	}

	const YAML::Node* rvec_u_node = camera_node.FindValue("r_vec_u");
	for (int i = 0; i < 3; ++i)
		(*rvec_u_node)[i] >> cam.r_vec_u[i];

	const YAML::Node* t_node = camera_node.FindValue("T");
	for (int i = 0; i < 3; ++i)
		 (*t_node)[i] >> cam.T[i];

	const YAML::Node* t_u_node = camera_node.FindValue("T_u");
	for (int i = 0; i < 3; ++i)
		 (*t_u_node)[i] >> cam.T_u[i];

	const YAML::Node* f_node = camera_node.FindValue("f");
	(*f_node)[0] >> cam.f[0];
	(*f_node)[1] >> cam.f[1];

	const YAML::Node* f_u_node = camera_node.FindValue("f_u");
	(*f_u_node)[0] >> cam.f_u[0];
	(*f_u_node)[1] >> cam.f_u[1];

	const YAML::Node* c_node = camera_node.FindValue("c");
	(*c_node)[0] >> cam.c[0];
	(*c_node)[1] >> cam.c[1];

	const YAML::Node* c_u_node = camera_node.FindValue("c_u");
	(*c_u_node)[0] >> cam.c_u[0];
	(*c_u_node)[1] >> cam.c_u[1];

	const YAML::Node* d_node = camera_node.FindValue("distortion");
	for (int i = 0; i < 4; ++i)                     //FIXME:  Allow number of distortion parameters to change
		(*d_node)[i] >> cam.dist_coeffs[i];

	const YAML::Node* d_u_node = camera_node.FindValue("distortion_uncertainty");
	for (int i = 0; i < 4; ++i)                     //FIXME:  Allow number of distortion parameters to change
		(*d_u_node)[i] >> cam.dist_coeffs_u[i];

	const YAML::Node* sensor_size_node = camera_node.FindValue("sensor_size");
	(*sensor_size_node)[0] >> cam.sensor_size[0];
	(*sensor_size_node)[1] >> cam.sensor_size[1];
	
	const YAML::Node* pixel_size_node = camera_node.FindValue("pixel_size");
	(*pixel_size_node)[0] >> cam.pixel_size[0];
	(*pixel_size_node)[1] >> cam.pixel_size[1];

	const YAML::Node* centroid_u_node = camera_node.FindValue("centriod_loc_uncertainty");
	(*centroid_u_node) >> cam.centriod_loc_uncertainty;

	const YAML::Node* e_node = camera_node.FindValue("avg_reprojection_error");
	(*e_node) >> cam.avg_reprojection_error;

	return camera_node;
}

YAML::Emitter& operator<< (YAML::Emitter& out, const CameraPair& camera_pair) {
	out <<YAML::BeginMap;
	out << YAML::Key << "Camera_A_id" << YAML::Value << camera_pair.cam_A.id;
	out << YAML::Key << "Camera_B_id" << YAML::Value << camera_pair.cam_B.id;
	out << YAML::Key << "FundamentalMatrix"<< YAML::Value << YAML::Flow;
	out << YAML::BeginSeq;
	for (int m = 0; m < 3; ++m)
		for (int n = 0; n < 3; ++n)
			out << camera_pair.F[m][n];
	out << YAML::EndSeq;
	out << YAML::Key << "reprojection_error" << YAML::Value << camera_pair.reprojection_error;
	out << YAML::Key << "epipolar_error" << YAML::Value << camera_pair.epipolar_error;
	out << YAML::EndMap;
	return out;
}

YAML::Emitter& operator<< (YAML::Emitter& out, const vector<CameraPair>& camera_pairs) {
	for (int i = 0; i < camera_pairs.size(); ++i)
		out << YAML::BeginDoc << camera_pairs[i] << YAML::EndDoc;

	return out;
}

lpt::CameraPair getCameraPairYAML(YAML::Node& doc, vector<lpt::Camera>& cameras) {
	int a_id;
	int b_id;
	doc["Camera_A_id"] >> a_id;
	doc["Camera_B_id"] >> b_id;
	CameraPair newpair(cameras[a_id], cameras[b_id]);
	const YAML::Node* F_node = doc.FindValue("FundamentalMatrix");
	for (int m = 0; m < 3; ++m){
		for (int n = 0; n < 3; ++n)
			(*F_node)[n + 3*m] >> newpair.F[m][n];
	}
	doc["reprojection_error"] >> newpair.reprojection_error;
	doc["epipolar_error"] >> newpair.epipolar_error;
	return newpair;
}

void readCamerasFile(const string filename, vector<lpt::Camera>& cameras) {
    cout << "Opening file to read in Camera data" << endl;
	cout << "\t" << filename << endl;
	ifstream fin( filename.c_str() );
	YAML::Parser parser(fin);
	YAML::Node doc;
	while (parser.GetNextDocument(doc)){
		lpt::Camera cameranew;
		doc >> cameranew;
		cameras.push_back(cameranew);
	}
	fin.close();
	cout << "Read in " << cameras.size() << " cameras" << endl;
}

void readCameraPairsFile( const string filename, vector<lpt::Camera>& cameras, vector<lpt::CameraPair>& camera_pairs) {
	ifstream fin( filename.c_str() );
	YAML::Parser parser(fin);
	YAML::Node doc;

	while (parser.GetNextDocument(doc)){
		camera_pairs.push_back( getCameraPairYAML(doc, cameras) );
	}
	fin.close();
}

void readImageFramesFile( const string filename, vector<lpt::ImageFrame>& frames) {
	ifstream fin( filename.c_str() );
	if ( fin.is_open() ) {
		YAML::Parser parser(fin);
		YAML::Node doc;
		int frame_index = 0;
		while (parser.GetNextDocument(doc)){
			lpt::ImageFrame newframe;
			frames.push_back(newframe);
			doc >> frames[frame_index];            //FIXME: Assumes frames are stored sequentially
			frame_index++;
		}
		fin.close();
	} else {
		cout << "File could not be opened: " << filename << endl;
	}
}

void generateCameraPairs(vector<lpt::Camera>& cameras, vector<lpt::CameraPair>& camera_pairs) {
	for (int ca = 0; ca < cameras.size(); ++ca) {
		for (int cb = ca + 1; cb < cameras.size(); ++cb) {
			lpt::CameraPair newpair(cameras[ca], cameras[cb]);
			camera_pairs.push_back( newpair );
		}
	}
	cout << "camera pairs size: " << camera_pairs.size() << endl;
}

void writeCamerasFile(const string filepath, const vector<lpt::Camera>& cameras) {
	cout << "Writing camera parameters to file" << endl;
	cout << "\t" << filepath << endl;
	YAML::Emitter cameras_output_buffer;
	cameras_output_buffer << cameras;
	ofstream fout( filepath.c_str() );
	fout << cameras_output_buffer.c_str();
	fout.close();
}

void writeCameraPairsFile(const string filepath, const vector<lpt::CameraPair>& camera_pairs) {
	cout << "Writing camera pairs parameters to file" << endl;
	cout << "\t" << filepath << endl;
	YAML::Emitter pairs_output_buffer;
    pairs_output_buffer << camera_pairs;
    
	ofstream fout( filepath.c_str() );
    fout << pairs_output_buffer.c_str();
    fout.close();
}

vector<Frame::Ptr> Input::frameinput (const char* str) {

  int f = 0;
  double temp = 0;
  string line;
  
  ifstream fin;
  fin.open(str);

  if(!fin){
    cerr << "Cannot find/open the input file \n" << str << "\n Program exiting (input.cpp)" << endl;
    exit(1);
  }
 
  vector<Frame::Ptr> frames;			
  while (getline(fin,line)){ 
    Frame::Ptr newFrm = Frame::create();
    stringstream ss(line);
          
    while (ss >> temp){               //input value in "line" into temp as a double
	 
    	Particle::Ptr newPart = Particle::create();
    	newPart->x = temp;           //set x value
	
    	ss >> temp;
    	newPart->y = temp;           //set y value
	
        ss >> temp;
        newPart->z = temp;           //set z value

        newPart->id = newPart->id = -1; // FIXME: drand48()*1000000000;
        newPart->frame_index = f;
        newFrm->particles.push_back(newPart);
      }
    
    //empty frame check
      if (newFrm->particles.empty())
         break;
    
    newFrm->frame_index = f;
   
    frames.push_back(newFrm);
    f++;
  }

 fin.close(); 

 return frames;
}

vector<Frame::Ptr> Input::frameinputPIVSTD (const char* str, int numframes) {
  char filename[100];
  vector<Frame::Ptr> frames;

  for (int f = 0; f < numframes; f++){

	  sprintf(filename, "%s%03d.dat",str,f);
	  printf("%s\n",filename);

	  string line;
	  ifstream fin;
  	  fin.open(filename);

  	  if(!fin){
  		  cerr << "Cannot find/open the input file \n" << filename << "\n Program exiting (input.cpp)" << endl;
  		  exit(1);
  	  }

  	  Frame::Ptr newFrm = Frame::create();

  	  while (getline(fin,line)){
  		  stringstream ss(line);
		  Particle::Ptr newPart = Particle::create();

		  ss >> newPart->id;
		  ss >> newPart->x;           //set x value
		  ss >> newPart->y;           //set y value
		  ss >> newPart->z;           //set z value

		  newPart->frame_index = f;
		  newFrm->particles.push_back(newPart);
	  }

  	  // empty frame check 
  	  if (newFrm->particles.empty())
              continue;

  	  newFrm->frame_index = f;

  	  frames.push_back(newFrm);
  	  fin.close();
  }

  return frames;
}

vector<Trajectory::Ptr> Input::frames2trajsPIVSTD(vector<Frame::Ptr> &frames){

	map<int, Trajectory::Ptr> trajs;

	for (int i = 0; i < frames.size(); i++){
		Frame::Ptr frame = frames[i];
		for (int p = 0; p < frame->particles.size(); p++){
			map<int,Trajectory::Ptr>::iterator it;
			Particle::Ptr P = frame->particles[p];
			it = trajs.find(P->id);
			if(it == trajs.end()){
				Trajectory::Ptr T = Trajectory::create();
				T->id = P->id;
				T->startframe = i;
				trajs[P->id] = T;
			}
			trajs[P->id]->particles.push_back(P);
		}
	}
	vector<Trajectory::Ptr> Vtrajs;
	for (map<int,Trajectory::Ptr>::iterator it = trajs.begin(); it != trajs.end(); ++it){
		Vtrajs.push_back(it->second);
	}
	return Vtrajs;
}


//*********************************************

void Input::readTrajectoryFile(string filename, vector<Trajectory3d_Ptr>& trajectories) {
	
  this->maxframes = 0;
  
  string line;
  ifstream fin( filename.c_str() );

  if(!fin){
    cerr << "Cannot find/open the input file \n" << filename << "\n Program exiting ( Input::trajinput(string filename) )" << endl;
    exit(1);
  }
  	
  Trajectory3d_Ptr newTraj = Trajectory3d::create();

  while (getline(fin,line)){
	  if (line == "" || line.length() < 10){
		  if (newTraj->objects.empty() == false){
			  newTraj->id = newTraj->objects[0]->id;//GoldTrajs.size()+1;
			  newTraj->startframe = newTraj->objects[0]->frame_index;//FIXME: uncomment
			  int fcount = newTraj->startframe + newTraj->objects.size();
			  if(fcount > this->maxframes){
				  this->maxframes = fcount;
			  }
			  trajectories.push_back(newTraj);
			  
			  newTraj = Trajectory3d::create();
		  }
	  }else{
		  Particle3d_Ptr newPart = Particle3d::create();
		  	stringstream ss(line);
			ss >> newPart->id;          //FIXME uncomment
			ss >> newPart->frame_index; //FIXME uncomment
			ss >> newPart->X[0];           //set x value
			ss >> newPart->X[1];           //set y value
			ss >> newPart->X[2];           //set z value
			newTraj->objects.push_back(newPart);
	  }
  }
 //Pickup last trajectory if no blank line at the end of file
 if (newTraj->objects.empty() == false)
      trajectories.push_back(newTraj);

 fin.close(); 

}

vector<Trajectory::Ptr> Input::trajinput (string filename) {

  this->maxframes = 0;
  
  vector<Trajectory::Ptr> GoldTrajs;
  
  string line;
  ifstream fin( filename.c_str() );

  if(!fin){
    cerr << "Cannot find/open the input file \n" << filename << "\n Program exiting ( Input::trajinput(string filename) )" << endl;
    exit(1);
  }
  	
  Trajectory::Ptr newTraj = createnewtraj();
  while (getline(fin,line)){
	  if (line == "" || line.length() < 10){
		  if (newTraj->particles.empty() == false){
			  newTraj->id = newTraj->particles[0]->id;//GoldTrajs.size()+1;
			  newTraj->startframe = newTraj->particles[0]->frame_index;//FIXME: uncomment
			  int fcount = newTraj->startframe + newTraj->particles.size();
			  if(fcount > this->maxframes){
				  this->maxframes = fcount;
			  }
			  GoldTrajs.push_back(newTraj);
			  
			  newTraj = createnewtraj();
		  }
	  }else{
		  Particle::Ptr newPart = createnewparticle(line);
		  //newPart->id = GoldTrajs.size();//FIXME:Comment out
		  //newPart->frame_index = newTraj->particles.size(); //FIXME: Comment Out
		  newTraj->particles.push_back(newPart);
	  }
  }
 //Pickup last trajectory if no blank line at the end of file
 if (newTraj->particles.empty() == false)
      GoldTrajs.push_back(newTraj);

 fin.close(); 
 
 return GoldTrajs;
}

void Input::convertTrajectoriesToFrames(vector<Trajectory3d_Ptr>& trajectories, vector<Frame3d_Ptr>& frames ) {

  for(int i = 0; i < maxframes; i++){
    	Frame3d_Ptr newFrm = Frame3d::create();
    	newFrm->frame_index = i;
    	frames.push_back(newFrm);
  }
     
  for(int j = 0; j < trajectories.size(); ++j){
	for(int p = 0; p < trajectories[j]->objects.size();p++){
	        Particle3d_Ptr P1 = trajectories[j]->objects[p];
		int index = P1->frame_index;
		frames[index]->objects.push_back(P1);
	}
  }
}

vector<Frame::Ptr> Input::trajs2frames(vector<Trajectory::Ptr>& GoldTrajs) {

  vector<Frame::Ptr> frames;			
  for(int i = 0; i < maxframes; i++){
    	Frame::Ptr newFrm = Frame::create();
    	newFrm->frame_index = i;
    	frames.push_back(newFrm);
  }
     
  for(int j = 0; j < GoldTrajs.size(); ++j){
	for(int p = 0; p < GoldTrajs[j]->particles.size();p++){
	        Particle::Ptr P1 = GoldTrajs[j]->particles[p];
		int index = P1->frame_index;
		frames[index]->particles.push_back(P1);
	}
  }

 return frames;
}

vector<Trajectory::Ptr> Input::resizetrajs(vector<Trajectory::Ptr>& trajs, int numtrajs, int numframes){ //FIXME:  This cannot handle trajectories of different sizes
	vector<Trajectory::Ptr> tnew;
	for (int i = 0; i < numtrajs; i++){
		tnew.push_back(trajs[i]);
		tnew[i]->particles.clear();
		for (int j = 0; j < numframes; j++){
			tnew[i]->particles.push_back(trajs[i]->particles[j]);
		}
	 }

#if(VERBOSE_SERIAL)
	printf("New Frames size = %d\n", tnew[0]->particles.size());
	printf("New Trajs size = %d \n", tnew.size());
#endif

	return tnew;
}

//*********Private Methods**************************
Trajectory::Ptr Input::createnewtraj(){
	Trajectory::Ptr newTraj = Trajectory::create();
	//newTraj->id = drand48()*1000000000;
   	newTraj->gap = 0;
   	//newTraj->startframe = 0;  //FIXME comment out
    return newTraj;
}

Particle::Ptr Input::createnewparticle(string line){
	stringstream ss(line);
	Particle::Ptr newPart = Particle::create();
	
	ss >> newPart->id;          //FIXME uncomment
	ss >> newPart->frame_index; //FIXME uncomment
	ss >> newPart->x;           //set x value
	ss >> newPart->y;           //set y value
	ss >> newPart->z;           //set z value
	
	return newPart;
}

void Output::fprintTrajs (const char* str, vector<Trajectory::Ptr>& trajs) {

  FILE* printFile = fopen (str,"w");
  
  if(!printFile){
    cerr << "Cannot find/open the output file \n" << str << "\n Program exiting (output.cpp)" << endl;
    exit(1);
  }
  
  for (unsigned int i = 0; i < trajs.size(); i++){
      for (unsigned int j = 0; j < trajs[i]->particles.size(); j++){
    	 Particle::Ptr P = trajs[i]->particles[j];
         fprintf(printFile, "%d\t%d\t%25.16e%25.16e%25.16e\n",P->id,P->frame_index,P->x , P->y , P->z );
      }
	fprintf(printFile, "\n");
  }
  fclose(printFile);

}
#if(EVALRESULT)
void Output::fprintTrajsDetail (const char* str, vector<Trajectory*> &trajs) {

  FILE* printFile = fopen (str,"w");

  if(!printFile){
    cerr << "Cannot find/open the output file \n" << str << "\n Program exiting (output.cpp)" << endl;
    exit(1);
  }

  for (unsigned int i = 0; i < trajs.size(); i++){
	  //fprintf(printFile,"%d\t%d\n",trajs[i]->id, trajs[i]->startframe);
      for (unsigned int j = 0; j < trajs[i]->particles.size(); j++){
    	 Particle* P = trajs[i]->particles[j];
         fprintf(printFile, "%d\t%d\t", P->id,P->frame_index);
         fprintf(printFile, "%f\t%f\t%f\t",P->cost, P->vel,P->accel);
         fprintf(printFile, "%d\t%25.16e%25.16e%25.16e\n",P->id,P->x, P->y, P->z);
      }
	fprintf(printFile, "\n");
  }
  fclose(printFile);

}

void Output::fprintTrajsStats (const char* str, vector<Trajectory*> &trajs) {

  FILE* printFile = fopen (str,"w");

  if(!printFile){
    cerr << "Cannot find/open the output file \n" << str << "\n Program exiting (output.cpp)" << endl;
    exit(1);
  }

  for (unsigned int i = 0; i < trajs.size(); i++){
	  fprintf(printFile,"%d\t%d\t%f\t%f\n",
			  trajs[i]->id, trajs[i]->particles.size(),
			  trajs[i]->S/trajs[i]->R, trajs[i]->Eff);/*,
			  trajs[i]->Rl, trajs[i]->R, trajs[i]->Rh);*/
  }
  fclose(printFile);

}


void Output::fprintFramesStats (const char* str, vector<Frame*> &frames) {

  FILE* printFile = fopen (str,"w");

  if(!printFile){
    cerr << "Cannot find/open the output file \n" << str << "\n Program exiting (output.cpp)" << endl;
    exit(1);
  }

  for (unsigned int i = 0; i < frames.size(); i++){
	  fprintf(printFile,"%d\t%d\t%f\n",frames[i]->frame_index, 
	  			frames[i]->particles.size(),frames[i]->S);
  }
  fclose(printFile);

}
#endif

} /* NAMESPACE_PT */
