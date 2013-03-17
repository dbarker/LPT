
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <fstream>
#include <iomanip>
#include <sstream>
#include "yaml-cpp/yaml.h" //TODO: make this an option?

#ifdef UNIX  // TODO: for LINUX only!! add proper function for windows timer
#include <sys/time.h>
#endif

namespace lpt {

using namespace std;

class Input;
class Output;

double get_clock();

ostream& operator<< (ostream& os, const Particle::Ptr particle);
ostream& operator<< (ostream& os, const vector<Particle::Ptr>& particles);
ostream& operator<< (ostream& os, const Match::Ptr match);
ostream& operator<< (ostream& os, const ImageFrame& frame);
ostream& operator<< (ostream& os, const Camera& cam);
ostream& operator<< (ostream& os, const CameraPair& pair);

YAML::Emitter& operator<< (YAML::Emitter& out, const Particle::Ptr particle);
YAML::Emitter& operator<< (YAML::Emitter& out, const vector<Particle::Ptr>& particles);
YAML::Node& operator>> (YAML::Node& node, vector<Particle*>& particles);

YAML::Emitter& operator<< (YAML::Emitter& out, const ImageFrame& frames);
YAML::Emitter& operator<< (YAML::Emitter& out, const vector<ImageFrame>& frames);
YAML::Node& operator>> (YAML::Node& node, ImageFrame& frame);

YAML::Emitter& operator<< (YAML::Emitter& out, const Frame* frames);
YAML::Emitter& operator<< (YAML::Emitter& out, const vector<Frame*>& frames);
YAML::Node& operator>> (YAML::Node& node, Frame& frame);

YAML::Emitter& operator<< (YAML::Emitter& out, const Camera& cam);
YAML::Emitter& operator<< (YAML::Emitter& out, const vector<Camera>& cameras);
YAML::Node& operator>> (YAML::Node& camera_node, Camera& cam);

YAML::Emitter& operator<< (YAML::Emitter& out, const CameraPair& camera_pair);
YAML::Emitter& operator<< (YAML::Emitter& out, const vector<CameraPair>& camera_pairs);
CameraPair getCameraPairYAML(YAML::Node& doc, vector<Camera>& cameras);

void readCamerasFile(const string filename, vector<lpt::Camera>& cameras);
void readCameraPairsFile(const string filename, vector<lpt::Camera>& cameras, vector<lpt::CameraPair>& camera_pairs);
void readImageFramesFile(const string filename, vector<lpt::ImageFrame>& frames);
void generateCameraPairs(vector<lpt::Camera>& cameras, vector<lpt::CameraPair>& camera_pairs);
void writeCamerasFile(const string filepath, const vector<lpt::Camera>& cameras);
void writeCameraPairsFile(const string filepath, const vector<lpt::CameraPair>& camera_pairs);

class Input {
 public:
  int maxframes;
  vector<Frame::Ptr> frameinput(const char* str);
  vector<Frame::Ptr> frameinputPIVSTD (const char* str, int numframes);
  vector<Trajectory::Ptr> trajinput(string filename);
  void readTrajectoryFile(string filename, vector<Trajectory3d_Ptr>& trajectories);
  void convertTrajectoriesToFrames(vector<Trajectory3d_Ptr>& trajectories, vector<Frame3d_Ptr>& frames );
  
  vector<Frame::Ptr> trajs2frames(vector<Trajectory::Ptr> &GoldTrajs);
  vector<Trajectory::Ptr> frames2trajsPIVSTD(vector<Frame::Ptr> &frames);
  vector<Trajectory::Ptr> resizetrajs(vector<Trajectory::Ptr> &trajs, int numtrajs, int numframes);

 private:
  Trajectory::Ptr createnewtraj();
  Particle::Ptr createnewparticle(string line);
};

class Output {
 public:
  int maxframes;
  void fprintTrajs(const char* str, vector<Trajectory::Ptr> &trajs);
#if(EVALRESULT)
  void fprintTrajsDetail(const char* str, vector<Trajectory*> &trajs);
  void fprintTrajsStats(const char* str, vector<Trajectory*> &trajs);
  void fprintFramesStats(const char* str, vector<Frame*> &frames);
#endif
};

} /* NAMESPACE_PT */

#endif /* UTILITIES_H_*/
