#ifndef OUTPUT_G
#define OUTPUT_G

#include "interface.h"
#include <charm++.h>
#include <list>
#include <map>
#include "Tracker.decl.h"

extern CProxy_Main mainProxy;
extern CProxy_Merger mergeProxy;
extern CProxy_FrameSet frameSetProxy;
extern int groups;
extern int frameset_size;
extern int blocksize;
extern char output_file[100];

class OutputG : public CBase_OutputG {
public:
  int procChareCount;
  int map_count, trajcount, trajsum;
  int finishcount, mytsize, mytstart, mytend;
  bool allIn, chareCountrec;
  map <int, bool > checkinMap;
  vector < map <int, int> > mapbuffer;

  vector< vector <int> > fullMap;

  vector< list< Trajectory > > builtTrajsLst;
    
  OutputG();
  OutputG(CkMigrateMessage*) {}
  void receivecounts(int);
  void recvMap(map <int, int> , int);
  void buildfullmap(); 
  void addTrajToList(Trajectory, int);
  void sendTrajRequest();
  void addShortTrajs(Trajectory);
  void addtocheckinMap(int);
  void finishmerge(int, int);
  void finished(int);
  void printTraj(Trajectory&);
  void fprintTrajs();
};



#endif /*OUTPUT_G*/
