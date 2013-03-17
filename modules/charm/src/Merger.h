#ifndef MERGER_H
#define MERGER_H

#include <vector>
#include <map>
#include <charm++.h>

#include "interface.h"
#include "Cost.h"
#include "Search.h"
#include "Tracker.decl.h"

extern CProxy_Main mainProxy;
extern CProxy_Merger mergeProxy;
extern CProxy_FrameSet frameSetProxy;
extern CProxy_OutputG outputProxy;

class Merger : public CBase_Merger {
public:
  
  CkVec<Trajectory> right_ends;
  CkVec<Trajectory> left_ends;
  int whoisright, whoisleft;
  int rec_count;

  map <int, int> TrajsMerged;

  Merger() : rec_count(0) {}
  Merger(CkMigrateMessage *m) {}

  void sendMerge(CkVec<Trajectory>, int);
  void mergeTrajs();
  void mergeTrajs4pt();
  void mergeTrajs4ptSearch();
  void finishMerge();
  void deallocateTrajs(CkVec<Trajectory> &trajs);
};

#endif /*MERGER_H*/

