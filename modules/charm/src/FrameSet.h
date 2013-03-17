#ifndef FRAME_SET_H
#define FRAME_SET_H

#include <pup.h>
#include <pup_stl.h>
#include "TrackParticles.h"
#include "Tracker.decl.h"

extern /*readonly*/ CProxy_Main mainProxy;
extern CProxy_Merger mergeProxy;
extern CProxy_FrameSet frameSetProxy;
extern CProxy_OutputG outputProxy;

extern int frameset_size;
extern int groups;


class FrameSet : public CBase_FrameSet {
public:
  vector<Frame*> frameset;
  vector<Trajectory*> FinalTrajs; //FIXME: tracker->FinalTrajs and this->FinalTrajs are not both needed! Pick one and eliminate the other
  CkVec<Trajectory> CkFinalTrajs;

  vector<Frame> actualframes;
  bool *sent;

  FrameSet() {}
  FrameSet(CkMigrateMessage*) {}

  void initialize(CkVec<Frame>);
  void merge();
  void sendShortTrajs();
  void output(CkVec<int>, CkVec<int>, int );
};

#endif /* FRAME_SET_H */
