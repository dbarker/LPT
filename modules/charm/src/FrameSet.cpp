#include "FrameSet.h"

void FrameSet::initialize(CkVec<Frame> Ckframeset) {
  for (int i = 0; i < Ckframeset.size(); i++) {
    this->actualframes.push_back(Ckframeset[i]);
  }

  for (int i = 0; i < this->actualframes.size(); i++) {  
    this->actualframes[i].particles.clear();
    for (int j = 0; j < this->actualframes[i]._particles.size(); j++) {
      this->actualframes[i].particles.push_back(&this->actualframes[i]._particles[j]);
    }
  }

  for (int i = 0; i < this->actualframes.size(); i++) {
    this->frameset.push_back(&this->actualframes[i]);
  }
  
  //*********Track Particles***********************************
  TrackParticles* tracker = new TrackParticles();
  tracker->element_index = thisIndex;
  
#if (VERBOSE_CHARM)
  CkPrintf("--FrameSet Element %d on PE# %d: Begin Tracking from frame_index = %d\n", thisIndex,CkMyPe(), this->frameset[0]->frame_index);
#endif

  this->FinalTrajs = tracker->track(this->frameset);
  
  delete tracker;
    
  this->frameset.clear();
  vector<Frame*>().swap(this->frameset);
#if (VERBOSE_CHARM)
  CkPrintf("--FrameSet Element %d on PE# %d: Finished Tracking -> FinalTraj = %d\n", thisIndex,CkMyPe(), this->FinalTrajs.size() );
#endif
  //***********************************************************
    
  sent = new bool[FinalTrajs.size()];

  for (int i = 0; i < FinalTrajs.size(); i++){
    CkFinalTrajs.push_back(*FinalTrajs[i]);
    sent[i] = 0;
  }
  
  contribute(CkCallback(CkIndex_Main::finishInit(), mainProxy));
  
  FinalTrajs.clear();
  vector<Trajectory*>().swap(FinalTrajs);
}

void FrameSet::merge() {
  int myIndex = thisIndex;
  //CkPrintf("Element %d: calling sendMerge size = %d\n", myIndex ,CkFinalTrajs.size());
    
  CkVec<Trajectory> right_ends, left_ends;
  
  for (int i = 0; i < CkFinalTrajs.size(); i++) {  
    Trajectory *R = (CkFinalTrajs[i].trim(true));  
    Trajectory *L = (CkFinalTrajs[i].trim(false));
    right_ends.push_back(*R);
    left_ends.push_back(*L);    
  }
  
  //CkPrintf("sendMerge myIndex = %d, right_ends size = %d\n", myIndex, right_ends.size());
  //CkPrintf("sendMerge myIndex = %d, left_ends size = %d\n", myIndex, left_ends.size());
  if (myIndex != groups -1)
  	mergeProxy[myIndex].sendMerge(right_ends, myIndex);
  if (myIndex != 0)  
	mergeProxy[myIndex-1].sendMerge(left_ends, myIndex);  

}

void FrameSet::sendShortTrajs(){
    int count = 0;
    for (int i = 0; i < CkFinalTrajs.size(); i++){
      if(sent[i] == 0){
      	if(CkFinalTrajs[i].particles.size() >= N_TRUE){
          outputProxy[CkMyPe()].addShortTrajs(CkFinalTrajs[i]);
          count++;
	}
      }
     }
     outputProxy[CkMyPe()].finished(count);
     
     CkFinalTrajs.clear();
         
     actualframes.clear();
     vector<Frame>().swap(actualframes);
     
     //CkPrintf("FrameSet Element %d: Sent %d Short Trajs to Output Element %d\n",thisIndex,count,CkMyPe());
}

void FrameSet::output(CkVec<int> localIDs, CkVec<int> globalIDs, int outIndex) {
  int trajidx, globID, trajcount = 0;
  for (int i = 0; i < localIDs.size(); i++){
      trajidx = localIDs[i];
      globID = globalIDs[i];
      CkFinalTrajs[trajidx].id = globID;
      outputProxy[outIndex].addTrajToList(CkFinalTrajs[trajidx], thisIndex);
      sent[trajidx] = 1;
      trajcount++;
   }
   outputProxy[outIndex].finishmerge(thisIndex,trajcount);
   
   //CkPrintf("FrameSet Element %d: Sent %d Full Trajs to Output Element %d\n",thisIndex,trajcount,outIndex);
}
