#include "interface.h"
#include "FrameSet.h"
#include "OutputG.h"
#include "Merger.h"
#include "Input.h"
#include <string.h>
#include <alloca.h>
#include <pup_stl.h>
#include <map>
#include "Tracker.decl.h"

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_Merger mergeProxy;
/*readonly*/ CProxy_FrameSet frameSetProxy;
/*readonly*/ CProxy_OutputG outputProxy;
/*readonly*/ int frameset_size;
/*readonly*/ int framesize;
/*readonly*/ int blocksize;
/*readonly*/ int groups;
/*readonly*/ char input_file[100];
/*readonly*/ char output_file[100];
/*readonly*/ double time1;

class rrmap : public CBase_rrmap {
public:
  int elementcount,iter;
  rrmap() {elementcount = 0;}

  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int elem = *(int *)idx.data();
    return procSelect(elem);
  }
  int procSelect(int elem){
    return elem % CkNumPes();
  }
  void sendcountstoOutput(int num){
     for (int i = 0; i < num; i++){
	if (CkMyPe() == procSelect(i))
	    elementcount++;
     }
     outputProxy[CkMyPe()].receivecounts(elementcount);
   }
};

class blockmap : public CBase_blockmap {
public:
  int elementcount;
  blockmap() {elementcount=0;}
  int registerArray(CkArrayIndexMax& numElements,CkArrayID aid) {
    return 0;
  } 
  int procNum(int arrayHdl, const CkArrayIndex &idx) {
    int elem = *(int *)idx.data();
    return procSelect(elem);
  }
  int procSelect(int elem){
    int penum = (elem/blocksize);
    if (penum >= CkNumPes() )
	penum = elem % CkNumPes();
    return penum;
  }
  void sendcountstoOutput(int num){
     for (int i = 0; i < num; i++){
	if (CkMyPe() == procSelect(i))
	    elementcount++;
     }
     outputProxy[CkMyPe()].receivecounts(elementcount);
   }
};

class Main : public Chare {
public:

  CProxy_FrameSet frameSet;
  CProxy_Merger merger;
  CProxy_OutputG output_group;
  vector<Frame*> frames;
  
 // usage is ./tracker <inputdatafile> <output file basename>
  Main(CkArgMsg* m) {

     if (m->argc > 1) {
    	 sprintf(input_file, "../../data/%s", m->argv[1]);
     } 
     else{
          sprintf(input_file,"../../data/scaling/256p512f");
     }
	
     if (m->argc > 2) {
      	  sprintf(output_file, "../output/%s", m->argv[2]);
     }
     else{
  	  sprintf(output_file,"../output/out");
     }

    delete m;
    frameset_size = FRAMESET_SIZE;
    	
    Input in;
    vector<Trajectory*> GoldTrajs = in.trajinput(input_file);
    this->frames = in.trajs2frames(GoldTrajs);
    
#if(TESTRESULT == 0)
    GoldTrajs.clear();
    vector<Trajectory*>().swap(GoldTrajs);
#endif
  

    framesize = this->frames.size();
    groups = framesize / frameset_size;
    blocksize = groups/CkNumPes();
    int mgroups = groups - 1;
    
    printsettings();
    
    CkPrintf("\n\n**** Initializing FrameSet, Merge, and Output Chares --> Begin Parallel Tracking\n");  
    
    mainProxy = thishandle;
    time1 = CmiWallTimer();

#if(BLOCKMAP)
    CkPrintf("Chare Mapping:  BLOCKED (min %d chares/proc)",blocksize);
    CProxy_blockmap map1 = CProxy_blockmap::ckNew();
    CProxy_blockmap map2 = CProxy_blockmap::ckNew(); 
#else
    CkPrintf("Chare Mapping:  ROUND ROBIN");
    CProxy_rrmap map1 = CProxy_rrmap::ckNew();
    CProxy_rrmap map2 = CProxy_rrmap::ckNew();   
#endif
    
    //***Create Frame Set Chares*****
    CkArrayOptions opts1(groups);
    opts1.setMap(map1);
    
    frameSet = CProxy_FrameSet::ckNew(opts1);
    frameSetProxy = frameSet;    
    
    //***Create Output Chares*****
    output_group = CProxy_OutputG::ckNew();
    outputProxy = output_group;

    map1.sendcountstoOutput(groups);
    
    //***Create Merge Chares*****
    CkArrayOptions opts2(mgroups);
    opts2.setMap(map2);

    merger = CProxy_Merger::ckNew(opts2);
    mergeProxy = merger;
    
    //***Initialize Frame set chares with frame data****  

    for (int i = 0; i < groups; i++) {
      CkVec<Frame> Ckframeset;
      int size2 = frameset_size;

      if (i == groups - 1)
	size2 = frames.size() - i * frameset_size;

      int x = 0;
      
      for (int j = i * frameset_size; j < i * frameset_size + size2; j++) {
	Ckframeset.push_back(*frames[j]);
      }
      
      frameSet[i].initialize(Ckframeset);
    }
    
  }

  void printsettings(){
    CkPrintf("\nRunning Charm++ VPTV Particle Tracking with the following settings:\n\n");
    CkPrintf("input file = %s \n", input_file);
    CkPrintf("output file base name = %s \n", output_file);

    if( CkNumPes() > groups ){
       CkPrintf("WARNING: ****Exiting Program**** The number of processors is greater than frame sets!\n");
       CkPrintf("Choose fewer processors or larger data set.\n #Pes = %d, #FrameSets = %d\n",CkNumPes(),groups);
       CkExit();
    }
    CkPrintf("Number of Particles = %d\n", this->frames[0]->particles.size());
    CkPrintf("Number of Frames = %d\n", this->frames.size());  
    CkPrintf("Number of Processors = %d\n",CkNumPes());
    CkPrintf("Number of Frame sets = %d\n",groups); 
    CkPrintf("Number of frames per set = %d\n\n",FRAMESET_SIZE);

    CkPrintf("BETA (Tracking Threshold) = %f\n",BETA);
    CkPrintf("R_MAX (Max search radius) = %f\n",R_MAX);
    CkPrintf("ALPHA (Search radius multiplier) = %f\n",ALPHA);
   
#if (FINITEDIFF)
    CkPrintf("SEARCH MODE = Finite Difference\n");
#else
    CkPrintf("SEARCH MODE = Regression\n");
#endif   
#if (TESTRESULT)
    CkPrintf("Trajectory Gold Test is ON\n");
#else
    CkPrintf("Trajectory Gold Test is OFF\n");
#endif
#if (OUTPUT)
    CkPrintf("Trajectory output to file is ON\n");
#else
    CkPrintf("Trajectory output to file is OFF\n");
#endif  
  }

  void finishInit() {
#if (VERBOSE_CHARM)
    CkPrintf("\n**** Finished Tracking --> Continuing to merge trajectories\n");
#endif  
    frameSet.merge();
    
    deallocateframes();
   
  }

  void startshortsend() {
#if (VERBOSE_CHARM)
    CkPrintf("\n**** Finished merge --> Continuing to send short trajectories\n");
#endif    
    frameSet.sendShortTrajs();
  }

  void stoptimer(){

   CkPrintf("\n**** Finished building all trajectories---\n");
   CkPrintf("total time = %f sec\n", CmiWallTimer() - time1);
#if (OUTPUT)
   CkPrintf("\n**** Printing trajectories to %d files---\n", CkNumPes());
#else
   CkPrintf("\n**** Printing trajectory stats from %d processors (No Files were written)---\n", CkNumPes());
#endif   
   output_group.fprintTrajs();
  }

  void finished() {
    CkPrintf("\n**** Charm++ VPTV Tracking Program Complete ****\n");
    CkExit();
  }
  
  
   void deallocateframes(){

   for (vector<Frame*>::const_iterator it = this->frames.begin(); it != this->frames.end(); ++it){
   	for (int p = 0; p < (*it)->particles.size(); p++){
    	     delete (*it)->particles[p];
    	}
    	delete (*it);
    }
    this->frames.clear();
    vector<Frame*>().swap(this->frames);
#if (VERBOSE_CHARM)     
    CkPrintf("Frames Deleted: capacity = %d, size = %d\n",this->frames.capacity(), this->frames.size()); 
#endif
    }

};

#include "Tracker.def.h"
