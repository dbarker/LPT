#include "OutputG.h"

OutputG::OutputG() {
  map_count = 0;
  finishcount = 0;
  procChareCount = 0;
  allIn = false;
  chareCountrec = false;
  trajcount = 0;
  trajsum = 0;
  int mgroups = groups - 1;
  map <int, int> empty;
  for(int i = 0; i < mgroups; i++)
  	mapbuffer.push_back(empty);
  
  //CkPrintf("Mapbuffer size = %d\n",mapbuffer.size());
}

void OutputG::receivecounts(int cnt){
     procChareCount = cnt;
#if (VERBOSE_CHARM)
     CkPrintf("--Output Element %d: contains %d FrameSet elements\n",CkMyPe(),procChareCount);
#endif
     chareCountrec = true;
     //call to finishmerge function to be sure procChareCount has been set prior to sending short trajs
     finishmerge(-1,0);
}

void OutputG::recvMap(map<int, int> MergedTrajs, int mergeindex) {
  map_count++;
  int mgroups = groups - 1;

  mapbuffer[mergeindex] = MergedTrajs;

  if (map_count >= mgroups){
#if (VERBOSE_CHARM)
   CkPrintf("--Output Element %d: %d/%d maps received --> building global trajectory map\n", CkMyPe(), mapbuffer.size() ,mgroups);
#endif
   buildfullmap();
  }
}

void OutputG::buildfullmap() {
  
  for(int i = 0; i < mapbuffer.size(); i++){

     map <int, int>::iterator itmap; 

     for(int j = 0; j < fullMap.size(); j++){  //fullMap should be empty when i = 0

         int TrajKey = fullMap[j][i];
         itmap = mapbuffer[i].find(TrajKey);   
         
         if (itmap != mapbuffer[i].end()){
           fullMap[j][i + 1] = itmap->second;   
           mapbuffer[i].erase(itmap); 
         }
     }
     //--Initialize fullMap or if key was not found above then add as new trajectory to fullMap
     for(map<int,int>::iterator itadd = mapbuffer[i].begin(); itadd != mapbuffer[i].end(); ++itadd){
        vector <int> TrajAdd(groups,-1);
    	TrajAdd[i] = itadd->first;
    	TrajAdd[i + 1] = itadd->second;

    	fullMap.push_back(TrajAdd);
     }   
  }
  
 mapbuffer.clear();
 vector< map <int,int> >().swap(mapbuffer);
 
 sendTrajRequest();

#if (VERBOSE_CHARM)
  CkPrintf("--Output Element %d: Global map construction complete -> total trajs = %d \n", CkMyPe(), fullMap.size());
#endif
fullMap.clear();
vector<vector <int> >().swap(fullMap);

}
void OutputG::sendTrajRequest(){
   
   int p = CkMyPe();
   int procs = CkNumPes();
   int Tsetsize = fullMap.size() / procs;
   
   if (Tsetsize >= 1){
    	mytsize = Tsetsize;
   	if (p == procs  - 1)
      		mytsize = fullMap.size() - p * Tsetsize;	
   }else{
   	if (p < fullMap.size())
        	mytsize = 1;
        else
        	mytsize = 0;        	
   }
   
   if (mytsize > 0){  
   	mytstart = p * Tsetsize;
   	mytend = p * Tsetsize + mytsize;
	   for(int i = 0; i < groups; i++){
		CkVec<int> localIDs;
		CkVec<int> globalIDs;

		for(int j = mytstart; j < mytend; j++){  //FIXME: Problem here?
		    int localid = fullMap[j][i];
		    
		    if (localid >= 0){  //if less than zero then no mapping made (remained initial value)
		     	localIDs.push_back(localid);
		     	globalIDs.push_back(j); 
		  }

		}
		if (localIDs.size() > 0){
		    addtocheckinMap(i);
		    frameSetProxy[i].output(localIDs, globalIDs, p);
	   	}
	  }
	  list<Trajectory> empty;
	  for(int i = 0; i < mytsize; i++)
	  	builtTrajsLst.push_back(empty);  
     }else{
     	contribute(CkCallback(CkIndex_Main::startshortsend(), mainProxy));
     }
}

void OutputG::addTrajToList(Trajectory traj, int fromIndex) {
  
  int listId = traj.id - mytstart;
  
  if (builtTrajsLst[listId].empty() == true){
        builtTrajsLst[listId].push_back(traj);    
  }else{
    
      for(list<Trajectory>::iterator it = builtTrajsLst[listId].begin(); it != builtTrajsLst[listId].end(); ++it){
        list<Trajectory>::iterator next = it; 
 	++next;
        if (traj.startframe < it->startframe && it == builtTrajsLst[listId].begin()){//FIXME:  This section needs dehacked.
	   builtTrajsLst[listId].insert(it,traj);
           break;
	}
        else if (next == builtTrajsLst[listId].end() ){
           if(traj.startframe < it->startframe)
           	builtTrajsLst[listId].insert(it,traj); 
	   else
		builtTrajsLst[listId].insert(next,traj);
           break;
        }
        else if (traj.startframe > it->startframe && traj.startframe < next->startframe){
           builtTrajsLst[listId].insert(next,traj);
           break;
        }
            
      }
  }
  trajcount++;
  //builtTrajs[traj.id].fixParts();

  //call to finishmerge function if all chares have checked in before trajs are finished processing
  if (allIn == true)
	finishmerge(-1, 0);
     
}

void OutputG::addShortTrajs(Trajectory traj) {
  list <Trajectory> shortie(1,traj);
  builtTrajsLst.push_back(shortie);
  trajcount++;
  if (allIn == true)
     finished(0);
}
  
void OutputG::finishmerge(int chareID, int count){
  if( chareID >= 0 ){ //chareID = -1 if called from addTrajToList()
    checkinMap[chareID] = true; //TODO: Add check to ensure this does not increment the map size
    trajsum += count;
    allIn = true;
    
    for(map<int,bool>::iterator mapit = checkinMap.begin(); mapit != checkinMap.end(); ++mapit){
       if (mapit->second == false){
	  allIn = false;
          break;
       }
    }
  }

  if(allIn == true && trajcount == trajsum && chareCountrec == true){ 
     //***reset checkin varables****
     trajsum = 0;
     trajcount = 0;
     checkinMap.clear();
     allIn = false;

     //***Callback to main proxy and start sending short trajectories*****
     contribute(CkCallback(CkIndex_Main::startshortsend(), mainProxy));
  }
}

void OutputG::addtocheckinMap(int chareID){
  checkinMap[chareID] = false;
}

void OutputG::finished(int count){
  if (allIn == false){
  	finishcount++; 
  	trajsum += count;
  }
  if(finishcount == procChareCount){
     allIn = true;
     if(trajcount == trajsum){
#if (VERBOSE_CHARM)     
       CkPrintf("--Ouput Element %d Finished: %d short trajectories received\n",CkMyPe(),trajcount);
#endif
       contribute(CkCallback(CkIndex_Main::stoptimer(), mainProxy));
     }
  }
}

void OutputG::printTraj(Trajectory& t) {
  for (int i = 0; i < t.particles.size(); i++)
    CkPrintf("(%.7e, %.7e, %.7e)\n", t.particles[i]->x, t.particles[i]->y, 
	     t.particles[i]->z);
}

void OutputG::fprintTrajs () {
#if(OUTPUT)
  char MyOutputFile[100];
  sprintf(MyOutputFile,"%s%d",output_file, CkMyPe());

  FILE* printFile = fopen (MyOutputFile,"w");
  
  if(!printFile){
    ckerr << "Cannot find/open the output file \n" << MyOutputFile << "\n Program exiting (OutputG.cpp)" << endl;
    CkExit();
  }
  
  for (unsigned int i = 0; i < builtTrajsLst.size(); i++){
    for (list<Trajectory>::iterator traj = builtTrajsLst[i].begin(); traj != builtTrajsLst[i].end(); ++traj){
      for (unsigned int k = 0; k < traj->_particles.size(); k++){
         /*if (count == 0)
         	fprintf(printFile, "%d\n", traj->startframe);*/         //TODO: Add startframe to output file

         fprintf(printFile, "%d\t%d\t%25.16e%25.16e%25.16e\n", traj->_particles[k].id, traj->_particles[k].frame_index,
              traj->_particles[k].x, traj->_particles[k].y, traj->_particles[k].z);
      }
      
    }
      fprintf(printFile, "\n");
  }
  fclose(printFile);
  CkPrintf("--Ouput Element %d Finished printing %d Trajectories to file: %s\n",CkMyPe(),builtTrajsLst.size(),MyOutputFile);
#else
  CkPrintf("--Ouput Element %d Finished with %d Trajectories --\n",CkMyPe(),builtTrajsLst.size());
#endif
  //---Reduction to END THE PROGRAM----
  contribute(CkCallback(CkIndex_Main::finished(), mainProxy));

}
