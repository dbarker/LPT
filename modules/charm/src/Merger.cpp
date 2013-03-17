#include "Merger.h"

void Merger::sendMerge(CkVec<Trajectory> trajs, int index) { 

	rec_count++;

	if (index == thisIndex){
		whoisleft = index;
		left_ends = trajs;

		for (int i = 0; i < left_ends.size(); i++)
			left_ends[i].fixParts();
	}
	else if (index > thisIndex){
		whoisright = index;
		right_ends = trajs;

		for (int i = 0; i < right_ends.size(); i++)
			right_ends[i].fixParts();
	}
	else{
		ckerr << "Error: Merge element "<< thisIndex <<" was sent a message by FrameSet element "<< index << endl;
	}

	if (rec_count >= 2){
#if (VERBOSE_CHARM)
		//CkPrintf("right_ends size = %d\n", right_ends.size());
		//CkPrintf("left_ends size = %d\n", left_ends.size());
		CkPrintf("--Merge Element %d: Recieved %d sets -> Merging now\n", thisIndex, rec_count);
#endif
		mergeTrajs4pt();
	}
}

void Merger::mergeTrajs() {
	Cost merge;
	double cost;

	vector<Particle*>::reverse_iterator P1;
	vector<Particle*>::iterator P2;

	//CkPrintf("Merge Element %d: inside mergeTrajs: Left Trajs size = %d, Right Trajs size= %d \n", thisIndex,left_ends.size(),right_ends.size());
	//CkPrintf("inside mergeTrajs: num trajs size T2 = %d \n", right_ends.size());

	for (int i = 0; i < left_ends.size(); i++) {

		Trajectory *T1 = &left_ends[i];
		int imergecount = 0;

		for (int j = 0; j < right_ends.size(); j++) {
			Trajectory *T2 = &right_ends[j];

			if (T1->particles.size() == 0 || T2->particles.size() == 0)
				ckerr << "Error: there is an empty traj in mergTrajs()"<< endl;

			cost = 0;
			P1 = T1->particles.rbegin();
			P2 = T2->particles.begin();

			if ((*P2)->frame_index - (*P1)->frame_index == 1){

				double diffx=((*P1)->x - (*P2)->x);
				double diffy=((*P1)->y - (*P2)->y);
				double diffz=((*P1)->z - (*P2)->z);
				double magXf = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

				if (magXf <= R_MAX) {
					cost = merge.cost(T1, T2);

					if (cost <= BETA){

						imergecount ++;

						if(imergecount == 1){
							TrajsMerged[i] = j;  //add merged trajs to map  :TODO try to impliment with a map <int,CkVec<int>> to allow more than 1 traj
							break;
						}
					}
				}
			}

		}
	}

	finishMerge();
}

void Merger::mergeTrajs4pt() {
	Cost merge;
	double cost;

	vector<Particle*>::reverse_iterator P1;
	vector<Particle*>::iterator P2;

	//CkPrintf("Merge Element %d: inside mergeTrajs: Left Trajs size = %d, Right Trajs size= %d \n", thisIndex,left_ends.size(),right_ends.size());
	//CkPrintf("inside mergeTrajs: num trajs size T2 = %d \n", right_ends.size());

	for (int i = 0; i < left_ends.size(); i++) {

		Trajectory *T1 = &left_ends[i];
		int imergecount = 0;

		for (int j = 0; j < right_ends.size(); j++) {
			Trajectory *T2 = &right_ends[j];

			if (T1->particles.size() == 0 || T2->particles.size() == 0)
				ckerr << "Error: there is an empty traj in mergTrajs()"<< endl;

			cost = 0;
			P1 = T1->particles.rbegin();
			P2 = T2->particles.begin();

			if ((*P2)->frame_index - (*P1)->frame_index == 1){

				double diffx=((*P1)->x - (*P2)->x);
				double diffy=((*P1)->y - (*P2)->y);
				double diffz=((*P1)->z - (*P2)->z);
				double magXf = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

				if (magXf <= R_MAX) {
					for (int u = 0; u < 4; u++){
						cost = merge.cost(T1, (*P2));
						if (cost <= BETA){
							T1->particles.push_back((*P2));
							++P2;
							if (u == 3){
								imergecount ++;
								if(imergecount == 1)
									TrajsMerged[i] = j;  //add merged trajs to map  :TODO try to impliment with a map <int,CkVec<int>> to allow more than 1 traj
							}
						}else{
							for (int z = u; z > 0; z--)
								T1->particles.pop_back();
							break;
						}
					}
				}
			}

		}
	}

	finishMerge();
}

void Merger::mergeTrajs4ptSearch() {
	Cost merge;
	Search S;
	double cost;

	vector<Particle*>::reverse_iterator P1;
	vector<Particle*>::reverse_iterator P2;
	vector<Particle*>::iterator P3c;
	//CkPrintf("Merge Element %d: inside mergeTrajs: Left Trajs size = %d, Right Trajs size= %d \n", thisIndex,left_ends.size(),right_ends.size());
	//CkPrintf("inside mergeTrajs: num trajs size T2 = %d \n", right_ends.size());

	for (int i = 0; i < left_ends.size(); i++) {

		Trajectory *T1 = &left_ends[i];
		int imergecount = 0;

		for (int j = 0; j < right_ends.size(); j++) {
			Trajectory *T2 = &right_ends[j];

			if (T1->particles.size() == 0 || T2->particles.size() == 0){
				ckerr << "Error: there is an empty traj in mergTrajs()"<< endl;
				break;
			}
			cost = 0;
			P2 = T1->particles.rbegin();
			P1 = T1->particles.rbegin() + 1;
			P3c = T2->particles.begin();

			if ( abs((*P3c)->frame_index - (*P2)->frame_index ) == 1){
				Particle P3;
				double ftime = (double)((*P3c)->frame_index)/(double)FRAMESPERSEC;
				double r = S.calcSearchRadius( (*(*P1)), (*(*P2)), ftime);
				int endindex = T1->particles.size() - 1;
				S.motion.extrapolatePosition(T1->particles, endindex, &P3);
				double d = S.calcDistance(P3, *(*P3c));

				if (d <= r) {
					for (int u = 0; u < 4; u++){
						cost = merge.cost(T1, (*P3c));
						if (cost <= BETA){
							T1->particles.push_back((*P3c));
							++P3c;
							if (u == 3){
								imergecount ++;
								if(imergecount == 1)
									TrajsMerged[i] = j;  //add merged trajs to map  :TODO try to impliment with a map <int,CkVec<int>> to allow more than 1 traj
							}
						}else{
							for (int z = u; z > 0; z--)
								T1->particles.pop_back();
							break;
						}
					}
				}

			}

		}
	}
	finishMerge();
}


void Merger::finishMerge() {
#if (VERBOSE_CHARM)
	CkPrintf("--Merge Element %d: Merge Finished with %d merges --> Broadcasting map to output chares\n", thisIndex, TrajsMerged.size());
#endif
	outputProxy.recvMap(TrajsMerged,thisIndex);

	deallocateTrajs(right_ends);
	deallocateTrajs(left_ends);
}

void Merger::deallocateTrajs(CkVec<Trajectory> &trajs){

	for (int i = 0; i < trajs.size(); i++){
		trajs[i].particles.clear();
		vector<Particle*>().swap((trajs[i].particles));

		trajs[i]._particles.clear();
		//delete &(trajs[i]);
	}
	trajs.clear();
}


