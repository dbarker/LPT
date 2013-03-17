
#include "Goldtest.h"

Goldtest::Goldtest()
{
  E_quality = 0.0;
  E_total = 0.0;
  N_l = 0;
  N_tot = 0;
  N_correct = 0;
}
void Goldtest::test(vector<Trajectory*> &Trajs, vector<Trajectory*> &TrajGold, int maxframes){

	printf("\n\n-----Running gold test------\n");

	for(int i = 0; i < TrajGold.size(); i++){
		vector<bool> temp(TrajGold[i]->particles.size(),0);
		TrajGold[i]->matches = temp;

	}


	printf("maxframes = %d\n", maxframes);

	for (int f = 0; f < maxframes - 1; f++){                  //FIXME Assumes that no gaps exist in the trajectorys which start at Traj->startframe
	   for (int i = 0; i < TrajGold.size(); i++){
		
		int sg = TrajGold[i]->particles[0]->frame_index;
		if (sg <= f   &&  f + 1 - sg < (int)(TrajGold[i]->particles.size()) ){
			Particle* P_g1 = TrajGold[i]->particles[f - sg];
			Particle* P_g2 = TrajGold[i]->particles[f + 1 - sg];
				
		   	for (int j = 0; j < Trajs.size(); j++){
				int st = Trajs[j]->particles[0]->frame_index;
				if (st <= f  &&  f + 1 - st < (int)(Trajs[j]->particles.size()) ){
					Particle* P_a1 = Trajs[j]->particles[f - st];
					Particle* P_a2 = Trajs[j]->particles[f + 1 - st];
		   	     
	 	            		if ( P_a1->id == P_g1->id  &&  P_a2->id == P_g2->id ){
	 	            			TrajGold[i]->matches[f - sg] = 1;
	 	            		}
		        	}
		    	}
	     	}
	   }
	}

	for (unsigned int j = 0; j < Trajs.size(); j++)
	    N_l += Trajs[j]->particles.size()-1; 
	   
	for (unsigned int i = 0; i < TrajGold.size(); i++){
	    if (TrajGold[i]->particles.size() >= N_TRUE){
	    	N_tot += TrajGold[i]->particles.size()-1;
	    	int temp = 0;
	    	for (unsigned int p = 0; p < TrajGold[i]->particles.size(); p++){ 
	    		temp += TrajGold[i]->matches[p];
	    	}
	      	N_correct += temp;
#if(EVALRESULT)	
	      	TrajGold[i]->Eff = (double)temp / (double)(TrajGold[i]->particles.size() - 1);
#endif
	    }
	}
#if(EVALRESULT)
	AnalyzeInput eval;
	Input in;
	Output out;
	char trajstat[100];
	sprintf(trajstat,"trajstats");
	char framestat[100];
	sprintf(framestat,"framestats");
	vector<Frame*> frames = in.trajs2frames(TrajGold);
	eval.evalframes(frames);
	eval.evaltrajresult(TrajGold,frames);
	
	out.fprintTrajsStats(trajstat, TrajGold);
	out.fprintFramesStats(framestat, frames);
	
	char results[100];
	sprintf(results,"testresults");
        
	FILE* printFile = fopen (results,"w");
	// Print out the gold trajectories with match indicators to show which trajs could be correctly reconstructed  
	if(printFile){  
	  for (unsigned int i = 0; i < TrajGold.size(); i++){
	    if(TrajGold[i]->particles.size() > 1){
	      for (unsigned int p = 0; p < TrajGold[i]->particles.size()-1; p++){
	    	Particle* PD = TrajGold[i]->particles[p];
		 fprintf(printFile, "%d\t%d\t%f\t%d\n", PD->id,PD->frame_index,(PD->ds)/(PD->dr),(int)TrajGold[i]->matches[p]);
		 /*if (TrajGold[i]->matches[p])
		 	fprintf(printFile, "\n");
		 else
		 	fprintf(printFile, "\t%f\t%f\n",TrajGold[i]->particles[p]->dr, TrajGold[i]->particles[p]->cost);*/
	      }
	      fprintf(printFile, "\n");
	    }
	  }  
	}
	fclose(printFile);
#endif	
	
	printf("Results:\n\n");

	printf("Total number of Trajectories Tracked:     = %d \n", Trajs.size());
	printf("Total number of Gold Trajectories:        = %d \n\n", TrajGold.size());

	printf("Average Trajectory Length tracked:        = %d frames\n", (int) ((double)N_l / (double)Trajs.size() + 1));
	printf("Average Trajectory Length (Goldtrajs)     = %d frames\n\n", (int) ((double)N_tot/(double)TrajGold.size() + 1));

	printf("Total number of links made:         Nlinked    = %d \n", N_l);
	printf("Total number of gold links:         Ntot  = %d \n", N_tot);
	printf("Total number of correct links:      Ncorrect = %d \n", N_correct);
 

	E_total = 100. * (double)N_correct / (double)N_tot;
	E_quality = 100. * (double)N_correct / (double)N_l;

	printf("\n");
	printf("Coverage Tracking Efficiency (Ncorrect / Ntot)    = %.3f percent\n", E_total);
	printf("Quality Tracking Efficiency (Ncorrect / Nlinked) = %.3f percent\n", E_quality); 
		  
}


