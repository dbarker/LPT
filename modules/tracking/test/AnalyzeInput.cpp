/*
 * AnalyzeInput.cpp
 *
 *  Created on: Nov 18, 2010
 *      Author: dbarker2
 */

#include "AnalyzeInput.h"

AnalyzeInput::AnalyzeInput() {
	// TODO Auto-generated constructor stub

}

double AnalyzeInput::distance(Particle* P1, Particle* P2){
	double diffx, diffy, diffz, d;

	diffx = (P1->x - P2->x);
    diffy = (P1->y - P2->y);
	diffz = (P1->z - P2->z);

	d = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
	return d;
}

#if(EVALRESULT)
void AnalyzeInput::evalframes(vector<Frame*> &frames){
//********Calculate particle spacing in frames********************
    this->DS = 0;
	for (int i = 0; i < frames.size(); i++){
		Frame* F = frames[i];
		F->S = 0;
		for (int p1 = 0; p1 < F->particles.size(); p1++){
			Particle* P1 = F->particles[p1];
			double dmin = 1E10;
			for (int p2 = 0; p2 < F->particles.size(); p2++){
				Particle* P2 = F->particles[p2];
				if (P1 != P2){
					double d = distance(P1, P2);
					if (d < dmin){
						dmin = d;
					}
				}
			}
			P1->ds = dmin;
			F->minmap[P1->id] = dmin;
			F->S += dmin;
		}
		F->S = F->S / (double)(F->particles.size());
		this->DS+= F->S;
	}

	this->DS = this->DS / frames.size();
}

void AnalyzeInput::evaltrajs(vector<Trajectory*> &trajs){
 //******Calculate Cost, Acceleration, and Distance traveled per time step*******************************
	Cost evalcost;
	this->DR = 0;
	this->C = 0;

	for (int i = 0; i < trajs.size(); i++){
		trajs[i]->R = 0;
		trajs[i]->C = 0;
		double Cmin = 100;
		double Cmax = 0;
		double Rmin = 1E9;
		double Rmax = 0;
		for (int j = 0; j < trajs[i]->particles.size()-1; j++){
			Particle* P = trajs[i]->particles[j];
			Particle* Pn = trajs[i]->particles[j+1];
			P->dr = distance(P,Pn);
			if(P->dr > Rmax)
				Rmax = P->dr;

			if(P->dr < Rmin)
				Rmin = P->dr;

            		trajs[i]->R += P->dr;

            		Pn->vel = P->dr * FRAMESPERSEC;
			if (j>1){
				Particle* Pm = trajs[i]->particles[j-1];
				P->accel = abs(Pn->vel - Pm->vel) * FRAMESPERSEC/2.0;
			}

			if (j >= COST_TRAJ_SIZE - 1){
				int s = (j) - (COST_TRAJ_SIZE -1);
				Trajectory tr;
				for (int x = s; x <= j; x++){
					tr.particles.push_back(trajs[i]->particles[x]);
				}
				Pn->cost = evalcost.cost(&tr, Pn);
			}else{
				Pn->cost = 0.0;
			}

			if(Pn->cost > Cmax)
				Cmax = Pn->cost;

			if(Pn->cost < Cmin && Pn->cost != 0)
				Cmin = Pn->cost;

			trajs[i]->C += Pn->cost;
		}
		trajs[i]->Cl = Cmin;
		trajs[i]->Ch = Cmax;
		trajs[i]->C = trajs[i]->C / ((double) trajs[i]->particles.size());
		this->C += trajs[i]->C;

		trajs[i]->Rl = Rmin;
		trajs[i]->Rh = Rmax;
		trajs[i]->R = trajs[i]->R / ((double) trajs[i]->particles.size());
		this->DR += trajs[i]->R;
	}
    	this->C = this->C / trajs.size();
	this->DR = this->DR / trajs.size();
}

void AnalyzeInput::evaltrajresult(vector<Trajectory*> &trajs,vector<Frame*> &frames){
 //******Calculate %tracked, Distance to closest particle per frame, and Distance traveled per time step*******************************
	Cost evalcost;
	this->DR = 0;
	this->C = 0;

	for (int i = 0; i < trajs.size(); i++){
		trajs[i]->R = 0;
		trajs[i]->S = 0;
		double Rmin = 1E9;
		double Rmax = 0;
		double Smin = 1E9;
		double Smax = 0;
		for (int j = 0; j < trajs[i]->particles.size()-1; j++){
			Particle* P = trajs[i]->particles[j];
			Particle* Pn = trajs[i]->particles[j+1];
			
			int f = P->frame_index;
			map<int,double> tmap = frames[f]->minmap;
			map<int,double>::iterator itmap;
			
			itmap = tmap.find(P->id);
			       
         		if (itmap != tmap.end()){
				if (itmap->second < Smin)
					Smin = itmap->second;
				if (itmap->second > Smax)
					Smax = itmap->second;
					
				trajs[i]->S += itmap->second;
			}
								
			
			P->dr = distance(P,Pn);
			if(P->dr > Rmax)
				Rmax = P->dr;

			if(P->dr < Rmin)
				Rmin = P->dr;

            		trajs[i]->R += P->dr;
            			
		}
		trajs[i]->Sl = Smin;
		trajs[i]->Sh = Smax;
		trajs[i]->S = trajs[i]->S / ((double) trajs[i]->particles.size());
		
		trajs[i]->Rl = Rmin;
		trajs[i]->Rh = Rmax;
		trajs[i]->R = trajs[i]->R / ((double) trajs[i]->particles.size());
	}
    	
}

#endif
AnalyzeInput::~AnalyzeInput() {
	// TODO Auto-generated destructor stub
}
