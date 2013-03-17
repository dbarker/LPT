#ifndef INTERFACE_H
#define INTERFACE_H

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <utility>

#if(CHARM) 
#include <charm++.h>
#endif

//SEARCH and COST parameters
#define ALPHA 0.8f  	  //0.8 for pivSTD352       //0.8 for CFD 1000p4096f100fps         //1.0 for scaling data sets
#define BETA  0.03f      //0.13 for pivSTD352      //0.03  for CFD 1000p4096f100fps       //0.01 for scaling data sets 
#define R_MAX 0.042f       //0.13 for pivSTD352      //0.042 for CFD 1000p4096f100fps	   //50.0f for scaling data sets
#define FINITEDIFF 1      //1 enables finite difference based search, 0 for regression based search

#define DIM 3             //Dimensions to track FIXME if 2d tracking is desired then Search.cpp, Regression.cpp and Cost.cpp need modification
#define FRAMESPERSEC 100

#define N_TRUE 6         //sets the min length of a trajectory to remain active
#define N_GAP 0		 //sets the max number of consectutive gaps in a trajectory

#define NUMLINKS 10   //sets the max number of candidates that can be found per trajectory..needs to be at least 12 for simdata2000
#define SETSIZE 4        //sets the max number of previous points used in candiate search 

//Charm++ settings
#define FRAMESET_SIZE 64 //sets the number of frames per chare in charm++ ***MINIMUM OF AT LEAST FRAMESET_SIZE == MERGE_COST_SIZE
#define BLOCKMAP 0       //1 for Blocked chare distribution and 0 for round robin 

//Output settings
#define VERBOSE_SERIAL 0  //0 for no output on serial tracking functions
#define VERBOSE_CHARM 0   //0 for no output on charm++ chare functions
#define OUTPUT 1          //0 for no trajectory output to file after program completion
#define EVALRESULT 0      //0 during timing runs to not include extra result evaluation data structures

//COST parameters
#define COST_TRAJ_SIZE 4  //Number of frames required to evaluate the cost of adding a candiate to a trajectory
#define MERGE_COST_SIZE 8  //Number of points used in merge evaluation (must be divisible by 2)

using namespace std;

class Particle {
public:
  double x, y, z;
  int id, frame_index;
#if(EVALRESULT)
  double vel, accel, cost, dr, ds;
#endif

#if(CHARM) 
  void pup(PUP::er &p) {
    p | x;
    p | y;
    p | z;
    p | frame_index;
    p | id;
  }
#endif
};

class Frame {
public:
  vector<Particle*> particles; 
  int frame_index;
  double time;
#if(EVALRESULT)
  map<int,double> minmap;
  double S;
#endif
  /*//SMAP
  //map for spatial mapping of particles
  map<pair<pair<int,int>,int>, vector<Particle*> > bins;
  vector<vector<vector<vector<Particle*> > > > vbins;
  */
#if(CHARM)    
  CkVec<Particle> _particles;
  void pup(PUP::er &p) {
    if (p.isSizing()) {
      _particles.clear();
           
      for (int i = 0; i < particles.size(); i++)
	_particles.push_back(*particles[i]);

   }

    p | frame_index;
    p | time;
    p | _particles;

    if (p.isUnpacking()) {
      particles.clear();
     
      for (int i = 0; i < _particles.size(); i++)
        particles.push_back(&_particles[i]);

    }
  }
#endif
};

class Trajectory {
public:
  vector<Particle*> particles;
  int id, gap, startframe;
  vector<bool> matches;
#if(EVALRESULT)
  double R,Rh,Rl,C,Ch,Cl,Eff,S,Sl,Sh;  
#endif

#if(CHARM) 
  CkVec<Particle> _particles;

  Trajectory *trim(bool end) {
    Trajectory *t = new Trajectory();
    t->id = this->id;
    t->gap = this->gap;
    t->startframe = this->startframe;

    // assuming that particles it at least MERGE_COST_SIZE/2 of size
    if (end) {
      for (int i = particles.size() - MERGE_COST_SIZE/2; i < particles.size(); i++)
	t->particles.push_back(this->particles[i]);
    } else {
      for (int i = 0; i < MERGE_COST_SIZE/2; i++)
	t->particles.push_back(this->particles[i]);
    }
    
    return t;
  }

  // TODO: JL push this into the copy constructor
  // TODO:  is this fixParts() needed? or is this taken care of in the pup routine below?

  void fixParts() {
    particles.clear();
    for (int i = 0; i < _particles.size(); i++) {
      particles.push_back(&_particles[i]);
    }
  }

  void pup(PUP::er &p) { 
    if (p.isSizing()) {
      _particles.clear();

      for (int i = 0; i < particles.size(); i++)
	_particles.push_back(*particles[i]);
    }

    p | id;
    p | gap;
    p | _particles;
    p | startframe;

    if (p.isUnpacking()) {
      particles.clear();

      for (int i = 0; i < _particles.size(); i++)
	particles.push_back(&_particles[i]);
    }
  }
#endif
};


#endif
