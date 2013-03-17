/*
 * AnalyzeInput.h
 *
 *  Created on: Nov 18, 2010
 *      Author: dbarker2
 */

#ifndef ANALYZEINPUT_H_
#define ANALYZEINPUT_H_
#include "interface.h"
#include "Cost.h"

class AnalyzeInput {
public:
    double C, DR, DS;
    map<int, Trajectory*> trajmap;
	AnalyzeInput();
	double distance(Particle* P1, Particle* P2);
	void evalframes(vector<Frame*> &frames);
	void evaltrajs(vector<Trajectory*> &trajs);
	void evaltrajresult(vector<Trajectory*> &trajs,vector<Frame*> &frames);
	virtual ~AnalyzeInput();
};

#endif /* ANALYZEINPUT_H_ */
