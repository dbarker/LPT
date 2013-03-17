#ifndef GOLDTEST_H
#define GOLDTEST_H
#include "interface.h"
#if(EVALRESULT)
  #include "AnalyzeInput.h"
  #include "Output.h"
  #include "Input.h"
#endif
class Goldtest{
 public:
  double E_quality;
  double E_total;
  int N_l;
  int N_tot;
  int N_correct;
  Goldtest();
  void test(vector<Trajectory*> &Trajs, vector<Trajectory*> &TrajGold, int maxframes);
 private:

};
#endif

