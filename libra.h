#pragma once
#include "NuFuncs.h"
#include <random>

class libra{
 public:

  NuFuncs* dfN1 = nullptr;
  NuFuncs* dfN4 = nullptr;
  NuFuncs* dfN17 = nullptr;
  NuFuncs* dfI1 = nullptr;
  NuFuncs* dfI4 = nullptr;
  NuFuncs* dfI17 = nullptr;
  
  
  ROOT::Math::Interpolator* oneND = nullptr;
  TGraph2D* twoND = nullptr;

  ROOT::Math::Interpolator* oneID = nullptr;
  TGraph2D* twoID = nullptr;

  double T23min;
  double T23max;
  
  double DCPmin;
  double DCPmax;
  
  double DMAminN;
  double DMAmaxN;
  double DMAminI;
  double DMAmaxI;
  
  double T13min;
  double T13max;
  
  double T12min;
  double T12max;
  
  double DMSmin;
  double DMSmax;
  
  double sDM221 = 1e-6;
  double sDM231 = 3e-6;
  double sT23 = 0.0095;
  double sT12 = 0.005;
  double sT13 = 0.0015;
  double sDCP = 2.5;
  double currentChi = 0,propChi = 0, U =0, check = 0;

  double TmpDCP = 0;
  double add = 180;
  double ipart = 0;
  double* ptr = &ipart;


  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<double> udist;

  libra();
  
  struct MCChain {
  public:
    double current[8];
    double proposal[8];
    unsigned long long int acceptCount, totCount;
    
    MCChain(double* currptr)
      : acceptCount(1), totCount(1)
    {
      for (int i = 0; i < 8; ++i) {
	current[i] = currptr[i];
      }
    }

    double getRatio(){
      if(totCount == 0) return 0;
      return static_cast<double>(acceptCount)/totCount;
    }
  };
   
  double getChiSquare(double,double,double,double,double,double);
  void proposalFunc(MCChain*);
  double MH(MCChain*);


  std::vector<MCChain*> MCS;

  std::vector<double> chis;

  double Uchain = 0;
  double Ucheck = 0;

  void initOnce();
  void MonteCarlo();
  
  
  ~libra();
  
  

};
