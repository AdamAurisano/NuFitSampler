#include "libra.h"

static constexpr double E = 1e-5;

libra::libra(){


  dfN1 = new NuFuncs(0,"N1");
  dfN4 = new NuFuncs(0,"N4");
  dfN17 = new NuFuncs(0,"N17");

  dfI1 = new NuFuncs(1,"I1");
  dfI4 = new NuFuncs(1,"I4");
  dfI17 = new NuFuncs(1,"I17");
  
  dfN1->putIntoDF("normalOrder.txt",1);
  dfN4->putIntoDF("normalOrder.txt",4);
  dfN17->putIntoDF("normalOrder.txt",17);
    
  oneND = dfN17->oneDI();
  twoND = dfN4->twoDI();
  dfN1->kDTree();
  
  dfI1->putIntoDF("invertedOrder.txt",1);
  dfI4->putIntoDF("invertedOrder.txt",4);
  dfI17->putIntoDF("invertedOrder.txt",17);
    
  oneID = dfI17->oneDI();
  twoID = dfI4->twoDI();
  dfI1->kDTree();

  int threeSizeN = dfN1->xData.size();
  int twoSize = dfN4->xData.size();
  int oneSize = dfN17->xData.size();

  int threeSizeI = dfI1->xData.size();
    
  T23min = dfN1->xData[0];
  T23max = dfN1->xData[threeSizeN-1];
  
  DCPmin = -180;
  DCPmax = 180;
  
  DMAminN = dfN1->yData[0];
  DMAmaxN = dfN1->yData[threeSizeN-1];
  DMAminI = dfI1->yData[0];
  DMAmaxI = dfI1->yData[threeSizeI-1];
  

  T13min = dfN17->xData[0];
  T13max = dfN17->xData[oneSize-1];
  
  T12min = dfN4->xData[0];
  T12max = dfN4->xData[twoSize-1];
  
  DMSmin = dfN4->yData[0];
  DMSmax = dfN4->yData[twoSize-1];

  std::mt19937 genTmp(rd());
  gen = genTmp;
  std::uniform_real_distribution<double> udistTmp(0,1);
  udist = udistTmp;
  
}

double snap(double val, double min, double max){
  if(std::abs(val-min) <= E) return min;
  if(std::abs(max-val) <= E) return max;
  return val;
}


double libra::getChiSquare(double dm221,double dm231,double t12,double t13,double t23,double dcpval){ //assumed order

  double dm231S = dm231 * 1000;
  double dm232S = (dm231 - dm221)*1000;
  double dm221S = log10(dm221);
  
  bool dm221Flag = ((dm221S > DMSmin - E) && (dm221S < DMSmax + E));
  
  bool dm231Flag = ((dm231S > DMAminN - E) && (dm231S < DMAmaxN + E)); // && (dm231S > 0)
  bool dm232Flag = ((dm232S > DMAminI - E) && (dm232S < DMAmaxI + E)); //&& (dm232S < 0)
  
  bool t12Flag = ((t12 > T12min - E) && (t12 < T12max + E));
  bool t13Flag = ((t13 > T13min - E) && (t13 < T13max + E));
  bool t23Flag = ((t23 > T23min - E) && (t23 < T23max + E));

  bool dcpFlag = ((dcpval > DCPmin - E) && (dcpval < DCPmax + E));

  bool NFlag = (t23Flag && dm231Flag && dcpFlag && t12Flag && dm221Flag && t13Flag);
  bool IFlag = (t23Flag && dm232Flag && dcpFlag && t12Flag && dm221Flag && t13Flag);
  
  if(NFlag){

    dm231S=snap(dm231S,DMAminN,DMAmaxN);
    t23=snap(t23,T23min,T23max);
    t12=snap(t12,T12min,T12max);
    t13=snap(t13,T13min,T13max);
    dm221S=snap(dm221S,DMSmin,DMSmax);
    
    double dpthree[3] = {t23,dm231S,dcpval};    
    return  minN + (dfN1->threeDI(dpthree) - minN) + (twoND->Interpolate(t12,dm221S)-minN) + (oneND->Eval(t13) - minN);
    
  }

  if(IFlag){

     dm231S=snap(dm231S,DMAminN,DMAmaxN);
    t23=snap(t23,T23min,T23max);
    t12=snap(t12,T12min,T12max);
    t13=snap(t13,T13min,T13max);
    dm221S=snap(dm221S,DMSmin,DMSmax);
    
    double dpthree[3] = {t23,dm232S,dcpval};
    //std::cout << dm232S << "\n";
    return minI + (dfI1->threeDI(dpthree) - minI) + (twoID->Interpolate(t12,dm221S)-minI) + (oneID->Eval(t13) - minI);
  }

  return -999999;

  
}

void libra::proposalFunc(MCChain* chain){
  //double dm221,double dm231,double t12,double t13,double t23,double dcpval

  //{7.5e-5,2.5e-3,0.3,0.022,0.5,0}
  
  std::normal_distribution<double> dm221R(chain->current[0],sDM221);
  std::normal_distribution<double> dm231R(chain->current[1],sDM231);
  std::normal_distribution<double> T23R(chain->current[4],sT23);
  std::normal_distribution<double> T12R(chain->current[2],sT12);
  std::normal_distribution<double> T13R(chain->current[3],sT13);
  std::normal_distribution<double> DCPR(chain->current[5],sDCP);
  
  chain->proposal[0] = dm221R(gen);
  chain->proposal[1] = dm231R(gen);
  chain->proposal[4] = T23R(gen);
  chain->proposal[2] = T12R(gen);
  chain->proposal[3] = T13R(gen);

  TmpDCP = DCPR(gen);
  double add = 180;
  if(TmpDCP < 0) add = -1*add;
  modf((TmpDCP+add)/360.0, ptr);
  chain->proposal[5] = TmpDCP - (360*ipart);

  //std::cout << TmpDCP << " " << prop[5] << "\n";
  
}

double libra::MH(MCChain* chain){
  chain->totCount++;
  
  currentChi = getChiSquare(chain->current[0],chain->current[1],chain->current[2],chain->current[3],chain->current[4],chain->current[5]);
  propChi = getChiSquare(chain->proposal[0],chain->proposal[1],chain->proposal[2],chain->proposal[3],chain->proposal[4],chain->proposal[5]);

  if(propChi < 0) return currentChi;
  
  check = std::exp(-propChi + currentChi);
  U = udist(gen);

  if(U < check){
    
    for(size_t i = 0; i < 6; i++){
      chain->current[i] = chain->proposal[i];
    }
    chain->acceptCount++;


    return propChi;
  }
  //std::cout << "check\n";
  return currentChi; 
}

void libra::initOnce(){
  //MCS.resize(8);
  chis.resize(8);


  //{7.5e-5,2.5e-3,0.3,0.022,0.5,0}
  
  double posNOp90list[6] = {7.5e-5, 2.5e-3, 0.30, 0.022, 0.54, 90};
  double negNOp90list[6] = {7.5e-5, 2.5e-3, 0.30, 0.022, 0.46, 90};
  double posIOp90list[6] = {7.5e-5, -2.5e-3, 0.30, 0.022,0.54, 90};
  double negIOp90list[6] = {7.5e-5, -2.5e-3, 0.30, 0.022,0.46,90};
  double posNOn90list[6] = {7.5e-5, 2.5e-3, 0.30, 0.022,0.54,-90};
  double negNOn90list[6] = {7.5e-5, 2.5e-3, 0.30, 0.022,0.46,-90};
  double posIOn90list[6] = {7.5e-5, -2.5e-3, 0.30, 0.022,0.54,-90};
  double negIOn90list[6] = {7.5e-5, -2.5e-3, 0.30, 0.022,0.46,-90};

  //std::cout << getChiSquare(posNOp90list[0],posNOp90list[1],posNOp90list[2],posNOp90list[3],posNOp90list[4],posNOp90list[5]);

  /*MCChain posNOp90(posNOp90list),negNOp90(negNOp90list),posIOp90(posIOp90list),
    negIOp90(negIOp90list),posNOn90(posNOn90list),negNOn90(negNOn90list),
    posIOn90(posIOn90list),negIOn90(negIOn90list);*/

  MCS.push_back(new MCChain(posNOp90list));
  MCS.push_back(new MCChain(negNOp90list));
  MCS.push_back(new MCChain(posIOp90list));
  MCS.push_back(new MCChain(negIOp90list));
  MCS.push_back(new MCChain(posNOn90list));
  MCS.push_back(new MCChain(negNOn90list));
  MCS.push_back(new MCChain(posIOn90list));
  MCS.push_back(new MCChain(negIOn90list));

}

void libra::MonteCarlo(){
  
  //std::cout << MCS.size() << "\n";
  
  for(size_t i = 0; i < 8; i++){
    proposalFunc(MCS[i]);
    chis[i] = MH(MCS[i]);
  }

  //std::cout << "check1\n";


  if(chis[0] < 0 && chis[1] < 0 && chis[2] < 0 && chis[3] < 0 && chis[4] < 0 && chis[5] < 0 && chis[6] < 0 && chis[7] < 0) {
    std::cout << chis[0] << chis[1] << chis[2] << chis[3] << chis[4] << chis[5] << chis[6] << chis[7] << "\n";
    return;
  }

  //std::cout << "check2\n";

  for(size_t i = 7; i > 0; i--){
    if(chis[i] < 0 || chis[0] < 0) continue;

    Ucheck = std::exp(-chis[i] + chis[0]);
    Uchain = udist(gen);

    if(Uchain < Ucheck){
      std::swap(MCS[0],MCS[i]);      
    }
  }

  //std::cout << "check3\n";

  return;
  
}

libra::~libra(){
    delete oneID;
    delete twoID;
    delete oneND;
    delete twoND;
    delete dfN1;
    delete dfN4;
    delete dfN17;
    delete dfI1;
    delete dfI4;
    delete dfI17;

    for(size_t i = 0; i < MCS.size(); i++){
      delete MCS[i];
    }
}


