#include <iostream>
#include <fstream>
#include <iomanip>
#include <Math/Interpolator.h>
#include <Math/InterpolationTypes.h>
#include <TMultiDimFit.h>
#include <TH3.h>
#include "NuFuncs.h"
#include <TCanvas.h>
#include <TGraph2D.h>
#include <algorithm>
#include <cmath>

NuFuncs::NuFuncs() {}

void NuFuncs::putIntoDF(std::string file_addy, int secnum) {
  std::ifstream datafile(file_addy);
  if (!datafile) {
    std::cerr << "Error opening file" << std::endl;
    return;
  }
  
  std::string line;
  int section = 0;
  std::string sub = "projection:";
  while (getline(datafile, line)) {
    if (line.empty()) continue;
    
    size_t start = line.find_first_not_of(" ");
    size_t end = line.find_last_not_of(" ");
    if (start != std::string::npos && end != std::string::npos) {
      line = line.substr(start, end - start + 1);
    }
    size_t pos;
    while ((pos = line.find("  ")) != std::string::npos) {
      line.erase(pos, 1);
    }
    
    int col_length;
    if (line[0] == '#' && line[1] != '#') {
      section++;
      
      if (section == secnum) {
	int pos = line.find(sub);
	if (pos != std::string::npos) {
	  line = line.substr(1, pos - 1);
	  size_t col_pos = 0;
	  size_t prev_col_pos = 0;
	  column_names.clear();
	  while ((col_pos = line.find("/", prev_col_pos)) != std::string::npos) {
	    column_names.push_back(line.substr(prev_col_pos, col_pos - prev_col_pos));
	    prev_col_pos = col_pos + 1;
 
	  }
	  column_names.push_back(line.substr(prev_col_pos));
	  col_length = column_names.size();

	  //std::cout << column_names[0] << column_names[1] << std::endl;
	}
      }
      
      continue;
    }
    
    if (section == secnum) {
      
      std::cout << std::setprecision(8);
      std::vector<float> tokens;
      tokens.reserve(col_length);
      
      size_t prev_pos = 0, pos = 0;
      while ((pos = line.find(" ", prev_pos)) != std::string::npos) {
	tokens.push_back(std::stof(line.substr(prev_pos, pos - prev_pos)));
	prev_pos = pos + 1;
      }
      tokens.push_back(std::stof(line.substr(prev_pos)));
      //std::cout << "Check1\n";

     
      switch (col_length) {
      case 3:
	//std::cout << "Case3\n";
	xData.push_back(tokens[0]);
	yData.push_back(tokens[1]);
	zData.push_back(tokens[2]);
	chiSq.push_back(tokens[3]);
	break;
	
      case 2:
	//std::cout << "Case2\n";
	xData.push_back(tokens[0]);
	yData.push_back(tokens[1]);
	chiSq.push_back(tokens[2]);
	break;
	
      case 1:
	//std::cout << "Case1\n";
	xData.push_back(tokens[0]);
	chiSq.push_back(tokens[1]);
	break;
	
      default:
	break;
      }
    }
  }
}

void NuFuncs::remDCP(){
  std::vector<double> xNew, yNew, zNew, chiNew;
  
  for (size_t i = 0; i < zData.size(); ++i) {
    if (zData[i] != 180) {
      xNew.push_back(xData[i]);
      yNew.push_back(yData[i]);
      zNew.push_back(zData[i]);
      chiNew.push_back(chiSq[i]);
    }
  }
  
  xData = std::move(xNew);
  yData = std::move(yNew);
  zData = std::move(zNew);
  chiSq = std::move(chiNew);
}

void NuFuncs::multiDimFit(){
  int nVars = 2;
  int nData = xData.size();

  TMultiDimFit* fit = new TMultiDimFit(nVars, TMultiDimFit::kMonomials,"v");
  int mPowers[nVars] = {1000,1000};
  fit->SetMaxPowers(mPowers);

  fit->SetMaxFunctions(5000000);
  fit->SetMaxStudy(1000000);
  fit->SetMaxTerms(200000);

  fit->SetPowerLimit(2000);
  fit->SetMinAngle(10);
  fit->SetMaxAngle(10);

  fit->SetMinRelativeError(0.0000001);

  double dataArray[nVars];
  
  for(int i = 0; i < nData; i++){
    dataArray[0] = xData[i];
    dataArray[1] = yData[i];
    //dataArray[2] = zData[i];
    //std::cout << "check3\n";
    fit->AddRow(dataArray,chiSq[i]);
  }

  for(int i = 0; i < nData; i++){
    dataArray[0] = xData[i];
    dataArray[1] = yData[i];
    //dataArray[2] = zData[i];
    std::cout << 100*(fit->Eval(dataArray) - chiSq[i])/chiSq[i] << "%\n";
  }
  //std::cout << "check1\n";

  delete fit;
}

void NuFuncs::triInterp(){
  size_t dataSize = xData.size();
  std::cout << "check-2\n";
  double* xEd = new double[dataSize+1];
  double* yEd = new double[dataSize+1];
  double* zEd = new double[dataSize+1];

  //std::sort(xData.begin(),xData.end());
  //std::sort(yData.begin(),yData.end());
  //std::sort(zData.begin(),zData.end());
  
  //std::cout << "check-1\n";
  
  //return;
  
  xEd[0] = 1.5*(xData[0]) - xData[1]/2;
  xEd[dataSize+1] = 1.5*(xData[dataSize]) - xData[dataSize-1]/2;
  yEd[0] = 1.5*(yData[0]) - yData[1]/2;
  yEd[dataSize+1] = 1.5*(yData[dataSize]) - yData[dataSize-1]/2;
  zEd[0] = 1.5*(zData[0]) - zData[1]/2;
  zEd[dataSize+1] = 1.5*(zData[dataSize]) - zData[dataSize-1]/2;

  //std::cout << "check0\n";


  for(size_t i = 0; i < dataSize-1; i++){
    //std::cout << xEd[i] << " " << yEd[i] << " " << zEd[i] << "\n";
    xEd[i+1] = (xData[i] + xData[i+1])/2;
    yEd[i+1] = (yData[i] + yData[i+1])/2;
    zEd[i+1] = (zData[i] + zData[i+1])/2;

  }
  return;

  //std::cout << "check1\n";

  TH3D* hist3D = new TH3D("Histogram","Hist for triInterp",dataSize,xEd,dataSize,yEd,dataSize,zEd);

  //std::cout << "check2\n";
}

void NuFuncs::otDI(){

  ROOT::Math::Interpolator inter(xData,chiSq,ROOT::Math::Interpolation::Type::kCSPLINE);
  /*
  double chiT[xData.size()];

  for(int i=0;i<xData.size();i++){
    chiT[i] = inter.Eval(xData[i]);
    }*/

}

void NuFuncs::twoDInterp(){
  TGraph2D* G = new TGraph2D(xData.size(), &(xData[0]),&(yData[0]),&(chiSq[0]));

  /*for(int i = 0; i< xData.size();i++){
    double err = 100*(G->Interpolate(xData[i],yData[i])-chiSq[i])/chiSq[i];
    if(err != 0) std::cout<< err << "\n";
    }*/
}

void NuFuncs::scale(){

  remDCP();

  int size = xData.size();
    
  datax = new double[size];
  datay = new double[size];
  dataz = new double[size];

 
  /*auto min_x = min_element(xData.begin(),xData.end());
  auto max_x = max_element(xData.begin(),xData.end());
  auto min_y = min_element(yData.begin(),yData.end());
  auto max_y = max_element(yData.begin(),yData.end());
  auto min_z = min_element(zData.begin(),zData.end());
  auto max_z = max_element(zData.begin(),zData.end());
  

  xmin = *min_x;
  ymin = *min_y;
  zmin = *min_z;
  xdist = *max_x - xmin;
  ydist = *max_y - ymin;
  zdist = *max_z - zmin;*/

  //std::cout << xmin << " " << ymin << " " << zmin << " " << xdist << " " << ydist << " " << zdist << "\n";

  for(int i = 0; i< xData.size(); i++){
    datax[i] = (xData[i]-xmin)/xdist;
    datay[i] = (yData[i]-ymin)/ydist;
    dataz[i] = (zData[i]-zmin)/zdist;
  }

}

void NuFuncs::scale(double* dptr){
  dptr[0] = (dptr[0]-xmin)/xdist;
  dptr[1] = (dptr[1]-ymin)/ydist;
  dptr[2] = (dptr[2]-zmin)/zdist;
}

double NuFuncs::IDW(double* dptr){

  int kNN = 3;
  float power = 12;

  std::vector<double> invD = {};
  invD.reserve(kNN);

  int* ind = new int[kNN];
  double* dist = new double[kNN];
  
  scale(dptr);
  
  tree->FindNearestNeighbors(dptr,kNN,ind,dist);

  if(dist[0] < 1e-2)return chiSq[ind[0]];

  //return chiSq[ind[0]];

  //double ratio = dist[kNN-1]/dist[0];
  //power = std::clamp(1.5*ratio , 4.0,12.0);
  
  double dSum = 0;
  double chi = 0;
  
  for(int i = 0;i < kNN;i++){
    invD[i] = 1/pow(dist[i],power);
    dSum += invD[i];
  }
  
  for(int i = 0; i < kNN; i++){
    invD[i] = invD[i]/dSum;
    chi += chiSq[ind[i]]*invD[i];
  }

  return chi;
  
}

void NuFuncs::kDTree(){
  scale();
  
  tree = new TKDTreeID(xData.size(),3,1);
  

  tree->SetData(0,datax);
  tree->SetData(1,datay);
  tree->SetData(2,dataz);
  tree->Build();
  tree->SetOwner(kTRUE);
  
}

NuFuncs::~NuFuncs(){}

