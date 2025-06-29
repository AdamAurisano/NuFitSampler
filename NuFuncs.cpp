#include "NuFuncs.h"

std::vector<double> getDMANor(){
  return {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.02,2.04,2.06,2.08,2.1,2.12,2.14,2.16,2.18,2.2,2.2,2.21,2.21,2.22,2.22,2.23,2.23,2.24,2.24,2.25,2.26,2.26,2.27,2.27,2.28,2.28,2.29,2.29,2.3,2.3,2.31,2.31,2.32,2.32,2.33,2.33,2.34,2.34,2.35,2.35,2.36,2.36,2.37,2.37,2.38,2.38,2.38,2.39,2.39,2.4,2.4,2.41,2.41,2.42,2.42,2.43,2.43,2.44,2.44,2.45,2.45,2.46,2.46,2.47,2.47,2.48,2.48,2.49,2.49,2.5,2.51,2.51,2.52,2.52,2.53,2.53,2.54,2.54,2.55,2.55,2.56,2.56,2.57,2.57,2.58,2.58,2.59,2.59,2.6,2.6,2.61,2.61,2.62,2.62,2.62,2.63,2.63,2.64,2.64,2.65,2.65,2.66,2.66,2.67,2.67,2.68,2.68,2.69,2.69,2.7,2.72,2.74,2.76,2.78,2.8,2.82,2.84,2.86,2.88,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7};
}

std::vector<double> getDMAInv(){
  return {-7,-6.8,-6.6,-6.4,-6.2,-6,-5.8,-5.6,-5.4,-5.2,-5,-4.8,-4.6,-4.4,-4.2,-4,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3,-2.9,-2.88,-2.86,-2.84,-2.82,-2.8,-2.78,-2.76,-2.74,-2.72,-2.7,-2.695,-2.69,-2.685,-2.68,-2.675,-2.67,-2.665,-2.66,-2.655,-2.65,-2.645,-2.64,-2.635,-2.63,-2.625,-2.62,-2.615,-2.61,-2.605,-2.6,-2.595,-2.59,-2.585,-2.58,-2.575,-2.57,-2.565,-2.56,-2.555,-2.55,-2.545,-2.54,-2.535,-2.53,-2.525,-2.52,-2.515,-2.51,-2.505,-2.5,-2.495,-2.49,-2.485,-2.48,-2.475,-2.47,-2.465,-2.46,-2.455,-2.45,-2.445,-2.44,-2.435,-2.43,-2.425,-2.42,-2.415,-2.41,-2.405,-2.4,-2.395,-2.39,-2.385,-2.38,-2.375,-2.37,-2.365,-2.36,-2.355,-2.35,-2.345,-2.34,-2.335,-2.33,-2.325,-2.32,-2.315,-2.31,-2.305,-2.3,-2.295,-2.29,-2.285,-2.28,-2.275,-2.27,-2.265,-2.26,-2.255,-2.25,-2.245,-2.24,-2.235,-2.23,-2.225,-2.22,-2.215,-2.21,-2.205,-2.2,-2.18,-2.16,-2.14,-2.12,-2.1,-2.08,-2.06,-2.04,-2.02,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2};
  
}

NuFuncs::NuFuncs(bool n,std::string s) :

  dma(n ? getDMAInv() : getDMANor()),

  DCPsize(dcp.size()),
  DMAmin(n ? -7 : 0.2),

  tree(nullptr),
  dataT23(nullptr),
  dataDMA(nullptr),
  dataDCP(nullptr),
  ID(s)

  {}

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
      
      size_t prev_pos = 0;
      size_t pos = 0;
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

ROOT::Math::Interpolator* NuFuncs::oneDI(){

  ROOT::Math::Interpolator* inter = new ROOT::Math::Interpolator(xData,chiSq,ROOT::Math::Interpolation::Type::kCSPLINE);
  //double chiT = inter.Eval(xData[i]);

  return inter;

   //ADD BOUND CHECK

}

TGraph2D* NuFuncs::twoDI(){

  std::string name = "Graph2D_" + this->ID;
  
  TGraph2D* G = new TGraph2D(name.c_str(),"2D_Fit",xData.size(), &(xData[0]),&(yData[0]),&(chiSq[0]));
  //std::cout << G->Interpolate(0.170,-6);
  //double err = G->Interpolate(xData[i],yData[i]);
  return G;

   //ADD BOUND CHECK
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

void NuFuncs::scale(){

  remDCP();

  int size = xData.size();
    
  dataT23 = new double[size];
  dataDMA = new double[size];
  dataDCP = new double[size];


  for(int i = 0; i< xData.size(); i++){
    dataT23[i] = (xData[i]-T23min)/T23dist;
    dataDMA[i] = (yData[i]-DMAmin)/DMAdist;
    dataDCP[i] = (zData[i]-DCPmin)/DCPdist;
  }

}

void NuFuncs::scale(double* dptr){
  dptr[0] = (dptr[0]-T23min)/T23dist;
  dptr[1] = (dptr[1]-DMAmin)/DMAdist;
  dptr[2] = (dptr[2]-DCPmin)/DCPdist;
}

void NuFuncs::kDTree(){
  scale();
  
  tree = new TKDTreeID(xData.size(),3,1);
  

  tree->SetData(0,dataT23);
  tree->SetData(1,dataDMA);
  tree->SetData(2,dataDCP);
  tree->Build();
  tree->SetOwner(kTRUE);
  
}

double NuFuncs::threeDI(double* dptr){
  //if((dptr[0] < t23[0]) || (dptr[1] < dma[0])) return -100000;
  //ADDED CHECK TO LIBRA
  if(std::abs(dptr[2] - 180) <= 1e-5) dptr[2] = -180;


  //assumes data alr refined CALL KDTREE BEFORE THIS
  auto lbx = lower_bound(t23.begin(),t23.end(),dptr[0]);
  auto lby = lower_bound(dma.begin(),dma.end(),dptr[1]);
  auto lbz = lower_bound(dcp.begin(),dcp.end(),dptr[2]);
  
  int Dindex = lbz-dcp.begin();
  

  bool xTrue = (*lbx - dptr[0] <= 1e-5);
  bool yTrue = (*lby - dptr[1] <= 1e-5);
  bool zTrue = (std::abs(read(Dindex) - dptr[2]) <= 1e-5);

  int flag = (xTrue << 2) | (yTrue << 1) | zTrue;
    
  int ind[1] = {0};
  double dist[1] = {0};
  std::vector<double> gridChi = {};
  gridChi.reserve(8);
  std::vector<std::vector<double>> pointVec = {};
  pointVec.reserve(8);
  double point[3] = {0,0,0};
  
  double ChiValue = 0;
  double chiX1, chiX2,chiX3,chiX4,chiY1,chiY2 = 0;

  double x = 0, x0=0, x1=0, y=0, y0=0, y1=0, z=0, z0=0, z1 = 0;

  double xd=0, yd=0, zd = 0;
  
  switch(flag){
  case 7: //exact point
    point[0] = dptr[0];
    point[1] = dptr[1];
    point[2] = read(Dindex);
    scale(point);
    tree->FindNearestNeighbors(point,1,ind,dist);
    ChiValue = chiSq[ind[0]];
    
    break;
  case 6: // x & y
    
    for(int i = 0; i < 2; i++){
      point[0] = dptr[0];
      point[1] = dptr[1];
      point[2] = read(Dindex - i);
      pointVec.emplace_back(point, point + 3);
      scale(point);
      tree->FindNearestNeighbors(point,1,ind,dist);
      gridChi.push_back(chiSq[ind[0]]);
    }
    
    x = dptr[2];
    y0 = gridChi[1];
    y1 = gridChi[0];
    x0 = pointVec[1][2];
    x1 = pointVec[0][2];

    ChiValue = y0*(x1-x)/(x1-x0) + y1*(x-x0)/(x1-x0);
    
    break;
  case 5: // x & z
        
    for(int i = 0; i < 2; i++){
      point[0] = dptr[0];
      point[1] = *(lby-i);
      point[2] = dptr[2];
      pointVec.emplace_back(point, point + 3);
      scale(point);
      tree->FindNearestNeighbors(point,1,ind,dist);
      gridChi.push_back(chiSq[ind[0]]);
    }

    x = dptr[1];
    y0 = gridChi[1];
    y1 = gridChi[0];
    x0 = pointVec[1][1];
    x1 = pointVec[0][1];

    ChiValue = y0*(x1-x)/(x1-x0) + y1*(x-x0)/(x1-x0);

    break;
  case 4: // x
    
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
	point[0] = dptr[0];
	point[1] = *(lby - i);
	point[2] = read(Dindex - j);
	
	pointVec.emplace_back(point, point + 3);
	scale(point);
	tree->FindNearestNeighbors(point,1,ind,dist);
	gridChi.push_back(chiSq[ind[0]]);
      }
    }

    x = dptr[1];
    y = dptr[2];
    x0 = pointVec[3][1];
    x1 = pointVec[0][1];
    y0 = pointVec[3][2];
    y1 = pointVec[0][2];

    
    chiX1 = gridChi[3]*(x1-x)/(x1-x0) + gridChi[1]*(x-x0)/(x1-x0);
    chiX2 = gridChi[2]*(x1-x)/(x1-x0) + gridChi[0]*(x-x0)/(x1-x0);

    ChiValue = chiX1*(y1-y)/(y1-y0) + chiX2*(y-y0)/(y1-y0);    
    
    break;
  case 3: // y & z
    
    for(int i = 0; i < 2; i++){
      point[0] = *(lbx-i);
      point[1] = dptr[1];
      point[2] = dptr[2];
      pointVec.emplace_back(point, point + 3);
      scale(point);
      tree->FindNearestNeighbors(point,1,ind,dist);
      gridChi.push_back(chiSq[ind[0]]);
    }

    x = dptr[0];
    y0 = gridChi[1];
    y1 = gridChi[0];
    x0 = pointVec[1][0];
    x1 = pointVec[0][0];

    ChiValue = y0*(x1-x)/(x1-x0) + y1*(x-x0)/(x1-x0);
    
    break;
  case 2: // y
    
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
	point[0] = *(lbx - i);
	point[1] = dptr[1];
	point[2] = read(Dindex - j);
	pointVec.emplace_back(point, point + 3);
	scale(point);
	tree->FindNearestNeighbors(point,1,ind,dist);
	gridChi.push_back(chiSq[ind[0]]);
      }
    }

    x = dptr[0];
    y = dptr[2];
    x0 = pointVec[3][0];
    x1 = pointVec[0][0];
    y0 = pointVec[3][2];
    y1 = pointVec[0][2];

    
    chiX1 = gridChi[3]*(x1-x)/(x1-x0) + gridChi[1]*(x-x0)/(x1-x0);
    chiX2 = gridChi[2]*(x1-x)/(x1-x0) + gridChi[0]*(x-x0)/(x1-x0);

    ChiValue = chiX1*(y1-y)/(y1-y0) + chiX2*(y-y0)/(y1-y0);
    
    break;
  case 1: // z
    
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
	point[0] = *(lbx - i);
	point[1] = *(lby - j);
	point[2] = dptr[2];
	pointVec.emplace_back(point, point + 3);
	scale(point);
	tree->FindNearestNeighbors(point,1,ind,dist);
	gridChi.push_back(chiSq[ind[0]]);
      }
    }

    x = dptr[0];
    y = dptr[1];
    x0 = pointVec[3][0];
    x1 = pointVec[0][0];
    y0 = pointVec[3][1];
    y1 = pointVec[0][1];

    
    chiX1 = gridChi[3]*(x1-x)/(x1-x0) + gridChi[1]*(x-x0)/(x1-x0);
    chiX2 = gridChi[2]*(x1-x)/(x1-x0) + gridChi[0]*(x-x0)/(x1-x0);

    ChiValue = chiX1*(y1-y)/(y1-y0) + chiX2*(y-y0)/(y1-y0);
    
    break;
  case 0: // triLinear
        
    for(int i = 0; i<2;i++){
      for(int j = 0; j<2;j++){
	for(int k = 0; k<2; k++){
	  point[0] = *(lbx-i);
	  point[1] = *(lby-j);
	  point[2] = read(Dindex-k);
	  pointVec.emplace_back(point, point + 3);
	  scale(point);
	  tree->FindNearestNeighbors(point,1,ind,dist);
	  gridChi.push_back(chiSq[ind[0]]);
	}
	
      }
      
    }
    x = dptr[0];
    y = dptr[1];
    z = dptr[2];
    x0 = pointVec[7][0];
    y0 = pointVec[7][1];
    z0 = pointVec[7][2];    
    x1 = pointVec[0][0];
    y1 = pointVec[0][1];
    z1 = pointVec[0][2];

    xd = (x-x0)/(x1-x0);
    yd = (y-y0)/(y1-y0);
    zd = (z-z0)/(z1-z0);

    chiX1 = gridChi[7]*(1-xd) + gridChi[3]*xd;
    chiX2 = gridChi[6]*(1-xd) + gridChi[2]*xd;
    chiX3 = gridChi[5]*(1-xd) + gridChi[1]*xd;
    chiX4 = gridChi[4]*(1-xd) + gridChi[0]*xd;

    chiY1 = chiX1*(1-yd) + chiX3*yd;
    chiY2 = chiX2*(1-yd) + chiX4*yd;

    ChiValue = chiY1*(1-zd) + chiY2*zd;
  
    break;

  }
  
  return ChiValue;
  
}

NuFuncs::~NuFuncs(){
  if(tree) delete tree;
}
