#include "NuFuncs.h"
#include <iostream>
#include <chrono>
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TH2D.h"
#include <TGraph.h>
#include <Math/Interpolator.h>
#include <Math/InterpolationTypes.h>
#include "TFile.h"




int main(int argc, char** argv){

  auto start = std::chrono::high_resolution_clock::now();
  
  NuFuncs df1;
  NuFuncs df2;
  df1.putIntoDF("normalOrder.txt",1);
  df2.putIntoDF("normalOrder.txt",7);
  TGraph2D* G = new TGraph2D(df2.xData.size(), &(df2.xData[0]),&(df2.yData[0]),&(df2.chiSq[0]));
  df1.kDTree();

  //std::cout << "check 1\n"; 

  int nbinx = 25;
  int nbiny = 25;
  double xlow = 0.42;
  double xhigh = 0.6;
  double ylow = 2.45;
  double yhigh = 2.62;

  TH2D* T2Dplot = new TH2D("T2DPlot","T23 vs DMA normorder 2D",nbinx,xlow,xhigh,nbiny,ylow,yhigh);
  TH2D* T3Dplot = new TH2D("T3DPlot","T23 vs DMA normorder 3D",nbinx,xlow,xhigh,nbiny,ylow,yhigh);
  
  double* xbinList = new double[nbinx];
  double* ybinList = new double[nbiny];  

  for(int i = 0; i < nbinx; i++){
    xbinList[i] = T2Dplot->GetXaxis()->GetBinCenter(i+1);
  }
  for(int i = 0; i < nbiny; i++){
    ybinList[i] = T2Dplot->GetYaxis()->GetBinCenter(i+1);
  }

  //std::cout << "check 2\n"; 

  float spacing = 5;
  
  for(int i = 0; i < nbinx; i++){
    for(int j = 0; j < nbiny; j++){

      double val = G->Interpolate(xbinList[i],ybinList[j]);
      T2Dplot->SetBinContent(i + 1, j + 1, val);
      
      double DCP = -180;
      double dp[3] = {xbinList[i],ybinList[j],0};
      double minChi = 1000000;
      double ChiSq = 0;
      
      while(DCP <= 180){
	dp[2] = DCP;
	ChiSq = df1.IDW(dp);
	if(minChi > ChiSq ) minChi = ChiSq;
	
	DCP += spacing;
      }
      //std::cout << minChi << " ";
      T3Dplot->SetBinContent(i + 1, j + 1, minChi);
      std::cout <<"[" <<minChi << " " << i << " " << j << "] ";
    }
    std::cout<< "\n\n";
    }
  
  
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  
  std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
  
  delete df1.tree;
  delete[] xbinList;
  delete[] ybinList;

  //std::cout << "check 3\n"; 
  T2Dplot->GetXaxis()->SetTitle("T23");
  T2Dplot->GetYaxis()->SetTitle("DMA");
  

  T3Dplot->GetXaxis()->SetTitle("T23");
  T3Dplot->GetYaxis()->SetTitle("DMA");
  


  TFile* file2D = new TFile("2DOutput.root","recreate");
  T2Dplot->Write();
  file2D->Close();

  TFile* file3D = new TFile("3DOutput.root","recreate");
  T3Dplot->Write();
  file3D->Close();

  
  T3Dplot->Draw("COLZ");

  
  TApplication app("app", &argc, argv);

  Double_t w = 600;
  Double_t h = 600;
  
  TCanvas* c2D = new TCanvas("Graph2D","Fit2D",w,h);
  c2D->SetWindowSize(w + (w - c2D->GetWw()), h + (h - c2D->GetWh()));
  T2Dplot->Draw("COLZ");
  c2D->Update();

  TCanvas* c3D = new TCanvas("Graph3D","Fit3D",w,h);
  c3D->SetWindowSize(w + (w - c3D->GetWw()), h + (h - c3D->GetWh()));
  T3Dplot->Draw("COLZ");
  c3D->Update();
  
  app.Run();
  
  delete T2Dplot;
  delete T3Dplot;
     
  return 0;
}
