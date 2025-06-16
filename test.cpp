#include "NuFuncs.h"
#include "libra.h"
#include <iostream>
#include <chrono>
#include "TFile.h"
#include "TGraph.h"
#include <numeric>
#include "TH2D.h"
#include "TTree.h"

int main(int argc, char** argv){
  auto start = std::chrono::high_resolution_clock::now();

  static constexpr double optimal = 0.234; 
  
  libra obj1;

  obj1.initOnce();

  unsigned long long int totIters = (2 << 15);

  std::cout << totIters << "\n";

  std::vector<double> currdm221, currdm231, currt23, currt12, currt13, currdcp;
  currdm221.reserve(totIters);
  currdm231.reserve(totIters);
  currt23.reserve(totIters);
  currt12.reserve(totIters);
  currt13.reserve(totIters);
  currdcp.reserve(totIters);

  double dcpVal = 0;

  for(unsigned long long int iterations = 0; iterations < totIters ; iterations++){
    

    //std::cout << iterations << "\n";

    if(iterations % 100 == 0) {
      currdm221.push_back(obj1.MCS[0]->current[0]);
      currdm231.push_back(obj1.MCS[0]->current[1]);
      currt12.push_back(obj1.MCS[0]->current[2]);
      currt13.push_back(obj1.MCS[0]->current[3]);
      currt23.push_back(obj1.MCS[0]->current[4]);
      
      dcpVal = obj1.MCS[0]->current[5];
      dcpVal = (dcpVal < 0) ? dcpVal+360 : dcpVal;
      currdcp.push_back(dcpVal);
    }

    if(iterations % 10000 == 0){
      std::cout << (static_cast<double>(iterations)/totIters)*100 << "%\n";
    }

    
    
    obj1.MonteCarlo();
    
  }
  
  
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Execution time: " << duration.count() << " ms" << std::endl;

  int nbinx = 100;
  int nbiny = 100;
  
  double T23low = 0.3;
  double T23high = 0.7;

  double DCPlow = 0;
  double DCPhigh = 360;

  double T12low = 0.26;
  double T12high = 0.36;

  double DM221low = 6e-5;
  double DM221high = 9e-5;

  double T13low = 0.018;
  double T13high = 0.026;

  double DM232low = -2.7e-3;
  double DM232high = -2.3e-3;

  double DM231low = 2.3e-3;
  double DM231high = 2.7e-3;

  TH2D* T23DCP = new TH2D("T23DCP","T23 vs DCP",nbinx,T23low,T23high,nbiny,DCPlow,DCPhigh);
  TH2D* T23DM232 = new TH2D("T23DM232","T23 vs DM232 (inv)",nbinx,T23low,T23high,nbiny,DM232low,DM232high);
  TH2D* T23DM231 = new TH2D("T23DM231","T23 vs DM231 (norm)",nbinx,T23low,T23high,nbiny,DM231low,DM231high);

   TH2D* T13DCP = new TH2D("T13DCP","T13 vs DCP",nbinx,T13low,T13high,nbiny,DCPlow,DCPhigh);
   TH2D* T13DM232 = new TH2D("T13DM232","T13 vs DM232 (inv)",nbinx,T13low,T13high,nbiny,DM232low,DM232high);
  TH2D* T13DM231 = new TH2D("T13DM231","T13 vs DM231 (norm)",nbinx,T13low,T13high,nbiny,DM231low,DM231high);

  TH2D* T12DM221 = new TH2D("T12DM221","T12 vs DM221",nbinx,T12low,T12high,nbiny,DM221low,DM221high);
  TH2D* T13DM221 = new TH2D("T13DM221","T13 vs DM221",nbinx,T13low,T13high,nbiny,DM221low,DM221high);

  size_t iterSize = currdm221.size();

  for(int i = 0; i < iterSize; i++){
    T23DCP->Fill(currt23[i],currdcp[i],1.0/iterSize);
    T23DM232->Fill(currt23[i],currdm221[i]-currdm231[i],1.0/iterSize);
    T23DM231->Fill(currt23[i],currdm231[i],1.0/iterSize);

    T13DCP->Fill(currt13[i],currdcp[i],1.0/iterSize);
    T13DM232->Fill(currt13[i],currdm221[i]-currdm231[i],1.0/iterSize);
    T13DM231->Fill(currt13[i],currdm231[i],1.0/iterSize);

    T12DM221->Fill(currt12[i],currdm221[i],1.0/iterSize);
    T13DM221->Fill(currt13[i],currdm221[i],1.0/iterSize);

  }

  TFile* histFile = new TFile("HistFile.root","recreate");
  T23DCP->Write();
  T23DM232->Write();
  T23DM231->Write();
  T13DCP->Write();
  T13DM232->Write();
  T13DM231->Write();
  T12DM221->Write();
  T13DM221->Write();
  histFile->Close();

  TFile* vecFile = new TFile("vectorFile.root","recreate");
  TTree* tree = new TTree("vecTree", "Tree with all Vectors");

  tree->Branch("dm221",&currdm221);
  tree->Branch("dm231",&currdm231);
  tree->Branch("t23",&currt23);
  tree->Branch("t12",&currt12);
  tree->Branch("t13",&currt13);
  tree->Branch("dcp",&currdcp);
  tree->Fill();
  tree->Write();
  vecFile->Close();
  

  //std::vector<double> currdm221, currdm231, currt23, currt12, currt13, currdcp;


  TCanvas* c1 = new TCanvas("c1","graph",1500,800);
  T23DCP->Draw("COLZ");
  c1->SaveAs("t23_v_dcp.pdf");

  T23DM232->Draw("COLZ");
  c1->SaveAs("t23_v_dm232.pdf");

  T23DM231->Draw("COLZ");
  c1->SaveAs("t23_v_dm231.pdf");

  T13DCP->Draw("COLZ");
  c1->SaveAs("t13_v_dcp.pdf");

  T13DM232->Draw("COLZ");
  c1->SaveAs("t13_v_dm232.pdf");

  T13DM231->Draw("COLZ");
  c1->SaveAs("t13_v_dm231.pdf");

  T12DM221->Draw("COLZ");
  c1->SaveAs("t12_v_dm221.pdf");

  T13DM221->Draw("COLZ");
  c1->SaveAs("t13_v_dm221.pdf");

  



  delete T23DCP;
  delete T23DM232;
  delete T23DM231;
  delete T13DCP;
  delete T13DM232;
  delete T13DM231;
  delete T12DM221;
  delete T13DM221;
  delete c1;
  
  
  

  
  return 0;
}
