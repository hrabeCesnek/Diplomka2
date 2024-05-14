#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
//#include "TTree.h"
#include <TCanvas.h>
#include "TMultiGraph.h"
#include <iostream>
#include <fstream>
#include <numeric>
#include <dirent.h>
#include <sys/types.h>
#include <TLegend.h>
#include <math.h>
#include <TLatex.h>

#include<algorithm>

#include<iterator>

#include <sstream>
#include <locale>
#include <iomanip>

#define skip_number 12
#define baseline 10
#define FWMH 12
#define z_max 5
#define w 7
#define eps 6 
#define dev_weight 2.4 
using namespace std;


//function for loading the spectrum datafile --------------------------------------------------------------------------------------------------------------


void loadDataFile(char *name, char *restrict, vector<double> * x, vector<double> * y)
{

  double x_single = 0;
  double y_single = 0;
  string str;
  ifstream f(name);



  for(int i = 0; i < (skip_number + baseline);i++){
    getline(f, str);
  }

  while (getline(f, str))
  {

    sscanf(str.c_str(), restrict, &y_single);




    std::cout << x_single +baseline << '\n';
    std::cout << y_single << '\n';
    x->push_back(x_single + baseline);
    y->push_back(y_single);
    x_single++;
  }
  x->erase(x->end() - 15, x->end());
  y->erase(y->end() - 15, y->end());


  //clear the ending:


}


//peak finding in spectra-------------------------------------------------------------------------------

void peakFind(vector<double> * x, vector<double> * y,vector<double> * peak_positions){
  


  const int m = 3;

  int i_start = 1;
  int i_max = y->size()- 1;

 
  const int j_max = y->size();

  vector<double> channels_altered;
  vector<double> second_derivations;
  vector<double> standart_deviation;
 

  vector<vector<int>> coef_old;
  vector<vector<int>> coef_new;
  vector<double> peak_positions_diff;
  vector<int> coef_temp;
  vector<double> s_diffs;
  int coef_sum;
  string str;



//coeficient calculation

  for(int i = i_start; i< i_max;i++){
    for(int j = 0; j < j_max;j++){

      if(abs(j-i)>=2){
        coef_temp.push_back(0);

      }
      else if(abs(j-i)==1){
        coef_temp.push_back(1);

      }
      else{
        coef_temp.push_back(-2);

      }


    }
    coef_old.push_back(coef_temp);
    coef_temp.clear();




  }


  


  for(int z = 1; z <z_max+1;z++){
   

    
    
    for(int i = (m); i< coef_old.size() - (m-1);i++){

     
     for(int j = 0; j < j_max;j++){

      coef_sum = 0;
      for(int l = i - m;l< i + m;l++){
        
        coef_sum = coef_sum + coef_old[l][j];
        
      }
     
      coef_temp.push_back(coef_sum);

    }
    
    coef_new.push_back(coef_temp);
    coef_temp.clear();

    }

    for(int i = 0; i < coef_old.size(); i++){
    coef_old.at(i).clear();
    
    }
    coef_old.clear();



    for(int i = 0; i < coef_new.size(); i++){
      coef_old.push_back(coef_new.at(i));
      coef_new.at(i).clear();
    
    }
    coef_new.clear();

  }






  //second difference calculation
  
  for(int i = 0; i< coef_old.size();i++){
    double sum_diff = 0;
    double sum_devi = 0;
    for(int j = 0; j < j_max;j++){
      sum_diff = sum_diff + (coef_old[i][j] * (y->at(j)));
      sum_devi = sum_devi + (pow(coef_old[i][j],2)*(y->at(j)));
  
 

    }
  second_derivations.push_back(sum_diff);
  standart_deviation.push_back(sqrt(sum_devi));
  

  }

  cout << "změna velikosti" << endl;
  cout << ((y->size() - second_derivations.size())) << endl;

  for(int i = ((y->size() - second_derivations.size())-4); i< second_derivations.size();i++){
      channels_altered.push_back(i);
  }

















  auto c2 = new TCanvas("c2","sds",400,10,4500,2500);
  
  TMultiGraph *mg = new TMultiGraph();

  TGraph *deriv = new TGraph((Int_t)channels_altered.size(),channels_altered.data(),second_derivations.data());
  

  TGraph *deviation = new TGraph((Int_t)channels_altered.size(),channels_altered.data(),standart_deviation.data());


  mg->Add(deriv);
  mg->Add(deviation);

  mg->SetTitle(";Channel; Second difference");

  deriv->Draw("AC");
  deviation->Draw("AC");

  deriv->SetMarkerColor(kBlue);
  deriv->SetMarkerStyle(20);
  deriv->SetMarkerSize(4);
  deriv->SetLineColor(kBlue);
  deriv->SetLineWidth(10);
  deviation->SetMarkerColor(kRed);
  deviation->SetMarkerStyle(20);
  deviation->SetMarkerSize(4);
  deviation->SetLineColor(kRed);
  deviation->SetLineWidth(10);


  mg->SetMinimum(-300000.);
  mg->SetMaximum(200000.);



  mg->Draw("AC");
  mg->GetXaxis()->SetLimits(15,380);


  

  //standart deviation calculation


  // secondary conditions:


  double sigma_first = 2.335 * FWMH;
  double sigma_second;
  double f;
  double f_z = 1;

  for(int z = 1; z <z_max+1/*5*/;z++){
    f = sqrt(0.5+0.5*sqrt(1 + (0.33333)*(pow(w/sigma_first,2))));
    f_z = f_z*f;
    sigma_second = f * sigma_first;
    sigma_first = sigma_second;
    
    //cout << f << endl;

  }
  //cout << f_z << endl;
  int n_1 = round(f_z*0.849*FWMH + 0.5);
  cout << "n_1 = " << n_1 << endl;
  







  cout << "channels are equal to " << second_derivations.size() << endl;
  //peak searching
  for(int i = 1; i< second_derivations.size()-1;i++){

    int i_1_3 = -1;
    int i_1_4 = -1;
    int i_2_3 = -1;
    int i_2_4 = -1;
    int i_3 = -1;
    int i_5 = -1;
    int lp = i;

    if(abs(second_derivations.at(i))>(dev_weight*standart_deviation.at(i)) && second_derivations.at(i) < 0 && (second_derivations.at(i) - second_derivations.at(i-1) < 0) && (second_derivations.at(i) - second_derivations.at(i+1) < 0)){
      
      //points in defined intervals calculation
      cout << "peak found at position " << channels_altered.at(i) << endl;
      cout << "deriv value " << second_derivations.at(i) << endl;
      cout << "deviation value " << standart_deviation.at(i) << endl;


      //scanning backwards


      for(int j = lp; j > 1;j--){



        
        if(second_derivations.at(j-1) >= 0){
          i_3 = channels_altered.at(j);
          lp = j;
          break;
        }
      }
      if(i_3 == -1){
        continue;
      }
      for(int j = lp; j > 1;j--){

        if(second_derivations.at(j) >= standart_deviation.at(j)){
          i_2_3 = channels_altered.at(j);
          lp = j;
          break;
        }
      }
      if(i_2_3 == -1){
        continue;
      }

      for(int j = lp; j > 1;j--){

        if(second_derivations.at(j-1) <= standart_deviation.at(j-1)){
          i_1_3 = channels_altered.at(j);
          lp = j;
          break;
        }
      }
       /*if(i_1_3 == -1){
          continue;
      }*/ //we are na kraji

      lp = i;
      //scanning forward

      for(int j = lp; j < second_derivations.size()-1;j++){

        
        if(second_derivations.at(j+1) >= 0){
          i_5 = channels_altered.at(j);
          lp = j;
          break;
        }
      }
      if(i_5 == -1){
        continue;
      }
      for(int j = lp;  j < second_derivations.size()-1;j++){

        if(second_derivations.at(j) >= standart_deviation.at(j)){
          i_1_4 = channels_altered.at(j);
          lp = j;
          break;
        }
      }
       if(i_1_4 == -1){
        continue;
      }

      for(int j = lp;  j < second_derivations.size()-1;j++){

        if(second_derivations.at(j+1) <= standart_deviation.at(j+1)){
          i_2_4 = channels_altered.at(j);
          lp = j;
          break;
        }
      }
       if(i_2_4 == -1){
        continue;
      }

    

    int n_2 = round(abs((((standart_deviation.at(i))/second_derivations.at(i))*(0.5)*(n_1 + eps)) + 0.5));
   
    int n_3 = round(abs(((n_1 - eps)*( 1 - ((-2)*(standart_deviation.at(i)/second_derivations.at(i)))) )+ 0.5 ));
    

    if(n_2 == 0){
      n_2= 1;
    }
   



  cout << "i_3 =" << i_3 << endl;
    cout << "i_2_3 =" << i_2_3 << endl;
    cout << "i_1_3 = " << i_1_3 << endl;
    cout << "i_5 =" << i_5 << endl;  

    
    cout << "i_1_4 =" << i_1_4 << endl;
    
    cout << "i_2_4 =" << i_2_4 << endl;
    
  
    

    cout << "n1 condition: " << i_5 - i_3 +1 << endl;
    cout << "n2_3 condition: " << i_3 - i_2_3 - 1 << endl;
    cout << "n2_4 condition: " << i_1_4 - i_5  - 1 << endl;
    cout << "n3_4 condition: " << i_2_4 - i_1_4 + 1 << endl;
    

   
    cout << "n_1: " << n_1 << endl;
    cout << "n_2: " << n_2 << endl;
    cout << "n_3: " << n_3 << endl;

    bool first_condition = (i_5 - i_3 +1 >= n_1 - eps) && (i_5 - i_3 +1 <= n_1 + eps);
    bool second_condition = (i_3 - i_2_3 - 1 <= n_2) || (i_1_4 - i_5  - 1 <= n_2);
    bool third_condition =  (i_2_3 - i_1_3 + 1 >= n_3) || (i_2_4 - i_1_4 + 1 >= n_3) ;
    cout << first_condition << endl;
    cout << second_condition <<  endl;
    cout << third_condition  << endl;
    //determination if the found peak is energy peak not edge
    if(first_condition /*&& second_condition*/ && third_condition){  
      peak_positions_diff.push_back(channels_altered.at(i));
      peak_positions->push_back(channels_altered.at(i));
      cout << channels_altered.at(i) << endl;
      
    }
    }


  }

  TLine l2(150,0.425,150,1000);
  l2.SetLineColor(kGreen);
  l2.SetLineWidth(6);
  TLatex t1(19.8,0.7,"#scale[0.65]{paper inserted}");

  for(int i = 0; i < (Int_t)peak_positions_diff.size();i++ ){
  cout << (peak_positions_diff.at(i) * (14.4/peak_positions_diff.at(0))) << '\n';

  l2.DrawLine(peak_positions_diff.at(i),-300000.,peak_positions_diff.at(i),200000.);

  str = to_string((peak_positions_diff.at(i) * (14.4/peak_positions_diff.at(0))));

  str = str.substr(0, str.find(".")+2);
  str = "#scale[0.35]{ " + str + " keV}";
 

}

  TLegend *leg = new TLegend(0.7,0.6,0.9,1);
  leg->SetFillColor(0);
  leg->AddEntry(deriv,"Second difference");
  leg->AddEntry(deviation,"Standart deviation of second difference");
  leg->Draw();

  c2->Modified();

  c2->Print("SecondDerivGraph.png");

}





int main(int argc, char *argv[]) {


  vector<double>  channels_OPF430_FULL;
  vector<double>  counts_OPF430_FULL;

  vector<double>  channels_OPF430_Cu;
  vector<double>  counts_OPF430_Cu;

  vector<double>  channels_OPF430_BACK;
  vector<double>  counts_OPF430_BACK;

  vector<double>  channels_OPF430_Al;
  vector<double>  counts_OPF430_Al;

  vector<double>  channels_OPF430_Pb;
  vector<double>  counts_OPF430_Pb;


  vector<double>  channels_S14_FULL;
  vector<double>  counts_S14_FULL;

  vector<double>  channels_S14_Cu;
  vector<double>  counts_S14_Cu;

  vector<double>  channels_S14_BACK;
  vector<double>  counts_S14_BACK;

  vector<double>  channels_S14_Al;
  vector<double>  counts_S14_Al;

  vector<double>  channels_S14_Pb;
  vector<double>  counts_S14_Pb;

  vector<double>  channels_BPW34_FULL;
  vector<double>  counts_BPW34_FULL;

  vector<double>  channels_BPW34_Cu;
  vector<double>  counts_BPW34_Cu;

  vector<double>  channels_BPW34_BACK;
  vector<double>  counts_BPW34_BACK;

  vector<double>  channels_BPW34_Al;
  vector<double>  counts_BPW34_Al;

  vector<double>  channels_BPW34_Pb;
  vector<double>  counts_BPW34_Pb;






  vector<double> peak_positions_OPF430;
  vector<double> peak_positions_BPW34;
  vector<double> peak_positions_S14;






  string str;

  
  auto c_OPF430 = new TCanvas("c_OPF430","spectres",400,10,4000,2500);
  auto c_BPW34 = new TCanvas("c_BPW34","spectres",400,10,4000,2500);
  auto c_S14 = new TCanvas("c_S14","spectres",400,10,4000,2500);
  

  




//-------------------------------------------loading files------------------------------------------------------------------------

//OPF430
loadDataFile("../ORTEC_MCA/Finaltests/OPF430Full1.Spe","%lf",&channels_OPF430_FULL,&counts_OPF430_FULL);
loadDataFile("../ORTEC_MCA/Finaltests/OPF430Cu1.Spe","%lf",&channels_OPF430_Cu,&counts_OPF430_Cu);
loadDataFile("../ORTEC_MCA/Finaltests/OPF430BACK1.Spe","%lf",&channels_OPF430_BACK,&counts_OPF430_BACK);
loadDataFile("../ORTEC_MCA/Finaltests/OPF430Al1.Spe","%lf",&channels_OPF430_Al,&counts_OPF430_Al);
loadDataFile("../ORTEC_MCA/Finaltests/OPF430Pb1.Spe","%lf",&channels_OPF430_Pb,&counts_OPF430_Pb);
//BPW34
loadDataFile("../ORTEC_MCA/Finaltests/BPW34FULL1.Spe","%lf",&channels_BPW34_FULL,&counts_BPW34_FULL);
loadDataFile("../ORTEC_MCA/Finaltests/BPW34Cu1.Spe","%lf",&channels_BPW34_Cu,&counts_BPW34_Cu);
loadDataFile("../ORTEC_MCA/Finaltests/BPW34BACK1.Spe","%lf",&channels_BPW34_BACK,&counts_BPW34_BACK);
loadDataFile("../ORTEC_MCA/Finaltests/BPW34Al1.Spe","%lf",&channels_BPW34_Al,&counts_BPW34_Al);
loadDataFile("../ORTEC_MCA/Finaltests/BPW34Pb1.Spe","%lf",&channels_BPW34_Pb,&counts_BPW34_Pb);
//S14605
loadDataFile("../ORTEC_MCA/Finaltests/S14FULL1.Spe","%lf",&channels_S14_FULL,&counts_S14_FULL);
loadDataFile("../ORTEC_MCA/Finaltests/S14Cu1.Spe","%lf",&channels_S14_Cu,&counts_S14_Cu);
loadDataFile("../ORTEC_MCA/Finaltests/S14BACK1.Spe","%lf",&channels_S14_BACK,&counts_S14_BACK);
loadDataFile("../ORTEC_MCA/Finaltests/S14Al1.Spe","%lf",&channels_S14_Al,&counts_S14_Al);
loadDataFile("../ORTEC_MCA/Finaltests/S14Pb1.Spe","%lf",&channels_S14_Pb,&counts_S14_Pb);

//---------------------------------------peak finding---------------------------------------------------------------------------------------------------------------------------------------

peakFind(&channels_BPW34_FULL,&counts_BPW34_FULL,&peak_positions_BPW34);
peakFind(&channels_S14_FULL,&counts_S14_FULL,&peak_positions_S14);
peakFind(&channels_OPF430_FULL,&counts_OPF430_FULL,&peak_positions_OPF430);


//--------------------------------------------Graph init-----------------------------------------------------------------------------------------------------------------------------------


TMultiGraph *multi_BPW34 = new TMultiGraph();
TMultiGraph *multi_OPF430 = new TMultiGraph();
TMultiGraph *multi_S14 = new TMultiGraph();

TGraph *BPW34_FULL = new TGraph((Int_t)channels_BPW34_FULL.size(),channels_BPW34_FULL.data(),counts_BPW34_FULL.data());
TGraph *BPW34_Cu = new TGraph((Int_t)channels_BPW34_Cu.size(),channels_BPW34_Cu.data(),counts_BPW34_Cu.data());
TGraph *BPW34_Al = new TGraph((Int_t)channels_BPW34_Al.size(),channels_BPW34_Al.data(),counts_BPW34_Al.data());
TGraph *BPW34_Pb = new TGraph((Int_t)channels_BPW34_Pb.size(),channels_BPW34_Pb.data(),counts_BPW34_Pb.data());
TGraph *BPW34_BACK = new TGraph((Int_t)channels_BPW34_BACK.size(),channels_BPW34_BACK.data(),counts_BPW34_BACK.data());

BPW34_FULL->SetMarkerColor(kOrange+7);
BPW34_FULL->SetMarkerStyle(20);
BPW34_FULL->SetMarkerSize(4);
BPW34_FULL->SetLineColor(kOrange+7);
BPW34_FULL->SetLineWidth(10);

BPW34_Cu->SetMarkerColor(kBlue);
BPW34_Cu->SetMarkerStyle(20);
BPW34_Cu->SetMarkerSize(4);
BPW34_Cu->SetLineColor(kBlue);
BPW34_Cu->SetLineWidth(10);

BPW34_Al->SetMarkerColor(kViolet);
BPW34_Al->SetMarkerStyle(20);
BPW34_Al->SetMarkerSize(4);
BPW34_Al->SetLineColor(kViolet);
BPW34_Al->SetLineWidth(10);

BPW34_Pb->SetMarkerColor(kYellow+3);
BPW34_Pb->SetMarkerStyle(20);
BPW34_Pb->SetMarkerSize(4);
BPW34_Pb->SetLineColor(kYellow+3);
BPW34_Pb->SetLineWidth(10);

BPW34_BACK->SetMarkerColor(kYellow);
BPW34_BACK->SetMarkerStyle(20);
BPW34_BACK->SetMarkerSize(4);
BPW34_BACK->SetLineColor(kYellow);
BPW34_BACK->SetLineWidth(10);



TGraph *S14_FULL = new TGraph((Int_t)channels_S14_FULL.size(),channels_S14_FULL.data(),counts_S14_FULL.data());
TGraph *S14_Cu = new TGraph((Int_t)channels_S14_Cu.size(),channels_S14_Cu.data(),counts_S14_Cu.data());
TGraph *S14_Al = new TGraph((Int_t)channels_S14_Al.size(),channels_S14_Al.data(),counts_S14_Al.data());
TGraph *S14_Pb = new TGraph((Int_t)channels_S14_Pb.size(),channels_S14_Pb.data(),counts_S14_Pb.data());
TGraph *S14_BACK = new TGraph((Int_t)channels_S14_BACK.size(),channels_S14_BACK.data(),counts_S14_BACK.data());

S14_FULL->SetMarkerColor(kOrange+7);
S14_FULL->SetMarkerStyle(20);
S14_FULL->SetMarkerSize(4);
S14_FULL->SetLineColor(kOrange+7);
S14_FULL->SetLineWidth(10);

S14_Cu->SetMarkerColor(kBlue);
S14_Cu->SetMarkerStyle(20);
S14_Cu->SetMarkerSize(4);
S14_Cu->SetLineColor(kBlue);
S14_Cu->SetLineWidth(10);

S14_Al->SetMarkerColor(kViolet);
S14_Al->SetMarkerStyle(20);
S14_Al->SetMarkerSize(4);
S14_Al->SetLineColor(kViolet);
S14_Al->SetLineWidth(10);

S14_Pb->SetMarkerColor(kYellow+3);
S14_Pb->SetMarkerStyle(20);
S14_Pb->SetMarkerSize(4);
S14_Pb->SetLineColor(kYellow+3);
S14_Pb->SetLineWidth(10);

S14_BACK->SetMarkerColor(kYellow);
S14_BACK->SetMarkerStyle(20);
S14_BACK->SetMarkerSize(4);
S14_BACK->SetLineColor(kYellow);
S14_BACK->SetLineWidth(10);



TGraph *OPF430_FULL = new TGraph((Int_t)channels_OPF430_FULL.size(),channels_BPW34_FULL.data(),counts_OPF430_FULL.data());
TGraph *OPF430_Cu = new TGraph((Int_t)channels_OPF430_Cu.size(),channels_BPW34_Cu.data(),counts_OPF430_Cu.data());
TGraph *OPF430_Al = new TGraph((Int_t)channels_OPF430_Al.size(),channels_BPW34_Al.data(),counts_OPF430_Al.data());
TGraph *OPF430_Pb = new TGraph((Int_t)channels_OPF430_Pb.size(),channels_BPW34_Pb.data(),counts_OPF430_Pb.data());
TGraph *OPF430_BACK = new TGraph((Int_t)channels_OPF430_BACK.size(),channels_BPW34_BACK.data(),counts_OPF430_BACK.data());

OPF430_FULL->SetMarkerColor(kOrange+7);
OPF430_FULL->SetMarkerStyle(20);
OPF430_FULL->SetMarkerSize(4);
OPF430_FULL->SetLineColor(kOrange+7);
OPF430_FULL->SetLineWidth(10);

OPF430_Cu->SetMarkerColor(kBlue);
OPF430_Cu->SetMarkerStyle(20);
OPF430_Cu->SetMarkerSize(4);
OPF430_Cu->SetLineColor(kBlue);
OPF430_Cu->SetLineWidth(10);

OPF430_Al->SetMarkerColor(kViolet);
OPF430_Al->SetMarkerStyle(20);
OPF430_Al->SetMarkerSize(4);
OPF430_Al->SetLineColor(kViolet);
OPF430_Al->SetLineWidth(10);

OPF430_Pb->SetMarkerColor(kYellow+3);
OPF430_Pb->SetMarkerStyle(20);
OPF430_Pb->SetMarkerSize(4);
OPF430_Pb->SetLineColor(kYellow+3);
OPF430_Pb->SetLineWidth(10);

OPF430_BACK->SetMarkerColor(kYellow);
OPF430_BACK->SetMarkerStyle(20);
OPF430_BACK->SetMarkerSize(4);
OPF430_BACK->SetLineColor(kYellow);
OPF430_BACK->SetLineWidth(10);






multi_BPW34->Add(BPW34_FULL);
multi_BPW34->Add(BPW34_Cu);
multi_BPW34->Add(BPW34_Al);
multi_BPW34->Add(BPW34_Pb);
multi_BPW34->Add(BPW34_BACK);
multi_BPW34->SetTitle(";Channel; Counts");
multi_BPW34->SetMinimum(0.);
multi_BPW34->SetMaximum(1200.);




multi_S14->Add(S14_FULL);
multi_S14->Add(S14_Cu);
multi_S14->Add(S14_Al);
multi_S14->Add(S14_Pb);
multi_S14->Add(S14_BACK);
multi_S14->SetTitle(";Channel; Counts");
multi_S14->SetMinimum(0.);
multi_S14->SetMaximum(29000.);




multi_OPF430->Add(OPF430_FULL);
multi_OPF430->Add(OPF430_Cu);
multi_OPF430->Add(OPF430_Al);
multi_OPF430->Add(OPF430_Pb);
multi_OPF430->Add(OPF430_BACK);
multi_OPF430->SetTitle(";Channel; Counts");
multi_OPF430->SetMinimum(0.);
multi_OPF430->SetMaximum(350.);



TLine l3(peak_positions_BPW34[0]/14.4 * 39.5,0.1,peak_positions_BPW34[0]/14.4 * 39.5,1000);
l3.SetLineColor(kRed);
l3.SetLineWidth(6);
TLatex t3(19.8,0.7,"");

TLine l2(150,0.425,150,1000);
l2.SetLineColor(kGreen);
l2.SetLineWidth(6);
TLatex t1(19.8,0.7,"");


// ------------------------------draving spectra for every photodiode--------------------------------------------------------------------------------
c_BPW34->cd();

multi_BPW34->Draw("AC");
multi_BPW34->GetXaxis()->SetLimits(16,550);



l3.DrawLine(peak_positions_BPW34[0]/14.4 * 39.5,0.1,peak_positions_BPW34[0]/14.4 * 39.5,1200);
t3.DrawLatex(peak_positions_BPW34[0]/14.4 * 39.5 +10,500, "#scale[0.8]{39.5 keV}");




for(int i = 0; i < (Int_t)peak_positions_BPW34.size();i++ ){

  l2.DrawLine(peak_positions_BPW34[i],0.1,peak_positions_BPW34[i],1200.);

}

for(int i = 0; i < (Int_t)peak_positions_BPW34.size();i++ ){
  cout << (peak_positions_BPW34[i] * (14.4/peak_positions_BPW34[0])) << '\n';

  

  str = to_string((peak_positions_BPW34[i] * (14.4/peak_positions_BPW34[0])));

  str = str.substr(0, str.find(".")+2);
  str = "#scale[0.8]{ " + str + " keV}";
  t1.DrawLatex(peak_positions_BPW34[i] - 1,1000 - 200*i, str.c_str());


}


c_BPW34->Modified();


TLegend *leg = new TLegend(0.7,0.6,0.9,1);//(0.85,0.7,0.75,0.9);
leg->SetFillColor(0);
leg->AddEntry(BPW34_FULL,"No filter");
leg->AddEntry(BPW34_Cu,"Cu filter");
leg->AddEntry(BPW34_Al,"Al filter");
leg->AddEntry(BPW34_Pb,"Pb filter");
leg->AddEntry(BPW34_BACK,"Background");

leg->Draw();


c_BPW34->Print("BPW34GammaTest.png");




c_S14->cd();

multi_S14->Draw("AC");
multi_S14->GetXaxis()->SetLimits(16,600);


l3.DrawLine(peak_positions_S14[0]/14.4 * 39.5,0.1,peak_positions_S14[0]/14.4 * 39.5,29000.);
t3.DrawLatex(peak_positions_S14[0]/14.4 * 39.5 +5,15000, "#scale[0.8]{39.5 keV}");


for(int i = 0; i < (Int_t)peak_positions_S14.size();i++ ){

  l2.DrawLine(peak_positions_S14[i],0.1,peak_positions_S14[i],29000.);

}

for(int i = 0; i < (Int_t)peak_positions_S14.size();i++ ){
  cout << (peak_positions_S14[i] * (14.4/peak_positions_S14[0])) << '\n';

  

  str = to_string((peak_positions_S14[i] * (14.4/peak_positions_S14[0])));

  str = str.substr(0, str.find(".")+2);
  str = "#scale[0.8]{ " + str + " keV}";
  t1.DrawLatex(peak_positions_S14[i] - 1,24000 - 1500*i, str.c_str());

  //l2.DrawLine(peak_positions_S14[i],0.1,peak_positions_S14[i],29000.);

}






c_S14->Modified();


TLegend *leg_S14 = new TLegend(0.7,0.6,0.9,1);//(0.85,0.7,0.75,0.9);
leg_S14->SetFillColor(0);
leg_S14->AddEntry(S14_FULL,"No filter");
leg_S14->AddEntry(S14_Cu,"Cu filter");
leg_S14->AddEntry(S14_Al,"Al filter");
leg_S14->AddEntry(S14_Pb,"Pb filter");
leg_S14->AddEntry(S14_BACK,"Background");

leg_S14->Draw();


c_S14->Print("S14GammaTest.png");





c_OPF430->cd();

multi_OPF430->Draw("AC");
multi_OPF430->GetXaxis()->SetLimits(16,550);


l3.DrawLine(peak_positions_OPF430[0]/14.4 * 39.5,0.1,peak_positions_OPF430[0]/14.4 * 39.5,350.);
t3.DrawLatex(peak_positions_OPF430[0]/14.4 * 39.5 + 5,175, "#scale[0.8]{39.5 keV}");




for(int i = 0; i < (Int_t)peak_positions_OPF430.size();i++ ){

  l2.DrawLine(peak_positions_OPF430[i],0.1,peak_positions_OPF430[i],350.);


}





for(int i = 0; i < (Int_t)peak_positions_OPF430.size();i++ ){
  cout << (peak_positions_OPF430[i] * (14.4/peak_positions_OPF430[0])) << '\n';

  

  str = to_string((peak_positions_OPF430[i] * (14.4/peak_positions_OPF430[0])));

  str = str.substr(0, str.find(".")+2);
  str = "#scale[0.8]{ " + str + " keV}";
  t1.DrawLatex(peak_positions_OPF430[i] - 1,280 - 50*i, str.c_str());

  //l2.DrawLine(peak_positions_OPF430[i],0.1,peak_positions_OPF430[i],350.);


}




c_OPF430->Modified();


TLegend *leg_OPF430 = new TLegend(0.7,0.6,0.9,1);//(0.85,0.7,0.75,0.9);
leg_OPF430->SetFillColor(0);
leg_OPF430->AddEntry(OPF430_FULL,"No filter");
leg_OPF430->AddEntry(OPF430_Cu,"Cu filter");
leg_OPF430->AddEntry(OPF430_Al,"Al filter");
leg_OPF430->AddEntry(OPF430_Pb,"Pb filter");
leg_OPF430->AddEntry(OPF430_BACK,"Background");

leg_OPF430->Draw();


c_OPF430->Print("OPF430GammaTest.png");

return 0;








}
