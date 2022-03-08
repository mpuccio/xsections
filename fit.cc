#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>

const int color[]={kRed+2,kRed,kOrange+7,kOrange,kYellow,kSpring,kSpring+2,kGreen+3,kGreen+1,kTeal,kTeal+3,kCyan,kBlue,kBlue+3,kGreen+3,kGreen+1,kTeal,kTeal+3,kCyan,kBlue,kBlue+3,kAzure,kAzure+2,kViolet,kViolet+1};
const char *inFileDataNames[]={"antip","xsections"};
const char *outFileNames[]={"fitxsections_antip","fitxsections_p"};
const char *g4directory[]={"antiproton","proton"};
const double range[]={4.5,3.};
const int a_location_in_string[]={2,1};
const char *mom_name[2][2]={{"",""},{"","1"}};
const double leg_pos[2][4]={{0.575188,0.547826,0.906015,0.845217},{0.414787,0.653913,0.890977,0.845217}};
const int lowExcludeDataMatt=5;
const int upExcludeDataMatt=14;
const double mom_frame_range[2][2]={{0.,0.},{0.,1000.}};

class TG2TF : public TNamed {
  public:
    TG2TF();
    TG2TF(const char* name){this->SetName(name);};
    void SetInputGraph (TGraph* g_in) {g_input = (TGraph*)g_in;}
    double Eval (double *x, double *p) {
      double xx = x[0];
      return p[0]*g_input->Eval(xx); // add scaling parameter
    }
  private:
    TGraph* g_input;
};

const int kNTrials=10000;

void fit(const bool matter=true,const bool mom=false){
  gStyle->SetOptStat(0);
  TFile inFileData(Form("%s.root",inFileDataNames[matter]));
  TFile inFileG4("g4xsection.root");
  TFile outFile(Form("%s%s.root",outFileNames[matter],mom_name[mom][matter]),"recreate");
  outFile.cd();
  int nKeys = inFileData.GetNkeys();
  auto lKeys = inFileData.GetListOfKeys();
  TGraphErrors **datasets=new TGraphErrors*[nKeys];
  TF1 **models=new TF1*[nKeys];
  double ndof=0;
  for (int iK=0; iK<nKeys; ++iK) {
    if (iK>lowExcludeDataMatt&&iK<upExcludeDataMatt&&matter)continue;
    TString dataName{lKeys->At(iK)->GetName()};
    auto modelName=(TObjString*)dataName.Tokenize("-")->At(0);
    if (!matter||mom){
      auto tmp_a=(TObjString*)dataName.Tokenize("-")->At(a_location_in_string[matter]);
      modelName->SetString(Form("%s-%s_mom%s",modelName->String().Data(),tmp_a->String().Data(),mom_name[mom][matter]));
    }
    datasets[iK]=(TGraphErrors*)inFileData.Get(dataName.Data());
    datasets[iK]->SetName(dataName);
    if (!matter){
      for (int iP=0; iP<datasets[iK]->GetN(); ++iP){
        datasets[iK]->SetPointY(iP,datasets[iK]->GetPointY(iP)*0.001);
        datasets[iK]->SetPointError(iP,0,datasets[iK]->GetErrorY(iP)*0.001);
      }
    }
    else if (matter&&mom){
      for (int iP=0; iP<datasets[iK]->GetN(); ++iP){
        double tmp_E=datasets[iK]->GetPointX(iP);
        datasets[iK]->SetPointX(iP,TMath::Sqrt(tmp_E*(tmp_E+2*938.272)));
      }
    }
    auto tmp_model=(TGraph*)inFileG4.Get(Form("%s/%s",g4directory[matter],modelName->String().Data()));
    auto tg2tf=new TG2TF(modelName->String());
    tg2tf->SetInputGraph(tmp_model);
    models[iK]=new TF1(modelName->String(),tg2tf,&TG2TF::Eval,100,5000+mom_frame_range[mom][matter],1);
    ndof+=datasets[iK]->GetN();
  }
  TGraph chiSqProfile(kNTrials);
  for (int iTrial=0; iTrial<kNTrials; ++iTrial){
    double chi2=0;
    double tmp_par=iTrial/5000.-0.;
    for (int iK=0; iK<nKeys; ++iK){
      if (iK>lowExcludeDataMatt&&iK<upExcludeDataMatt&&matter)continue;
      models[iK]->SetParameter(0,tmp_par);
      for (int iP=0; iP<datasets[iK]->GetN(); ++iP){
        double model_y=models[iK]->Eval(datasets[iK]->GetPointX(iP));
        double point_y=datasets[iK]->GetPointY(iP);
        double point_y_err=datasets[iK]->GetErrorY(iP);
        chi2+=(point_y-model_y)*(point_y-model_y)/point_y_err/point_y_err;
      }
    }
    chiSqProfile.SetPoint(iTrial,tmp_par,chi2);
  }

  double minChiSq=TMath::MinElement(kNTrials,chiSqProfile.GetY());
  double minPar=-999.;
  double minParUpErr=-999.;
  double minParLowErr=-999.;
  for (int iP=0;iP<kNTrials;++iP){
    if ((chiSqProfile.GetPointY(iP)-1.e-9)<minChiSq){
      minPar=chiSqProfile.GetPointX(iP);
      break;
    }
  }
  for (int iP=kNTrials-1;iP>-1;--iP){
    if ((chiSqProfile.GetPointY(iP)-1.e-9)<(minChiSq+1)){
      minParUpErr=chiSqProfile.GetPointX(iP);
      break;
    }
  }
  for (int iP=0;iP<kNTrials;++iP){
    if ((chiSqProfile.GetPointY(iP)-1.e-9)<(minChiSq+1.)){
      minParLowErr=chiSqProfile.GetPointX(iP);
      break;
    }
  }
  std::cout<<"minPar = "<<minPar<<" +/- "<<0.5*(minParUpErr-minParLowErr)<<", chi2NDF = "<<minChiSq<<"/"<<ndof-1<<std::endl;
  TCanvas c("cFit","cFit");
  TH1D frame("frame"," ",1,0,5000+mom_frame_range[mom][matter]);
  frame.GetYaxis()->SetRangeUser(0.,range[matter]);
  frame.GetYaxis()->SetTitle("#sigma_{inel} (barn)");
  (!matter||mom) ? frame.GetXaxis()->SetTitle("#it{p} (MeV/#it{c})") : frame.GetXaxis()->SetTitle("T (MeV)");
  frame.Draw("pe");
  TLegend l(leg_pos[matter][0],leg_pos[matter][1],leg_pos[matter][2],leg_pos[matter][3]);
  for (int iK=0; iK<nKeys;++iK){
    if (iK>lowExcludeDataMatt&&iK<upExcludeDataMatt&&matter)continue;
    datasets[iK]->SetLineColor(color[iK]);
    datasets[iK]->SetMarkerColor(color[iK]);
    datasets[iK]->SetMarkerStyle(20);
    datasets[iK]->SetMarkerSize(0.9);
    models[iK]->SetLineColor(color[iK]);
    datasets[iK]->Draw("pesame");
    models[iK]->SetParameter(0,minPar);
    models[iK]->Draw("same");
    TString string_name(datasets[iK]->GetName());
    auto tmp_z=(TObjString*)string_name.Tokenize("-")->At(0);
    auto tmp_a=(TObjString*)string_name.Tokenize("-")->At(a_location_in_string[matter]);
    l.AddEntry(datasets[iK],Form("Z=%s, A=%s",tmp_z->String().Data(),tmp_a->String().Data()));
  }
  if(matter)l.SetNColumns(3);
  l.Draw("same");
  c.Write();
  TCanvas cChi2("cChi2","cChi2");
  chiSqProfile.SetLineColor(kBlue);
  TLine lH(minPar-3*(minPar-minParLowErr),minChiSq+1,minParUpErr,minChiSq+1);
  TLine lHlow(minPar-3*(minPar-minParLowErr),minChiSq,minPar,minChiSq);
  TLine lVsx(minParLowErr,minChiSq-1,minParLowErr,minChiSq+1);
  TLine lVc(minPar,minChiSq-1,minPar,minChiSq);
  TLine lVdx(minParUpErr,minChiSq-1,minParUpErr,minChiSq+1);
  lH.SetLineStyle(kDashed);
  lHlow.SetLineStyle(kDashed);
  lVsx.SetLineStyle(kDashed);
  lVc.SetLineStyle(kDashed);
  lVdx.SetLineStyle(kDashed);
  lH.SetLineColor(kGray+2);
  lHlow.SetLineColor(kGray+2);
  lVsx.SetLineColor(kGray+2);
  lVc.SetLineColor(kGray+2);
  lVdx.SetLineColor(kGray+2);
  TH1D frameChi2("frameChi2"," ",1,minPar-3*(minPar-minParLowErr),minPar+3*(minParUpErr-minPar));
  frameChi2.GetYaxis()->SetTitle("#chi^{2}");
  frameChi2.GetXaxis()->SetTitle("G4 scaling factor");
  frameChi2.Draw("pe");
  frameChi2.GetYaxis()->SetRangeUser(minChiSq-1,minChiSq+10);
  chiSqProfile.Draw("lsame");
  lH.Draw("same");
  lHlow.Draw("same");
  lVsx.Draw("same");
  lVc.Draw("same");
  lVdx.Draw("same");
  cChi2.Write();
  outFile.Close();
}