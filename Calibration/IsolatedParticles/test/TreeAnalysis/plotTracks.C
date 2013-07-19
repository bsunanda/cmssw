#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"


#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

/*
 void plotEnergy(std::string fname, std::string HLT, int var, int ien, int eta, bool logy, int pos)
 Plots energy response distribution measured energy/track momentum for tracks
 of given momentum range in a eta window
   fname = name of the i/p root file 
   HLT   = type of HLT used (to be given in the figure legend
   var   = type of energy reponse with values between 0:5
           E_{7x7}/p, H_{3x3}/p, (E_{7x7}+H_{3x3})/p,
	   E_{11x11}/p, H_{5x5}/p, (E_{11x11}+H_{5x5})/p
   ien   = Track momentum range with values between 0:9
           1:2,2:3,3:4,4:5,5:6,6:7,7:9,9:11,11:15,15:20
   eta   = calorimeter cell where track will reach 0:3
           ieta (HCAL) values 1:7, 7-13, 13:17, 17:23
   logy  = If y-axis scale shuld by linear/logarithmic
   pos   = position of the statistics boxes

 void plotTrack(std::string fname, std::string HLT, int var, bool logy, int pos)
 Plots kinematic propeties of the track 
   fname = name of the i/p root file 
   HLT   = type of HLT used (to be given in the figure legend
   var   = kinematic variable 0:3 --> p, pt, eta, phi
   logy  = If y-axis scale shuld by linear/logarithmic
   pos   = position of the statistics boxes

 void plotIsolation(std::string fname, std::string HLT, int var, bool logy, int pos)
 Plots variables used for deciding on track isolation
   fname = name of the i/p root file 
   HLT   = type of HLT used (to be given in the figure legend
   var   = isolation variable 0:3 --> Charge isolation energy,
	   Neutral isolation energy, Energy in smaller cone,
	   Energy in larger cone
   logy  = If y-axis scale shuld by linear/logarithmic
   pos   = position of the statistics boxes

 void plotHLT(std::string fname, std::string HLT, int run, bool logy, int pos)
 Plots HLT accept information for a given run or a summary
   fname = name of the i/p root file 
   HLT   = type of HLT used (to be given in the figure legend
   run   = run number; if <=0 the overall summary
   logy  = If y-axis scale shuld by linear/logarithmic
   pos   = position of the statistics boxes
 */

int         styles[7]  = {20, 24, 31, 29, 21, 22, 23};
int         colors[7]  = {1, 8, 2, 4, 7, 38, 3};
std::string names[7]   = {"All", "Quality", "okEcal", "EcalCharIso", 
			  "HcalCharIso", "EcalNeutIso", "HcalNeutIso"};
std::string namefull[7]= {"All tracks", "Good quality tracks", 
			  "Tracks reaching ECAL", "Charge isolation in ECAL", 
			  "Charge isolation in HCAL", "Isolated in ECAL", 
			  "Isolated in HCAL"};
std::string nameEta[4] = {"i#eta 1:6","i#eta 7:12","i#eta 13:16","i#eta 17:22"};
std::string namePV[5]  = {"all PV","PV 1:1","PV 2:2","PV 3:5","PV > 5"};
std::string varname[4] = {"p", "pt", "eta", "phi"};
std::string vartitle[4]= {"p (GeV/c)", "p_{T} (GeV/c)", "#eta", "#phi"};
std::string nameC[2]   = {"Ecal", "Hcal"};
std::string nameCF[2]  = {"ECAL", "HCAL"};
std::string varnameC[4]= {"maxNearP", "ediff", "ene1", "ene2"};
std::string vartitlC[4]= {"Charge isolation energy",
			  "Neutral isolation energy",
			  "Energy in smaller cone",
			  "Energy in larger cone"};
std::string varPs[10] = {"1:2","2:3","3:4","4:5","5:6","6:7","7:9",
			 "9:11","11:15","15:20"};
std::string varEta[4] = {"1:6", "7:12", "13:16", "17:23"};
std::string varEne[6] = {"E_{7x7}", "H_{3x3}", "(E_{7x7}+H_{3x3})",
			 "E_{11x11}", "H_{5x5}", "(E_{11x11}+H_{5x5})"};


void plotEMean(std::string fname="MinBias_TuneZ2star_8TeV_pythia6.root", std::string hlt="MinBias PyThia6 Z2Star", int var=0, int eta=0, int pv=0, bool data=true, int savePlot=-1);
TCanvas* plotEMeanDraw(std::vector<std::string> fnames, std::vector<std::string> hlts, int var, int eta, int pv=0, std::string dtype="Data");
TCanvas* plotEnergy(std::string fname="hlt.root", std::string HLT="All HLTs", int var=0, int ien=0, int eta=0, bool logy=true, int pos=0);
TCanvas* plotEnergyPV(std::string fnamee="hlt.root", std::string HLT="All HLTs", int var=0, int ien=0, int eta=0, bool logy=true, int pos=0);
TCanvas* plotTrack(std::string fname="hlt.root", std::string HLT="All HLTs", int var=0, bool logy=true, int pos=0);
TCanvas* plotIsolation(std::string fname="hlt.root", std::string HLT="All HLTs", int var=0, bool logy=true, int pos=0);
TCanvas* plotHLT(std::string fname="hlt.root", std::string HLT="All HLTs", int run=-1, bool logy=true, int pos=0);
TCanvas* plotHisto(char* cname, std::string HLT, TObjArray& histArr, std::vector<std::string>& labels, std::vector<int>& color, char* name, double ymx0, bool logy, int pos, double yloff, double yhoff, double xmax=-1);
void getHistStats(TH1D *h, int& entries, int& integral, double& mean, double& meanE, double& rms, double& rmsE, int& uflow, int& oflow);
TFitResultPtr getHistFitStats(TH1F *h, const char* formula, double xlow, double xup, unsigned int& nPar, double* par, double* epar);
void setHistAttr(TH1F *h, int icol, int lwid=1, int ltype=1);
double getWeightedMean (int NPT, int Start, std::vector<double>& mean, std::vector<double>& emean);

void plotAll(std::string fname="hlt.root", std::string HLT="All HLTs", int var=-1, int ien=-1, int eta=-1, bool logy=true, int pos=0, bool pv=false, int savePlot=-1) {
  int varmin(0), varmax(5), enemin(0), enemax(9), etamin(0), etamax(3);
  if (var >= 0) varmin = varmax = var;
  if (ien >= 0) enemin = enemax = ien;
  if (eta >= 0) etamin = etamax = eta;
  std::cout << "Var " << varmin << ":" << varmax << " Ene " << enemin << ":"
	    << enemax << " Eta " << etamin << ":" << etamax << std::endl;
  for (int var=varmin; var<=varmax; ++var) {
    for (int ene=enemin; ene<=enemax; ++ene) {
      for (int eta=etamin; eta<=etamax; ++eta) {
	TCanvas *c(0);
	if (pv) {
	  c = plotEnergyPV(fname, HLT, var, ene, eta, logy, pos);
	} else {
	  c = plotEnergy(fname, HLT, var, ene, eta, logy, pos);
	}
	if (c != 0 && savePlot >= 0 && savePlot < 2) {
	  std::string ext[2] = {"eps", "gif"};
	  char name[200];
	  sprintf (name, "%s.%s", c->GetName(), ext[savePlot].c_str());
	  c->Print(name);
	}
      }
    }
  }
}

void plotEMean(std::string fname, std::string hlt, int var, int eta, int pv, bool data, int savePlot) {

  std::string files[3]={"StudyHLT_PixelTrack_1PV.root","StudyHLT_ZeroBias_1PV.root","StudyHLT_1PV.root"};
//  std::string files[3]={"pixelTrack.root","zerobias.root","physics.root"};
  std::string hltx[3]={"Pixel Track HLT","Zero Bias HLT","All HLTs"};
  std::string filem[2]={"StudyHLT_MCPileup.root","StudyHLT_MCNopileup.root"};
  std::string typem[2]={"Pythia8 MinBias (Pileup)","Pythia8 MinBias (No Pileup)"};
  std::vector<std::string> fnames, hlts;
  std::string dtype("MC");
  if (data) dtype = "Data";
  if (fname == "") {
    if (data) {
      for (int i=0; i<3; ++i) {
	fnames.push_back(files[i]); hlts.push_back(hltx[i]);
      }
    } else {
      for (int i=0; i<2; ++i) {
	fnames.push_back(filem[i]); hlts.push_back(typem[i]);
      }
    }
  } else {
    fnames.push_back(fname); hlts.push_back(hlt);
  }
  int varmin(0), varmax(5), etamin(0), etamax(3), pvmin(0), pvmax(4);
  if (var >= 0) varmin = varmax = var;
  if (eta >= 0) etamin = etamax = eta;
  if (pv  >= 0) pvmin  = pvmax  = pv;
  std::cout << "Var " << varmin << ":" << varmax << " Eta " << etamin << ":" 
	    << etamax << std::endl;
  for (int var=varmin; var<=varmax; ++var) {
    for (int eta=etamin; eta<=etamax; ++eta) {
      for (int pv=pvmin; pv<=pvmax; ++pv) {
	TCanvas* c = plotEMeanDraw(fnames, hlts, var, eta, pv, dtype);
	if (c != 0 && savePlot >= 0 && savePlot < 2) {
	  std::string ext[2] = {"eps", "gif"};
	  char name[200];
	  sprintf (name, "%s.%s", c->GetName(), ext[savePlot].c_str());
	  c->Print(name);
	}
      }
    }
  }
}

TCanvas* plotEMeanDraw(std::vector<std::string> fnames, std::vector<std::string> hlts, int var, int eta, int pv, std::string dtype) {

  std::vector<TGraphAsymmErrors*> graphs;
  TLegend*  legend = new TLegend(0.60, 0.35, 0.90, 0.5);
  legend->SetBorderSize(1); legend->SetFillColor(kWhite);
  legend->SetMargin(0.4);
  for (unsigned k=0; k<fnames.size(); ++k) {
    TFile *file = TFile::Open(fnames[k].c_str());
    const int NPT=10;
    double mom[NPT]={1.5,2.5,3.5,4.5,5.5,6.5,8.0,10.0,13.0,17.5};
    double dmom[NPT]={0.5,0.5,0.5,0.5,0.5,0.5,1.0,1.0,2.0,2.5};
    double mean[NPT], dmean[NPT];
    for (int i=0; i<NPT; ++i) {
      char  name[100];
      sprintf (name, "h_energy_%d_%d_%d_%d", pv+3, i, eta, var);
      TH1D *histo = (TH1D*) file->FindObjectAny(name);
      if (histo) {
	mean[i] = histo->GetMean();
	dmean[i]= histo->GetMeanError();
      } else {
	mean[i] = -100.;
	dmean[i]= 0;
      }
    }
    /*
    std::cout << "Get mean for " << NPT << " points" << std::endl;
    for (i=0; i<NPT; ++i) 
      std::cout << "[" << i << "]" << " Momentum " << mom[i] << " +- "
		<< dmom[i] << " Mean " << mean[i] << " +- " << dmean[i] 
		<< std::endl;
    */
    TGraphAsymmErrors *graph = new TGraphAsymmErrors(NPT, mom, mean, dmom,dmom, dmean,dmean);
    graph->SetMarkerStyle(styles[k]);
    graph->SetMarkerColor(colors[k]);
    graph->SetMarkerSize(1.6);
    graph->SetLineColor(colors[k]);
    graphs.push_back(graph);
    //    TLegend* legend;
    legend->AddEntry(graph, hlts[k].c_str(), "lp");
    /*
    std::cout << "Complete " << hlts[k] << std::endl;
    */
    file->Close();
  }

  char cname[100], name[200];
  sprintf (cname, "c_%s_%d_%d_%s", varEne[var].c_str(), eta, pv, dtype.c_str());  
  gStyle->SetCanvasBorderMode(0); gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);    gStyle->SetFillColor(kWhite);
  gStyle->SetOptTitle(kFALSE);    gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  TCanvas *canvas = new TCanvas(cname, cname, 500, 400);
  gStyle->SetOptStat(0);     gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.025);
  gPad->SetBottomMargin(0.20);
  TH1F *vFrame = canvas->DrawFrame(0.0, 0.01, 50.0, 0.5);
  vFrame->GetYaxis()->SetRangeUser(0.0,1.5);
  vFrame->GetXaxis()->SetLabelSize(0.06);
  vFrame->GetYaxis()->SetLabelSize(0.05);
  vFrame->GetXaxis()->SetTitleSize(0.06);
  vFrame->GetYaxis()->SetTitleSize(0.06);
  vFrame->GetYaxis()->SetTitleOffset(0.9);
  vFrame->GetXaxis()->SetRangeUser(1.0,20.0);
  sprintf (name, "<%s/p_{Track}>", varEne[var].c_str());  
  vFrame->GetYaxis()->SetTitle(name);  
  sprintf (name, "p_{Track} (GeV/c)");
  vFrame->GetXaxis()->SetTitle(name);
  for (unsigned int ii=0; ii<graphs.size(); ++ii) graphs[ii]->Draw("P");
  legend->Draw();
  TLine *line = new TLine(1.0, 1.0, 20.0, 1.0);
  line->SetLineStyle(2); line->SetLineWidth(2); line->SetLineColor(kRed);
  line->Draw();
  TPaveText* text = new TPaveText(0.5, 0.85, 0.90, 0.90, "brNDC");
  sprintf (name, "(%s, %s)", nameEta[eta].c_str(), namePV[pv].c_str());
  std::cout << "Name " << name << " |" << std::endl;
  text->AddText(name);
  text->Draw("same");
  return canvas;
}

TCanvas* plotEnergy(std::string fname, std::string HLT, int var, int ien, int eta, bool logy, int pos) {

  TFile *file = TFile::Open(fname.c_str());
  char                     name[100];
  TObjArray                histArr;
  std::vector<std::string> labels;
  std::vector<int>         color;
  double                   ymx0(0);
  for (int i=0; i<4; ++i) {
    sprintf (name, "h_energy_%d_%d_%d_%d", i, ien, eta, var);
    TH1D *histo = (TH1D*) file->FindObjectAny(name);
    if (histo) {
      histArr.AddLast(histo);
      sprintf (name, "p=%s, #eta=%s %s", varPs[ien].c_str(), varEta[eta].c_str(), namefull[i+3].c_str());
      labels.push_back(name);
      color.push_back(colors[i]);
      int ibin = histo->GetMaximumBin();
      if (histo->GetBinContent(ibin) > ymx0) ymx0 = histo->GetBinContent(ibin);
    }
  }
  if (histArr.GetEntries()>0) {
    char cname[50];
    sprintf (cname, "c_%s_%d_%d", varEne[var].c_str(), ien, eta);  
    sprintf ( name, "%s/p", varEne[var].c_str());
    return plotHisto(cname, HLT, histArr, labels, color, name, ymx0, logy, pos, 0.10, 0.05, 2.5);
  } else {
    return 0;
  }
}

TCanvas* plotEnergyPV(std::string fname, std::string HLT, int var, int ien, int eta, bool logy, int pos) {

  const int nPVBin=4;
  int    pvBins[nPVBin+1] = {1, 2, 3, 5, 100};
  TFile *file = TFile::Open(fname.c_str());
  char                     name[100];
  TObjArray                histArr;
  std::vector<std::string> labels;
  std::vector<int>         color;
  double                   ymx0(0);
  for (int i=4; i<nPVBin+4; ++i) {
    sprintf (name, "h_energy_%d_%d_%d_%d", i, ien, eta, var);
    TH1D *histo = (TH1D*) file->FindObjectAny(name);
    if (histo) {
      histArr.AddLast(histo);
      sprintf (name, "p=%s, #eta=%s, PV=%d:%d (%s)", varPs[ien].c_str(), varEta[eta].c_str(), pvBins[i-4], pvBins[i-3]-1, namefull[6].c_str());
      labels.push_back(name);
      color.push_back(colors[i-4]);
      int ibin = histo->GetMaximumBin();
      if (histo->GetBinContent(ibin) > ymx0) ymx0 = histo->GetBinContent(ibin);
    }
  }
  if (histArr.GetEntries()>0) {
    char cname[50];
    sprintf (cname, "c_%s_%d_%d", varEne[var].c_str(), ien, eta);  
    sprintf ( name, "%s/p", varEne[var].c_str());
    plotHisto(cname, HLT, histArr, labels, color, name, ymx0, logy, pos, 0.10, 0.05, 2.5);
  }
}

TCanvas* plotTrack(std::string fname, std::string HLT, int var, bool logy, int pos) {

  TFile *file = TFile::Open(fname.c_str());
  char                     name[100];
  TObjArray                histArr;
  std::vector<std::string> labels;
  std::vector<int>         color;
  double                   ymx0(0);
  for (int i=0; i<7; ++i) {
    sprintf (name, "h_%s_%s", varname[var].c_str(), names[i].c_str());
    TH1D *histo = (TH1D*) file->FindObjectAny(name);
    if (histo) {
      histArr.AddLast(histo);
      labels.push_back(namefull[i]);
      color.push_back(colors[i]);
      int ibin = histo->GetMaximumBin();
      if (histo->GetBinContent(ibin) > ymx0) ymx0 = histo->GetBinContent(ibin);
    }
  }
  if (histArr.GetEntries()>0) {
    char cname[50];
    sprintf (cname, "c_%s", varname[var].c_str());  
    sprintf ( name, "%s", vartitle[var].c_str());
    return plotHisto(cname, HLT, histArr, labels, color, name, ymx0, logy, pos, 0.10, 0.05, -1);
  } else {
    return 0;
  }
}

TCanvas* plotIsolation(std::string fname, std::string HLT, int var, bool logy, int pos) {

  TFile *file = TFile::Open(fname.c_str());
  char                     name[100];
  TObjArray                histArr;
  std::vector<std::string> labels;
  std::vector<int>         color;
  double                   ymx0(0);
  for (int i=0; i<2; ++i) {
    sprintf (name, "h_%s_%s", varnameC[var].c_str(), nameC[i].c_str());
    TH1D *histo = (TH1D*) file->FindObjectAny(name);
    if (histo) {
      histArr.AddLast(histo);
      labels.push_back(nameCF[i]);
      color.push_back(colors[i]);
      int ibin = histo->GetMaximumBin();
      if (histo->GetBinContent(ibin) > ymx0) ymx0 = histo->GetBinContent(ibin);
    }
  }
  if (histArr.GetEntries()>0) {
    char cname[50];
    sprintf (cname, "c_%s", varnameC[var].c_str());  
    sprintf ( name, "%s (GeV)", vartitlC[var].c_str());
    return plotHisto(cname, HLT, histArr, labels, color, name, ymx0, logy, pos, 0.10, 0.05, -1);
  } else {
    return 0;
  }
}

TCanvas* plotHLT(std::string fname, std::string HLT, int run, bool logy, int pos) {

  TFile *file = TFile::Open(fname.c_str());
  char                     name[100];
  TObjArray                histArr;
  std::vector<std::string> labels;
  std::vector<int>         color;
  double                   ymx0(0);
  if (run > 0) sprintf (name, "h_HLTAccepts_%d", run);
  else         sprintf (name, "h_HLTAccept");
  TH1D *histo = (TH1D*) file->FindObjectAny(name);
  if (histo) {
    histArr.AddLast(histo);
    labels.push_back(HLT);
    color.push_back(colors[3]);
    int ibin = histo->GetMaximumBin();
    ymx0 = histo->GetBinContent(ibin);
  }
  if (histArr.GetEntries()>0) {
    char cname[50], hname[50];
    if (run > 0) {sprintf (cname,"c_HLT_%d",run); sprintf (name,"Run %d",run);
    }  else      {sprintf (cname,"c_HLTs");       sprintf (name, "All runs");}
    sprintf (hname, "HLT");
    return plotHisto(cname, "HLT", histArr, labels, color, hname, ymx0, logy, pos, 0.25, 0.01, -1);
  } else {
    return 0;
  }
}

TCanvas* plotHisto(char* cname, std::string HLT, TObjArray& histArr, std::vector<std::string>& labels, std::vector<int>& color, char* name, double ymx0, bool logy, int pos, double yloff, double yhoff, double xmax) { 

  int nentry = histArr.GetEntries();
  double ymax = 10.;
  for (int i=0; i<10; ++i) { 
    if (ymx0 < ymax) break;
    ymax *= 10.;
  }
  double ystep = ymax*0.1;
  for (int i=0; i<9; ++i) {
    if (ymax-ystep < ymx0) break;
    ymax -= ystep;
  }
  double ymin(0);
  if (logy) ymin = 1;
  std::cout << "ymin " << ymin << " ymax " << ymx0 << ":" << ymax  << std::endl;

  gStyle->SetCanvasBorderMode(0); gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);    gStyle->SetFillColor(kWhite);
  gStyle->SetOptTitle(kFALSE);    gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  TCanvas *canvas = new TCanvas(cname, cname, 500, 500);
  gStyle->SetOptStat(1110);  gPad->SetTopMargin(yhoff);
  gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.025);
  gPad->SetBottomMargin(yloff);
  if (logy) canvas->SetLogy();
  TLegend  *legend(0);
  TPaveText *text(0);
  double height = 0.08;
  double dx = (nentry > 2) ? 0.50 : 0.30;
  double dy = (nentry > 2) ? 0.12 : 0.04;
  double xmin1 = (pos > 1) ? 0.375 : 0.75-dx;
  double xmin2 = (pos > 1) ? 0.15 : 0.75;
  double ymin2 = (pos%2 == 0) ? (0.96-yhoff-nentry*height) : (yloff+0.025+nentry*height);
  std::cout << nentry << " " << pos << " " << dx << " " << dy << " " << xmin1 << " " << xmin2 << std::endl;
  if (pos%2 == 0) {
    legend  = new TLegend(xmin1, 1.0-yhoff-dy, xmin1+dx, 1.0-yhoff);   
    text    = new TPaveText(xmin2, ymin2, xmin2+0.225, ymin2+0.035, "brNDC");
  } else {
    legend  = new TLegend(xmin1, yloff+0.02, xmin1+dx, yloff+0.02+dy);
    text    = new TPaveText(xmin2, ymin2, xmin2+0.225, ymin2+0.035, "brNDC");
  }
  legend->SetBorderSize(1); legend->SetFillColor(kWhite);
  text->AddText(HLT.c_str());
  THStack *Hs      = new THStack("hs2"," ");
  for (int i=0; i<nentry; i++) {
    TH1D *h =  (TH1D*)histArr[i];
    h->SetLineColor(color[i]);
    h->SetLineWidth(2);
    h->SetMarkerSize(0.8);
    h->GetYaxis()->SetRangeUser(ymin,ymax);
    if (xmax > 0) h->GetXaxis()->SetRangeUser(0,xmax);
    Hs->Add(h, "hist sames");
    legend->AddEntry(h,labels[i].c_str(),"l");
  }
  Hs->Draw("nostack");
  canvas->Update();
  Hs->GetHistogram()->GetXaxis()->SetTitle(name);
  Hs->GetHistogram()->GetYaxis()->SetTitle("Tracks");
  Hs->GetHistogram()->GetXaxis()->SetLabelSize(0.035);
  Hs->GetHistogram()->GetYaxis()->SetTitleOffset(1.6);
  if (xmax > 0) Hs->GetHistogram()->GetXaxis()->SetRangeUser(0,xmax);
  canvas->Modified();
  
  canvas->Update();
  for (int i=0; i<nentry; i++) {
    TH1D *h =  (TH1D*)histArr[i];
    if (h != NULL) {
      TPaveStats* st1 = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
      if (st1 != NULL) {
	if (pos%2 == 0) {
	  st1->SetY1NDC(1.0-yhoff-(i+1)*height); st1->SetY2NDC(1.0-yhoff-i*height);
	} else {
	  st1->SetY1NDC(yloff+0.02+i*height); st1->SetY2NDC(yloff+0.02+(i+1)*height);
	}
	if (pos > 1) {
	  st1->SetX1NDC(0.15); st1->SetX2NDC(.375);
	} else {
	  st1->SetX1NDC(0.75); st1->SetX2NDC(.975);
	}
	st1->SetTextColor(colors[i]);
      }
    }
  }
  legend->Draw("");
  text->Draw("same");
  
  return canvas;
}

void getHistStats(TH1D *h, int& entries, int& integral, double& mean, double& meanE, double& rms, double& rmsE, int& uflow, int& oflow) {
  entries  = h->GetEntries();
  integral = h->Integral();
  mean     = h->GetMean();
  meanE    = h->GetMeanError();
  rms      = h->GetRMS();
  rmsE     = h->GetRMSError();
  uflow    = h->GetBinContent(0);
  oflow    = h->GetBinContent( h->GetNbinsX()+1 );
}

TFitResultPtr getHistFitStats(TH1F *h, const char* formula, double xlow, 
			      double xup, unsigned int& nPar, double* par, 
			      double* epar) {

  TFitResultPtr fit =  h->Fit(formula, "+qRB0", "", xlow, xup);
  nPar = fit->NPar();
  const std::vector<double> errors = fit->Errors();
  for (unsigned int i=0; i<nPar; i++) {
    par[i]  = fit->Value(i);
    epar[i] = errors[i];
  }
  return fit;
}

void setHistAttr(TH1F *h, int icol, int lwid, int ltype) {
  h->SetLineColor(icol);
  h->SetLineStyle(ltype);
  h->SetLineWidth(lwid);
  TF1 *f = h->GetFunction("gaus");
  if (!f->IsZombie()) {
    f->SetLineColor(icol);
    f->SetLineStyle(2);
  }
}

double getWeightedMean (int NPT, int Start, std::vector<double>& mean, std::vector<double>& emean) {
  double sumDen=0, sumNum=0;
  for (int i=Start; i<NPT; i++){
    if (mean[i]==0.0 || emean[i]==0.0) { 
      sumNum += 0;
      sumDen += 0;
    } else { 
      sumNum += mean[i]/(emean[i]*emean[i]);
      sumDen += 1.0/(emean[i]*emean[i]);
    }
  }
  double WeightedMean = sumNum/sumDen;
  return WeightedMean;
}

