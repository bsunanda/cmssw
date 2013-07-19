#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int styles[6] = {20, 24, 31, 29, 21, 22};

int colors[6] = {2, 8, 1, 4, 41, 38};

void plotScale(char* det="HB", int var=0, std::string filename="run2012a.txt", char* title="Run 2012A");
std::vector<std::string> readFileNames(std::string filename);
TH1I* getStack(char* dirname, char* histName, std::vector<std::string>& files, bool accumulate);
void plotTrigger(int var=0, int plots=-1, std::string fname="AllRuns.root", char *gTitle="Run 2012 A,B,C",  bool logy=true, int pos=0);
void plotEta(int var=0, int icut=0, int jcut=0, int kcut=0, 
	     std::string fname="AllRuns.root", bool logy=true, int pos=0);
TCanvas* plotHisto(char* cname, char* globaltitle, TObjArray& histArr, std::vector<std::string>& labels, std::vector<int>& color, char* xtitle, char* ytitle, double ymx0, bool logy, bool largeLegend, int pos, double yloff, double yhoff);

void plotScale(char* det, int var, std::string filename, char* title) {

  std::string varName[]  = {"PreL1", "PreHLT", "Pre"};
  std::string varTitle[] = {"L1", "HLT", "Overall"};

  std::vector<std::string> files = readFileNames(filename);
  std::cout << files.size() << " input files to be used" << std::endl;
  char histname[20], dirname[20];
  sprintf(histname, "h_%svsRN", varName[var].c_str());
  sprintf(dirname,  "IsoTrig%s", det);
  std::cout << "Directory " << dirname << " histogram " << histname << setd::endl;
  TH1I* hist = getStack(dirname, histname, files, false);

  gStyle->SetCanvasBorderMode(0); gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);    gStyle->SetFillColor(kWhite);
  gStyle->SetOptTitle(kFALSE);    gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0); gStyle->SetOptStat(10);
  char name[100];
  sprintf (name, "c_%s", histname);
  std::cout << "Canvas " << name << " for " << hist << std::endl;
  TCanvas *canvas = new TCanvas(name, name, 500, 500);
  gPad->SetTopMargin(0.065);  gPad->SetBottomMargin(1.0);
  gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.02);
  hist->GetXaxis()->SetTitle("Run number");
  sprintf (name, "%s prescale factor", varTitle[var].c_str());
  hist->GetYaxis()->SetTitle(name); hist->GetYaxis()->SetTitleOffset(1.6);
  hist->Draw();
  TPaveText *text1 = new TPaveText(0.84, 0.850, 0.975, 0.885, "brNDC");
  text1->SetFillColor(kWhite); text1->AddText(title); text1->Draw("same");
  TPaveText *text2 = new TPaveText(0.92, 0.815, 0.975, 0.850, "brNDC");
  text2->SetFillColor(kWhite); text2->AddText(det);   text2->Draw("same");
}

std::vector<std::string> readFileNames(std::string filename) {  

  std::vector<std::string> files;
  std::string buffer;
  ifstream infile(filename.c_str());
  std::cout << "Reading from " << filename << std::endl;
  if (!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << filename
	      << "' for the list of input files" << std::endl;
  } else {
    while(1) {
      infile >> buffer;
      if (!infile.good()) break;
      files.push_back(buffer);
      //      std::cout << "File[" << files.size() << "] = " << buffer << std::endl;
    }
  }
  return files;
}

TH1I* getStack(char* dirname, char* histName, std::vector<std::string>& files, 
	       bool accumulate) {

  int    nbin(0);
  double lowx(0), highx(0), xstep(0);
  for (unsigned int ifile=0; ifile<files.size(); ++ifile) {
    TFile *file = TFile::Open(files[ifile].c_str());
    if (file != 0) {
      TDirectory *dir = (TDirectory*)file->FindObjectAny(dirname);
      if (dir) {
	TH1I *histin = (TH1I*)dir->FindObjectAny(histName);
	if (histin) {
	  nbin = histin->GetNbinsX();
	  lowx  = histin->GetBinLowEdge(1);
	  xstep = histin->GetBinWidth(1);
	  highx = lowx + nbin;
	  std::cout << nbin << " bins with step size " << xstep << " between " << lowx << ":" << highx << std::endl;
	  file->Close();
	  break;
	}
      }
      file->Close();
    }
  }

  char  name[100];
  sprintf (name, "%s_1", histName);
  TH1I* hist = new TH1I(name, "Modified", nbin, lowx, highx);
  std::cout << "Histogram " << name << " booked with " << nbin << " bins between " << lowx << ":" << highx << std::endl;

  for (unsigned int ifile=0; ifile<files.size(); ++ifile) {
    //    std::cout << "Opening file " << files[ifile] << std::endl;
    TFile *file = TFile::Open(files[ifile].c_str());
    if (file != 0) {
      TDirectory *dir = (TDirectory*)file->FindObjectAny(dirname);
      if (dir) {
	TH1I *histin = (TH1I*)dir->FindObjectAny(histName);
	//      std::cout << "Finds histogram " << histName << " at " << histin << std::endl;
	if (histin) {
	  int nbin = histin->GetNbinsX();
	  for (int i=1; i<=nbin; ++i) {
	    if (histin->GetBinContent(i) > 0) {
	      if (hist->GetBinContent(i) == 0) {
		hist->SetBinContent(i, histin->GetBinContent(i));
		//		std::cout << "Insert bin " << i << " from file " << ifile << " with " << histin->GetBinContent(i) << std::endl;
	      } else if (accumulate) {
		//		std::cout << "Modify bin " << i << " with " << hist->GetBinContent(i) << " using from file " << ifile << " to " << (hist->GetBinContent(i)+histin->GetBinContent(i)) << std::endl;
		hist->SetBinContent(i, (hist->GetBinContent(i)+histin->GetBinContent(i)));
	      }
	    }
	  }
	}
      }
      file->Close();
    }
  }
  return hist;
}

void plotTrigger(int var, int plots, std::string fname, char *gTitle, 
		 bool logy, int pos) {

  std::string varName[]  = {"HLT", "HLTAcceptsvsRN", "PreL1", "PreHLT", "Pre"};
  std::string varTitlx[] = {"HLT Trigger bit", "Run Number", 
			    "L1 Pre-scale Factor", "HLT Pre-scale Factor",
			    "Overall Pre-scale Factor"};
  std::string varTitly[] = {"Events", "HLT passed events", "Events", "Events", 
			    "Events"};
  std::string labelx[]   = {"HB", "HE"};
  TFile *file = TFile::Open(fname.c_str());
  if (file) {
    char                     name[100], histname[20], dirname[20];
    TObjArray                histArr;
    std::vector<std::string> labels;
    std::vector<int>         color;
    double                   ymx0(0);
    sprintf(histname, "h_%s", varName[var].c_str());
    int                      plotFirst(0), plotLast(1);
    if (plots >= 0) plotFirst = plotLast = plots;
    for (int i=plotFirst; i<=plotLast; ++i) {
      sprintf(dirname,  "IsoTrig%s", labelx[i].c_str());
      TDirectory *dir = (TDirectory*)file->FindObjectAny(dirname);
      if (dir) {
	TH1I *histo = (TH1I*)dir->FindObjectAny(histname);
	if (histo) {
	  histArr.AddLast(histo);
	  labels.push_back(labelx[i]);
	  color.push_back(colors[i]);
	  int ibin = histo->GetMaximumBin();
	  if (histo->GetBinContent(ibin) > ymx0) ymx0 = histo->GetBinContent(ibin);
	}
      }
    }
    if (histArr.GetEntries()>0) {
      char cname[50];
      sprintf (cname, "c_%s", varName[var].c_str());  
      plotHisto(cname, gTitle, histArr, labels, color, varTitlx[var].c_str(),
		varTitly[var].c_str(), ymx0, logy, false, pos, 0.10, 0.05);
    }
  } else {
    std::cout << "Cannot open file " << fname << std::endl;
  }
}

void plotEta(int var, int icut, int jcut, int kcut, std::string fname, 
	     bool logy, int pos) {

  std::string varNam1[]  = {"etaCalibTracks", "etaMipTracks"};
  std::string varNam2[]  = {"HLTMatched", "HLTNotMatched"};
  std::string varTitlx1[]= {"Totally isolated", "Charge isolated"}; 
  std::string varTitlx2[]= {"matched", "not matched"}; 
  std::string varTitlx3[]= {"all", "away from L1"};
  std::string varTitlx4[]= {"20:30", "30:40", "40:60", "60:80", "80:120"};
  std::string labelx[]   = {"HB", "HE"};
  TFile *file = TFile::Open(fname.c_str());
  std::cout << "Opens " << fname << " at " << file << std::endl;
  if (file) {
    char                     name[100], histname[100], dirname[20], gtitle[200];
    TObjArray                histArr;
    std::vector<std::string> labels;
    std::vector<int>         color;
    double                   ymx0(0);
    sprintf(histname, "h_%s%sCut%dLim%d", varNam1[var].c_str(), varNam2[icut].c_str(), jcut, kcut);
    sprintf(gtitle, "%s %s %s tracks (p=%s)",varTitlx1[var].c_str(), 
	    varTitlx3[jcut].c_str(), varTitlx2[icut].c_str(), 
	    varTitlx4[kcut].c_str());
    std::cout << "Title: " << gtitle << "\nHistogram: " << histname << std::endl;
    for (int i=0; i<2; ++i) {
      sprintf(dirname,  "IsoTrig%s", labelx[i].c_str());
      TDirectory *dir = (TDirectory*)file->FindObjectAny(dirname);
      std::cout << "Directory : " << dirname << " at " << dir << std::endl;
      if (dir) {
	TH1D *histo = (TH1D*)dir->FindObjectAny(histname);
	if (histo) {
	  histArr.AddLast(histo);
	  labels.push_back(labelx[i]);
	  color.push_back(colors[i]);
	  int ibin = histo->GetMaximumBin();
	  if (histo->GetBinContent(ibin) > ymx0) ymx0 = histo->GetBinContent(ibin);
	}
      }
    }
    if (histArr.GetEntries()>0) {
      char cname[100];
      sprintf (cname, "c_%s%sCut%dLim%d", varNam1[var].c_str(), varNam2[icut].c_str(), jcut, kcut);
      std::cout << "Canvas: " << cname << std::endl;
      plotHisto(cname, gtitle, histArr, labels, color, "#eta", "Tracks",
		ymx0, logy, true, pos, 0.10, 0.05);
    }
  }
}

TCanvas* plotHisto(char* cname, char* globaltitle, TObjArray& histArr, 
		   std::vector<std::string>& labels, std::vector<int>& color,
		   char* xtitle, char* ytitle, double ymx0, bool logy, 
		   bool largeLegend, int pos, double yloff, double yhoff) { 
  
  int nentry = histArr.GetEntries();
  double ymax = 10.;
  for (int i=0; i<10; ++i) { 
    if (ymx0 < ymax) break;
    ymax *= 10.;
  }
  double ystep = ymax*0.1;
  if (ymax > 10) {
    for (int i=0; i<9; ++i) {
      if (ymax-ystep < ymx0) break;
      ymax -= ystep;
    }
  }
  double ymin(0);
  if (logy) ymin = 0.5;
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
  double dx    = 0.20;
  double dy    = 0.08;
  double dxleg = (largeLegend) ? 0.59 : 0.45;
  double xmin1 = (pos > 1) ? 0.375 : 0.75-dx;
  double xmin2 = (pos > 1) ? 0.375 : 0.75-dxleg;
  if (largeLegend) 
  std::cout << largeLegend << " " << pos << " " << dx << " " << dy << " " << xmin1 << " " << xmin2 << std::endl;
  if (pos%2 == 0) {
    legend  = new TLegend(xmin1, 1.0-yhoff-dy, xmin1+dx, 1.0-yhoff);   
    text    = new TPaveText(xmin2, 0.95-yhoff-dy, xmin2+dxleg, 0.99-yhoff-dy, "brNDC");
  } else {
    legend  = new TLegend(xmin1, yloff+0.02, xmin1+dx, yloff+0.02+dy);
    text    = new TPaveText(xmin2, yloff+0.03+dy, xmin2+dxleg, yloff+0.07+dy, "brNDC");
  }
  legend->SetBorderSize(1); legend->SetFillColor(kWhite);
  text->AddText(globaltitle);
  THStack *Hs      = new THStack("hs2"," ");
  for (int i=0; i<nentry; i++) {
    TH1 *h =  (TH1*)histArr[i];
    h->SetLineColor(color[i]);
    h->SetLineWidth(2);
    h->SetMarkerSize(0.8);
    h->GetYaxis()->SetRangeUser(ymin,ymax);
    Hs->Add(h, "hist sames");
    legend->AddEntry(h,labels[i].c_str(),"l");
  }
  Hs->Draw("nostack");
  canvas->Update();
  Hs->GetHistogram()->GetXaxis()->SetTitle(xtitle);
  Hs->GetHistogram()->GetYaxis()->SetTitle(ytitle);
  Hs->GetHistogram()->GetXaxis()->SetLabelSize(0.035);
  Hs->GetHistogram()->GetYaxis()->SetTitleOffset(1.6);
  canvas->Modified();
  
  double height =0.08;
  canvas->Update();
  for (int i=0; i<histArr.GetEntries(); i++) {
    TH1 *h =  (TH1*)histArr[i];
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
	st1->SetTextColor(color[i]);
      }
    }
  }
  legend->Draw("");
  text->Draw("same");
  
  return canvas;
}
