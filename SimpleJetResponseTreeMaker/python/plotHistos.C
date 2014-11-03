{

  using namespace std;
gStyle->SetCanvasBorderMode(0);
cout<<endl;  
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(500); //Height of canvas
  gStyle->SetCanvasDefW(800); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);

  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

  // For the histo:
  // gStyle->SetHistFillColor(1);
  // gStyle->SetHistFillStyle(0);
  //cut out by jim//  gStyle->SetHistLineColor(1);
  //cut out by jim//  gStyle->SetHistLineStyle(0);
  //cut out by jim//  
  gStyle->SetHistLineWidth(2);
  // gStyle->SetLegoInnerR(Float_t rad = 0.5);
  // gStyle->SetNumberContours(Int_t number = 20);

  //cut out by jim// gStyle->SetEndErrorSize(2);
  //gStyle->SetErrorMarker(20);
  //cut out by jim// gStyle->SetErrorX(0.);

  //gStyle->SetMarkerStyle(20);

  //For the fit/function:
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

  //For the date:
  gStyle->SetOptDate(0);
  // gStyle->SetDateX(Float_t x = 0.01);
  // gStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(111111); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);
  // gStyle->SetStatStyle(Style_t style = 1001);
  // gStyle->SetStatX(Float_t x = 0);
  // gStyle->SetStatY(Float_t y = 0);

  // Margins:
  gStyle->SetPadTopMargin(0.14);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.17);//0.14);//0.16);
  gStyle->SetPadRightMargin(0.13);//0.3);//0.16);

  // For the Global title:

  //  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  //gStyle->SetTitleFontSize(0.05);
    gStyle->SetTitleFontSize(0.04);

  // gStyle->SetTitleH(0); // Set the height of the title box
  gStyle->SetTitleW(0); // Set the width of the title box
  gStyle->SetTitleX(0.5); // Set the position of the title box
  gStyle->SetTitleY(0.94); // Set the position of the title box
  // gStyle->SetTitleStyle(Style_t style = 1001);
  // gStyle->SetTitleBorderSize(2);

  // For the axis titles:

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  //gStyle->SetTitleSize(0.06, "XYZ");
  //gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // gStyle->SetTitleYSize(Float_t size = 0.02);
  //gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleXOffset(1.2);
  //   gStyle->SetTitleYOffset(1.5);
  gStyle->SetTitleOffset(1.5, "Y"); // Another way to set the Offset

  // For the axis labels:

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "YZ");
  gStyle->SetLabelOffset(0.012, "X");
  //gStyle->SetLabelSize(0.05, "XYZ");
  ///gStyle->SetLabelSize(0.035, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");

  // For the axis:

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

  // Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);

  gStyle->SetPalette(1,0);



  //gStyle->SetLegendFont(42);
  gStyle->SetOptStat(0);
  gROOT->UseCurrentStyle();
  gROOT->ForceStyle();


  gROOT->Reset();
  gROOT->ProcessLine(".L /Users/jdolen/Dropbox/Code/MyAnalysisRootFunctions_NEW.C");

  gROOT->Reset();
 gROOT->ForceStyle(); 
 //gROOT->SetBatch(); 

  //Define my color table
  const Int_t NRGBs = 7;
  const Int_t NCont = 555;
  //Each column is a color{ violet, blue, light blue, green, yellow, orange,  red }
  Double_t stops[NRGBs] = {   0.00, 0.04,       0.15,  0.25,   0.41,   0.74, 1.00 };
  Double_t red[NRGBs]   = {   0.60, 0.00,       0.00,  0.00,   0.87,   1.00, 0.51 };
  Double_t green[NRGBs] = {   0.00, 0.00,       0.81,  1.00,   1.00,   0.20, 0.00 };
  Double_t blue[NRGBs]  = {   1.00, 0.51,       1.00,  0.00,   0.12,   0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  Int_t nb = 555;


  #include <algorithm> 


  TCanvas *c1236= new TCanvas("c1236","",200,10,800,700);



   TFile *F1   = new TFile("OutTreeQCDflat.root"); 
   //QCDvPFjetMass              = DrawHistogramFromTreeAndFindMedian(QCD, "JetTree", 100, 0, 500, 160, 250, "PFjetMass" );


   TFile *F2   = new TFile("TGraphs.root"); 

    // gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300->SetTitle(";GenJet p_{T} (GeV/c);<p_{T}(gen)/p_{T}(reco)>");

  // gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20     
  // gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25     
  // gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300    
  // gr_MEAN_PT_RESPONSE_VS_GEN_PT_central
  // gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap 
  // gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward
  // gr_PT_RESOLUTION_VS_GEN_PT_central 
  // gr_PT_RESOLUTION_VS_GEN_PT_endcap  
  // gr_PT_RESOLUTION_VS_GEN_PT_forward 




  TGraph *gr_MEAN_PT_RESPONSE_VS_GEN_PT_central = F2->Get("gr_MEAN_PT_RESPONSE_VS_GEN_PT_central");
  TGraph *gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  = F2->Get("gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap");
  TGraph *gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward = F2->Get("gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward");

  gr_MEAN_PT_RESPONSE_VS_GEN_PT_central->SetTitle(";GenJet p_{T} (GeV/c);<p_{T}^{gen}/p_{T}^{reco}>");

  gr_MEAN_PT_RESPONSE_VS_GEN_PT_central ->SetLineColor(4);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  ->SetLineColor(3);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward ->SetLineColor(2);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_central ->SetMarkerColor(4);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  ->SetMarkerColor(3);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward ->SetMarkerColor(2);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_central ->SetMarkerStyle(20);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  ->SetMarkerStyle(21);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward ->SetMarkerStyle(22);
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_central ->Draw("Ap");
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  ->Draw("p");
  gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward ->Draw("p");

  leg = new TLegend(0.19,0.67,0.45,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(gr_MEAN_PT_RESPONSE_VS_GEN_PT_central ,"| #eta | < 1.3","LP");
  leg->AddEntry(gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  ,"1.3 < | #eta | < 3","LP");
  leg->AddEntry(gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward ,"3 < | #eta | < 5","LP");
  leg->Draw("same");

  TPaveText *textbox1 = new TPaveText(0.55,0.55,0.89,0.65,"NDC");
  textbox1->SetFillColor(0);
  textbox1->SetLineColor(0);
  TText *line1 = textbox1->AddText("QCD Flat - 72X");
  TText *line2 = textbox1->AddText("AK4 PF Jets");
  TText *line3 = textbox1->AddText("AK5 Gen Jets");
  TText *line4 = textbox1->AddText("CSA14 AK4 Corrections");
  line1->SetTextColor(1);
  line2->SetTextColor(1);
  line3->SetTextColor(1);
  line4->SetTextColor(1);
  line1->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)
  line2->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)
  line3->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)
  line4->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)

  //textbox1->Draw("same");

  gPad->RedrawAxis();

  c1236 ->SetGrid();
  c1236 ->SetLogx();
  c1236 ->SaveAs("plots_pretty/MEAN_PT_RESPONSE_VS_GEN_PT.pdf")
  c1236 ->SetLogx(0);


//------------------------------------------------------------


  TGraph *gr_PT_RESOLUTION_VS_GEN_PT_central = F2->Get("gr_PT_RESOLUTION_VS_GEN_PT_central");
  TGraph *gr_PT_RESOLUTION_VS_GEN_PT_endcap  = F2->Get("gr_PT_RESOLUTION_VS_GEN_PT_endcap");
  TGraph *gr_PT_RESOLUTION_VS_GEN_PT_forward = F2->Get("gr_PT_RESOLUTION_VS_GEN_PT_forward");

  gr_PT_RESOLUTION_VS_GEN_PT_central->SetTitle(";GenJet p_{T} (GeV/c);#sigma(p_{T}^{gen}/p_{T}^{reco})/<p_{T}^{gen}/p_{T}^{reco}>");

  gr_PT_RESOLUTION_VS_GEN_PT_central ->SetLineColor(4);
  gr_PT_RESOLUTION_VS_GEN_PT_endcap  ->SetLineColor(3);
  gr_PT_RESOLUTION_VS_GEN_PT_forward ->SetLineColor(2);
  gr_PT_RESOLUTION_VS_GEN_PT_central ->SetMarkerColor(4);
  gr_PT_RESOLUTION_VS_GEN_PT_endcap  ->SetMarkerColor(3);
  gr_PT_RESOLUTION_VS_GEN_PT_forward ->SetMarkerColor(2);
  gr_PT_RESOLUTION_VS_GEN_PT_central ->SetMarkerStyle(20);
  gr_PT_RESOLUTION_VS_GEN_PT_endcap  ->SetMarkerStyle(21);
  gr_PT_RESOLUTION_VS_GEN_PT_forward ->SetMarkerStyle(22);
  gr_PT_RESOLUTION_VS_GEN_PT_central ->Draw("Ap");
  gr_PT_RESOLUTION_VS_GEN_PT_endcap  ->Draw("p");
  gr_PT_RESOLUTION_VS_GEN_PT_forward ->Draw("p");

  leg = new TLegend(0.55,0.67,0.89,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(gr_PT_RESOLUTION_VS_GEN_PT_central ,"| #eta | < 1.3","LP");
  leg->AddEntry(gr_PT_RESOLUTION_VS_GEN_PT_endcap  ,"1.3 < | #eta | < 3","LP");
  leg->AddEntry(gr_PT_RESOLUTION_VS_GEN_PT_forward ,"3 < | #eta | < 5","LP");
  leg->Draw("same");

  TPaveText *textbox1 = new TPaveText(0.55,0.55,0.89,0.65,"NDC");
  textbox1->SetFillColor(0);
  textbox1->SetLineColor(0);
  TText *line1 = textbox1->AddText("QCD Flat - 72X");
  TText *line2 = textbox1->AddText("AK4 PF Jets");
  TText *line3 = textbox1->AddText("AK5 Gen Jets");
  TText *line4 = textbox1->AddText("CSA14 AK4 Corrections");
  line1->SetTextColor(1);
  line2->SetTextColor(1);
  line3->SetTextColor(1);
  line4->SetTextColor(1);
  line1->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)
  line2->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)
  line3->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)
  line4->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)

  textbox1->Draw("same");

  gPad->RedrawAxis();
  c1236 ->SaveAs("plots_pretty/PT_RESOLUTION_VS_GEN_PT.pdf")


 }
