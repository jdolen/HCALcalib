/*
mkdir plots
mkdir plots/ResolutionPlots_PFJet_NoCorr
mkdir plots/ResolutionPlots_PFJet_Corr
mkdir plots/ResolutionPlots_CaloJet_NoCorr
mkdir plots/ResolutionPlots_CaloJet_Corr
.L plotResolution.C
run()
*/

#include <iostream>
#include <TTree.h>
#include <TH1F.h>
#include <sstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>

void plot_graphs(int, int, int, int ); 
void plot_hists(string, string, string, string, TFile*, TFile*, TFile *,TFile *  ); 
void run()
{

	using namespace std;
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetCanvasDefH(500); //Height of canvas
	gStyle->SetCanvasDefW(800); //Width of canvas
	gStyle->SetCanvasDefX(0);   //POsition on screen
	gStyle->SetCanvasDefY(0);

	// For the Pad:
	gStyle->SetPadBorderMode(0);
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
	gStyle->SetHistLineWidth(2);

	//For the fit/function:
	gStyle->SetOptFit(1);
	gStyle->SetFitFormat("5.4g");
	gStyle->SetFuncColor(2);
	gStyle->SetFuncStyle(1);
	gStyle->SetFuncWidth(1);

	//For the date:
	gStyle->SetOptDate(0);

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

	// Margins:
	gStyle->SetPadTopMargin(0.14);
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetPadLeftMargin(0.17);//0.14);//0.16);
	gStyle->SetPadRightMargin(0.13);//0.3);//0.16);

	// For the Global title:
	gStyle->SetTitleFont(42);
	gStyle->SetTitleColor(1);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleFontSize(0.04);
	gStyle->SetTitleW(0); // Set the width of the title box
	gStyle->SetTitleX(0.5); // Set the position of the title box
	gStyle->SetTitleY(0.94); // Set the position of the title box


	// For the axis titles:
	gStyle->SetTitleColor(1, "XYZ");
	gStyle->SetTitleFont(42, "XYZ");
	gStyle->SetTitleSize(0.05, "XYZ");
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetTitleOffset(1.5, "Y"); // Another way to set the Offset

	// For the axis labels:
	gStyle->SetLabelColor(1, "XYZ");
	gStyle->SetLabelFont(42, "XYZ");
	gStyle->SetLabelOffset(0.007, "YZ");
	gStyle->SetLabelOffset(0.012, "X");
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

	gStyle->SetOptStat(0);
	gROOT->UseCurrentStyle();
	gROOT->ForceStyle();

	gROOT->Reset();
	 gROOT->ForceStyle(); 
	gROOT->SetBatch(); 


	for (int jettype = 0; jettype<=1; jettype++){
		for (int corr = 0; corr<=1; corr++){
			for (int var = 0; var<=2; var++){
				for (int etaregion = 0; etaregion<=1; etaregion++){
					plot_graphs(jettype, corr, var, etaregion ); 
				}
			}
		}
	}
} 

// jettype = 0 "PFJet", 1 "CaloJet"
// corr = 0 (nocorr), 1 (default), 2( rhoArea), 3
// etaregion = 0 (central), 1 (endcap)
void plot_graphs(int jettype, int corr, int var, int etaregion ) 
{
	cout<<"plot_graphs"<<endl;
	string sjet;
	if ( jettype==0 ) sjet = "PFJet";
	if ( jettype==1 ) sjet = "CaloJet";

	string scorr;
	if ( corr==0 ) scorr = "NoCorr";
	if ( corr==1 ) scorr = "Corr";
	if ( corr==2 ) scorr = "CorrPtRhoArea";
	if ( corr==3 ) scorr = "CorrPtRhocentralArea";

	string scorrLabel;
	if ( corr==0 ) scorrLabel = "No PU subtract";
	if ( corr==1 ) scorrLabel = "L1 PU subtract";
	if ( corr==2 ) scorrLabel = "#rho area subtract";
	if ( corr==3 ) scorrLabel = "central calo     #rho area subtract";

	string svar;
	if (var==0) svar = "RESPONSE";
	if (var==1) svar = "RESOLUTION";
	if (var==2) svar = "WIDTH";

	string svar2;
	if (var==0) svar2 = "ABSOLUTE_RESOLUTION_WIDTH";
	if (var==1) svar2 = "ABSOLUTE_RESOLUTION_MEAN";
	if (var==2) svar2 = "ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN";

	string svar3;
	if (var==0) svar3 = "ABSOLUTE_NORM_RESOLUTION_WIDTH";
	if (var==1) svar3 = "ABSOLUTE_NORM_RESOLUTION_MEAN";
	if (var==2) svar3 = "ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN";

	string seta;
	if (etaregion==0) seta = "central";
	if (etaregion==1) seta = "endcap";

	string setaregion;
	if (etaregion==0) setaregion = "| #eta | < 1.3";
	if (etaregion==1) setaregion = "1.3< | #eta | < 3";

	if (var==1 &&  jettype==1 && corr==1 )  cout<<"THIS IS THE ONE THAT IS GOING BAD"<<endl;

	string file_pre = "ResolutionPlots_" + sjet + "_" + scorr + "_";
	string graph_name = "gr_" + svar + "_" + seta;
	string graph_name2 = "gr_" + svar2 + "_" + seta;
	string graph_name3 = "gr_" + svar3 + "_" + seta;

	string sfolder = "ResolutionPlots_" + sjet + "_" + scorr + "/";
	string save_name = "plots/"+sfolder+graph_name+".pdf";

	cout<<save_name<<endl;

	string file_config0               = file_pre + "tt25_Hcal2_config0.root"            ;
	string file_config1               = file_pre + "tt25_Hcal2_config1.root"            ;        
	string file_config2               = file_pre + "tt25_Hcal2_config2correct.root"     ;        

	TFile * F_config0                 = new TFile( file_config0              .c_str() );
	TFile * F_config1                 = new TFile( file_config1              .c_str() );
	TFile * F_config2                 = new TFile( file_config2              .c_str() );
	
	TGraphErrors * gr_config0                = F_config0               ->Get( graph_name.c_str() );
	TGraphErrors * gr_config1                = F_config1               ->Get( graph_name.c_str() );
	TGraphErrors * gr_config2                = F_config2               ->Get( graph_name.c_str() );


	//--------------------------------------------------------------------------------
	// All canvas stuff
	TCanvas *c1236= new TCanvas("c1236","",200,10,800,700);

	TPaveText *textbox2 = new TPaveText(0.14,0.86,0.79,0.91,"NDC");
	textbox2->SetFillStyle(3000);
	textbox2->SetFillColor(0);
	textbox2->SetLineColor(0);
	TText *line2a = textbox2->AddText("CMS Simulation, #sqrt{s} = 13 TeV (25ns, <#mu>=35)");
	line2a->SetTextColor(1);
	line2a->SetTextAlign(12); //first number = horizontal alignment (1 left, 2 center, 3 right). second number =vertical alignment (1 top, 2 center, 3 bottom)
	line2a->SetTextFont(43);
	line2a->SetTextSizePixels(28); //40

	TPaveText *textbox3 = new TPaveText(0.37,0.67,0.57,0.82,"NDC");
	//TPaveText *textbox3 = new TPaveText(0.57,0.52,0.83,0.67,"NDC");
	textbox3->SetFillColor(0);
	textbox3->SetLineColor(0);
	TText *line1 = textbox3->AddText("t#bar{t} - 73X");
	string tb = "AK4 "+sjet;
	TText *line2 = textbox3->AddText(tb.c_str());
	TText *line3 = textbox3->AddText(scorrLabel.c_str());
	TText *line4 = textbox3->AddText(setaregion.c_str());
	line1->SetTextFont(43);
	line2->SetTextFont(43);
	line3->SetTextFont(43);
	line4->SetTextFont(43);
	line1->SetTextSizePixels(23);
	line2->SetTextSizePixels(23);
	line3->SetTextSizePixels(23);
	line4->SetTextSizePixels(23);
	line1->SetTextColor(1);
	line2->SetTextColor(1);
	line3->SetTextColor(1);
	line4->SetTextColor(1);
	line1->SetTextAlign(12); 
	line3->SetTextAlign(12); 
	line2->SetTextAlign(12); 
	line4->SetTextAlign(12); 
	//textbox3->Draw("same");
	//textbox2->Draw("same");


	TPaveText *textbox3b = new TPaveText(0.37,0.27,0.57,0.4,"NDC");
	//TPaveText *textbox3 = new TPaveText(0.57,0.52,0.83,0.67,"NDC");
	textbox3b->SetFillColor(0);
	textbox3b->SetLineColor(0);
	TText *line1 = textbox3b->AddText("t#bar{t} - 73X");
	string tb = "AK4 "+sjet;
	TText *line2 = textbox3b->AddText(tb.c_str());
	TText *line3 = textbox3b->AddText(scorrLabel.c_str());
	TText *line4 = textbox3b->AddText(setaregion.c_str());
	line1->SetTextFont(43);
	line2->SetTextFont(43);
	line3->SetTextFont(43);
	line4->SetTextFont(43);
	line1->SetTextSizePixels(23);
	line2->SetTextSizePixels(23);
	line3->SetTextSizePixels(23);
	line4->SetTextSizePixels(23);
	line1->SetTextColor(1);
	line2->SetTextColor(1);
	line3->SetTextColor(1);
	line4->SetTextColor(1);
	line1->SetTextAlign(12); 
	line3->SetTextAlign(12); 
	line2->SetTextAlign(12); 
	line4->SetTextAlign(12); 
	//textbox3->Draw("same");
	//textbox2->Draw("same");

	//--------------------------------------------------------------------------------
	// Colors and such
	gr_config0              ->SetLineColor(1); 	
	gr_config1              ->SetLineColor(2); 
	gr_config2              ->SetLineColor(kOrange-3);

	gr_config0              ->SetMarkerColor(1);
	gr_config1              ->SetMarkerColor(2);
	gr_config2              ->SetMarkerColor(kOrange-3);

	gr_config0              ->SetMarkerStyle(20);
	gr_config1              ->SetMarkerStyle(21);
	gr_config2              ->SetMarkerStyle(22);
	gr_config2              ->SetMarkerSize(1.3);

	gr_config0              ->SetLineWidth(2);
	gr_config1              ->SetLineWidth(2);
	gr_config2              ->SetLineWidth(2);


// 	if (var==0 && corr == 0)
// 	{
// 		double down = 0.5;
// 		double up = 1.5;
// 		gr_config0              ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config1              ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config1p5            ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_low_plus     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_old_L1       ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_10fC    ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_3fC     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_5fC     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_5fC_P105->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_vlow_plus    ->GetYaxis()->SetRangeUser(down,up);
// 	}

// 	if (var==0 && corr > 0)
// 	{
// 		double down = 0.0;
// 		double up = 1.2;
// 		gr_config0              ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config1              ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config1p5            ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_low_plus     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_old_L1       ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_10fC    ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_3fC     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_5fC     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_5fC_P105->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_vlow_plus    ->GetYaxis()->SetRangeUser(down,up);
// 	}

// 	if (var==1 && corr == 0 && etaregion==0)
// 	{
// 		double down = 0.05;
// 		double up = 0.35;
// 		gr_config0              ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config1              ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config1p5            ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_low_plus     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_old_L1       ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_10fC    ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_3fC     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_5fC     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_5fC_P105->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_vlow_plus    ->GetYaxis()->SetRangeUser(down,up);
// 	}
// 	if (var==2 && corr == 0 && etaregion==0)
// 	{
// 		double down = 0.05;
// 		double up = 0.35;
// 		gr_config0              ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config1              ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config1p5            ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_low_plus     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_old_L1       ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_10fC    ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_3fC     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_5fC     ->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_plus_5fC_P105->GetYaxis()->SetRangeUser(down,up);
// 		gr_config2_vlow_plus    ->GetYaxis()->SetRangeUser(down,up);
// 	}
// 	// gr_config0              ->Draw("ap");
// 	// gr_config1              ->Draw("p");
// 	// gr_config1p5            ->Draw("p");
// 	// gr_config2_low_plus     ->Draw("p");
// 	// gr_config2_old_L1       ->Draw("p");
// 	// gr_config2_plus_10fC    ->Draw("p");
// 	// gr_config2_plus_3fC     ->Draw("p");
// 	// gr_config2_plus_5fC     ->Draw("p");
// 	// gr_config2_plus_5fC_P105->Draw("p");
// 	// gr_config2_vlow_plus    ->Draw("p");
// 	// c1236->SaveAs(save_name.c_str());
// cout<<"debug4a"<<endl;
// 	gr_config1->Draw("ap");
	if (var==1 &&  jettype==1 && corr==1 ) {	string test = "plots/"+sfolder+"_config2.pdf";  	gr_config2->Draw("ap"); c1236->SaveAs(test.c_str());}

	if (var==1 &&  jettype==1 && corr==1 ) {	string test = "plots/"+sfolder+"_config2.pdf";  	gr_config2->Draw("ap"); c1236->SaveAs(test.c_str());}
	if (var==1 &&  jettype==1 && corr==1 ) {	string test = "plots/"+sfolder+"_config0.pdf";  	gr_config0->Draw("ap"); c1236->SaveAs(test.c_str());}
	if (var==1 &&  jettype==1 && corr==1 ) {	string test = "plots/"+sfolder+"_config1.pdf";  	gr_config1->Draw("ap"); c1236->SaveAs(test.c_str());}


	//-----------------------------------------------------------------------------------
	// MoneyPlot
	gr_config0->Draw("ap");
	gr_config1->Draw("p");
	gr_config2->Draw("p");

	leg = new TLegend(0.57,0.62,0.83,0.84);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextFont(43);
	leg->SetTextSizePixels(32);
	leg->SetMargin(0.25);
	leg->AddEntry(gr_config0               , "Method 0", "LP");
	leg->AddEntry(gr_config1               , "Method 1", "LP");
	leg->AddEntry(gr_config2               , "Method 2", "LP");
	leg->Draw("same");

	textbox2->Draw("same");
	textbox3->Draw("same");

	string save_name_money = "plots/"+sfolder+graph_name+"_money.pdf";
			cout<<"save_name_money " <<save_name_money<<endl;

	string save_name_money_log = "plots/"+sfolder+graph_name+"_money_log.pdf";
	c1236 ->SetGrid();
	gPad->RedrawAxis();
	c1236->SaveAs(save_name_money.c_str());

	c1236->SetLogx(1);
	c1236->SaveAs(save_name_money_log.c_str());
	c1236->SetLogx(0);

	
	//---------------------------------------------------------------------------------------------------------------
	// Absolute response

	cout<<"graph_name2 "<<graph_name2<<endl;
	TGraphErrors * gr_abs_config0                = F_config0               ->Get( graph_name2.c_str() );
	TGraphErrors * gr_abs_config1                = F_config1               ->Get( graph_name2.c_str() );
	TGraphErrors * gr_abs_config2                = F_config2               ->Get( graph_name2.c_str() );

	// if (corr == 0 && etaregion==0)
	// {
	// 	double down = -50;
	// 	double up = 5;
	// 	gr_abs_config0              ->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config1              ->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config1p5            ->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config2_low_plus     ->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config2_old_L1       ->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config2_plus_10fC    ->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config2_plus_3fC     ->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config2_plus_5fC     ->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config2_plus_5fC_P105->GetYaxis()->SetRangeUser(down,up);
	// 	gr_abs_config2_vlow_plus    ->GetYaxis()->SetRangeUser(down,up);
	// }

	//--------------------------------------------------------------------------------
	// Colors and such
	gr_abs_config0              ->SetLineColor(1);
	gr_abs_config1              ->SetLineColor(2);
	gr_abs_config2              ->SetLineColor(kOrange-3);

	gr_abs_config0              ->SetMarkerColor(1);
	gr_abs_config1              ->SetMarkerColor(2);
	gr_abs_config2              ->SetMarkerColor(kOrange-3);


	gr_abs_config0              ->SetMarkerStyle(20);
	gr_abs_config1              ->SetMarkerStyle(21);
	gr_abs_config2              ->SetMarkerStyle(22);
	gr_abs_config2              ->SetMarkerSize(1.3);

	gr_abs_config0              ->SetLineWidth(2);
	gr_abs_config1              ->SetLineWidth(2);
	gr_abs_config2              ->SetLineWidth(2);


	//-----------------------------------------------------------------------------------
	// MoneyPlot
	gr_abs_config0->Draw("ap");
	gr_abs_config1->Draw("p");
	gr_abs_config2->Draw("p");
	
	
	leg = new TLegend(0.57,0.62,0.83,0.84);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextFont(43);
	leg->SetTextSizePixels(32);
	leg->SetMargin(0.25);
	leg->AddEntry(gr_abs_config0               , "Method 0", "LP");
	leg->AddEntry(gr_abs_config1               , "Method 1", "LP");
	leg->AddEntry(gr_abs_config2               , "Method 2", "LP");
	leg->Draw("same");

	textbox2->Draw("same");
	textbox3b->Draw("same");

	string save_name_money = "plots/"+sfolder+graph_name2+"_money.pdf";
	string save_name_money_log = "plots/"+sfolder+graph_name2+"_money_log.pdf";
	c1236 ->SetGrid();
	gPad->RedrawAxis();
	c1236->SaveAs(save_name_money.c_str());
	c1236->SetLogx(1);
	c1236->SaveAs(save_name_money_log.c_str());
	c1236->SetLogx(0);


	//---------------------------------------------------------------------------------------------------------------
	// Absolute response norm

	cout<<"graph_name3 "<<graph_name3<<endl;
	TGraphErrors * gr_abs_norm_config0                = F_config0               ->Get( graph_name3.c_str() );
	TGraphErrors * gr_abs_norm_config1                = F_config1               ->Get( graph_name3.c_str() );
	TGraphErrors * gr_abs_norm_config2                = F_config2               ->Get( graph_name3.c_str() );
	

	//--------------------------------------------------------------------------------
	// Colors and such
	gr_abs_norm_config0              ->SetLineColor(1);
	gr_abs_norm_config1              ->SetLineColor(2);
	gr_abs_norm_config2             ->SetLineColor(kOrange-3);

	gr_abs_norm_config0              ->SetMarkerColor(1);
	gr_abs_norm_config1              ->SetMarkerColor(2);
	gr_abs_norm_config2              ->SetMarkerColor(kOrange-3);


	gr_abs_norm_config0              ->SetMarkerStyle(20);
	gr_abs_norm_config1              ->SetMarkerStyle(21);
	gr_abs_norm_config2->SetMarkerStyle(22);
	gr_abs_norm_config2 ->SetMarkerSize(1.3);



	gr_abs_norm_config0              ->SetLineWidth(2);
	gr_abs_norm_config1              ->SetLineWidth(2);
	gr_abs_norm_config2    ->SetLineWidth(2);


	//-----------------------------------------------------------------------------------
	// MoneyPlot
	gr_abs_norm_config0->Draw("ap");
	gr_abs_norm_config1->Draw("p");
	gr_abs_norm_config2->Draw("p");
	
	leg = new TLegend(0.57,0.62,0.83,0.84);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextFont(43);
	leg->SetTextSizePixels(31);
	leg->SetMargin(0.19);
	leg->AddEntry(gr_abs_norm_config0               , "Method 0", "LP");
	leg->AddEntry(gr_abs_norm_config1               , "Method 1", "LP");
	leg->AddEntry(gr_abs_norm_config2      , "Method 2", "LP");
	leg->Draw("same");
	textbox2->Draw("same");
	textbox3->Draw("same");

	string save_name_money = "plots/"+sfolder+graph_name3+"_money.pdf";
	string save_name_money_log = "plots/"+sfolder+graph_name3+"_money_log.pdf";
	c1236 ->SetGrid();
	gPad->RedrawAxis();
	c1236->SaveAs(save_name_money.c_str());
	c1236->SetLogx(1);
	c1236->SaveAs(save_name_money_log.c_str());
	c1236->SetLogx(0);


	//---------------------------------------------------------------------------------------------------------------
	// plot hists
	// string title  = "Gaussian fit of jet p_{T} response";// ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}(reco)/p_{T}(gen);Number of jets";
	// string title2  = "Gaussian fit of p_{T}^{reco}/p_{T}^{gen}";// ;Number of jets";
	// plot_hists("RESPONSE", sfolder, title, title2,  F_config0, F_config1, F_config2);
	// string title  = "Gaussian fit of jet p_{T} response";// ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}(reco)/p_{T}(gen);Number of jets";
	// string title2  = "Gaussian fit of p_{T}^{reco}-p_{T}^{gen}";// ;Number of jets";
	// plot_hists("ABSOLUTE_RESOLUTION", sfolder, title, title2,  F_config0, F_config1, F_config2);
	// string title  = "Gaussian fit of jet p_{T} response";// ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}(reco)/p_{T}(gen);Number of jets";
	// string title2  = "Gaussian fit of (p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}";// ;Number of jets";
	// plot_hists("ABSOLUTE_RESOLUTION_NORM", sfolder, title, title2,  F_config0, F_config1, F_config2);

	
	

}
// plot_money(TGraph * g1,TGraph * g2,TGraph * g3,TGraph * g4,)
// {
// 	g1->Draw("ap");
// 	g2->Draw("p");
// 	if ( jettype==0 ) g3->Draw("p");
// 	g4->Draw("p");
	
// 	leg = new TLegend(0.57,0.62,0.83,0.84);
// 	leg->SetBorderSize(0);
// 	leg->SetFillColor(0);
// 	leg->SetTextFont(43);
// 	leg->SetTextSizePixels(32);
// 	leg->SetMargin(0.25);
// 	leg->AddEntry(g1               , "Method 0", "LP");
// 	leg->AddEntry(g2               , "Method 1", "LP");
// 	if ( jettype==0 ) leg->AddEntry(g3             , "Method 1.5", "LP");
// 	leg->AddEntry(g4      , "Method 2", "LP");
// 	leg->Draw("same");

// 	textbox2->Draw("same");
// 	textbox3->Draw("same");

// 	string save_name_money = "plots/"+sfolder+graph_name+"_money.pdf";
// 	string save_name_money_log = "plots/"+sfolder+graph_name+"_money_log.pdf";
// 	c1236 ->SetGrid();
// 	gPad->RedrawAxis();
// 	c1236->SaveAs(save_name_money.c_str());
// 	c1236->SetLogx(1);
// 	c1236->SaveAs(save_name_money_log.c_str());
// 	c1236->SetLogx(0);

// }
plot_hists(string name, string sfolder, string title1, string title2, TFile * F1, TFile * F2, TFile * F3){
	const int nbins_eta = 3;
	const int nbins_pt = 9;
	double etabins[nbins_eta]  =  {0.0, 1.3, 3};//, 5};
	double ptbins[nbins_pt]    =  {10., 20., 40., 60., 80., 100., 150., 200., 300. };
	string histName_etabins[nbins_eta] = {"0", "1p3", "3p0"};//, "5p0"};
	string histName_ptbins[nbins_pt]   = {"10", "20", "40", "60", "80", "100", "150", "200", "300"};
	string histTitle_etabins[nbins_eta] = {"0", "1.3", "3.0"};//, "5.0"};
	string histTitle_ptbins[nbins_pt]   = {"10", "20", "40", "60", "80", "100", "150", "200", "300" };
	
	for (int etabin =0; etabin< nbins_eta-1; etabin++)
  	{
		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
		{		  			
	  		string histoname   = name+"_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
			cout<<histoname<<endl;
			TH1F * h0              = F1              ->Get( histoname.c_str() ); 
			TH1F * h1              = F2             ->Get( histoname.c_str() );
			TH1F * h2              = F3             ->Get( histoname.c_str() );

			TH1F * myh0              = h0->Clone(); 
			TH1F * myh1              = h1->Clone();
			TH1F * myh2              = h2->Clone();

			// h0->Rebin(8);
			// h1->Rebin(8);
			// h2->Rebin(8);

			h0->SetLineColor(1);
			h1->SetLineColor(2);
			h2->SetLineColor(kOrange-3);

			// myh0->SetLineColor(1);
			// myh1->SetLineColor(2);
			// myh2->SetLineColor(kOrange-3);

			// myh0->Rebin(8);
			// myh1->Rebin(8);
			// myh2->Rebin(8);
			// 			myh2->Draw("");

			// myh0->Draw("same");
			// myh1->Draw("same");
			// myh2->Draw("same");

			// string save_name = "plots/"+sfolder+histoname+"_rebin.pdf";
			// c1236->SaveAs(save_name.c_str());

			TF1 * fit0 = h0->GetFunction("Gaussian"); cout<<"here1"<<endl;
			TF1 * fit1 = h1->GetFunction("Gaussian"); cout<<"here2"<<endl;
			TF1 * fit2 = h2->GetFunction("Gaussian"); cout<<"here3"<<endl;


			double fit0_mean   = fit0->GetParameter(1);
			double fit0_width  = fit0->GetParameter(2);
			double fit1_mean   = fit1->GetParameter(1);
			double fit1_width  = fit1->GetParameter(2);
			double fit2_mean   = fit2->GetParameter(1);
			double fit2_width  = fit2->GetParameter(2);


			stringstream t_fit0_mean  ;
			stringstream t_fit0_width ;
			stringstream t_fit1_mean  ;
			stringstream t_fit1_width ;
			stringstream t_fit2_mean  ;
			stringstream t_fit2_width ;
			t_fit0_mean  << fit0_mean ;
			t_fit0_width << fit0_width;
			t_fit1_mean  << fit1_mean ;
			t_fit1_width << fit1_width;
			t_fit2_mean  << fit2_mean ;
			t_fit2_width << fit2_width;
			string s_fit0_mean  = t_fit0_mean .str();
			string s_fit0_width = t_fit0_width.str();
			string s_fit1_mean  = t_fit1_mean .str();
			string s_fit1_width = t_fit1_width.str();
			string s_fit2_mean  = t_fit2_mean .str();
			string s_fit2_width = t_fit2_width.str();
			string text_fit0_mean    =  "Method 0 - mean "+s_fit0_mean;
			string text_fit1_mean    =  "Method 1 - mean "+s_fit1_mean;
			string text_fit2_mean    =  "Method 2 - mean  "+s_fit2_mean;
			string text_fit0_width   =  "Method 0 - width "+s_fit0_width;
			string text_fit1_width   =  "Method 1 - width "+s_fit1_width;
			string text_fit2_width   =  "Method 2 - width "+s_fit2_width;
		

			TPaveText *textbox = new TPaveText(0.35,0.16,0.6,0.36,"NDC");
			textbox->SetFillColor(0);
			textbox->SetLineColor(0);
			TText *line1 = textbox->AddText(text_fit0_mean.c_str());
			TText *line2 = textbox->AddText(text_fit0_width.c_str());
			TText *line3 = textbox->AddText(text_fit1_mean.c_str());
			TText *line4 = textbox->AddText(text_fit1_width.c_str());
			TText *line5 = textbox->AddText(text_fit2_mean.c_str());
			TText *line6 = textbox->AddText(text_fit2_width.c_str());
			line1->SetTextFont(43);
			line2->SetTextFont(43);
			line3->SetTextFont(43);
			line4->SetTextFont(43);
			line5->SetTextFont(43);
			line6->SetTextFont(43);
			line1->SetTextSizePixels(20);
			line2->SetTextSizePixels(20);
			line3->SetTextSizePixels(20);
			line4->SetTextSizePixels(20);
			line5->SetTextSizePixels(20);
			line6->SetTextSizePixels(20);
			line1->SetTextColor(1);
			line2->SetTextColor(1);
			line3->SetTextColor(2);
			line4->SetTextColor(2);
			line5->SetTextColor(kOrange-3);
			line6->SetTextColor(kOrange-3);
			line1->SetTextAlign(12); 
			line3->SetTextAlign(12); 
			line2->SetTextAlign(12); 
			line4->SetTextAlign(12); 
			line5->SetTextAlign(12); 
			line6->SetTextAlign(12); 

			TH1F *htemp0 = fit0->GetHistogram();
			TH1F *htemp1 = fit1->GetHistogram();
			TH1F *htemp2 = fit2->GetHistogram();
			// fit0->SetLineColor(1);
			// fit1->SetLineColor(2);
			// fit2->SetLineColor(kOrange-3);
			// h0->Scale(h0->Integral()/htemp0->Integral())
			// h1->Scale(h0->Integral()/htemp0->Integral())
			// h2->Scale(h0->Integral()/htemp0->Integral())
			// htemp0->Scale(htemp0->Integral()/htemp0->Integral())
			// htemp1->Scale(htemp1->Integral()/htemp0->Integral())
			// htemp2->Scale(htemp2->Integral()/htemp0->Integral())

			cout<<htemp1->FindFirstBinAbove(5.0)<<endl;
			cout<<htemp1->FindFirstBinAbove(20.0)<<endl;
			cout<<htemp1->FindLastBinAbove(5.0)<<endl;
			cout<<htemp1->FindLastBinAbove(20.0)<<endl;
			cout<<htemp1->GetMaximumBin()<<endl;

			// h0->GetXaxis()->SetRange(htemp1->FindFirstBinAbove(10.0,1), htemp1->GetMaximumBin()*4 );
			// h1->GetXaxis()->SetRange(htemp1->FindFirstBinAbove(10.0,1), htemp1->GetMaximumBin()*4 );
			// h2->GetXaxis()->SetRange(htemp1->FindFirstBinAbove(10.0,1), htemp1->GetMaximumBin()*4 );

			h2->Draw("HIST");
			h0->Draw("HISTsame");
			h1->Draw("HISTsame");
			h2->Draw("HISTsame");
			htemp0->Draw("HISTsame");
			htemp1->Draw("HISTsame");
			htemp2->Draw("HISTsame");
			string save_name = "plots/"+sfolder+histoname+".pdf";
			c1236->SaveAs(save_name.c_str());
			
			string histotitle  ="";
			//histotitle =  title1 + " ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );"+title2+";Normalized";
			histotitle =  histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c;"+title2;

			htemp0->SetTitle(histotitle.c_str());
			htemp1->SetTitle(histotitle.c_str());
			htemp2->SetTitle(histotitle.c_str());
			// htemp0->GetXaxis()->SetRangeUser(0,4);
			// htemp1->GetXaxis()->SetRangeUser(0,4);
			// htemp2->GetXaxis()->SetRangeUser(0,4);
			htemp0->DrawNormalized("");
			htemp1->DrawNormalized("same");
			htemp2->DrawNormalized("same");
			textbox->Draw("same");
			string save_name = "plots/"+sfolder+histoname+"_fit.pdf";
			c1236->SaveAs(save_name.c_str());

			delete textbox;
			delete fit0;
			delete fit1;
			delete fit2;
			delete h0;
			delete h1;
			delete h2;

		}
	}
}
// 	cout<<"ABSOLUTE_RESOLUTION_ETA_"<<endl;
// 	for (int etabin =0; etabin< nbins_eta-1; etabin++)
//   	{
// 		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
// 		{		  			
// 	  		string histoname   = "ABSOLUTE_RESOLUTION_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
// 			cout<<histoname<<endl;
// 			TH1F * hh0              = F_config0                           ->Get( histoname.c_str() );
// 			TH1F * hh1              = F_config1                           ->Get( histoname.c_str() );
// 			TH1F * hh2              = F_config2_plus_5fC_P105             ->Get( histoname.c_str() );

// 			hh0->Rebin(8);
// 			hh1->Rebin(8);
// 			hh2->Rebin(8);


// 			hh0->SetLineColor(1);
// 			hh1->SetLineColor(2);
// 			hh2->SetLineColor(kOrange-3);
// 			hh2->DrawNormalized("");
// 			hh0->DrawNormalized("same");
// 			hh1->DrawNormalized("same");
// 			hh2->DrawNormalized("same");
// 			string save_name = "plots/"+sfolder+histoname+".pdf";
// 			c1236->SaveAs(save_name.c_str());
// 			delete hh0;
// 			delete hh1;
// 			delete hh2;
// 		}
// 	}
// 	cout<<"ABSOLUTE_RESOLUTION_NORM_ETA_"<<endl;

// 	for (int etabin =0; etabin< nbins_eta-1; etabin++)
//   	{
// 		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
// 		{		  			
// 	  		string histoname   = "ABSOLUTE_RESOLUTION_NORM_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
// 			cout<<histoname<<endl;

// 			// for some reason new hist name needed to avoid Error: Incorrect assignment to h0, wrong type 'TObject*' plotPretty.C:706:

// 			TH1F * hhh0              = F_config0                           ->Get( histoname.c_str() );
// 			TH1F * hhh1              = F_config1                           ->Get( histoname.c_str() );
// 			TH1F * hhh2              = F_config2_plus_5fC_P105             ->Get( histoname.c_str() );

// 			hhh0 ->Rebin(8);
// 			hhh1 ->Rebin(8);
// 			hhh2 ->Rebin(8);

// 			hhh0->SetLineColor(1);
// 			hhh1->SetLineColor(2);
// 			hhh2->SetLineColor(kOrange-3);
// 			hhh2->DrawNormalized("");
// 			hhh0->DrawNormalized("same");
// 			hhh1->DrawNormalized("same");
// 			hhh2->DrawNormalized("same");



// 			string save_name = "plots/"+sfolder+histoname+".pdf";
// 			c1236->SaveAs(save_name.c_str());
// 			delete hhh0;
// 			delete hhh1;
// 			delete hhh2;
// 		}
// 	}


// }