{

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


	string file_pre = "ResolutionPlots_PFJet_Corr_";
	string sjet = "PFJet";
	string scorr = "Corr";
	string scorrLabel = "L1 PU subtract";
	string graph_name = "gr_RESOLUTION_central";
	string setaregion = "| #eta | < 1.3";

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

	//-----------------------------------------------------------------------------------
	// Plot
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

	string save_name = "relative_resolution.pdf";
	string save_name_log = "relative_resolution_log.pdf";
	c1236 ->SetGrid();
	gPad->RedrawAxis();
	c1236->SaveAs(save_name.c_str());

	c1236->SetLogx(1);
	c1236->SaveAs(save_name_log.c_str());
	c1236->SetLogx(0);


}