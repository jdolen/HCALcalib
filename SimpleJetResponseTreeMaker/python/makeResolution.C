// .L makeResolution.C+
// run()

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
#include <TColor.h>

void make_resolution(string, Color_t);
void make_pfjet_resolution(string, int, Color_t);
void make_calojet_resolution(string, int, Color_t);

void run()
{
	make_resolution("tt25_Hcal2_config0.root", 1);
	make_resolution("tt25_Hcal2_config1.root", 2);
	make_resolution("tt25_Hcal2_config2correct.root", kOrange-3);
	// make_resolution("tt25_cfg0.root", 1);
	// make_resolution("tt25_cfg1.root", 2);
	// make_resolution("tt25_cfg1p5_L1.root", 3);
	// make_resolution("tt25_config2_low_plus.root", 4);
	// make_resolution("tt25_config2_old_L1.root", 6);
	// make_resolution("tt25_config2_plus_10fC.root", 7);
	// make_resolution("tt25_config2_plus_3fC.root", 8);
	// make_resolution("tt25_config2_plus_5fC.root", 9);
	// make_resolution("tt25_config2_plus_5fC_P105.root", kOrange-3);
	// make_resolution("tt25_config2_vlow_plus.root", 10);

}
void make_resolution(string input_file, Color_t color)
{
	for (int i=0; i<2; i++){
		make_pfjet_resolution(input_file, i, color) ;
	}
	for (int i=0; i<2; i++){
		make_calojet_resolution(input_file, i, color) ;
	}
}
void make_pfjet_resolution(string input_file, int numerator, Color_t color) {

	gROOT->SetBatch(); 
  	TCanvas *c1236= new TCanvas("c1236","",200,10,800,700);

	TFile *Out;
	string output_file;
	if (numerator==0) output_file = "ResolutionPlots_PFJet_NoCorr_"+input_file;
	if (numerator==1) output_file = "ResolutionPlots_PFJet_Corr_"+input_file;

	cout<<output_file<<endl;
	Out = new TFile(output_file.c_str() ,"RECREATE");
	TFile *F1   = new TFile(input_file.c_str() );

	// Get Tree entries
	Float_t PFJet_Pt               ;
	Float_t PFJet_CorrPt           ;
	Float_t PFJet_CorrPtRhoArea    ;
	Float_t PFJet_Eta              ;
	Float_t PFJet_Phi              ;
	Float_t PFJet_MatchedGenJet_Pt              ;
	Float_t PFJet_MatchedGenJet_Eta             ;
	Float_t PFJet_MatchedGenJet_Phi             ;

	Float_t PFJet_Nconst;
	Float_t PFJet_chargedMultiplicity          ;
	Float_t PFJet_chargedEmEnergyFraction      ;
	Float_t PFJet_chargedHadronEnergyFraction  ;
	Float_t PFJet_chargedMuEnergyFraction      ;
	Float_t PFJet_neutralEmEnergyFraction      ;
	Float_t PFJet_neutralHadronEnergyFraction  ;

	TTree *T1    = (TTree*)  F1     ->Get("tree/PFJetTree");
	T1->SetBranchAddress("PFJet_Pt"                          , & PFJet_Pt                          );
	T1->SetBranchAddress("PFJet_CorrPt"                      , & PFJet_CorrPt                      );
	T1->SetBranchAddress("PFJet_CorrPtRhoArea"               , & PFJet_CorrPtRhoArea               );
	T1->SetBranchAddress("PFJet_Eta"                         , & PFJet_Eta                         );
	T1->SetBranchAddress("PFJet_Phi"                         , & PFJet_Phi                         );
	T1->SetBranchAddress("PFJet_MatchedGenJet_Pt"            , & PFJet_MatchedGenJet_Pt            );
	T1->SetBranchAddress("PFJet_MatchedGenJet_Eta"           , & PFJet_MatchedGenJet_Eta           );
	T1->SetBranchAddress("PFJet_MatchedGenJet_Phi"           , & PFJet_MatchedGenJet_Phi           );
	T1->SetBranchAddress("PFJet_Nconst"                      , & PFJet_Nconst                      );
	T1->SetBranchAddress("PFJet_chargedMultiplicity"         , & PFJet_chargedMultiplicity         );
	T1->SetBranchAddress("PFJet_chargedEmEnergyFraction"     , & PFJet_chargedEmEnergyFraction     );
	T1->SetBranchAddress("PFJet_chargedHadronEnergyFraction" , & PFJet_chargedHadronEnergyFraction );
	T1->SetBranchAddress("PFJet_chargedMuEnergyFraction"     , & PFJet_chargedMuEnergyFraction     );
	T1->SetBranchAddress("PFJet_neutralEmEnergyFraction"     , & PFJet_neutralEmEnergyFraction     );
	T1->SetBranchAddress("PFJet_neutralHadronEnergyFraction" , & PFJet_neutralHadronEnergyFraction );


	double entries = T1->GetEntries();
	cout<<"entries = "<< entries <<endl;

	 // // Define genJet pt an eta bins
	static const int nbins_eta = 3;
	static const int nbins_pt = 16;
	double etabins[nbins_eta]           = {0.0, 1.3, 3};
	double ptbins[nbins_pt]             = {10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100., 125., 150., 200., 300.};
	string histName_etabins[nbins_eta]  = {"0", "1p3", "3p0"};
	string histName_ptbins[nbins_pt]    = {"010", "015", "020", "025", "030", "040", "050", "060", "070", "080", "090", "100", "125", "150", "200", "300"};
	string histTitle_etabins[nbins_eta] = {"0", "1.3", "3.0"};
	string histTitle_ptbins[nbins_pt]   = {"10", "15", "20", "25", "30", "40", "50", "60", "70", "80", "90", "100", "125", "150", "200", "300"};

	// static const int nbins_eta = 3;
	// static const int nbins_pt = 9;
	// double etabins[nbins_eta]  =  {0.0, 1.3, 3};//, 5};
	// double ptbins[nbins_pt]    =  {10., 20., 40., 60., 80., 100., 150., 200., 300. };
	// string histName_etabins[nbins_eta] = {"0", "1p3", "3p0"};//, "5p0"};
	// string histName_ptbins[nbins_pt]   = {"10", "20", "40", "60", "80", "100", "150", "200", "300"};
	// string histTitle_etabins[nbins_eta] = {"0", "1.3", "3.0"};//, "5.0"};
	// string histTitle_ptbins[nbins_pt]   = {"10", "20", "40", "60", "80", "100", "150", "200", "300" };
	



	// Arrays for TGraph x and y positions and uncertainties

	Float_t ptbins_x  [nbins_pt-1] = {};
	Float_t ptbins_ex [nbins_pt-1] = {};

	// Relataive resolution 
	Float_t response_mean_y    [nbins_eta-1][nbins_pt-1] = {};
	Float_t response_mean_ey   [nbins_eta-1][nbins_pt-1] = {};

	Float_t response_width_y   [nbins_eta-1][nbins_pt-1] = {};
	Float_t response_width_ey  [nbins_eta-1][nbins_pt-1] = {};
	
	Float_t resolution_y       [nbins_eta-1][nbins_pt-1] = {};
	Float_t resolution_ey      [nbins_eta-1][nbins_pt-1] = {};
 
	// Absolute resolution
	Float_t absolute_resolution_mean_y              [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_resolution_mean_ey             [nbins_eta-1][nbins_pt-1] = {};
 
	Float_t absolute_resolution_width_y             [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_resolution_width_ey            [nbins_eta-1][nbins_pt-1] = {};

	Float_t absolute_resolution_width_over_mean_y   [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_resolution_width_over_mean_ey  [nbins_eta-1][nbins_pt-1] = {};

	// Absoulte normalized by response
	Float_t absolute_norm_resolution_width_y             [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_norm_resolution_width_ey            [nbins_eta-1][nbins_pt-1] = {};

	Float_t absolute_norm_resolution_mean_y        [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_norm_resolution_mean_ey       [nbins_eta-1][nbins_pt-1] = {};

	Float_t absolute_norm_resolution_width_over_mean_y   [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_norm_resolution_width_over_mean_ey  [nbins_eta-1][nbins_pt-1] = {};

	// response histogram array
	TH1F * PT_RESPONSE[nbins_eta-1][nbins_pt-1];
	TH1F * ABSOLUTE_RESOLUTION[nbins_eta-1][nbins_pt-1];
	TH1F * ABSOLUTE_RESOLUTION_NORM[nbins_eta-1][nbins_pt-1];
	// For each eta and pt bin loop over the jet tree, fill the response histogram, and then fit
	for (int etabin =0; etabin< nbins_eta-1; etabin++)
  	{
		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
		{		  			
	  		string histoname   = "RESPONSE_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
			string histotitle  = "Jet p_{T} response ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}^{reco}/p_{T}^{gen};Number of jets";
			PT_RESPONSE[etabin][ptbin]  = new TH1F(histoname.c_str() ,histotitle.c_str(),400,0,4); 
			cout<<histoname <<endl;

			string histoname2   = "ABSOLUTE_RESOLUTION_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
			string histotitle2  = "Jet p_{T} absolute resolution ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}^{reco}-p_{T}^{gen};Number of jets";
			ABSOLUTE_RESOLUTION[etabin][ptbin]  = new TH1F(histoname2.c_str() ,histotitle2.c_str(),400,-50,60); 

			string histoname3   = "ABSOLUTE_RESOLUTION_NORM_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
			string histotitle3  = "Jet p_{T} absolute response ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );(p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen};Number of jets";
			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]  = new TH1F(histoname3.c_str() ,histotitle3.c_str(),400,-2,4); 


			for (int i=0; i<entries; i++ )
		 	{  		    
		 		T1->GetEntry(i);

		 		// all-eta PF Jet ID
		 		if ( PFJet_neutralHadronEnergyFraction > 0.99 ) continue;
		 		if ( PFJet_neutralEmEnergyFraction     > 0.99 ) continue;
		 		if ( PFJet_Nconst                     <= 1    ) continue;
		 		if ( PFJet_chargedMuEnergyFraction     > 0.8  ) continue; 

		 		// central PF Jet ID
		 		if ( fabs(PFJet_Eta)<2.4 && PFJet_chargedEmEnergyFraction      > 0.99 ) continue;
		 		if ( fabs(PFJet_Eta)<2.4 && PFJet_chargedMultiplicity         <= 0    ) continue;
		 		if ( fabs(PFJet_Eta)<2.4 && PFJet_chargedHadronEnergyFraction <= 0    ) continue;
		 		if ( fabs(PFJet_Eta)<2.4 && PFJet_chargedEmEnergyFraction      > 0.99 ) continue;
     
     			// Fill PF Jet Response histogram for the selected eta and pt bin
		  		if ( fabs(PFJet_MatchedGenJet_Eta)>etabins[etabin] && fabs(PFJet_MatchedGenJet_Eta)<etabins[etabin+1] && PFJet_MatchedGenJet_Pt > ptbins[ptbin] && PFJet_MatchedGenJet_Pt < ptbins[ptbin+1] )
		  		{		 		
	  				//if (PFJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(PFJet_CorrPt/PFJet_MatchedGenJet_Pt);
	  				if (numerator ==0 && PFJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(PFJet_Pt/PFJet_MatchedGenJet_Pt);
	  				if (numerator ==1 && PFJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(PFJet_CorrPt/PFJet_MatchedGenJet_Pt);
	  				if (numerator ==2 && PFJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(PFJet_CorrPtRhoArea/PFJet_MatchedGenJet_Pt);
	  				if (numerator ==0 && PFJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( PFJet_Pt-PFJet_MatchedGenJet_Pt                     );
	  				if (numerator ==1 && PFJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( PFJet_CorrPt-PFJet_MatchedGenJet_Pt                 );
	  				if (numerator ==2 && PFJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( PFJet_CorrPtRhoArea-PFJet_MatchedGenJet_Pt          );
	  				if (numerator ==0 && PFJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( PFJet_Pt-PFJet_MatchedGenJet_Pt            ) / PFJet_MatchedGenJet_Pt          );
	  				if (numerator ==1 && PFJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( PFJet_CorrPt-PFJet_MatchedGenJet_Pt        ) / PFJet_MatchedGenJet_Pt          );
	  				if (numerator ==2 && PFJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( PFJet_CorrPtRhoArea-PFJet_MatchedGenJet_Pt ) / PFJet_MatchedGenJet_Pt          );
	  			}
		 	}

			ptbins_x  [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
			ptbins_ex [ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;

		 	//----------------------------------------
		 	// Fit response (fit twice, first to find rough mean and sigma, second with correct range)
			TF1 *gfit = new TF1("Gaussian","gaus",0.4,3); 
			gfit->SetLineColor(color);
			gfit->SetLineWidth(2.5);

			PT_RESPONSE[etabin][ptbin]->Fit(gfit,"RQ"); // Fit histogram in range 
			double amp    = gfit->GetParameter(0);
			double eamp   = gfit->GetParError(0); 
			double mean   = gfit->GetParameter(1);
			double emean  = gfit->GetParError(1); 
			double width  = gfit->GetParameter(2);
			double ewidth = gfit->GetParError(2); 

			TF1 *gfit0 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
			gfit0->SetLineColor(color);
			gfit0->SetLineWidth(2.5);

			PT_RESPONSE[etabin][ptbin]->Fit(gfit0,"RQ"); // Fit histogram in range  R
			amp    = gfit0->GetParameter(0);
			eamp   = gfit0->GetParError(0); 
			mean   = gfit0->GetParameter(1);
			emean  = gfit0->GetParError(1); 
			width  = gfit0->GetParameter(2);
			ewidth = gfit0->GetParError(2); 

			response_mean_y [etabin][ptbin] = mean;
			response_mean_ey[etabin][ptbin] = emean;
			
			response_width_y [etabin][ptbin] = width;
			response_width_ey[etabin][ptbin] = ewidth ;
			
			resolution_y [etabin][ptbin] = width/mean;
			resolution_ey[etabin][ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			
		 	//----------------------------------------
			// Fit Absolute Resolution (fit twice, first to find rough mean and sigma, second with correct range)
			TF1 *gfit2 = new TF1("Gaussian","gaus");//,0.4,3); 
			gfit2->SetLineColor(color);
			gfit2->SetLineWidth(2.5);

			ABSOLUTE_RESOLUTION[etabin][ptbin]->Fit(gfit2,"Q"); // Fit histogram in range  R
			amp    = gfit2->GetParameter(0);
			eamp   = gfit2->GetParError(0); 
			mean   = gfit2->GetParameter(1);
			emean  = gfit2->GetParError(1); 
			width  = gfit2->GetParameter(2);
			ewidth = gfit2->GetParError(2); 

			TF1 *gfit3 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
			gfit3->SetLineColor(color);
			gfit3->SetLineWidth(2.5);

			ABSOLUTE_RESOLUTION[etabin][ptbin]->Fit(gfit3,"RQ"); // Fit histogram in range  R
			amp    = gfit3->GetParameter(0);
			eamp   = gfit3->GetParError(0); 
			mean   = gfit3->GetParameter(1);
			emean  = gfit3->GetParError(1); 
			width  = gfit3->GetParameter(2);
			ewidth = gfit3->GetParError(2); 

			absolute_resolution_mean_y [etabin][ptbin] = 1+mean;
			absolute_resolution_mean_ey[etabin][ptbin] = emean;

			absolute_resolution_width_over_mean_y [etabin][ptbin] = width/(1+mean);
			absolute_resolution_width_over_mean_ey[etabin][ptbin] = width/(1+mean)*sqrt( (emean/(1+mean))*(emean/(1+mean)) + (ewidth/width)*(ewidth/width)   );
			
			absolute_resolution_width_y [etabin][ptbin] = width;
			absolute_resolution_width_ey[etabin][ptbin] = ewidth ;
			
		 	//----------------------------------------
			// Fit Absolute Resolution (fit twice, first to find rough mean and sigma, second with correct range)

			TF1 *gfit4 = new TF1("Gaussian","gaus");//,0.4,3); 
			gfit4->SetLineColor(color);
			gfit4->SetLineWidth(2.5);

			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fit(gfit4,"Q"); // Fit histogram in range  R
			amp    = gfit4->GetParameter(0);
			eamp   = gfit4->GetParError(0); 
			mean   = gfit4->GetParameter(1);
			emean  = gfit4->GetParError(1); 
			width  = gfit4->GetParameter(2);
			ewidth = gfit4->GetParError(2); 

			TF1 *gfit5 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
			gfit5->SetLineColor(color);
			gfit5->SetLineWidth(2.5);

			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fit(gfit5,"RQ"); // Fit histogram in range  R
			amp    = gfit5->GetParameter(0);
			eamp   = gfit5->GetParError(0); 
			mean   = gfit5->GetParameter(1);
			emean  = gfit5->GetParError(1); 
			width  = gfit5->GetParameter(2);
			ewidth = gfit5->GetParError(2); 

			// central response array
			absolute_norm_resolution_mean_y [etabin][ptbin] = 1+mean;
			absolute_norm_resolution_mean_ey[etabin][ptbin] = emean;

   			// central resoabsolute_norm_lution array
			absolute_norm_resolution_width_over_mean_y [etabin][ptbin] = width/(1+mean);
			absolute_norm_resolution_width_over_mean_ey[etabin][ptbin] = width/(1+mean)*sqrt( (emean/(1+mean))*(emean/(1+mean)) + (ewidth/width)*(ewidth/width)   );
			
			// central widtabsolute_norm_h array
			absolute_norm_resolution_width_y [etabin][ptbin] = width;
			absolute_norm_resolution_width_ey[etabin][ptbin] = ewidth ;
	  	}
	}
	// Relative resolution
	TGraphErrors * gr_RESPONSE_central = new TGraphErrors(nbins_pt-1,ptbins_x,response_mean_y[0],ptbins_ex,response_mean_ey[0]);
   	TGraphErrors * gr_RESPONSE_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,response_mean_y[1],ptbins_ex,response_mean_ey[1]);
	string title_response = ";GenJet p_{T} (GeV/c);<p_{T}^{pfjet}/p_{T}^{gen}>";
   	gr_RESPONSE_central ->SetTitle( title_response.c_str() );
   	gr_RESPONSE_endcap  ->SetTitle( title_response.c_str() );

   	TGraphErrors * gr_RESOLUTION_central = new TGraphErrors(nbins_pt-1,ptbins_x,resolution_y[0],ptbins_ex,resolution_ey[0]);
   	TGraphErrors * gr_RESOLUTION_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,resolution_y[1],ptbins_ex,resolution_ey[1]);
	string title_resolution = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{pfjet}/p_{T}^{gen})/<p_{T}^{pfjet}/p_{T}^{gen}>";
   	gr_RESOLUTION_central->SetTitle( title_resolution.c_str() );
   	gr_RESOLUTION_endcap ->SetTitle( title_resolution.c_str() );

	TGraphErrors * gr_WIDTH_central = new TGraphErrors(nbins_pt-1,ptbins_x,response_width_y[0],ptbins_ex,response_width_ey[0]);
   	TGraphErrors * gr_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,response_width_y[1],ptbins_ex,response_width_ey[1]);
	string title_width = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{pfjet}/p_{T}^{gen})";
   	gr_WIDTH_central->SetTitle( title_width.c_str() );
   	gr_WIDTH_endcap ->SetTitle( title_width.c_str() );

	// Absolute resolution
	TGraphErrors * gr_ABSOLUTE_RESOLUTION_MEAN_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_mean_y[0],ptbins_ex,absolute_resolution_mean_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_RESOLUTION_MEAN_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_mean_y[1],ptbins_ex,absolute_resolution_mean_ey[1]);
	title_response = ";GenJet p_{T} (GeV/c);1+<p_{T}^{pfjet}-p_{T}^{gen}>";
   	gr_ABSOLUTE_RESOLUTION_MEAN_central ->SetTitle( title_response.c_str() );
   	gr_ABSOLUTE_RESOLUTION_MEAN_endcap  ->SetTitle( title_response.c_str() );

   	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_width_over_mean_y[0],ptbins_ex,absolute_resolution_width_over_mean_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_width_over_mean_y[1],ptbins_ex,absolute_resolution_width_over_mean_ey[1]);
	title_resolution = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{pfjet}-p_{T}^{gen})/(1+<p_{T}^{pfjet}-p_{T}^{gen}>)";
   	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central->SetTitle( title_resolution.c_str() );
   	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap ->SetTitle( title_resolution.c_str() );

	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_width_y[0],ptbins_ex,absolute_resolution_width_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_width_y[1],ptbins_ex,absolute_resolution_width_ey[1]);
	title_width = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{pfjet}-p_{T}^{gen})";
   	gr_ABSOLUTE_RESOLUTION_WIDTH_central->SetTitle( title_width.c_str() );
   	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap ->SetTitle( title_width.c_str() );





	// Make Tgraphs
	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_mean_y[0],ptbins_ex,absolute_norm_resolution_mean_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_mean_y[1],ptbins_ex,absolute_norm_resolution_mean_ey[1] );
	title_response = ";GenJet p_{T} (GeV/c);1+ < (p_{T}^{pfjet}-p_{T}^{gen}) / p_{T}^{gen}>";
   	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central ->SetTitle( title_response.c_str() );
   	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap  ->SetTitle( title_response.c_str() );

   	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_width_over_mean_y[0],ptbins_ex,absolute_norm_resolution_width_over_mean_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_width_over_mean_y[1],ptbins_ex,absolute_norm_resolution_width_over_mean_ey[1]);
	title_resolution = ";GenJet p_{T} (GeV/c);#sigma( (p_{T}^{pfjet}-p_{T}^{gen}) / p_{T}^{gen} )/(1+< (p_{T}^{pfjet}-p_{T}^{gen})/p_{T}^{gen}>)";
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central->SetTitle( title_resolution.c_str() );
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap ->SetTitle( title_resolution.c_str() );

	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_width_y[0],ptbins_ex,absolute_norm_resolution_width_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_width_y[1],ptbins_ex,absolute_norm_resolution_width_ey[1]);
	title_width = ";GenJet p_{T} (GeV/c);#sigma( (p_{T}^{pfjet}-p_{T}^{gen})/p_{T}^{gen} )";
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central->SetTitle( title_width.c_str() );
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap ->SetTitle( title_width.c_str() );




 //   	// Save TGraphs to output file	
	Out->cd();
	gr_RESPONSE_central      -> SetName("gr_RESPONSE_central"  );
	gr_RESPONSE_endcap       -> SetName("gr_RESPONSE_endcap"   );
	gr_RESOLUTION_central    -> SetName("gr_RESOLUTION_central");
	gr_RESOLUTION_endcap     -> SetName("gr_RESOLUTION_endcap" );
    gr_WIDTH_central         -> SetName("gr_WIDTH_central");
	gr_WIDTH_endcap          -> SetName("gr_WIDTH_endcap" );
   	gr_ABSOLUTE_RESOLUTION_WIDTH_central                -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_central"  );
	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap                 -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_endcap"   );
   	gr_ABSOLUTE_RESOLUTION_MEAN_central                 -> SetName("gr_ABSOLUTE_RESOLUTION_MEAN_central"  );
	gr_ABSOLUTE_RESOLUTION_MEAN_endcap                  -> SetName("gr_ABSOLUTE_RESOLUTION_MEAN_endcap"   );
   	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central      -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central"  );
	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap       -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap"   );

	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central                -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central"  );
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap                 -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap"   );
   	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central                 -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central"  );
	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap                  -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap"   );
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central      -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central"  );
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap       -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap"   );

	gr_RESPONSE_central   ->Write();
	gr_RESPONSE_endcap    ->Write();
	gr_RESOLUTION_central ->Write();
	gr_RESOLUTION_endcap  ->Write();
	gr_WIDTH_central      ->Write();
	gr_WIDTH_endcap       ->Write();
	
	gr_ABSOLUTE_RESOLUTION_WIDTH_central                ->Write();
	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap                 ->Write();
   	gr_ABSOLUTE_RESOLUTION_MEAN_central                 ->Write();
	gr_ABSOLUTE_RESOLUTION_MEAN_endcap                  ->Write();
   	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central      ->Write();
	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap       ->Write();
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central           ->Write();
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap            ->Write();
   	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central            ->Write();
	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap             ->Write();
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central ->Write();
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap  ->Write();
	
	
	for (int etabin =0; etabin< nbins_eta-1; etabin++)
  	{
		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
		{		  			
			PT_RESPONSE[etabin][ptbin] ->Write();
			ABSOLUTE_RESOLUTION[etabin][ptbin] ->Write();
			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin] ->Write();
		}
	}

	Out->ls();
	Out->Write();
 }



void make_calojet_resolution(string input_file, int numerator, Color_t color) {

	gROOT->SetBatch(); 
  	TCanvas *c1236= new TCanvas("c1236","",200,10,800,700);

	TFile *Out;
	string output_file;
	if (numerator==0) output_file = "ResolutionPlots_CaloJet_NoCorr_"+input_file;
	if (numerator==1) output_file = "ResolutionPlots_CaloJet_Corr_"+input_file;

	cout<<output_file<<endl;
	Out = new TFile(output_file.c_str() ,"RECREATE");
	TFile *F1   = new TFile(input_file.c_str() );

	// Get Tree entries
	Float_t CaloJet_Pt               ;
	Float_t CaloJet_CorrPt           ;
	Float_t CaloJet_CorrPtRhoArea    ;
	Float_t CaloJet_Eta              ;
	Float_t CaloJet_Phi              ;
	Float_t CaloJet_MatchedGenJet_Pt              ;
	Float_t CaloJet_MatchedGenJet_Eta             ;
	Float_t CaloJet_MatchedGenJet_Phi             ;


	TTree *T1    = (TTree*)  F1     ->Get("tree/CaloJetTree");
	T1->SetBranchAddress("CaloJet_Pt"                          , & CaloJet_Pt                          );
	T1->SetBranchAddress("CaloJet_CorrPt"                      , & CaloJet_CorrPt                      );
	T1->SetBranchAddress("CaloJet_CorrPtRhoArea"               , & CaloJet_CorrPtRhoArea               );
	T1->SetBranchAddress("CaloJet_Eta"                         , & CaloJet_Eta                         );
	T1->SetBranchAddress("CaloJet_Phi"                         , & CaloJet_Phi                         );
	T1->SetBranchAddress("CaloJet_MatchedGenJet_Pt"            , & CaloJet_MatchedGenJet_Pt            );
	T1->SetBranchAddress("CaloJet_MatchedGenJet_Eta"           , & CaloJet_MatchedGenJet_Eta           );
	T1->SetBranchAddress("CaloJet_MatchedGenJet_Phi"           , & CaloJet_MatchedGenJet_Phi           );


	double entries = T1->GetEntries();
	cout<<"entries = "<< entries <<endl;

	 // // Define genJet pt an eta bins
	static const int nbins_eta = 3;
	static const int nbins_pt = 16;
	double etabins[nbins_eta]           = {0.0, 1.3, 3};
	double ptbins[nbins_pt]             = {30., 40., 50., 60., 70., 80., 90., 100., 125., 150., 200., 300.}; //{10., 15., 20., 25., 
	string histName_etabins[nbins_eta]  = {"0", "1p3", "3p0"};
	string histName_ptbins[nbins_pt]    = { "030", "040", "050", "060", "070", "080", "090", "100", "125", "150", "200", "300"}; //"010", "015", "020", "025",
	string histTitle_etabins[nbins_eta] = {"0", "1.3", "3.0"};
	string histTitle_ptbins[nbins_pt]   = { "30", "40", "50", "60", "70", "80", "90", "100", "125", "150", "200", "300"}; //"10", "15", "20", "25",

	// static const int nbins_eta = 3;
	// static const int nbins_pt = 9;
	// double etabins[nbins_eta]  =  {0.0, 1.3, 3};//, 5};
	// double ptbins[nbins_pt]    =  {10., 20., 40., 60., 80., 100., 150., 200., 300. };
	// string histName_etabins[nbins_eta] = {"0", "1p3", "3p0"};//, "5p0"};
	// string histName_ptbins[nbins_pt]   = {"10", "20", "40", "60", "80", "100", "150", "200", "300"};
	// string histTitle_etabins[nbins_eta] = {"0", "1.3", "3.0"};//, "5.0"};
	// string histTitle_ptbins[nbins_pt]   = {"10", "20", "40", "60", "80", "100", "150", "200", "300" };
	



	// Arrays for TGraph x and y positions and uncertainties

	Float_t ptbins_x  [nbins_pt-1] = {};
	Float_t ptbins_ex [nbins_pt-1] = {};

	// Relataive resolution 
	Float_t response_mean_y    [nbins_eta-1][nbins_pt-1] = {};
	Float_t response_mean_ey   [nbins_eta-1][nbins_pt-1] = {};

	Float_t response_width_y   [nbins_eta-1][nbins_pt-1] = {};
	Float_t response_width_ey  [nbins_eta-1][nbins_pt-1] = {};
	
	Float_t resolution_y       [nbins_eta-1][nbins_pt-1] = {};
	Float_t resolution_ey      [nbins_eta-1][nbins_pt-1] = {};
 
	// Absolute resolution
	Float_t absolute_resolution_mean_y              [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_resolution_mean_ey             [nbins_eta-1][nbins_pt-1] = {};
 
	Float_t absolute_resolution_width_y             [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_resolution_width_ey            [nbins_eta-1][nbins_pt-1] = {};

	Float_t absolute_resolution_width_over_mean_y   [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_resolution_width_over_mean_ey  [nbins_eta-1][nbins_pt-1] = {};

	// Absoulte normalized by response
	Float_t absolute_norm_resolution_width_y             [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_norm_resolution_width_ey            [nbins_eta-1][nbins_pt-1] = {};

	Float_t absolute_norm_resolution_mean_y        [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_norm_resolution_mean_ey       [nbins_eta-1][nbins_pt-1] = {};

	Float_t absolute_norm_resolution_width_over_mean_y   [nbins_eta-1][nbins_pt-1] = {};
	Float_t absolute_norm_resolution_width_over_mean_ey  [nbins_eta-1][nbins_pt-1] = {};

	// response histogram array
	TH1F * PT_RESPONSE[nbins_eta-1][nbins_pt-1];
	TH1F * ABSOLUTE_RESOLUTION[nbins_eta-1][nbins_pt-1];
	TH1F * ABSOLUTE_RESOLUTION_NORM[nbins_eta-1][nbins_pt-1];
	// For each eta and pt bin loop over the jet tree, fill the response histogram, and then fit
	for (int etabin =0; etabin< nbins_eta-1; etabin++)
  	{
		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
		{		  			
	  		string histoname   = "RESPONSE_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
			string histotitle  = "Jet p_{T} response ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}^{reco}/p_{T}^{gen};Number of jets";
			PT_RESPONSE[etabin][ptbin]  = new TH1F(histoname.c_str() ,histotitle.c_str(),400,0,4); 
			cout<<histoname <<endl;

			string histoname2   = "ABSOLUTE_RESOLUTION_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
			string histotitle2  = "Jet p_{T} absolute resolution ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}^{reco}-p_{T}^{gen};Number of jets";
			ABSOLUTE_RESOLUTION[etabin][ptbin]  = new TH1F(histoname2.c_str() ,histotitle2.c_str(),400,-50,60); 

			string histoname3   = "ABSOLUTE_RESOLUTION_NORM_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
			string histotitle3  = "Jet p_{T} absolute response ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );(p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen};Number of jets";
			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]  = new TH1F(histoname3.c_str() ,histotitle3.c_str(),400,-2,4); 


			for (int i=0; i<entries; i++ )
		 	{  		    
		 		T1->GetEntry(i);

		 		// // all-eta PF Jet ID
		 		// if ( CaloJet_neutralHadronEnergyFraction > 0.99 ) continue;
		 		// if ( CaloJet_neutralEmEnergyFraction     > 0.99 ) continue;
		 		// if ( CaloJet_Nconst                     <= 1    ) continue;
		 		// if ( CaloJet_chargedMuEnergyFraction     > 0.8  ) continue; 

		 		// // central PF Jet ID
		 		// if ( fabs(CaloJet_Eta)<2.4 && CaloJet_chargedEmEnergyFraction      > 0.99 ) continue;
		 		// if ( fabs(CaloJet_Eta)<2.4 && CaloJet_chargedMultiplicity         <= 0    ) continue;
		 		// if ( fabs(CaloJet_Eta)<2.4 && CaloJet_chargedHadronEnergyFraction <= 0    ) continue;
		 		// if ( fabs(CaloJet_Eta)<2.4 && CaloJet_chargedEmEnergyFraction      > 0.99 ) continue;
     
     			// Fill PF Jet Response histogram for the selected eta and pt bin
		  		if ( fabs(CaloJet_MatchedGenJet_Eta)>etabins[etabin] && fabs(CaloJet_MatchedGenJet_Eta)<etabins[etabin+1] && CaloJet_MatchedGenJet_Pt > ptbins[ptbin] && CaloJet_MatchedGenJet_Pt < ptbins[ptbin+1] )
		  		{		 		
	  				//if (CaloJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(CaloJet_CorrPt/CaloJet_MatchedGenJet_Pt);
	  				if (numerator ==0 && CaloJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(CaloJet_Pt/CaloJet_MatchedGenJet_Pt);
	  				if (numerator ==1 && CaloJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(CaloJet_CorrPt/CaloJet_MatchedGenJet_Pt);
	  				if (numerator ==2 && CaloJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(CaloJet_CorrPtRhoArea/CaloJet_MatchedGenJet_Pt);
	  				if (numerator ==0 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( CaloJet_Pt-CaloJet_MatchedGenJet_Pt                     );
	  				if (numerator ==1 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( CaloJet_CorrPt-CaloJet_MatchedGenJet_Pt                 );
	  				if (numerator ==2 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( CaloJet_CorrPtRhoArea-CaloJet_MatchedGenJet_Pt          );
	  				if (numerator ==0 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( CaloJet_Pt-CaloJet_MatchedGenJet_Pt            ) / CaloJet_MatchedGenJet_Pt          );
	  				if (numerator ==1 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( CaloJet_CorrPt-CaloJet_MatchedGenJet_Pt        ) / CaloJet_MatchedGenJet_Pt          );
	  				if (numerator ==2 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( CaloJet_CorrPtRhoArea-CaloJet_MatchedGenJet_Pt ) / CaloJet_MatchedGenJet_Pt          );
	  			}
		 	}

			ptbins_x  [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
			ptbins_ex [ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;

		 	//----------------------------------------
		 	// Fit response (fit twice, first to find rough mean and sigma, second with correct range)
			TF1 *gfit = new TF1("Gaussian","gaus",0.4,3); 
			gfit->SetLineColor(color);
			gfit->SetLineWidth(2.5);

			PT_RESPONSE[etabin][ptbin]->Fit(gfit,"RQ"); // Fit histogram in range 
			double amp    = gfit->GetParameter(0);
			double eamp   = gfit->GetParError(0); 
			double mean   = gfit->GetParameter(1);
			double emean  = gfit->GetParError(1); 
			double width  = gfit->GetParameter(2);
			double ewidth = gfit->GetParError(2); 

			TF1 *gfit0 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
			gfit0->SetLineColor(color);
			gfit0->SetLineWidth(2.5);

			PT_RESPONSE[etabin][ptbin]->Fit(gfit0,"RQ"); // Fit histogram in range  R
			amp    = gfit0->GetParameter(0);
			eamp   = gfit0->GetParError(0); 
			mean   = gfit0->GetParameter(1);
			emean  = gfit0->GetParError(1); 
			width  = gfit0->GetParameter(2);
			ewidth = gfit0->GetParError(2); 

			response_mean_y [etabin][ptbin] = mean;
			response_mean_ey[etabin][ptbin] = emean;
			
			response_width_y [etabin][ptbin] = width;
			response_width_ey[etabin][ptbin] = ewidth ;
			
			resolution_y [etabin][ptbin] = width/mean;
			resolution_ey[etabin][ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			
		 	//----------------------------------------
			// Fit Absolute Resolution (fit twice, first to find rough mean and sigma, second with correct range)
			TF1 *gfit2 = new TF1("Gaussian","gaus");//,0.4,3); 
			gfit2->SetLineColor(color);
			gfit2->SetLineWidth(2.5);

			ABSOLUTE_RESOLUTION[etabin][ptbin]->Fit(gfit2,"Q"); // Fit histogram in range  R
			amp    = gfit2->GetParameter(0);
			eamp   = gfit2->GetParError(0); 
			mean   = gfit2->GetParameter(1);
			emean  = gfit2->GetParError(1); 
			width  = gfit2->GetParameter(2);
			ewidth = gfit2->GetParError(2); 

			TF1 *gfit3 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
			gfit3->SetLineColor(color);
			gfit3->SetLineWidth(2.5);

			ABSOLUTE_RESOLUTION[etabin][ptbin]->Fit(gfit3,"RQ"); // Fit histogram in range  R
			amp    = gfit3->GetParameter(0);
			eamp   = gfit3->GetParError(0); 
			mean   = gfit3->GetParameter(1);
			emean  = gfit3->GetParError(1); 
			width  = gfit3->GetParameter(2);
			ewidth = gfit3->GetParError(2); 

			absolute_resolution_mean_y [etabin][ptbin] = 1+mean;
			absolute_resolution_mean_ey[etabin][ptbin] = emean;

			absolute_resolution_width_over_mean_y [etabin][ptbin] = width/(1+mean);
			absolute_resolution_width_over_mean_ey[etabin][ptbin] = width/(1+mean)*sqrt( (emean/(1+mean))*(emean/(1+mean)) + (ewidth/width)*(ewidth/width)   );
			
			absolute_resolution_width_y [etabin][ptbin] = width;
			absolute_resolution_width_ey[etabin][ptbin] = ewidth ;
			
		 	//----------------------------------------
			// Fit Absolute Resolution (fit twice, first to find rough mean and sigma, second with correct range)

			TF1 *gfit4 = new TF1("Gaussian","gaus");//,0.4,3); 
			gfit4->SetLineColor(color);
			gfit4->SetLineWidth(2.5);

			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fit(gfit4,"Q"); // Fit histogram in range  R
			amp    = gfit4->GetParameter(0);
			eamp   = gfit4->GetParError(0); 
			mean   = gfit4->GetParameter(1);
			emean  = gfit4->GetParError(1); 
			width  = gfit4->GetParameter(2);
			ewidth = gfit4->GetParError(2); 

			TF1 *gfit5 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
			gfit5->SetLineColor(color);
			gfit5->SetLineWidth(2.5);

			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fit(gfit5,"RQ"); // Fit histogram in range  R
			amp    = gfit5->GetParameter(0);
			eamp   = gfit5->GetParError(0); 
			mean   = gfit5->GetParameter(1);
			emean  = gfit5->GetParError(1); 
			width  = gfit5->GetParameter(2);
			ewidth = gfit5->GetParError(2); 

			// central response array
			absolute_norm_resolution_mean_y [etabin][ptbin] = 1+mean;
			absolute_norm_resolution_mean_ey[etabin][ptbin] = emean;

   			// central resoabsolute_norm_lution array
			absolute_norm_resolution_width_over_mean_y [etabin][ptbin] = width/(1+mean);
			absolute_norm_resolution_width_over_mean_ey[etabin][ptbin] = width/(1+mean)*sqrt( (emean/(1+mean))*(emean/(1+mean)) + (ewidth/width)*(ewidth/width)   );
			
			// central widtabsolute_norm_h array
			absolute_norm_resolution_width_y [etabin][ptbin] = width;
			absolute_norm_resolution_width_ey[etabin][ptbin] = ewidth ;
	  	}
	}
	// Relative resolution
	TGraphErrors * gr_RESPONSE_central = new TGraphErrors(nbins_pt-1,ptbins_x,response_mean_y[0],ptbins_ex,response_mean_ey[0]);
   	TGraphErrors * gr_RESPONSE_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,response_mean_y[1],ptbins_ex,response_mean_ey[1]);
	string title_response = ";GenJet p_{T} (GeV/c);<p_{T}^{calojet}/p_{T}^{gen}>";
   	gr_RESPONSE_central ->SetTitle( title_response.c_str() );
   	gr_RESPONSE_endcap  ->SetTitle( title_response.c_str() );

   	TGraphErrors * gr_RESOLUTION_central = new TGraphErrors(nbins_pt-1,ptbins_x,resolution_y[0],ptbins_ex,resolution_ey[0]);
   	TGraphErrors * gr_RESOLUTION_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,resolution_y[1],ptbins_ex,resolution_ey[1]);
	string title_resolution = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{calojet}/p_{T}^{gen})/<p_{T}^{calojet}/p_{T}^{gen}>";
   	gr_RESOLUTION_central->SetTitle( title_resolution.c_str() );
   	gr_RESOLUTION_endcap ->SetTitle( title_resolution.c_str() );

	TGraphErrors * gr_WIDTH_central = new TGraphErrors(nbins_pt-1,ptbins_x,response_width_y[0],ptbins_ex,response_width_ey[0]);
   	TGraphErrors * gr_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,response_width_y[1],ptbins_ex,response_width_ey[1]);
	string title_width = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{calojet}/p_{T}^{gen})";
   	gr_WIDTH_central->SetTitle( title_width.c_str() );
   	gr_WIDTH_endcap ->SetTitle( title_width.c_str() );

	// Absolute resolution
	TGraphErrors * gr_ABSOLUTE_RESOLUTION_MEAN_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_mean_y[0],ptbins_ex,absolute_resolution_mean_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_RESOLUTION_MEAN_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_mean_y[1],ptbins_ex,absolute_resolution_mean_ey[1]);
	title_response = ";GenJet p_{T} (GeV/c);1+<p_{T}^{calojet}-p_{T}^{gen}>";
   	gr_ABSOLUTE_RESOLUTION_MEAN_central ->SetTitle( title_response.c_str() );
   	gr_ABSOLUTE_RESOLUTION_MEAN_endcap  ->SetTitle( title_response.c_str() );

   	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_width_over_mean_y[0],ptbins_ex,absolute_resolution_width_over_mean_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_width_over_mean_y[1],ptbins_ex,absolute_resolution_width_over_mean_ey[1]);
	title_resolution = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{calojet}-p_{T}^{gen})/(1+<p_{T}^{calojet}-p_{T}^{gen}>)";
   	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central->SetTitle( title_resolution.c_str() );
   	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap ->SetTitle( title_resolution.c_str() );

	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_width_y[0],ptbins_ex,absolute_resolution_width_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_resolution_width_y[1],ptbins_ex,absolute_resolution_width_ey[1]);
	title_width = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{calojet}-p_{T}^{gen})";
   	gr_ABSOLUTE_RESOLUTION_WIDTH_central->SetTitle( title_width.c_str() );
   	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap ->SetTitle( title_width.c_str() );





	// Make Tgraphs
	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_mean_y[0],ptbins_ex,absolute_norm_resolution_mean_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_mean_y[1],ptbins_ex,absolute_norm_resolution_mean_ey[1] );
	title_response = ";GenJet p_{T} (GeV/c);1+ < (p_{T}^{calojet}-p_{T}^{gen}) / p_{T}^{gen}>";
   	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central ->SetTitle( title_response.c_str() );
   	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap  ->SetTitle( title_response.c_str() );

   	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_width_over_mean_y[0],ptbins_ex,absolute_norm_resolution_width_over_mean_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_width_over_mean_y[1],ptbins_ex,absolute_norm_resolution_width_over_mean_ey[1]);
	title_resolution = ";GenJet p_{T} (GeV/c);#sigma( (p_{T}^{calojet}-p_{T}^{gen}) / p_{T}^{gen} )/(1+< (p_{T}^{calojet}-p_{T}^{gen})/p_{T}^{gen}>)";
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central->SetTitle( title_resolution.c_str() );
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap ->SetTitle( title_resolution.c_str() );

	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_width_y[0],ptbins_ex,absolute_norm_resolution_width_ey[0]);
   	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,ptbins_x,absolute_norm_resolution_width_y[1],ptbins_ex,absolute_norm_resolution_width_ey[1]);
	title_width = ";GenJet p_{T} (GeV/c);#sigma( (p_{T}^{calojet}-p_{T}^{gen})/p_{T}^{gen} )";
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central->SetTitle( title_width.c_str() );
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap ->SetTitle( title_width.c_str() );




 //   	// Save TGraphs to output file	
	Out->cd();
	gr_RESPONSE_central      -> SetName("gr_RESPONSE_central"  );
	gr_RESPONSE_endcap       -> SetName("gr_RESPONSE_endcap"   );
	gr_RESOLUTION_central    -> SetName("gr_RESOLUTION_central");
	gr_RESOLUTION_endcap     -> SetName("gr_RESOLUTION_endcap" );
    gr_WIDTH_central         -> SetName("gr_WIDTH_central");
	gr_WIDTH_endcap          -> SetName("gr_WIDTH_endcap" );
   	gr_ABSOLUTE_RESOLUTION_WIDTH_central                -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_central"  );
	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap                 -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_endcap"   );
   	gr_ABSOLUTE_RESOLUTION_MEAN_central                 -> SetName("gr_ABSOLUTE_RESOLUTION_MEAN_central"  );
	gr_ABSOLUTE_RESOLUTION_MEAN_endcap                  -> SetName("gr_ABSOLUTE_RESOLUTION_MEAN_endcap"   );
   	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central      -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central"  );
	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap       -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap"   );

	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central                -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central"  );
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap                 -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap"   );
   	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central                 -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central"  );
	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap                  -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap"   );
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central      -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central"  );
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap       -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap"   );

	gr_RESPONSE_central   ->Write();
	gr_RESPONSE_endcap    ->Write();
	gr_RESOLUTION_central ->Write();
	gr_RESOLUTION_endcap  ->Write();
	gr_WIDTH_central      ->Write();
	gr_WIDTH_endcap       ->Write();
	
	gr_ABSOLUTE_RESOLUTION_WIDTH_central                ->Write();
	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap                 ->Write();
   	gr_ABSOLUTE_RESOLUTION_MEAN_central                 ->Write();
	gr_ABSOLUTE_RESOLUTION_MEAN_endcap                  ->Write();
   	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central      ->Write();
	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap       ->Write();
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central           ->Write();
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap            ->Write();
   	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central            ->Write();
	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap             ->Write();
   	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central ->Write();
	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap  ->Write();
	
	
	for (int etabin =0; etabin< nbins_eta-1; etabin++)
  	{
		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
		{		  			
			PT_RESPONSE[etabin][ptbin] ->Write();
			ABSOLUTE_RESOLUTION[etabin][ptbin] ->Write();
			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin] ->Write();
		}
	}

	Out->ls();
	Out->Write();
 }

// void make_calojet_resolution(string input_file, int numerator, Color_t color) {
// 	gROOT->SetBatch(); 
//   	TCanvas *c1236= new TCanvas("c1236","",200,10,800,700);



// 	TFile *Out;

// 	string output_file = "";
// 	if (numerator==0) output_file = "ResolutionPlots_CaloJet_NoCorr_"+input_file;
// 	if (numerator==1) output_file = "ResolutionPlots_CaloJet_Corr_"+input_file;
// 	if (numerator==2) output_file = "ResolutionPlots_CaloJet_CorrPtRhoArea_"+input_file;
// 	if (numerator==3) output_file = "ResolutionPlots_CaloJet_CorrPtRhocentralArea_"+input_file;


// 	cout<<output_file<<endl;
// 	Out = new TFile(output_file.c_str() ,"RECREATE");
// 	TFile *F1   = new TFile(input_file.c_str() );

// 	// Get Tree entries
// 	Float_t CaloJet_Pt        ;
// 	Float_t CaloJet_CorrPt              ;
// 	Float_t CaloJet_CorrPtRhoArea       ;
// 	Float_t CaloJet_CorrPtRhocentralArea;

// 	Float_t CaloJet_Eta       ;
// 	Float_t CaloJet_Phi       ;
// 	Float_t CaloJet_MatchedGenJet_Pt       ;
// 	Float_t CaloJet_MatchedGenJet_Eta      ;
// 	Float_t CaloJet_MatchedGenJet_Phi      ;

// 	TTree *T1    = (TTree*)  F1     ->Get("tree/CaloJetTree");
// 	T1->SetBranchAddress("CaloJet_Pt"                        , & CaloJet_Pt                        );
// 	T1->SetBranchAddress("CaloJet_CorrPt"                    , & CaloJet_CorrPt                    );
// 	T1->SetBranchAddress("CaloJet_CorrPtRhoArea"             , & CaloJet_CorrPtRhoArea             );
// 	T1->SetBranchAddress("CaloJet_CorrPtRhocentralArea"      , & CaloJet_CorrPtRhocentralArea      );
// 	T1->SetBranchAddress("CaloJet_Eta"                       , & CaloJet_Eta                       );
// 	T1->SetBranchAddress("CaloJet_Phi"                       , & CaloJet_Phi                       );
// 	T1->SetBranchAddress("CaloJet_MatchedGenJet_Pt"          , & CaloJet_MatchedGenJet_Pt          );
// 	T1->SetBranchAddress("CaloJet_MatchedGenJet_Eta"         , & CaloJet_MatchedGenJet_Eta         );
// 	T1->SetBranchAddress("CaloJet_MatchedGenJet_Phi"         , & CaloJet_MatchedGenJet_Phi         );

// 	double entries = T1->GetEntries();
// 	cout<<"entries = "<< entries <<endl;


// 	// Define genJet pt an eta bins
// 	static const int nbins_eta = 4;
// 	static const int nbins_pt = 16;
// 	double etabins[nbins_eta]           = {0.0, 1.3, 3, 5};
// 	double ptbins[nbins_pt]             = {10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 90., 100., 125., 150., 200., 300.};
// 	string histName_etabins[nbins_eta]  = {"0", "1p3", "3p0", "5p0"};
// 	string histName_ptbins[nbins_pt]    = {"010", "015", "020", "025", "030", "040", "050", "060", "070", "080", "090", "100", "125", "150", "200", "300"};
// 	string histTitle_etabins[nbins_eta] = {"0", "1.3", "3.0", "5.0"};
// 	string histTitle_ptbins[nbins_pt]   = {"10", "15", "20", "25", "30", "40", "50", "60", "70", "80", "90", "100", "125", "150", "200", "300"};



// 	// static const int nbins_eta = 3;
// 	// static const int nbins_pt = 9;
// 	// double etabins[nbins_eta]  =  {0.0, 1.3, 3};//, 5};
// 	// double ptbins[nbins_pt]    =  {10., 20., 40., 60., 80., 100., 150., 200., 300. };
// 	// string histName_etabins[nbins_eta] = {"0", "1p3", "3p0"};//, "5p0"};
// 	// string histName_ptbins[nbins_pt]   = {"10", "20", "40", "60", "80", "100", "150", "200", "300"};
// 	// string histTitle_etabins[nbins_eta] = {"0", "1.3", "3.0"};//, "5.0"};
// 	// string histTitle_ptbins[nbins_pt]   = {"10", "20", "40", "60", "80", "100", "150", "200", "300" };
	

// 	// Arrays for TGraph x and y positions and uncertainties
// 	Float_t resolution_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t resolution_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t resolution_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t resolution_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t response_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t response_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t response_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t response_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t width_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t width_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t width_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t width_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_width_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_width_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_width_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_mean_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_mean_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_mean_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_mean_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_width_over_mean_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_width_over_mean_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_resolution_width_over_mean_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_resolution_width_over_mean_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};






// 	Float_t absolute_norm_resolution_width_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_norm_resolution_width_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_norm_resolution_width_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_norm_resolution_mean_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_norm_resolution_mean_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_norm_resolution_mean_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_mean_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_central_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_central_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_central_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_central_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_endcap_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_endcap_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_endcap_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_endcap_ey[nbins_eta-1][nbins_pt-1] = {};

// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_forward_x [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_forward_y [nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_forward_ex[nbins_eta-1][nbins_pt-1] = {};
// 	Float_t absolute_norm_resolution_width_over_mean_vs_pt_forward_ey[nbins_eta-1][nbins_pt-1] = {};


// 	// response histogram array
// 	TH1F * PT_RESPONSE[nbins_eta-1][nbins_eta-1][nbins_pt-1];
// 	TH1F * ABSOLUTE_RESOLUTION[nbins_eta-1][nbins_eta-1][nbins_pt-1];
// 	TH1F * ABSOLUTE_RESOLUTION_NORM[nbins_eta-1][nbins_eta-1][nbins_pt-1];

// 	// For each eta and pt bin loop over the jet tree, fill the response histogram, and then fit
// 	for (int etabin =0; etabin< nbins_eta-1; etabin++)
//   	{
// 		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
// 		{		  			
// 	  		string histoname   = "RESPONSE_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
// 			string histotitle  = "Jet p_{T} response ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}^{reco}/p_{T}^{gen};Number of jets";
// 			PT_RESPONSE[etabin][ptbin]  = new TH1F(histoname.c_str() ,histotitle.c_str(),400,0,4); 
			
// 			string histoname2   = "ABSOLUTE_RESOLUTION_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
// 			string histotitle2  = "Jet p_{T} absolute resolution ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );p_{T}^{reco}-p_{T}^{gen};Number of jets";
// 			ABSOLUTE_RESOLUTION[etabin][ptbin]  = new TH1F(histoname2.c_str() ,histotitle2.c_str(),500,-100,100); 

// 			string histoname3   = "ABSOLUTE_RESOLUTION_NORM_ETA_"+histName_etabins[etabin]+"_"+histName_etabins[etabin+1]+"_PT_"+histName_ptbins[ptbin]+"_"+histName_ptbins[ptbin+1];
// 			string histotitle3  = "Jet p_{T} absolute response ( "+histTitle_etabins[etabin]+" < | #eta | < "+histTitle_etabins[etabin+1]+" ,  "+histTitle_ptbins[ptbin]+" < p_{T} < "+histTitle_ptbins[ptbin+1]+" GeV/c );(p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen};Number of jets";
// 			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]  = new TH1F(histoname3.c_str() ,histotitle3.c_str(),500,-4,4);//100,-40,60); 

// 			cout<<histoname <<endl;
// 			int countn  =0;
// 			int countp  =0;
// 			for (int i=0; i<entries; i++ )
// 		 	{  		    
// 		 		T1->GetEntry(i);

		 		
//      			// Fill PF Jet Response histogram for the selected eta and pt bin
// 		  		if ( fabs(CaloJet_MatchedGenJet_Eta)>etabins[etabin] && fabs(CaloJet_MatchedGenJet_Eta)<etabins[etabin+1] && CaloJet_MatchedGenJet_Pt > ptbins[ptbin] && CaloJet_MatchedGenJet_Pt < ptbins[ptbin+1] )
// 		  		{		 		
// 	  				//if (GenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(PFJet_CorrPt/GenJet_Pt);
// 	  				if (numerator ==0 && CaloJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(CaloJet_Pt/CaloJet_MatchedGenJet_Pt);
// 	  				if (numerator ==1 && CaloJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(CaloJet_CorrPt/CaloJet_MatchedGenJet_Pt);
// 	  				if (numerator ==2 && CaloJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(CaloJet_CorrPtRhoArea/CaloJet_MatchedGenJet_Pt);
// 	  				if (numerator ==3 && CaloJet_MatchedGenJet_Pt!=0) PT_RESPONSE[etabin][ptbin]->Fill(CaloJet_CorrPtRhocentralArea/CaloJet_MatchedGenJet_Pt);
// 	  				if (numerator ==0 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( CaloJet_Pt-CaloJet_MatchedGenJet_Pt                     );
// 	  				if (numerator ==1 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( CaloJet_CorrPt-CaloJet_MatchedGenJet_Pt                 );
// 	  				if (numerator ==2 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( CaloJet_CorrPtRhoArea-CaloJet_MatchedGenJet_Pt          );
// 	  				if (numerator ==3 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION[etabin][ptbin]->Fill( CaloJet_CorrPtRhocentralArea-CaloJet_MatchedGenJet_Pt          );
// 	  				if (numerator ==0 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( CaloJet_Pt-CaloJet_MatchedGenJet_Pt            ) / CaloJet_MatchedGenJet_Pt          );
// 	  				if (numerator ==1 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( CaloJet_CorrPt-CaloJet_MatchedGenJet_Pt        ) / CaloJet_MatchedGenJet_Pt          );
// 	  				if (numerator ==2 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( CaloJet_CorrPtRhoArea-CaloJet_MatchedGenJet_Pt ) / CaloJet_MatchedGenJet_Pt          );
// 	  				if (numerator ==3 && CaloJet_MatchedGenJet_Pt!=0) ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fill( ( CaloJet_CorrPtRhocentralArea-CaloJet_MatchedGenJet_Pt ) / CaloJet_MatchedGenJet_Pt          );


// 	  			}
// 		 	}
// 		 	cout<<"countn "<<countn<<endl;
// 		 	cout<<"countp "<<countp<<endl;
// 		 	// Fit
// 			TF1 *gfit = new TF1("Gaussian","gaus",0,3); 
// 			gfit->SetLineColor(color);
// 			gfit->SetLineWidth(2.5);

// 			PT_RESPONSE[etabin][ptbin]->Fit(gfit,"RQ"); // Fit histogram in range 
// 			double amp    = gfit->GetParameter(0);
// 			double eamp   = gfit->GetParError(0); 
// 			double mean   = gfit->GetParameter(1);
// 			double emean  = gfit->GetParError(1); 
// 			double width  = gfit->GetParameter(2);
// 			double ewidth = gfit->GetParError(1); 


// 			cout<<"mean "<<mean<<" width "<<width<<endl;
// 			cout<<"mean-1.75*width "<<mean-1.75*width<<" mean+1.75*width "<<mean+1.75*width<<endl;
// 			TF1 *gfit0 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
// 			gfit0->SetLineColor(color);
// 			gfit0->SetLineWidth(2.5);

// 			PT_RESPONSE[etabin][ptbin]->Fit(gfit0,"RQ"); // Fit histogram in range  R
// 			amp    = gfit0->GetParameter(0);
// 			eamp   = gfit0->GetParError(0); 
// 			mean   = gfit0->GetParameter(1);
// 			emean  = gfit0->GetParError(1); 
// 			width  = gfit0->GetParameter(2);
// 			ewidth = gfit0->GetParError(1); 


// 			// central response array
// 			if (etabin==0) response_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) response_vs_pt_central_y [ptbin] = mean;
// 			if (etabin==0) response_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) response_vs_pt_central_ey[ptbin] = emean;

// 			// endcap response array
// 			if (etabin==1) response_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) response_vs_pt_endcap_y [ptbin] = mean;
// 			if (etabin==1) response_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) response_vs_pt_endcap_ey[ptbin] = emean;

// 			// forward response array
// 			if (etabin==2) response_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) response_vs_pt_forward_y [ptbin] = mean;
// 			if (etabin==2) response_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) response_vs_pt_forward_ey[ptbin] = emean;

//    			// central resolution array
// 			if (etabin==0) resolution_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) resolution_vs_pt_central_y [ptbin] = width/mean;
// 			if (etabin==0) resolution_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) resolution_vs_pt_central_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			
// 			// endcap resolution array
// 			if (etabin==1) resolution_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) resolution_vs_pt_endcap_y [ptbin] = width/mean;
// 			if (etabin==1) resolution_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) resolution_vs_pt_endcap_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			
// 			// forward resolution array
// 			if (etabin==2) resolution_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) resolution_vs_pt_forward_y [ptbin] = width/mean;
// 			if (etabin==2) resolution_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) resolution_vs_pt_forward_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );


// 			// central width array
// 			if (etabin==0) width_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) width_vs_pt_central_y [ptbin] = width;
// 			if (etabin==0) width_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) width_vs_pt_central_ey[ptbin] = ewidth ;
			
// 			// endcap width array
// 			if (etabin==1) width_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) width_vs_pt_endcap_y [ptbin] = width;
// 			if (etabin==1) width_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) width_vs_pt_endcap_ey[ptbin] = ewidth;
			
// 			// forward width array
// 			if (etabin==2) width_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) width_vs_pt_forward_y [ptbin] = width;
// 			if (etabin==2) width_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) width_vs_pt_forward_ey[ptbin] = ewidth;






// 			TF1 *gfit2 = new TF1("Gaussian","gaus");//,0.4,3); 
// 			gfit2->SetLineColor(color);
// 			gfit2->SetLineWidth(2.5);

// 			ABSOLUTE_RESOLUTION[etabin][ptbin]->Fit(gfit2,"Q"); // Fit histogram in range  R
// 			amp    = gfit2->GetParameter(0);
// 			eamp   = gfit2->GetParError(0); 
// 			mean   = gfit2->GetParameter(1);
// 			emean  = gfit2->GetParError(1); 
// 			width  = gfit2->GetParameter(2);
// 			ewidth = gfit2->GetParError(1); 

// 			TF1 *gfit3 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
// 			gfit3->SetLineColor(color);
// 			gfit3->SetLineWidth(2.5);

// 			ABSOLUTE_RESOLUTION[etabin][ptbin]->Fit(gfit3,"RQ"); // Fit histogram in range  R
// 			amp    = gfit3->GetParameter(0);
// 			eamp   = gfit3->GetParError(0); 
// 			mean   = gfit3->GetParameter(1);
// 			emean  = gfit3->GetParError(1); 
// 			width  = gfit3->GetParameter(2);
// 			ewidth = gfit3->GetParError(1); 

// 			// central response array
// 			if (etabin==0) absolute_resolution_mean_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) absolute_resolution_mean_vs_pt_central_y [ptbin] = 1+mean;
// 			if (etabin==0) absolute_resolution_mean_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) absolute_resolution_mean_vs_pt_central_ey[ptbin] = emean;

// 			// endcap respoabsolute_nse array
// 			if (etabin==1) absolute_resolution_mean_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) absolute_resolution_mean_vs_pt_endcap_y [ptbin] = 1+mean;
// 			if (etabin==1) absolute_resolution_mean_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) absolute_resolution_mean_vs_pt_endcap_ey[ptbin] = emean;

// 			// forward respabsolute_onse array
// 			if (etabin==2) absolute_resolution_mean_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) absolute_resolution_mean_vs_pt_forward_y [ptbin] = 1+mean;
// 			if (etabin==2) absolute_resolution_mean_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) absolute_resolution_mean_vs_pt_forward_ey[ptbin] = emean;

//    			// central resoabsolute_lution array
// 			if (etabin==0) absolute_resolution_width_over_mean_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) absolute_resolution_width_over_mean_vs_pt_central_y [ptbin] = width/(1+mean);
// 			if (etabin==0) absolute_resolution_width_over_mean_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) absolute_resolution_width_over_mean_vs_pt_central_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			
// 			// endcap resolabsolute_ution arra_width_over_meany
// 			if (etabin==1) absolute_resolution_width_over_mean_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) absolute_resolution_width_over_mean_vs_pt_endcap_y [ptbin] = width/(1+mean);
// 			if (etabin==1) absolute_resolution_width_over_mean_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) absolute_resolution_width_over_mean_vs_pt_endcap_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			
// 			// forward resoabsolute_lution arr_width_over_meanay
// 			if (etabin==2) absolute_resolution_width_over_mean_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) absolute_resolution_width_over_mean_vs_pt_forward_y [ptbin] = width/(1+mean);
// 			if (etabin==2) absolute_resolution_width_over_mean_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) absolute_resolution_width_over_mean_vs_pt_forward_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );


// 			// central widtabsolute_h array
// 			if (etabin==0) absolute_resolution_width_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) absolute_resolution_width_vs_pt_central_y [ptbin] = width;
// 			if (etabin==0) absolute_resolution_width_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) absolute_resolution_width_vs_pt_central_ey[ptbin] = ewidth ;
			
// 			// endcap widthabsolute_resolution_ array
// 			if (etabin==1) absolute_resolution_width_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) absolute_resolution_width_vs_pt_endcap_y [ptbin] = width;
// 			if (etabin==1) absolute_resolution_width_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) absolute_resolution_width_vs_pt_endcap_ey[ptbin] = ewidth;
			
// 			// forward widtabsolute_resolution_h array
// 			if (etabin==2) absolute_resolution_width_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) absolute_resolution_width_vs_pt_forward_y [ptbin] = width;
// 			if (etabin==2) absolute_resolution_width_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) absolute_resolution_width_vs_pt_forward_ey[ptbin] = ewidth;




//  			//PT_RESPONSE[etabin][ptbin]->Write();
// 		 	// string savename1 = "plots/"+histoname1+".pdf";
// 		 	// c1236->SaveAs(savename1.c_str());


// 			TF1 *gfit4 = new TF1("Gaussian","gaus");//,0.4,3); 
// 			gfit4->SetLineColor(color);
// 			gfit4->SetLineWidth(2.5);

// 			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fit(gfit4,"Q"); // Fit histogram in range  R
// 			amp    = gfit4->GetParameter(0);
// 			eamp   = gfit4->GetParError(0); 
// 			mean   = gfit4->GetParameter(1);
// 			emean  = gfit4->GetParError(1); 
// 			width  = gfit4->GetParameter(2);
// 			ewidth = gfit4->GetParError(1); 

		
// 			TF1 *gfit5 = new TF1("Gaussian","gaus",mean-1.75*width,mean+1.75*width); 
// 			gfit5->SetLineColor(color);
// 			gfit5->SetLineWidth(2.5);

// 			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin]->Fit(gfit5,"RQ"); // Fit histogram in range  R
// 			amp    = gfit5->GetParameter(0);
// 			eamp   = gfit5->GetParError(0); 
// 			mean   = gfit5->GetParameter(1);
// 			emean  = gfit5->GetParError(1); 
// 			width  = gfit5->GetParameter(2);
// 			ewidth = gfit5->GetParError(1); 


// 			// central response array
// 			if (etabin==0) absolute_norm_resolution_mean_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) absolute_norm_resolution_mean_vs_pt_central_y [ptbin] = 1+mean;
// 			if (etabin==0) absolute_norm_resolution_mean_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) absolute_norm_resolution_mean_vs_pt_central_ey[ptbin] = emean;

// 			// endcap respoabsolute_norm_nse array
// 			if (etabin==1) absolute_norm_resolution_mean_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) absolute_norm_resolution_mean_vs_pt_endcap_y [ptbin] = 1+mean;
// 			if (etabin==1) absolute_norm_resolution_mean_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) absolute_norm_resolution_mean_vs_pt_endcap_ey[ptbin] = emean;

// 			// forward respabsolute_norm_onse array
// 			if (etabin==2) absolute_norm_resolution_mean_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) absolute_norm_resolution_mean_vs_pt_forward_y [ptbin] = 1+mean;
// 			if (etabin==2) absolute_norm_resolution_mean_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) absolute_norm_resolution_mean_vs_pt_forward_ey[ptbin] = emean;

//    			// central resoabsolute_norm_lution array
// 			if (etabin==0) absolute_norm_resolution_width_over_mean_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) absolute_norm_resolution_width_over_mean_vs_pt_central_y [ptbin] = width/(1+mean);
// 			if (etabin==0) absolute_norm_resolution_width_over_mean_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) absolute_norm_resolution_width_over_mean_vs_pt_central_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			
// 			// endcap resolabsolute_norm_ution arra_width_over_meany
// 			if (etabin==1) absolute_norm_resolution_width_over_mean_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) absolute_norm_resolution_width_over_mean_vs_pt_endcap_y [ptbin] = width/(1+mean);
// 			if (etabin==1) absolute_norm_resolution_width_over_mean_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) absolute_norm_resolution_width_over_mean_vs_pt_endcap_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			
// 			// forward resoabsolute_norm_lution arr_width_over_meanay
// 			if (etabin==2) absolute_norm_resolution_width_over_mean_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) absolute_norm_resolution_width_over_mean_vs_pt_forward_y [ptbin] = width/(1+mean);
// 			if (etabin==2) absolute_norm_resolution_width_over_mean_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) absolute_norm_resolution_width_over_mean_vs_pt_forward_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );


// 			// central widtabsolute_norm_h array
// 			if (etabin==0) absolute_norm_resolution_width_vs_pt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==0) absolute_norm_resolution_width_vs_pt_central_y [ptbin] = width;
// 			if (etabin==0) absolute_norm_resolution_width_vs_pt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==0) absolute_norm_resolution_width_vs_pt_central_ey[ptbin] = ewidth ;
			
// 			// endcap widthabsolute_norm_resolution_ array
// 			if (etabin==1) absolute_norm_resolution_width_vs_pt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==1) absolute_norm_resolution_width_vs_pt_endcap_y [ptbin] = width;
// 			if (etabin==1) absolute_norm_resolution_width_vs_pt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==1) absolute_norm_resolution_width_vs_pt_endcap_ey[ptbin] = ewidth;
			
// 			// forward widtabsolute_norm_resolution_h array
// 			if (etabin==2) absolute_norm_resolution_width_vs_pt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
// 			if (etabin==2) absolute_norm_resolution_width_vs_pt_forward_y [ptbin] = width;
// 			if (etabin==2) absolute_norm_resolution_width_vs_pt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
// 			if (etabin==2) absolute_norm_resolution_width_vs_pt_forward_ey[ptbin] = ewidth;



// 	  	}
// 	}

// 	// Make Tgraphs
// 	TGraphErrors * gr_RESPONSE_central = new TGraphErrors(nbins_pt-1,response_vs_pt_central_x,response_vs_pt_central_y,response_vs_pt_central_ex,response_vs_pt_central_ey);
//    	TGraphErrors * gr_RESPONSE_endcap  = new TGraphErrors(nbins_pt-1,response_vs_pt_endcap_x ,response_vs_pt_endcap_y ,response_vs_pt_endcap_ex ,response_vs_pt_endcap_ey );
// 	TGraphErrors * gr_RESPONSE_forward = new TGraphErrors(nbins_pt-1,response_vs_pt_forward_x,response_vs_pt_forward_y,response_vs_pt_forward_ex,response_vs_pt_forward_ey);

// 	string title_response = ";GenJet p_{T} (GeV/c);<p_{T}^{calo}/p_{T}^{gen}>";
//    	gr_RESPONSE_central ->SetTitle( title_response.c_str() );
//    	gr_RESPONSE_endcap  ->SetTitle( title_response.c_str() );
//    	gr_RESPONSE_forward ->SetTitle( title_response.c_str() );

//    	TGraphErrors * gr_RESOLUTION_central = new TGraphErrors(nbins_pt-1,resolution_vs_pt_central_x,resolution_vs_pt_central_y,resolution_vs_pt_central_ex,resolution_vs_pt_central_ey);
//    	TGraphErrors * gr_RESOLUTION_endcap  = new TGraphErrors(nbins_pt-1,resolution_vs_pt_endcap_x ,resolution_vs_pt_endcap_y ,resolution_vs_pt_endcap_ex ,resolution_vs_pt_endcap_ey );
// 	TGraphErrors * gr_RESOLUTION_forward = new TGraphErrors(nbins_pt-1,resolution_vs_pt_forward_x,resolution_vs_pt_forward_y,resolution_vs_pt_forward_ex,resolution_vs_pt_forward_ey);

// 	string title_resolution = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{calo}/p_{T}^{gen})/<p_{T}^{calo}/p_{T}^{gen}>";
//    	gr_RESOLUTION_central->SetTitle( title_resolution.c_str() );
//    	gr_RESOLUTION_endcap ->SetTitle( title_resolution.c_str() );
//    	gr_RESOLUTION_forward->SetTitle( title_resolution.c_str() );

// 	TGraphErrors * gr_WIDTH_central = new TGraphErrors(nbins_pt-1,width_vs_pt_central_x,width_vs_pt_central_y,width_vs_pt_central_ex,width_vs_pt_central_ey);
//    	TGraphErrors * gr_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,width_vs_pt_endcap_x ,width_vs_pt_endcap_y ,width_vs_pt_endcap_ex ,width_vs_pt_endcap_ey );
// 	TGraphErrors * gr_WIDTH_forward = new TGraphErrors(nbins_pt-1,width_vs_pt_forward_x,width_vs_pt_forward_y,width_vs_pt_forward_ex,width_vs_pt_forward_ey);

// 	string title_width = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{calo}/p_{T}^{gen})";
//    	gr_WIDTH_central->SetTitle( title_width.c_str() );
//    	gr_WIDTH_endcap ->SetTitle( title_width.c_str() );
//    	gr_WIDTH_forward->SetTitle( title_width.c_str() );



// 	// Make Tgraphs
// 	TGraphErrors * gr_ABSOLUTE_RESOLUTION_MEAN_central = new TGraphErrors(nbins_pt-1,absolute_resolution_mean_vs_pt_central_x,absolute_resolution_mean_vs_pt_central_y,absolute_resolution_mean_vs_pt_central_ex,absolute_resolution_mean_vs_pt_central_ey);
//    	TGraphErrors * gr_ABSOLUTE_RESOLUTION_MEAN_endcap  = new TGraphErrors(nbins_pt-1,absolute_resolution_mean_vs_pt_endcap_x ,absolute_resolution_mean_vs_pt_endcap_y ,absolute_resolution_mean_vs_pt_endcap_ex ,absolute_resolution_mean_vs_pt_endcap_ey );
// 	TGraphErrors * gr_ABSOLUTE_RESOLUTION_MEAN_forward = new TGraphErrors(nbins_pt-1,absolute_resolution_mean_vs_pt_forward_x,absolute_resolution_mean_vs_pt_forward_y,absolute_resolution_mean_vs_pt_forward_ex,absolute_resolution_mean_vs_pt_forward_ey);

// 	title_response = ";GenJet p_{T} (GeV/c);1+<p_{T}^{calo}-p_{T}^{gen}>";
//    	gr_ABSOLUTE_RESOLUTION_MEAN_central ->SetTitle( title_response.c_str() );
//    	gr_ABSOLUTE_RESOLUTION_MEAN_endcap  ->SetTitle( title_response.c_str() );
//    	gr_ABSOLUTE_RESOLUTION_MEAN_forward ->SetTitle( title_response.c_str() );

//    	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central = new TGraphErrors(nbins_pt-1,absolute_resolution_width_over_mean_vs_pt_central_x,absolute_resolution_width_over_mean_vs_pt_central_y,absolute_resolution_width_over_mean_vs_pt_central_ex,absolute_resolution_width_over_mean_vs_pt_central_ey);
//    	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap  = new TGraphErrors(nbins_pt-1,absolute_resolution_width_over_mean_vs_pt_endcap_x ,absolute_resolution_width_over_mean_vs_pt_endcap_y ,absolute_resolution_width_over_mean_vs_pt_endcap_ex ,absolute_resolution_width_over_mean_vs_pt_endcap_ey );
// 	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_forward = new TGraphErrors(nbins_pt-1,absolute_resolution_width_over_mean_vs_pt_forward_x,absolute_resolution_width_over_mean_vs_pt_forward_y,absolute_resolution_width_over_mean_vs_pt_forward_ex,absolute_resolution_width_over_mean_vs_pt_forward_ey);

// 	title_resolution = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{calo}-p_{T}^{gen})/(1+<p_{T}^{calo}-p_{T}^{gen}>)";
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central->SetTitle( title_resolution.c_str() );
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap ->SetTitle( title_resolution.c_str() );
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_forward->SetTitle( title_resolution.c_str() );

// 	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_central = new TGraphErrors(nbins_pt-1,absolute_resolution_width_vs_pt_central_x,absolute_resolution_width_vs_pt_central_y,absolute_resolution_width_vs_pt_central_ex,absolute_resolution_width_vs_pt_central_ey);
//    	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,absolute_resolution_width_vs_pt_endcap_x ,absolute_resolution_width_vs_pt_endcap_y ,absolute_resolution_width_vs_pt_endcap_ex ,absolute_resolution_width_vs_pt_endcap_ey );
// 	TGraphErrors * gr_ABSOLUTE_RESOLUTION_WIDTH_forward = new TGraphErrors(nbins_pt-1,absolute_resolution_width_vs_pt_forward_x,absolute_resolution_width_vs_pt_forward_y,absolute_resolution_width_vs_pt_forward_ex,absolute_resolution_width_vs_pt_forward_ey);

// 	title_width = ";GenJet p_{T} (GeV/c);#sigma(p_{T}^{calo}-p_{T}^{gen})";
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_central->SetTitle( title_width.c_str() );
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap ->SetTitle( title_width.c_str() );
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_forward->SetTitle( title_width.c_str() );





// 	// Make Tgraphs
// 	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_mean_vs_pt_central_x,absolute_norm_resolution_mean_vs_pt_central_y,absolute_norm_resolution_mean_vs_pt_central_ex,absolute_norm_resolution_mean_vs_pt_central_ey);
//    	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap  = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_mean_vs_pt_endcap_x ,absolute_norm_resolution_mean_vs_pt_endcap_y ,absolute_norm_resolution_mean_vs_pt_endcap_ex ,absolute_norm_resolution_mean_vs_pt_endcap_ey );
// 	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_MEAN_forward = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_mean_vs_pt_forward_x,absolute_norm_resolution_mean_vs_pt_forward_y,absolute_norm_resolution_mean_vs_pt_forward_ex,absolute_norm_resolution_mean_vs_pt_forward_ey);

// 	title_response = ";GenJet p_{T} (GeV/c);1 + <(p_{T}^{calo}-p_{T}^{gen}) / p_{T}^{gen}>";
//    	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central ->SetTitle( title_response.c_str() );
//    	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap  ->SetTitle( title_response.c_str() );
//    	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_forward ->SetTitle( title_response.c_str() );

//    	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_width_over_mean_vs_pt_central_x,absolute_norm_resolution_width_over_mean_vs_pt_central_y,absolute_norm_resolution_width_over_mean_vs_pt_central_ex,absolute_norm_resolution_width_over_mean_vs_pt_central_ey);
//    	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap  = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_width_over_mean_vs_pt_endcap_x ,absolute_norm_resolution_width_over_mean_vs_pt_endcap_y ,absolute_norm_resolution_width_over_mean_vs_pt_endcap_ex ,absolute_norm_resolution_width_over_mean_vs_pt_endcap_ey );
// 	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_forward = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_width_over_mean_vs_pt_forward_x,absolute_norm_resolution_width_over_mean_vs_pt_forward_y,absolute_norm_resolution_width_over_mean_vs_pt_forward_ex,absolute_norm_resolution_width_over_mean_vs_pt_forward_ey);

// 	title_resolution = ";GenJet p_{T} (GeV/c);#sigma( (p_{T}^{calo}-p_{T}^{gen}) / p_{T}^{gen} )/(1+ < (p_{T}^{calo}-p_{T}^{gen})/p_{T}^{gen}>)";
//    	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central->SetTitle( title_resolution.c_str() );
//    	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap ->SetTitle( title_resolution.c_str() );
//    	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_forward->SetTitle( title_resolution.c_str() );

// 	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_width_vs_pt_central_x,absolute_norm_resolution_width_vs_pt_central_y,absolute_norm_resolution_width_vs_pt_central_ex,absolute_norm_resolution_width_vs_pt_central_ey);
//    	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap  = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_width_vs_pt_endcap_x ,absolute_norm_resolution_width_vs_pt_endcap_y ,absolute_norm_resolution_width_vs_pt_endcap_ex ,absolute_norm_resolution_width_vs_pt_endcap_ey );
// 	TGraphErrors * gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_forward = new TGraphErrors(nbins_pt-1,absolute_norm_resolution_width_vs_pt_forward_x,absolute_norm_resolution_width_vs_pt_forward_y,absolute_norm_resolution_width_vs_pt_forward_ex,absolute_norm_resolution_width_vs_pt_forward_ey);

// 	title_width = ";GenJet p_{T} (GeV/c);#sigma( (p_{T}^{calo}-p_{T}^{gen})/p_{T}^{gen} )";
//    	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central->SetTitle( title_width.c_str() );
//    	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap ->SetTitle( title_width.c_str() );
//    	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_forward->SetTitle( title_width.c_str() );

//    	// Save TGraphs to output file	
// 	Out->cd();
// 	gr_RESPONSE_central      -> SetName("gr_RESPONSE_central"  );
// 	gr_RESPONSE_endcap       -> SetName("gr_RESPONSE_endcap"   );
// 	gr_RESPONSE_forward      -> SetName("gr_RESPONSE_forward"  );
// 	gr_RESOLUTION_central    -> SetName("gr_RESOLUTION_central");
// 	gr_RESOLUTION_endcap     -> SetName("gr_RESOLUTION_endcap" );
// 	gr_RESOLUTION_forward    -> SetName("gr_RESOLUTION_forward");
//     gr_WIDTH_central         -> SetName("gr_WIDTH_central");
// 	gr_WIDTH_endcap          -> SetName("gr_WIDTH_endcap" );
// 	gr_WIDTH_forward         -> SetName("gr_WIDTH_forward");
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_central                -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_central"  );
// 	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap                 -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_endcap"   );
// 	gr_ABSOLUTE_RESOLUTION_WIDTH_forward                -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_forward"  );
//    	gr_ABSOLUTE_RESOLUTION_MEAN_central                 -> SetName("gr_ABSOLUTE_RESOLUTION_MEAN_central"  );
// 	gr_ABSOLUTE_RESOLUTION_MEAN_endcap                  -> SetName("gr_ABSOLUTE_RESOLUTION_MEAN_endcap"   );
// 	gr_ABSOLUTE_RESOLUTION_MEAN_forward                 -> SetName("gr_ABSOLUTE_RESOLUTION_MEAN_forward"  );
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central      -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central"  );
// 	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap       -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap"   );
// 	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_forward      -> SetName("gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_forward"  );

// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central                -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central"  );
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap                 -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap"   );
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_forward                -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_forward"  );
//    	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central                 -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central"  );
// 	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap                  -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap"   );
// 	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_forward                 -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_MEAN_forward"  );
//    	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central      -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central"  );
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap       -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap"   );
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_forward      -> SetName("gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_forward"  );

//    	gr_RESPONSE_central   ->Write();
// 	gr_RESPONSE_endcap    ->Write();
// 	gr_RESPONSE_forward   ->Write();
// 	gr_RESOLUTION_central ->Write();
// 	gr_RESOLUTION_endcap  ->Write();
// 	gr_RESOLUTION_forward ->Write();
// 	gr_WIDTH_central      ->Write();
// 	gr_WIDTH_endcap       ->Write();
// 	gr_WIDTH_forward      ->Write();

// 	gr_ABSOLUTE_RESOLUTION_WIDTH_central                ->Write();
// 	gr_ABSOLUTE_RESOLUTION_WIDTH_endcap                 ->Write();
// 	gr_ABSOLUTE_RESOLUTION_WIDTH_forward                ->Write();
//    	gr_ABSOLUTE_RESOLUTION_MEAN_central                 ->Write();
// 	gr_ABSOLUTE_RESOLUTION_MEAN_endcap                  ->Write();
// 	gr_ABSOLUTE_RESOLUTION_MEAN_forward                 ->Write();
//    	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_central      ->Write();
// 	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_endcap       ->Write();
// 	gr_ABSOLUTE_RESOLUTION_WIDTH_OVER_MEAN_forward      ->Write();
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_central           ->Write();
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_endcap            ->Write();
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_forward           ->Write();
//    	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_central            ->Write();
// 	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_endcap             ->Write();
// 	gr_ABSOLUTE_NORM_RESOLUTION_MEAN_forward            ->Write();
//    	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_central ->Write();
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_endcap  ->Write();
// 	gr_ABSOLUTE_NORM_RESOLUTION_WIDTH_OVER_MEAN_forward ->Write();
	

// 	for (int etabin =0; etabin< nbins_eta-1; etabin++)
//   	{
// 		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
// 		{		  			
// 			PT_RESPONSE[etabin][ptbin] ->Write();
// 			ABSOLUTE_RESOLUTION[etabin][ptbin] ->Write();
// 			ABSOLUTE_RESOLUTION_NORM[etabin][ptbin] ->Write();
// 		}
// 	}

// 	Out->ls();
// 	Out->Write();
//  }
