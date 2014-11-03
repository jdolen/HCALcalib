{
	TCanvas *c1236= new TCanvas("c1236","",200,10,1100,700);
  
  	gStyle->SetOptFit(1);

	TFile *F1   = new TFile("OutTreeQCDflat.root");
	Int_t Run         ;
	Int_t Lumi        ;
	Int_t Event       ;
	Int_t Nvtx        ;
	Int_t NvtxGood    ;
	Float_t Mass      ;
	Float_t Energy    ;
	Float_t Pt        ;
	Float_t Eta       ;
	Float_t Phi       ;
	Float_t GenMass   ;
	Float_t GenEnergy ;
	Float_t GenPt     ;
	Float_t GenEta    ;
	Float_t GenPhi    ;
	Float_t PFJet_Nconst;
	Float_t PFJet_chargedEmEnergy              ;
	Float_t PFJet_chargedEmEnergyFraction      ;
	Float_t PFJet_chargedHadronEnergy          ;
	Float_t PFJet_chargedHadronEnergyFraction  ;
	Float_t PFJet_chargedHadronMultiplicity    ;
	Float_t PFJet_chargedMuEnergy              ;
	Float_t PFJet_chargedMuEnergyFraction      ;
	Float_t PFJet_chargedMultiplicity          ;
	Float_t PFJet_electronEnergy               ;
	Float_t PFJet_electronEnergyFraction       ;
	Float_t PFJet_electronMultiplicity         ;
	Float_t PFJet_HFEMEnergy                   ;
	Float_t PFJet_HFEMEnergyFraction           ;
	Float_t PFJet_HFEMMultiplicity             ;
	Float_t PFJet_HFHadronEnergy               ;
	Float_t PFJet_HFHadronEnergyFraction       ;
	Float_t PFJet_HFHadronMultiplicity         ;
	Float_t PFJet_muonEnergy                   ;
	Float_t PFJet_muonEnergyFraction           ;
	Float_t PFJet_muonMultiplicity             ;
	Float_t PFJet_neutralEmEnergy              ;
	Float_t PFJet_neutralEmEnergyFraction      ;
	Float_t PFJet_neutralHadronEnergy          ;
	Float_t PFJet_neutralHadronEnergyFraction  ;
	Float_t PFJet_neutralHadronMultiplicity    ;
	Float_t PFJet_neutralMultiplicity          ;
	Float_t PFJet_photonEnergy                 ;
	Float_t PFJet_photonEnergyFraction         ;
	Float_t PFJet_photonMultiplicity           ;

	TTree *T1    = (TTree*)  F1     ->Get("tree/JetTree");

	T1->SetBranchAddress("Run"      ,  & Run        );
	T1->SetBranchAddress("Lumi"     ,  & Lumi       );
	T1->SetBranchAddress("Event"    ,  & Event      );
	T1->SetBranchAddress("Nvtx"     ,  & Nvtx       );
	T1->SetBranchAddress("NvtxGood" ,  & NvtxGood   );
	T1->SetBranchAddress("Mass"     ,  & Mass       );
	T1->SetBranchAddress("Energy"   ,  & Energy     );
	T1->SetBranchAddress("Pt"       ,  & Pt         );
	T1->SetBranchAddress("Eta"      ,  & Eta        );
	T1->SetBranchAddress("Phi"      ,  & Phi        );
	T1->SetBranchAddress("GenMass"  ,  & GenMass    );
	T1->SetBranchAddress("GenEnergy",  & GenEnergy  );
	T1->SetBranchAddress("GenPt"    ,  & GenPt      );
	T1->SetBranchAddress("GenEta"   ,  & GenEta     );
	T1->SetBranchAddress("GenPhi"   ,  & GenPhi     );


	T1->SetBranchAddress("PFJet_Nconst"                      , & PFJet_Nconst             );
	T1->SetBranchAddress("PFJet_chargedEmEnergy"             , & PFJet_chargedEmEnergy             );
	T1->SetBranchAddress("PFJet_chargedEmEnergyFraction"     , & PFJet_chargedEmEnergyFraction     );
	T1->SetBranchAddress("PFJet_chargedHadronEnergy"         , & PFJet_chargedHadronEnergy         );
	T1->SetBranchAddress("PFJet_chargedHadronEnergyFraction" , & PFJet_chargedHadronEnergyFraction );
	T1->SetBranchAddress("PFJet_chargedHadronMultiplicity"   , & PFJet_chargedHadronMultiplicity   );
	T1->SetBranchAddress("PFJet_chargedMuEnergy"             , & PFJet_chargedMuEnergy             );
	T1->SetBranchAddress("PFJet_chargedMuEnergyFraction"     , & PFJet_chargedMuEnergyFraction     );
	T1->SetBranchAddress("PFJet_chargedMultiplicity"         , & PFJet_chargedMultiplicity         );
	T1->SetBranchAddress("PFJet_electronEnergy"              , & PFJet_electronEnergy              );
	T1->SetBranchAddress("PFJet_electronEnergyFraction"      , & PFJet_electronEnergyFraction      );
	T1->SetBranchAddress("PFJet_electronMultiplicity"        , & PFJet_electronMultiplicity        );
	T1->SetBranchAddress("PFJet_HFEMEnergy"                  , & PFJet_HFEMEnergy                  );
	T1->SetBranchAddress("PFJet_HFEMEnergyFraction"          , & PFJet_HFEMEnergyFraction          );
	T1->SetBranchAddress("PFJet_HFEMMultiplicity"            , & PFJet_HFEMMultiplicity            );
	T1->SetBranchAddress("PFJet_HFHadronEnergy"              , & PFJet_HFHadronEnergy              );
	T1->SetBranchAddress("PFJet_HFHadronEnergyFraction"      , & PFJet_HFHadronEnergyFraction      );
	T1->SetBranchAddress("PFJet_HFHadronMultiplicity"        , & PFJet_HFHadronMultiplicity        );
	T1->SetBranchAddress("PFJet_muonEnergy"                  , & PFJet_muonEnergy                  );
	T1->SetBranchAddress("PFJet_muonEnergyFraction"          , & PFJet_muonEnergyFraction          );
	T1->SetBranchAddress("PFJet_muonMultiplicity"            , & PFJet_muonMultiplicity            );
	T1->SetBranchAddress("PFJet_neutralEmEnergy"             , & PFJet_neutralEmEnergy             );
	T1->SetBranchAddress("PFJet_neutralEmEnergyFraction"     , & PFJet_neutralEmEnergyFraction     );
	T1->SetBranchAddress("PFJet_neutralHadronEnergy"         , & PFJet_neutralHadronEnergy         );
	T1->SetBranchAddress("PFJet_neutralHadronEnergyFraction" , & PFJet_neutralHadronEnergyFraction );
	T1->SetBranchAddress("PFJet_neutralHadronMultiplicity"   , & PFJet_neutralHadronMultiplicity   );
	T1->SetBranchAddress("PFJet_neutralMultiplicity"         , & PFJet_neutralMultiplicity         );
	T1->SetBranchAddress("PFJet_photonEnergy"                , & PFJet_photonEnergy                );
	T1->SetBranchAddress("PFJet_photonEnergyFraction"        , & PFJet_photonEnergyFraction        );
	T1->SetBranchAddress("PFJet_photonMultiplicity"          , & PFJet_photonMultiplicity          );


	//TH1F * PT_RESPONSE               = new TH1F("RESPONSE"             ,"",100,0,4); 


	double entries = T1->GetEntries();
	cout<<"entries = "<< entries <<endl;

	// int nbins_eta = 5;
	// int nbins_pt = 17;
	// double etabins[nbins_eta] = {0.0, 1.4, 2.6, 3.2, 4.7};
	// double ptbins[nbins_pt]   = {20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.};
	// string setabins[nbins_eta] = {"0", "1p4", "2p4", "3p2", "4p6"};
	// string sptbins[nbins_pt]   = {"20", "25", "30", "35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85", "90", "95", "100"};
	// string tetabins[nbins_eta] = {"0", "1.4", "2.4", "3.2", "4.6"};
	// string tptbins[nbins_pt]   = {"20", "25", "30", "35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85", "90", "95", "100"};

	int nbins_eta = 4;
	int nbins_pt = 29;
	double etabins[nbins_eta] = {0.0, 1.3, 3, 5 };
	double ptbins[nbins_pt]   = {20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 125., 150., 200., 300., 400., 500., 700., 1000., 1500., 2000., 2500., 3000.,};
	string setabins[nbins_eta] = {"0", "1p3", "3p0", "5p0"};
	string sptbins[nbins_pt]   = {"020", "025", "030", "035", "040", "045", "050", "055", "060", "065", "070", "075", "080", "085", "090", "095", "100", "125", "150", "200", "300","400", "500", "700", "1000", "1500", "2000", "2500", "3000"};
	string tetabins[nbins_eta] = {"0", "1.3", "3.0", "5.0"};
	string tptbins[nbins_pt]   = {"20", "25", "30", "35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85", "90", "95", "100", "125", "150", "200", "300","400", "500", "700", "1000", "1500", "2000", "2500", "3000"};


	Float_t pt_resolution_vs_genpt_central_x [nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_central_y [nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_central_ex[nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_central_ey[nbins_pt-1] = {};

	Float_t pt_resolution_vs_genpt_endcap_x [nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_endcap_y [nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_endcap_ex[nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_endcap_ey[nbins_pt-1] = {};

	Float_t pt_resolution_vs_genpt_forward_x [nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_forward_y [nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_forward_ex[nbins_pt-1] = {};
	Float_t pt_resolution_vs_genpt_forward_ey[nbins_pt-1] = {};

	Float_t pt_response_mean_vs_genpt_central_x [nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_central_y [nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_central_ex[nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_central_ey[nbins_pt-1] = {};

	Float_t pt_response_mean_vs_genpt_endcap_x [nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_endcap_y [nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_endcap_ex[nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_endcap_ey[nbins_pt-1] = {};

	Float_t pt_response_mean_vs_genpt_forward_x [nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_forward_y [nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_forward_ex[nbins_pt-1] = {};
	Float_t pt_response_mean_vs_genpt_forward_ey[nbins_pt-1] = {};

	Float_t pt_response_mean_vs_geneta_pt20_x [nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt20_y [nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt20_ex[nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt20_ey[nbins_eta-1] = {};

	Float_t pt_response_mean_vs_geneta_pt25_x [nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt25_y [nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt25_ex[nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt25_ey[nbins_eta-1] = {};

	Float_t pt_response_mean_vs_geneta_pt300_x [nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt300_y [nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt300_ex[nbins_eta-1] = {};
	Float_t pt_response_mean_vs_geneta_pt300_ey[nbins_eta-1] = {};



	TH1F * PT_RESPONSE[nbins_eta-1][nbins_pt-1];

	for (int etabin =0; etabin< nbins_eta-1; etabin++)
  	{
		for (int ptbin =0; ptbin< nbins_pt-1; ptbin++)
		{		  			
			cout<<endl;
	  		cout<<etabin<<" "<<ptbin<<endl;

	  		string histoname1  = "PT_RESPONSE_ETA_"+setabins[etabin]+"_"+setabins[etabin+1]+"_PT_"+sptbins[ptbin]+"_"+sptbins[ptbin+1];
			string histotitle1  = "Jet p_{T} response ( "+tetabins[etabin]+" < | #eta | < "+tetabins[etabin+1]+" ,  "+sptbins[ptbin]+" < p_{T} < "+sptbins[ptbin+1]+" GeV/c );p_{T}(reco)/p_{T}(gen);Number of jets";
			cout<<histoname1 <<endl;//<<"       -       " <<histotitle1<<endl;
			PT_RESPONSE[etabin][ptbin]  = new TH1F(histoname1.c_str() ,histotitle1.c_str(),100,0,4); 

			for (int i=0; i<entries; i++ )
		 	{  		    
		 		T1->GetEntry(i);

				//  Two more cuts in the PF Jet ID definition are recommended:
				// Muon Fraction (muf) < 0.8 AND Charge Electromagnetic Fraction (elf) < 0.9

				// PF Jet ID 	Loose (Recommended) 	Medium 	Tight
				// Neutral Hadron Fraction 	< 0.99 	< 0.95 	< 0.90
				// Neutral EM Fraction 	< 0.99 	< 0.95 	< 0.90
				// Number of Constituents 	> 1 	> 1 	> 1
				// Muon Fraction 	< 0.8 	< 0.8 	< 0.8
				// Charge Electromagnetic Fraction 	< 0.9 	< 0.9 	< 0.9
				// And for η < 2.4 , η > -2.4 in addition apply
				// Charged Hadron Fraction 	> 0 	> 0 	> 0
				// Charged Multiplicity 	> 0 	> 0 	> 0
				// Charged EM Fraction 	< 0.99 	< 0.99 	< 0.99 

				double muf   = 0;
				double muef  = 0;
				double chf   = 0;

		 		if (PFJet_Nconst!=0) muf  = PFJet_muonMultiplicity / PFJet_Nconst;
		 		if (Energy!=0) muef = PFJet_muonEnergyFraction / Energy;
		 		if (PFJet_Nconst!=0) chf  = PFJet_chargedHadronMultiplicity / PFJet_Nconst;

		 		// all-eta PF ID
		 		if ( PFJet_neutralHadronEnergyFraction > 0.99 ) continue;
		 		if ( PFJet_neutralEmEnergyFraction     > 0.99 ) continue;
		 		if ( PFJet_Nconst                     <= 1    ) continue;
		 		if ( PFJet_chargedMuEnergyFraction     > 0.8  ) continue; 

		 		// central PF ID
		 		if ( fabs(Eta)<2.4 && PFJet_chargedEmEnergyFraction      > 0.99 ) continue;
		 		if ( fabs(Eta)<2.4 && PFJet_chargedMultiplicity         <= 0    ) continue;
		 		if ( fabs(Eta)<2.4 && PFJet_chargedHadronEnergyFraction <= 0    ) continue;
		 		if ( fabs(Eta)<2.4 && PFJet_chargedEmEnergyFraction      > 0.99 ) continue;

				// PFJet_chargedEmEnergy            
				// PFJet_chargedEmEnergyFraction    
				// PFJet_chargedHadronEnergy        
				// PFJet_chargedHadronEnergyFraction
				// PFJet_chargedHadronMultiplicity  
				// PFJet_chargedMuEnergy            
				// PFJet_chargedMuEnergyFraction    
				// PFJet_chargedMultiplicity        
				// PFJet_electronEnergy             
				// PFJet_electronEnergyFraction     
				// PFJet_electronMultiplicity       
				// PFJet_HFEMEnergy                 
				// PFJet_HFEMEnergyFraction         
				// PFJet_HFEMMultiplicity           
				// PFJet_HFHadronEnergy             
				// PFJet_HFHadronEnergyFraction     
				// PFJet_HFHadronMultiplicity       
				// PFJet_muonEnergy                 
				// PFJet_muonEnergyFraction         
				// PFJet_muonMultiplicity           
				// PFJet_neutralEmEnergy            
				// PFJet_neutralEmEnergyFraction    
				// PFJet_neutralHadronEnergy        
				// PFJet_neutralHadronEnergyFraction
				// PFJet_neutralHadronMultiplicity  
				// PFJet_neutralMultiplicity        
				// PFJet_photonEnergy               
				// PFJet_photonEnergyFraction       
				// PFJet_photonMultiplicity         

		  		if ( fabs(GenEta)>etabins[etabin] && fabs(GenEta)<etabins[etabin+1] && GenPt > ptbins[ptbin] && GenPt < ptbins[ptbin+1] )
		  		{
	  				if (GenPt!=0) PT_RESPONSE[etabin][ptbin]->Fill(Pt/GenPt);
	  			}
		 	}

		 

			TF1 *gfit = new TF1("Gaussian","gaus",0.6,1.4); // Create the fit function
			gfit->SetLineColor(2);
			gfit->SetLineWidth(2.5);

			PT_RESPONSE[etabin][ptbin]->Fit(gfit,"RQ"); // Fit histogram h"
			//*-*     gaus(0) is a substitute for [0]*exp(-0.5*((x-[1])/[2])**2)
			//double chi2ndf= gfit->GetChisquare()/gfit->GetNDF();
			double amp    = gfit->GetParameter(0);
			double eamp   = gfit->GetParError(0); 
			double mean   = gfit->GetParameter(1);
			double emean  = gfit->GetParError(1); 
			double width  = gfit->GetParameter(2);
			double ewidth = gfit->GetParError(1); 

			if (etabin==0) pt_response_mean_vs_genpt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
			if (etabin==0) pt_response_mean_vs_genpt_central_y [ptbin] = mean;
			if (etabin==0) pt_response_mean_vs_genpt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
			if (etabin==0) pt_response_mean_vs_genpt_central_ey[ptbin] = emean;

			if (etabin==1) pt_response_mean_vs_genpt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
			if (etabin==1) pt_response_mean_vs_genpt_endcap_y [ptbin] = mean;
			if (etabin==1) pt_response_mean_vs_genpt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
			if (etabin==1) pt_response_mean_vs_genpt_endcap_ey[ptbin] = emean;

			if (etabin==2) pt_response_mean_vs_genpt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
			if (etabin==2) pt_response_mean_vs_genpt_forward_y [ptbin] = mean;
			if (etabin==2) pt_response_mean_vs_genpt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
			if (etabin==2) pt_response_mean_vs_genpt_forward_ey[ptbin] = emean;

   
			if (etabin==0) pt_resolution_vs_genpt_central_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
			if (etabin==0) pt_resolution_vs_genpt_central_y [ptbin] = width/mean;
			if (etabin==0) pt_resolution_vs_genpt_central_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
			if (etabin==0) pt_resolution_vs_genpt_central_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			if (etabin==1) pt_resolution_vs_genpt_endcap_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
			if (etabin==1) pt_resolution_vs_genpt_endcap_y [ptbin] = width/mean;
			if (etabin==1) pt_resolution_vs_genpt_endcap_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
			if (etabin==1) pt_resolution_vs_genpt_endcap_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );
			if (etabin==2) pt_resolution_vs_genpt_forward_x [ptbin] = (ptbins[ptbin+1]+ptbins[ptbin])/2;
			if (etabin==2) pt_resolution_vs_genpt_forward_y [ptbin] = width/mean;
			if (etabin==2) pt_resolution_vs_genpt_forward_ex[ptbin] = (ptbins[ptbin+1]-ptbins[ptbin])/2;
			if (etabin==2) pt_resolution_vs_genpt_forward_ey[ptbin] = width/mean*sqrt( (emean/mean)*(emean/mean) + (ewidth/width)*(ewidth/width)   );

   			if (ptbin == 0)  { pt_response_mean_vs_geneta_pt20_x [etabin] = (etabins[etabin+1]+etabins[etabin])/2;; cout<<" etabins[etabin] "<<etabins[etabin]<<" mean "<<mean<<endl; }
   			if (ptbin == 0)  { pt_response_mean_vs_geneta_pt20_y [etabin] = mean;            }
   			if (ptbin == 0)  { pt_response_mean_vs_geneta_pt20_ex[etabin] = (etabins[etabin+1]-etabins[etabin])/2;               }
   			if (ptbin == 0)  { pt_response_mean_vs_geneta_pt20_ey[etabin] = emean;           }

   			if (ptbin == 1)  { pt_response_mean_vs_geneta_pt25_x [etabin] = (etabins[etabin+1]+etabins[etabin])/2; }
   			if (ptbin == 1)  { pt_response_mean_vs_geneta_pt25_y [etabin] = mean;            } 
   			if (ptbin == 1)  { pt_response_mean_vs_geneta_pt25_ex[etabin] = (etabins[etabin+1]-etabins[etabin])/2; ;               }
   			if (ptbin == 1)  { pt_response_mean_vs_geneta_pt25_ey[etabin] = emean;           }

   			if (ptbin == 20)  pt_response_mean_vs_geneta_pt300_x [etabin] = (etabins[etabin+1]+etabins[etabin])/2;
   			if (ptbin == 20)  pt_response_mean_vs_geneta_pt300_y [etabin] = mean;
   			if (ptbin == 20)  pt_response_mean_vs_geneta_pt300_ex[etabin] = (etabins[etabin+1]-etabins[etabin])/2; ;
   			if (ptbin == 20)  pt_response_mean_vs_geneta_pt300_ey[etabin] = emean;





 			PT_RESPONSE[etabin][ptbin]->Draw();
		 	string savename1 = "plots/"+histoname1+".pdf";
		 	c1236->SaveAs(savename1.c_str());


	  	}
	}
	
	TGraphErrors *
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_central = new TGraphErrors(nbins_pt-1,pt_response_mean_vs_genpt_central_x,pt_response_mean_vs_genpt_central_y,pt_response_mean_vs_genpt_central_ex,pt_response_mean_vs_genpt_central_ey);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_central->SetTitle(";GenJet p_{T} (GeV/c);<p_{T}(gen)/p_{T}(reco)>");
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_central->SetMarkerColor(4);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_central->SetMarkerStyle(21);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_central->Draw("AP");
   	c1236->SetGrid();
   	c1236->SaveAs("plots/MEAN_PT_RESPONSE_VS_GEN_PT_central.pdf");

	TGraphErrors *
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap = new TGraphErrors(nbins_pt-1,pt_response_mean_vs_genpt_endcap_x,pt_response_mean_vs_genpt_endcap_y,pt_response_mean_vs_genpt_endcap_ex,pt_response_mean_vs_genpt_endcap_ey);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap->SetTitle(";GenJet p_{T} (GeV/c);<p_{T}(gen)/p_{T}(reco)>");
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap->SetMarkerColor(3);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap->SetMarkerStyle(21);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap->Draw("AP");
   	c1236->SetGrid();
   	c1236->SaveAs("plots/MEAN_PT_RESPONSE_VS_GEN_PT_endcap.pdf");

	TGraphErrors *
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward = new TGraphErrors(nbins_pt-1,pt_response_mean_vs_genpt_forward_x,pt_response_mean_vs_genpt_forward_y,pt_response_mean_vs_genpt_forward_ex,pt_response_mean_vs_genpt_forward_ey);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward->SetTitle(";GenJet p_{T} (GeV/c);<p_{T}(gen)/p_{T}(reco)>");
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward->SetMarkerColor(2);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward->SetMarkerStyle(21);
   	gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward->Draw("AP");
   	c1236->SetGrid();
   	c1236->SaveAs("plots/MEAN_PT_RESPONSE_VS_GEN_PT_forward.pdf");

	gr_MEAN_PT_RESPONSE_VS_GEN_PT_central ->Draw("AP");  
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  ->Draw("P"); 
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward ->Draw("P"); 

	TAxis *axis = gr_MEAN_PT_RESPONSE_VS_GEN_PT_central->GetYaxis();
    axis->SetLimits(0.5,1.2);      

	c1236->SetLogx(); 
   	c1236->SaveAs("plots/MEAN_PT_RESPONSE_VS_GEN_PT_all.pdf");
	c1236->SetLogx(0); 

	TGraphErrors *
	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20 = new TGraphErrors(nbins_eta-1,pt_response_mean_vs_geneta_pt20_x,pt_response_mean_vs_geneta_pt20_y,pt_response_mean_vs_geneta_pt20_ex,pt_response_mean_vs_geneta_pt20_ey);
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20->SetTitle(";GenJet p_{T} (GeV/c);<p_{T}(gen)/p_{T}(reco)>");
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20->SetMarkerColor(1);
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20->SetMarkerStyle(21);
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20->Draw("AP");
   	c1236->SetGrid();
   	c1236->SaveAs("plots/MEAN_PT_RESPONSE_VS_GEN_ETA_pt20.pdf");

	TGraphErrors *
	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25 = new TGraphErrors(nbins_eta-1,pt_response_mean_vs_geneta_pt25_x,pt_response_mean_vs_geneta_pt25_y,pt_response_mean_vs_geneta_pt25_ex,pt_response_mean_vs_geneta_pt25_ey);
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25->SetTitle(";GenJet p_{T} (GeV/c);<p_{T}(gen)/p_{T}(reco)>");
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25->SetMarkerColor(2);
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25->SetMarkerStyle(21);
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25->Draw("AP");
   	c1236->SetGrid();
   	c1236->SaveAs("plots/MEAN_PT_RESPONSE_VS_GEN_ETA_pt25.pdf");

   	TGraphErrors *
	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300 = new TGraphErrors(nbins_eta-1,pt_response_mean_vs_geneta_pt300_x,pt_response_mean_vs_geneta_pt300_y,pt_response_mean_vs_geneta_pt300_ex,pt_response_mean_vs_geneta_pt300_ey);
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300->SetTitle(";GenJet p_{T} (GeV/c);<p_{T}(gen)/p_{T}(reco)>");
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300->SetMarkerColor(3);
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300->SetMarkerStyle(21);
   	c1236->SetGrid();
   	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300->Draw("AP");
   	c1236->SaveAs("plots/MEAN_PT_RESPONSE_VS_GEN_ETA_pt300.pdf");


   	TGraphErrors *
	gr_PT_RESOLUTION_VS_GEN_PT_central = new TGraphErrors(nbins_pt-1,pt_resolution_vs_genpt_central_x,pt_resolution_vs_genpt_central_y,pt_resolution_vs_genpt_central_ex,pt_resolution_vs_genpt_central_ey);
   	gr_PT_RESOLUTION_VS_GEN_PT_central->SetTitle(";GenJet p_{T} (GeV/c);#sigma(p_{T}(gen)/p_{T}(reco))/<p_{T}(gen)/p_{T}(reco)>");
   	gr_PT_RESOLUTION_VS_GEN_PT_central->SetMarkerColor(4);
   	gr_PT_RESOLUTION_VS_GEN_PT_central->SetMarkerStyle(21);
   	gr_PT_RESOLUTION_VS_GEN_PT_central->Draw("AP");
   	c1236->SetGrid();
   	c1236->SaveAs("plots/PT_RESOLUTION_VS_GEN_PT_central.pdf");

	TGraphErrors *
	gr_PT_RESOLUTION_VS_GEN_PT_endcap = new TGraphErrors(nbins_pt-1,pt_resolution_vs_genpt_endcap_x,pt_resolution_vs_genpt_endcap_y,pt_resolution_vs_genpt_endcap_ex,pt_resolution_vs_genpt_endcap_ey);
   	gr_PT_RESOLUTION_VS_GEN_PT_endcap->SetTitle(";GenJet p_{T} (GeV/c);#sigma(p_{T}(gen)/p_{T}(reco))/<p_{T}(gen)/p_{T}(reco)>");
   	gr_PT_RESOLUTION_VS_GEN_PT_endcap->SetMarkerColor(3);
   	gr_PT_RESOLUTION_VS_GEN_PT_endcap->SetMarkerStyle(21);
   	gr_PT_RESOLUTION_VS_GEN_PT_endcap->Draw("AP");
   	c1236->SetGrid();
   	c1236->SaveAs("plots/PT_RESOLUTION_VS_GEN_PT_endcap.pdf");

	TGraphErrors *
	gr_PT_RESOLUTION_VS_GEN_PT_forward = new TGraphErrors(nbins_pt-1,pt_resolution_vs_genpt_forward_x,pt_resolution_vs_genpt_forward_y,pt_resolution_vs_genpt_forward_ex,pt_resolution_vs_genpt_forward_ey);
   	gr_PT_RESOLUTION_VS_GEN_PT_forward->SetTitle(";GenJet p_{T} (GeV/c);#sigma(p_{T}(gen)/p_{T}(reco))/<p_{T}(gen)/p_{T}(reco)>");
   	gr_PT_RESOLUTION_VS_GEN_PT_forward->SetMarkerColor(2);
   	gr_PT_RESOLUTION_VS_GEN_PT_forward->SetMarkerStyle(21);
   	gr_PT_RESOLUTION_VS_GEN_PT_forward->Draw("AP");
   	c1236->SetGrid();
   	c1236->SaveAs("plots/PT_RESOLUTION_VS_GEN_PT_forward.pdf");

	gr_PT_RESOLUTION_VS_GEN_PT_central ->Draw("AP");  
	gr_PT_RESOLUTION_VS_GEN_PT_endcap  ->Draw("P"); 
	gr_PT_RESOLUTION_VS_GEN_PT_forward ->Draw("P"); 

	TAxis *axis = gr_PT_RESOLUTION_VS_GEN_PT_central->GetYaxis();
    axis->SetLimits(0.5,1.2);      

	c1236->SetLogx(); 
   	c1236->SaveAs("plots/PT_RESOLUTION_VS_GEN_PT_all.pdf");
	c1236->SetLogx(0); 




	TFile *Out;
	Out = new TFile("TGraphs.root","RECREATE");
	Out->cd();

	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20   -> SetName("gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20");
	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25   -> SetName("gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25");
	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300  -> SetName("gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300");
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_central -> SetName("gr_MEAN_PT_RESPONSE_VS_GEN_PT_central");
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  -> SetName("gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap");
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward -> SetName("gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward");
	gr_PT_RESOLUTION_VS_GEN_PT_central    -> SetName("gr_PT_RESOLUTION_VS_GEN_PT_central");
	gr_PT_RESOLUTION_VS_GEN_PT_endcap     -> SetName("gr_PT_RESOLUTION_VS_GEN_PT_endcap");
	gr_PT_RESOLUTION_VS_GEN_PT_forward    -> SetName("gr_PT_RESOLUTION_VS_GEN_PT_forward");

	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt20  ->Write();         
	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt25  ->Write();         
	gr_MEAN_PT_RESPONSE_VS_GEN_ETA_pt300 ->Write();          

	gr_MEAN_PT_RESPONSE_VS_GEN_PT_central ->Write();
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_endcap  ->Write();
	gr_MEAN_PT_RESPONSE_VS_GEN_PT_forward ->Write();

	gr_PT_RESOLUTION_VS_GEN_PT_central ->Write();
	gr_PT_RESOLUTION_VS_GEN_PT_endcap  ->Write();
	gr_PT_RESOLUTION_VS_GEN_PT_forward ->Write();

 
	Out->ls();
	Out->Write();

 }