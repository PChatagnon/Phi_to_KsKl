#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "bib/KsKl_analysis.h"
#include "bib/KsKl_Event.h"
#include "bib/InputParser.h"

#include "hipo4/reader.h"
#include "rcdb_reader.h"

// QADB header and namespace
#include "QADB.h"
using namespace QA;

#include <ctime> // time_t
#include <cstdio>
using namespace std;

#define ADDVAR(x, name, t, tree) tree->Branch(name, x, TString(name) + TString(t))



int analysis_Phi_KsKl()
{

	time_t begin, intermediate, end; // time_t is a datatype to store time values.

	time(&begin); // note time before execution

	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(111);
	gStyle->SetPalette(55);
	gStyle->SetLabelSize(.05, "xyz");
	gStyle->SetTitleSize(.05, "xyz");
	gStyle->SetTitleSize(.07, "t");
	gStyle->SetMarkerStyle(13);
	gStyle->SetOptFit(1);


	Int_t argc = gApplication->Argc();
	char **argv = gApplication->Argv();
	Input input(argc, argv);

	

	/////////Parse command line/////////////
	bool option = input.cmdOptionExists("-option");
	all_Gen_vector = input.cmdOptionExists("-all_Gen_vector");
	

	/////////////////////////////////////////

	if (input.cmdOptionExists("-energy"))
	{
		ebeam = std::stof(input.getCmdOption("-energy"));
	}

	cout<<"Energy of the beam: "<<ebeam<<endl;
	/////////End parse command line/////////////


	double nbrecEvent = 0;
	TString nameFiles = "";
	TString type = "REC";


	///////////////////////////////////////////
	// Setup the TTree output
	///////////////////////////////////////////
	TString output_file = (TString)(input.getCmdOption("-o")); // argv[4]);
	TFile *outFile = new TFile(Form("output_muCLAS12_%s.root", output_file.Data()), "recreate");
	
	TTree *outT = new TTree("tree", "tree");
	TTree *outT_Gen = new TTree("tree_Gen", "tree_Gen");


	///////////////////////////////////////////
	// Create REC tree
	///////////////////////////////////////////
	TLorentzVector tree_Electron, tree_Pi_plus, tree_Pi_minus, tree_Proton, tree_Missing;
	outT->Branch("Electron", "TLorentzVector", &tree_Electron);
	outT->Branch("PiPlus", "TLorentzVector", &tree_Pi_plus);
	outT->Branch("PiMinus", "TLorentzVector", &tree_Pi_minus);
	outT->Branch("Proton", "TLorentzVector", &tree_Proton);
	outT->Branch("Missing", "TLorentzVector", &tree_Missing);


	std::vector<TString> fvars = {
		"evt_num",
		"run",
		"analysis_stage",
		"topology",
		
		"Phi_Mass",
	    "Ks_Mass",
		"Kl_Mass",
	    "Q2",
        "t",
	    "W",

		"vt_elec",
		"vt_pi_plus",
		"vt_pi_minus",
		"vt_proton",

		"rec_e",
		"rec_pi_p",
		"rec_pi_m",
		"rec_p",

		"status_elec",
		"status_pi_plus",
		"status_pi_minus",
		"status_proton",



	};

	std::map<TString, Float_t> outVars;
	for (size_t i = 0; i < fvars.size(); i++)
	{
		outVars[fvars[i]] = 0.;
		ADDVAR(&(outVars[fvars[i]]), fvars[i], "/F", outT);
	}

	///////////////////////////////////////////
	// Create GEN tree
	///////////////////////////////////////////

	TString fvars_Gen[] = {
		"M_Gen", 
		"Q2_Gen",
		"t_Gen",
		};

	std::map<TString, Float_t> outVars_Gen;
	for (size_t i = 0; i < sizeof(fvars_Gen) / sizeof(TString); i++)
	{
		outVars_Gen[fvars_Gen[i]] = 0.;
		ADDVAR(&(outVars_Gen[fvars_Gen[i]]), fvars_Gen[i], "/F", outT_Gen);
	}

	// Add the 4vectors to the GEN tree
	TLorentzVector gen_Electron, gen_Pi_plus, gen_Pi_minus, gen_Proton;
	cout<<"Include all gen particles"<<endl;
	outT_Gen->Branch("gen_Electron", "TLorentzVector", &gen_Electron);
	outT_Gen->Branch("gen_Pi_plus", "TLorentzVector", &gen_Pi_plus);
	outT_Gen->Branch("gen_Pi_minus", "TLorentzVector", &gen_Pi_minus);
	outT_Gen->Branch("gen_Proton", "TLorentzVector", &gen_Proton);

	outT->Branch("gen_Electron", "TLorentzVector", &gen_Electron);
	outT->Branch("gen_Pi_plus", "TLorentzVector", &gen_Pi_plus);
	outT->Branch("gen_Pi_minus", "TLorentzVector", &gen_Pi_minus);
	outT->Branch("gen_Proton", "TLorentzVector", &gen_Proton);
	
	///////////////////////////////////////////
	


	////////////////////////////////////////////
	// Get file name
	////////////////////////////////////////////
	int nbf = 0;
	int nbEvent = 0;
	for (Int_t i = input.getCmdIndex("-f") + 2; i < input.getCmdIndex("-ef") + 1; i++)
	{
		nbf++;
		nameFiles = TString(argv[i]);

		////////////////////////////////////////////
		// hipo reader
		hipo::reader reader;
		hipo::dictionary factory;
		hipo::event hipo_event;
		////////////////////////////////////////////
		
		reader.open(nameFiles);
		reader.readDictionary(factory);
		
		hipo::bank EVENT(factory.getSchema("REC::Event"));
		hipo::bank PART(factory.getSchema("REC::Particle"));
		hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
		hipo::bank CHE(factory.getSchema("REC::Cherenkov"));
		hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
		hipo::bank RUN(factory.getSchema("RUN::config"));
		hipo::bank MCPART(factory.getSchema("MC::Particle"));
		hipo::bank MCEVENT(factory.getSchema("MC::Event"));
		hipo::bank TRACK(factory.getSchema("REC::Track"));
		hipo::bank TRAJ(factory.getSchema("REC::Traj"));

		outFile->cd();

		while (reader.next())
		{

			nbEvent++;
			if (nbEvent % 500000 == 0)
			{
				time(&intermediate);
				double intermediate_time = difftime(intermediate, begin);
				cout << nbEvent << " events processed in " << intermediate_time << "s" << "\n";
			}

			KsKl_Event ev;
			//muMCEvent MC_ev;

			int run = 0;
			int event_nb = 0;
			
			// Get banks
			reader.read(hipo_event);
			hipo_event.getStructure(MCPART);
			hipo_event.getStructure(MCEVENT);
			hipo_event.getStructure(RUN);
			hipo_event.getStructure(PART);
			hipo_event.getStructure(SCIN);
			hipo_event.getStructure(CHE);
			hipo_event.getStructure(CALO);
			hipo_event.getStructure(EVENT);
			hipo_event.getStructure(TRAJ);
			hipo_event.getStructure(TRACK);
            
			int np_input = PART.getRows();
			ev.Set_nb_part(np_input);
			

			//if(false){
			//	//Assign MC particle and compute kinematics
			//	MC_ev.Set_MC_Particles(MCEVENT, MCPART, isElSpectro, isGrape, isCoincidence, isCoincidence_Quasi, IsInelastic);
			//	MC_ev.Get_Kinematics();
			//	//Fill the branches
			//	outVars_Gen["M_Gen"] = MC_ev.M_Gen;
			//	outVars_Gen["Q2_Gen"] = MC_ev.Q2_Gen;
			//	outVars_Gen["t_Gen"] = MC_ev.t_Gen;
			//	gen_Electron = MC_ev.Electron;
			//	gen_mu_plus = MC_ev.mu_plus;
			//	gen_mu_minus = MC_ev.mu_minus;
			//	gen_Proton = MC_ev.Proton;
			//	//Fill the tree
            //	outT_Gen->Fill();
			//}
			
			///////////////////////////////////////////
			// Get Particles and cut on event topology
			///////////////////////////////////////////
			ev.Set_Particles(PART);

			///////////////////////////////////////////
			// Compute kinematics (before the topology cut using these kinematics)
			///////////////////////////////////////////
			ev.Get_Kinematics();

			///////////////////////////////////////////
			// Topology cuts and loose kinematic cuts
			///////////////////////////////////////////
			if (!ev.pass_topology_cut())
			{
				continue;
			}

			///////////////////////////////////////////
			// Associate detector responses and do EC cuts
			///////////////////////////////////////////
			ev.Associate_detector_resp(CHE, SCIN, CALO);
			//ev.Associate_DC_traj(TRAJ);
			//ev.Set_Nphe_HTCC();
			///////////////////////////////////////////

			


			///////////////////////////////////////////
			// Fill the branches of the tree
			///////////////////////////////////////////
			outVars["evt_num"] = nbEvent;
			outVars["run"] = ev.run;

			outVars["Phi_Mass"] = ev.Phi_Mass;
			outVars["Ks_Mass"] = ev.Ks_Mass;
			outVars["Kl_Mass"] = ev.Kl_Mass;
			outVars["Q2"] = ev.Q2;
			outVars["t"] = ev.t;
			outVars["W"] = ev.W;

			outVars["vt_elec"] = ev.Electron.vt;
			outVars["vt_pi_plus"] = ev.Pi_plus.vt;
			outVars["vt_pi_minus"] = ev.Pi_minus.vt;
			outVars["vt_proton"] = ev.Proton.vt;

			outVars["rec_e"] = ev.rec_e;
			outVars["rec_pi_p"] = ev.rec_pi_p;
			outVars["rec_pi_m"] = ev.rec_pi_m;
			outVars["rec_p"] = ev.rec_p;

			outVars["status_elec"] = ev.Electron.status;
			outVars["status_pi_plus"] = ev.Pi_plus.status;
			outVars["status_pi_minus"] = ev.Pi_minus.status;
			outVars["status_proton"] = ev.Proton.status;

			tree_Electron = ev.Electron.Vector;
			tree_Pi_plus = ev.Pi_plus.Vector;
			tree_Pi_minus = ev.Pi_minus.Vector;
			tree_Proton = ev.Proton.Vector;
			tree_Missing = ev.vMissing;
			
			///////////////////////////////////////////
			// Fill the tree
			///////////////////////////////////////////
			outT->Fill();
			
		}
	}

	///////////////////////////////////////////
	// Save REC and GEN tree + timing
	///////////////////////////////////////////
	

	outFile->cd();
	outT->Write();
	outFile->Write();
	outFile->Close();

	cout << "Tree written" << endl;
	cout << "nb of file " << nbf << "\n";
	cout << "nb of events " << nbEvent << "\n";

	// gROOT->ProcessLine(".q");

	time(&end); // note time after execution

	double difference = difftime(end, begin);
	printf("All this work done in only %.2lf seconds. Congratulations !\n", difference);

	gApplication->Terminate();

	return 0;
}
