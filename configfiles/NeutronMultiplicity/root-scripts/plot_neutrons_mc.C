void plot_neutrons_mc(const char* infile, const char* outfile, bool verbose = false){

	TFile *f = new TFile(infile,"READ");
	TTree *tTrigger = (TTree*) f->Get("phaseIITriggerTree");
	TTree *tMRD = (TTree*) f->Get("phaseIIMRDClusterTree");
	TTree *tPMT = (TTree*) f->Get("phaseIITankClusterTree");

	TFile *fout = new TFile(outfile,"RECREATE");

	int nentries_trigger = tTrigger->GetEntries();
	int nentries_mrd = tMRD->GetEntries();
	int nentries_pmt = tPMT->GetEntries();

	std::cout << "Entries Trigger/MRD/PMT: "<<nentries_trigger<<","<<nentries_mrd<<","<<nentries_pmt<<std::endl;

	int run_trigger;
	int run_mrd;
	int run_pmt;

	int ev_trigger;
	int ev_mrd;
	int ev_pmt;

	int beam_ok;
	int num_mrd_tracks;
        int has_tank;
	std::vector<double>* eloss = new std::vector<double>;
	std::vector<double>* tracklength = new std::vector<double>;
	std::vector<bool>* mrdstop = new std::vector<bool>;
	std::vector<double>* mrdstartx = new std::vector<double>;
	std::vector<double>* mrdstarty = new std::vector<double>;
	std::vector<double>* mrdstartz = new std::vector<double>;
	std::vector<double>* mrdstopx = new std::vector<double>;
	std::vector<double>* mrdstopy = new std::vector<double>;
	std::vector<double>* mrdstopz = new std::vector<double>;

	double cluster_time;
	double cluster_cb;
	double cluster_pe;
	int extended;

	int tankmrdcoinc;
	int noveto;
	double true_muone;
	int true_pdg;
	double true_vtxx;
	double true_vtxy;
	double true_vtxz;
	
	int true_neutrons;
	int true_protons;
	int true_mec;
	int true_dis;
	int true_res;
	int true_qel;
	int true_cc;
	double true_q2;
	double true_neutrino_energy;
	std::vector<double> *true_neut_cap_e = new std::vector<double>;
	std::vector<double> *true_neut_cap_time = new std::vector<double>;
	std::vector<double> *true_neut_cap_nucl = new std::vector<double>;
	std::vector<double> *true_neut_cap_vtxx = new std::vector<double>;
	std::vector<double> *true_neut_cap_vtxy = new std::vector<double>;
	std::vector<double> *true_neut_cap_vtxz = new std::vector<double>;
	std::vector<double> *true_neut_cap_gammas = new std::vector<double>;
	std::vector<int> *true_primary_pdgs = new std::vector<int>;
	int trueMultiRing;

	tTrigger->SetBranchAddress("runNumber",&run_trigger);
	tMRD->SetBranchAddress("runNumber",&run_mrd);
	tPMT->SetBranchAddress("runNumber",&run_pmt);
	tTrigger->SetBranchAddress("eventNumber",&ev_trigger);
	tMRD->SetBranchAddress("eventNumber",&ev_mrd);
	tPMT->SetBranchAddress("eventNumber",&ev_pmt);
	tTrigger->SetBranchAddress("beam_ok",&beam_ok);
	tTrigger->SetBranchAddress("numMRDTracks",&num_mrd_tracks);
	tTrigger->SetBranchAddress("MRDTrackLength",&tracklength);
	tTrigger->SetBranchAddress("MRDEnergyLoss",&eloss);
        tTrigger->SetBranchAddress("MRDTrackStartX",&mrdstartx);
        tTrigger->SetBranchAddress("MRDTrackStartY",&mrdstarty);
        tTrigger->SetBranchAddress("MRDTrackStartZ",&mrdstartz);
        tTrigger->SetBranchAddress("MRDTrackStopX",&mrdstopx);
        tTrigger->SetBranchAddress("MRDTrackStopY",&mrdstopy);
        tTrigger->SetBranchAddress("MRDTrackStopZ",&mrdstopz);
	tTrigger->SetBranchAddress("MRDStop",&mrdstop);
        tTrigger->SetBranchAddress("HasTank",&has_tank);
	tPMT->SetBranchAddress("clusterTime",&cluster_time);
	tPMT->SetBranchAddress("clusterPE",&cluster_pe);
	tPMT->SetBranchAddress("clusterChargeBalance",&cluster_cb);
	tPMT->SetBranchAddress("Extended",&extended);
	tTrigger->SetBranchAddress("TankMRDCoinc",&tankmrdcoinc);
	tTrigger->SetBranchAddress("NoVeto",&noveto);
	tTrigger->SetBranchAddress("trueMuonEnergy",&true_muone);
	tTrigger->SetBranchAddress("truePrimaryPdg",&true_pdg);
	tTrigger->SetBranchAddress("trueVtxX",&true_vtxx);
	tTrigger->SetBranchAddress("trueVtxY",&true_vtxy);
	tTrigger->SetBranchAddress("trueVtxZ",&true_vtxz);
	tTrigger->SetBranchAddress("truePrimaryPdgs",&true_primary_pdgs);
	tTrigger->SetBranchAddress("trueNeutCapVtxX",&true_neut_cap_vtxx);
	tTrigger->SetBranchAddress("trueNeutCapVtxY",&true_neut_cap_vtxy);
	tTrigger->SetBranchAddress("trueNeutCapVtxZ",&true_neut_cap_vtxz);
	tTrigger->SetBranchAddress("trueNeutCapNucleus",&true_neut_cap_nucl);
	tTrigger->SetBranchAddress("trueNeutCapTime",&true_neut_cap_time);
	tTrigger->SetBranchAddress("trueNeutCapGammas",&true_neut_cap_gammas);
	tTrigger->SetBranchAddress("trueNeutCapE",&true_neut_cap_e);
	tTrigger->SetBranchAddress("trueNeutrinoEnergy",&true_neutrino_energy);
	tTrigger->SetBranchAddress("trueQ2",&true_q2);
	tTrigger->SetBranchAddress("trueCC",&true_cc);
	tTrigger->SetBranchAddress("trueQEL",&true_qel);
	tTrigger->SetBranchAddress("trueRES",&true_res);
	tTrigger->SetBranchAddress("trueDIS",&true_dis);
	tTrigger->SetBranchAddress("trueMEC",&true_mec);
	tTrigger->SetBranchAddress("trueNeutrons",&true_neutrons);
	tTrigger->SetBranchAddress("trueProtons",&true_protons);
	tTrigger->SetBranchAddress("trueMultiRing",&trueMultiRing);

	TH1F *h_time_neutrons = new TH1F("h_time_neutrons","Cluster time beam neutrons",200,10000,70000);
	TH1F *h_time_neutrons_mrdstop = new TH1F("h_time_neutrons_mrdstop","Cluster time beam neutrons",200,10000,70000);
	TH1F *h_neutrons = new TH1F("h_neutrons","Number of neutrons in beam events",10,0,10);
	TH1F *h_neutrons_mrdstop = new TH1F("h_neutrons_mrdstop","Number of of neutrons in beam events (MRD stop)",10,0,10);
	TH1F *h_neutrons_mrdstop_fv = new TH1F("h_neutrons_mrdstop_fv","Number of of neutrons in beam events (MRD stop, FV)",10,0,10);
	TH2F *h_neutrons_energy = new TH2F("h_neutrons_energy","Neutron multiplicity vs muon energy",10,0,2000,20,0,20);
	TH2F *h_neutrons_energy_fv = new TH2F("h_neutrons_energy_fv","Neutron multiplicity vs muon energy (FV)",10,0,2000,20,0,20);
	TH2F *h_neutrons_energy_zoom = new TH2F("h_neutrons_energy_zoom","Neutron multiplicity vs muon energy",8,400,1200,20,0,20);
	TH2F *h_neutrons_energy_fv_zoom = new TH2F("h_neutrons_energy_fv_zoom","Neutron multiplicity vs muon energy (FV)",6,600,1200,20,0,20);
	TH2F *h_primneutrons_energy = new TH2F("h_primneutrons_energy","Primary Neutron multiplicity vs muon energy",10,0,2000,20,0,20);
	TH2F *h_primneutrons_energy_fv = new TH2F("h_primneutrons_energy_fv","Primary Neutron multiplicity vs muon energy (FV)",10,0,2000,20,0,20);
	TH2F *h_primneutrons_energy_zoom = new TH2F("h_primneutrons_energy_zoom","Primary Neutron multiplicity vs muon energy",8,400,1200,20,0,20);
	TH2F *h_primneutrons_energy_fv_zoom = new TH2F("h_primneutrons_energy_fv_zoom","Primary Neutron multiplicity vs muon energy (FV)",6,600,1200,20,0,20);
	TH2F *h_totalneutrons_energy = new TH2F("h_totalneutrons_energy","Total Neutron multiplicity vs muon energy",10,0,2000,20,0,20);
	TH2F *h_totalneutrons_energy_fv = new TH2F("h_totalneutrons_energy_fv","Total Neutron multiplicity vs muon energy (FV)",10,0,2000,20,0,20);
	TH2F *h_totalneutrons_energy_zoom = new TH2F("h_totalneutrons_energy_zoom","Total Neutron multiplicity vs muon energy",8,400,1200,20,0,20);
	TH2F *h_totalneutrons_energy_fv_zoom = new TH2F("h_totalneutrons_energy_fv_zoom","Total Neutron multiplicity vs muon energy (FV)",6,600,1200,20,0,20);
	TH2F *h_pmtvolneutrons_energy = new TH2F("h_pmtvolneutrons_energy","PMTVol Neutron multiplicity vs muon energy",10,0,2000,20,0,20);
	TH2F *h_pmtvolneutrons_energy_fv = new TH2F("h_pmtvolneutrons_energy_fv","PMTVol Neutron multiplicity vs muon energy (FV)",10,0,2000,20,0,20);
	TH2F *h_pmtvolneutrons_energy_zoom = new TH2F("h_pmtvolneutrons_energy_zoom","PMTVol Neutron multiplicity vs muon energy",8,400,1200,20,0,20);
	TH2F *h_pmtvolneutrons_energy_fv_zoom = new TH2F("h_pmtvolneutrons_energy_fv_zoom","PMTVol Neutron multiplicity vs muon energy (FV)",6,600,1200,20,0,20);
	TH2F *h_neutrons_costheta = new TH2F("h_neutrons_costheta","Neutron multiplicity vs muon angle",6,0.7,1.0,20,0,20);
        TH2F *h_neutrons_costheta_fv = new TH2F("h_neutrons_costheta_fv","Neutron multiplicity vs muon angle (FV)",6,0.7,1.0,20,0,20);
	TH2F *h_primneutrons_costheta = new TH2F("h_primneutrons_costheta","Primary Neutron multiplicity vs muon angle",6,0.7,1.0,20,0,20);
        TH2F *h_primneutrons_costheta_fv = new TH2F("h_primvolneutrons_costheta_fv","Primary Neutron multiplicity vs muon angle (FV)",6,0.7,1.0,20,0,20);
	TH2F *h_totalneutrons_costheta = new TH2F("h_totalneutrons_costheta","Total Neutron multiplicity vs muon angle",6,0.7,1.0,20,0,20);
        TH2F *h_totalneutrons_costheta_fv = new TH2F("h_totalneutrons_costheta_fv","Total Neutron multiplicity vs muon angle (FV)",6,0.7,1.0,20,0,20);
	TH2F *h_pmtvolneutrons_costheta = new TH2F("h_pmtvolneutrons_costheta","PMTVol Neutron multiplicity vs muon angle",6,0.7,1.0,20,0,20);
        TH2F *h_pmtvolneutrons_costheta_fv = new TH2F("h_pmtvolneutrons_costheta_fv","PMTVol Neutron multiplicity vs muon angle (FV)",6,0.7,1.0,20,0,20);

	TH1F *h_muon_energy = new TH1F("h_muon_energy","Muon energy distribution",100,0,2000);
	TH1F *h_muon_energy_fv = new TH1F("h_muon_energy_fv","Muon energy distribution (FV)",100,0,2000);
	TH2F *h_muon_vtx_yz = new TH2F("h_muon_vtx_yz","Muon vertex (tank) Y-Z",200,-3,3,200,-3,3);
	TH2F *h_muon_vtx_xz = new TH2F("h_muon_vtx_xz","Muon vertex (tank) X-Z",200,-2.5,2.5,200,-2.5,2.5);
 	TH1F *h_muon_costheta = new TH1F("h_muon_costheta","Muon cos(#theta) distribution",100,-1,1);
        TH1F *h_muon_costheta_fv = new TH1F("h_muon_costheta_fv","Muon cos(#theta) distribution (FV)",100,-1,1);
        TH1F *h_muon_vtx_x = new TH1F("h_muon_vtx_x","Muon vertex (x)",200,-3,3);
        TH1F *h_muon_vtx_y = new TH1F("h_muon_vtx_y","Muon vertex (y)",200,-3,3);
        TH1F *h_muon_vtx_z = new TH1F("h_muon_vtx_z","Muon vertex (z)",200,-3,3);

	int nPrimNeut;
	int nPrimProt;
	int nCaptures;
	int nCapturesPMTVol;
	double Emu;
	double Enu;
	double Q2;
	double vtxX;
	double vtxY;
	double vtxZ;
	std::vector<double> *clusterCB = new std::vector<double>;
	std::vector<double> *clusterTime = new std::vector<double>;
	std::vector<double> *clusterPE = new std::vector<double>;
	std::vector<double> *nCandCB = new std::vector<double>;
	std::vector<double> *nCandTime = new std::vector<double>;
	std::vector<double> *nCandPE = new std::vector<double>;
	double mrdEnergyLoss;
	std::vector<double> *neutVtxX = new std::vector<double>;
	std::vector<double> *neutVtxY = new std::vector<double>;
	std::vector<double> *neutVtxZ = new std::vector<double>;
	std::vector<double> *neutCapNucl = new std::vector<double>;
	std::vector<double> *neutCapTime = new std::vector<double>;
	std::vector<double> *neutCapETotal = new std::vector<double>;
	std::vector<double> *neutCapNGamma = new std::vector<double>;
	std::vector<int> *neutCapPrimary = new std::vector<int>;
	int isCC;
	int isQEL;
	int isDIS;
	int isRES;
	int isCOH;
	int isMEC;
	int tankMRDCoinc;
	int multiRing;
	int clusters;
	int nCandidates;
	double RecoVtxX;
	double RecoVtxY;
	double RecoVtxZ;
	int FV;
	double cosTheta;
	std::vector<int> *primaryPdgs = new std::vector<int>;

	TTree *tout = new TTree("neutron_tree","Neutron tree");
	tout->Branch("TankMRDCoinc",&tankMRDCoinc);
	tout->Branch("nPrimNeut",&nPrimNeut);
	tout->Branch("nPrimProt",&nPrimProt);
	tout->Branch("nCaptures",&nCaptures);
	tout->Branch("nCapturesPMTVol",&nCapturesPMTVol);
	tout->Branch("nCandidates",&nCandidates);
	tout->Branch("Emu",&Emu);
	tout->Branch("Enu",&Enu);
	tout->Branch("Q2",&Q2);
	tout->Branch("vtxX",&vtxX);
	tout->Branch("vtxY",&vtxY);
	tout->Branch("vtxZ",&vtxZ);
	tout->Branch("RecoVtxX",&RecoVtxX);
	tout->Branch("RecoVtxY",&RecoVtxY);
	tout->Branch("RecoVtxZ",&RecoVtxZ);
        tout->Branch("FV",&FV);
        tout->Branch("cosTheta",&cosTheta);
	tout->Branch("Clusters",&clusters);
	tout->Branch("ClusterCB",&clusterCB);
	tout->Branch("ClusterTime",&clusterTime);
	tout->Branch("ClusterPE",&clusterPE);
	tout->Branch("nCandCB",&nCandCB);
	tout->Branch("nCandTime",&nCandTime);
	tout->Branch("nCandPE",&nCandPE);
	tout->Branch("mrdEnergyLoss",&mrdEnergyLoss);	
	tout->Branch("neutVtxX",&neutVtxX);
	tout->Branch("neutVtxY",&neutVtxY);
	tout->Branch("neutVtxZ",&neutVtxZ);
	tout->Branch("neutCapNucl",&neutCapNucl);
	tout->Branch("neutCapTime",&neutCapTime);
	tout->Branch("neutCapETotal",&neutCapETotal);
	tout->Branch("neutCapNGamma",&neutCapNGamma);
	tout->Branch("neutCapPrimary",&neutCapPrimary);
	tout->Branch("isCC",&isCC);
	tout->Branch("isQEL",&isQEL);
	tout->Branch("isDIS",&isDIS);
	tout->Branch("isRES",&isRES);
	tout->Branch("isMEC",&isMEC);
	tout->Branch("isCOH",&isCOH);
	tout->Branch("multiRing",&multiRing);
	tout->Branch("primaryPdgs",&primaryPdgs);

	int i_pmt=0;
	for (int i_trig=0; i_trig < nentries_trigger; i_trig++){
		std::cout <<"Load trig entry "<<i_trig<<"/"<<nentries_trigger<<std::endl;
		tTrigger->GetEntry(i_trig);
		if (!noveto) continue;
		if (!tankmrdcoinc) continue;
		bool found_coinc=false;
		bool is_extended = false;
		int n_neutrons=0;
		double max_pe = 0;
		std::vector<double> n_times;
		Emu = 0.;
		Enu = 0.;
		Q2 = 0.;
		nPrimNeut = 0;
		nPrimProt = 0;
		nCaptures = 0;
		nCapturesPMTVol = 0;
		mrdEnergyLoss = 0;
		clusterCB->clear();
		clusterTime->clear();
		clusterPE->clear();
		nCandCB->clear();
		nCandTime->clear();
		nCandPE->clear();
		primaryPdgs->clear();
		neutVtxX->clear();
		neutVtxY->clear();
		neutVtxZ->clear();
		neutCapNucl->clear();
		neutCapETotal->clear();
		neutCapNGamma->clear();
		isCC = 0;
		isQEL = 0;
		isDIS = 0;
		isRES = 0;
		isMEC = 0;
		isCOH = 0;
		multiRing = 0;
		clusters = 0;
		nCandidates = 0;
		RecoVtxX = -999;
                RecoVtxY = -999;
                RecoVtxZ = -999;
                cosTheta = -999;
                FV = 0;

		bool check_above=false;
                bool check_below = false;
                bool run_switch=false;

		while (i_pmt < nentries_pmt){
			tPMT->GetEntry(i_pmt);
			if (verbose){
				std::cout <<"i_pmt: "<<i_pmt<<"run (PMT): "<<run_pmt<<", run (Trigger): "<<run_trigger<<std::endl;
				std::cout <<"i_pmt: "<<i_pmt<<"event (PMT): "<<ev_pmt<<", event (Trigger): "<<ev_trigger<<std::endl;
			}
			if (run_trigger == run_pmt){
				if (ev_trigger == ev_pmt){
					if (cluster_time < 2000 && cluster_pe > max_pe) max_pe = cluster_pe;
					if (cluster_time > 10000 && extended > 0 && beam_ok == 1 && cluster_pe < 120 && cluster_cb < 0.4){
					if (cluster_cb <= (1.-cluster_pe/150.)*0.5) {

						n_times.push_back(cluster_time);
						n_neutrons++;
						nCandCB->push_back(cluster_cb);
						nCandTime->push_back(cluster_time);
						nCandPE->push_back(cluster_pe);
						}
					}
					i_pmt++;
					found_coinc = true;
					if (extended) is_extended = true;
					clusters++;
					clusterCB->push_back(cluster_cb);
					clusterTime->push_back(cluster_time);
					clusterPE->push_back(cluster_pe);
					
				} else if (found_coinc) break;
				/*else if (ev_pmt < ev_trigger) i_pmt++;
				else if (ev_pmt > ev_trigger) {
					//i_pmt--;	//Data case
					break;		//MC case
				}*/
				else if (ev_pmt < ev_trigger) {i_pmt++; check_below=true;}
                                else if (ev_pmt > ev_trigger && i_pmt > 0 && !run_switch) {i_pmt--; check_above=true;}
                                if (run_switch && ev_pmt > ev_trigger) break;
                                if (i_pmt ==0) break;
                                if (check_below && check_above) break;
			}
			else if (run_pmt < run_trigger) {i_pmt++; run_switch = true;}
			else if (i_pmt > 0) break;
		}
		if (found_coinc && is_extended && num_mrd_tracks == 1){
			h_neutrons->Fill(n_neutrons);
			for (int i_t = 0; i_t < (int) n_times.size(); i_t++){
                                h_time_neutrons->Fill(n_times.at(i_t));
                        }
			if (mrdstop->at(0)){
				for (int i_t = 0; i_t < (int) n_times.size(); i_t++){
                                	h_time_neutrons_mrdstop->Fill(n_times.at(i_t));
                        	}
				h_neutrons_mrdstop->Fill(n_neutrons);

				//Energy reconstruction
				double muon_eloss = eloss->at(0);
				double mrdTrackLength = tracklength->at(0);
				double reco_muon_energy = muon_eloss + max_pe*0.08534;				
				//Direction information (from MRD)
				double startx = mrdstartx->at(0);
				double starty = mrdstarty->at(0);
				double startz = mrdstartz->at(0);
				double stopx = mrdstopx->at(0);
				double stopy = mrdstopy->at(0);
				double stopz = mrdstopz->at(0);
				double diffx = stopx-startx;
				double diffy = stopy-starty;
				double diffz = stopz-startz;
				double startz_c = startz-1.681;
				bool hit_tank = false;
				double a = pow(diffx,2)+pow(diffz,2);
				double b = -2*diffx*startx-2*diffz*startz_c;
				double c = pow(startx,2)+pow(startz_c,2)-1.0*1.0;
				double t1 = (-b+sqrt(b*b-4*a*c))/(2*a);
				double t2 = (-b-sqrt(b*b-4*a*c))/(2*a);
				double t = 0;
				if (t1 < 0) t = t2;
				else if (t2 < 0) t = t1;
				else t = (t1 < t2)? t1 : t2;
				double exitx = startx - t*diffx;
				double exity = starty - t*diffy;
				double exitz = startz - t*diffz;
				double dirx = diffx/(sqrt(diffx*diffx+diffy*diffy+diffz*diffz));
				double diry = diffy/(sqrt(diffx*diffx+diffy*diffy+diffz*diffz));
				double dirz = diffz/(sqrt(diffx*diffx+diffy*diffy+diffz*diffz));
				cosTheta = dirz;
				
				double dist_mrd = sqrt(pow(startx-stopx,2)+pow(starty-stopy,2)+pow(startz-stopz,2));
                                double dist_air = sqrt(pow(startx-exitx,2)+pow(starty-exity,2)+pow(startz-exitz,2));
                                double dist_air_x = sqrt(pow(startx-exitx,2));
                                double dist_air_y = sqrt(pow(starty-exity,2));
                                double dist_air_z = sqrt(pow(startz-exitz,2));

                                double ap = pow(diffx,2)+pow(diffz,2);
                                double bp = -2*diffx*startx-2*diffz*startz_c;
                                double cp = pow(startx,2)+pow(startz_c,2)-1.5*1.5;
                                double t1p = (-bp+sqrt(bp*bp-4*ap*cp))/(2*ap);
                                double t2p = (-bp-sqrt(bp*bp-4*ap*cp))/(2*ap);
                                std::cout <<"t1, t2: "<<t1p<<","<<t2p<<std::endl;
                                double tp = 0;
                                if (t1p < 0) tp = t2p;
                                else if (t2p < 0) tp = t1p;
                                else tp = (t1p < t2p)? t1p : t2p;
                                std::cout <<"t: "<<tp<<std::endl;

                                double exitxp = startx - tp*diffx;
                                double exityp = starty - tp*diffy;
                                double exitzp = startz - tp*diffz;

                                double diff_exit_x = exitxp - exitx;
                                double diff_exit_y = exityp - exity;
                                double diff_exit_z = exitzp - exitz;
                                double diff_pmtvol_tank = sqrt(pow(diff_exit_x,2)+pow(diff_exit_y,2)+pow(diff_exit_z,2));
                                std::cout <<"old energy: "<<reco_muon_energy<<std::endl;
                                if (!std::isnan(diff_pmtvol_tank)){
                                        reco_muon_energy += 200*(diff_pmtvol_tank);
                                }
                                std::cout <<"new energy: "<<reco_muon_energy<<std::endl;
				
				//Add offset (potentially from slight bias in MRD energy estimate)
                                reco_muon_energy += 87.3;

				h_neutrons_energy->Fill(reco_muon_energy,n_neutrons);
                                h_neutrons_energy_zoom->Fill(reco_muon_energy,n_neutrons);
                                h_muon_energy->Fill(reco_muon_energy);
                                h_neutrons_costheta->Fill(cosTheta,n_neutrons);
                                h_muon_costheta->Fill(cosTheta);

				//Vertex reconstruction
				double reco_vtxx = exitx - max_pe/9/2./100.*dirx;
                                double reco_vtxy = exity - max_pe/9/2./100.*diry;
                                double reco_vtxz = exitz - max_pe/9/2./100.*dirz;

                                h_muon_vtx_x->Fill(reco_vtxx);
                                h_muon_vtx_y->Fill(reco_vtxy);
                                h_muon_vtx_z->Fill(reco_vtxz);
                                h_muon_vtx_yz->Fill(reco_vtxz,reco_vtxy);
                                h_muon_vtx_xz->Fill(reco_vtxz,reco_vtxx);

				FV = 0;

                                if (sqrt(reco_vtxx*reco_vtxx+(reco_vtxz-1.681)*(reco_vtxz-1.681))<1.0 && fabs(reco_vtxy)<0.5 && ((reco_vtxz-1.681) < 0.)) {
                                        h_neutrons_energy_fv->Fill(reco_muon_energy,n_neutrons);
                                        h_neutrons_energy_fv_zoom->Fill(reco_muon_energy,n_neutrons);
                                        h_muon_energy_fv->Fill(reco_muon_energy);
                                        h_muon_costheta_fv->Fill(cosTheta);
                                        h_neutrons_mrdstop_fv->Fill(n_neutrons);
                                        h_neutrons_costheta_fv->Fill(cosTheta,n_neutrons);
                                        FV = 1;
                                }

				//MC version of vtx x/y/z
				double vtxx = true_vtxx/100.;
				double vtxy = true_vtxy/100.;
				double vtxz = true_vtxz/100.;
				vtxX = vtxx;
				vtxY = vtxy;
				vtxZ = vtxz;
				Emu = reco_muon_energy;
				RecoVtxX = reco_vtxx;
                                RecoVtxY = reco_vtxy;
                                RecoVtxZ = reco_vtxz;
                                mrdEnergyLoss = muon_eloss;
				tankMRDCoinc = tankmrdcoinc;
				nPrimNeut = true_neutrons;
				isMEC = true_mec;
				isDIS = true_dis;
				isRES = true_res;
				isQEL = true_qel;
				isCC = true_cc;
				Q2 = -true_q2;	//For some reason Q2 in Genie is saved with negative sign
				nPrimProt = true_protons;
				multiRing = trueMultiRing;
				nCaptures = true_neut_cap_vtxx->size();
				nCandidates = n_neutrons;
				nCapturesPMTVol = 0;
				for (int i_neut= 0; i_neut < (int) true_neut_cap_vtxx->size(); i_neut++){
					neutVtxX->push_back(true_neut_cap_vtxx->at(i_neut));
					neutVtxY->push_back(true_neut_cap_vtxy->at(i_neut));
					neutVtxZ->push_back(true_neut_cap_vtxz->at(i_neut));
					neutCapNucl->push_back(true_neut_cap_nucl->at(i_neut));
					neutCapTime->push_back(true_neut_cap_time->at(i_neut));
					neutCapETotal->push_back(true_neut_cap_e->at(i_neut));
					neutCapNGamma->push_back(true_neut_cap_gammas->at(i_neut));
					//PMT volume cut
					double zcorr = true_neut_cap_vtxz->at(i_neut)/100.-1.681;
					double ycorr = true_neut_cap_vtxy->at(i_neut)/100.+0.144;
					double xcorr = true_neut_cap_vtxx->at(i_neut)/100.;
					double rad = sqrt(xcorr*xcorr+zcorr*zcorr);
					if (rad < 1.0 && fabs(ycorr)<1.5) nCapturesPMTVol++;
				}
				h_primneutrons_energy->Fill(reco_muon_energy,nPrimNeut);
                                h_primneutrons_energy_zoom->Fill(reco_muon_energy,nPrimNeut);
				h_totalneutrons_energy->Fill(reco_muon_energy,nCaptures);
                                h_totalneutrons_energy_zoom->Fill(reco_muon_energy,nCaptures);
				h_pmtvolneutrons_energy->Fill(reco_muon_energy,nCapturesPMTVol);
                                h_pmtvolneutrons_energy_zoom->Fill(reco_muon_energy,nCapturesPMTVol);
                                h_primneutrons_costheta->Fill(cosTheta,nPrimNeut);
                                h_totalneutrons_costheta->Fill(cosTheta,nCaptures);
                                h_pmtvolneutrons_costheta->Fill(cosTheta,nCapturesPMTVol);
				if (FV){
					h_primneutrons_energy_fv->Fill(reco_muon_energy,nPrimNeut);
                                	h_primneutrons_energy_fv_zoom->Fill(reco_muon_energy,nPrimNeut);
					h_totalneutrons_energy_fv->Fill(reco_muon_energy,nCaptures);
                                	h_totalneutrons_energy_fv_zoom->Fill(reco_muon_energy,nCaptures);
					h_pmtvolneutrons_energy_fv->Fill(reco_muon_energy,nCapturesPMTVol);
                                	h_pmtvolneutrons_energy_fv_zoom->Fill(reco_muon_energy,nCapturesPMTVol);
                                	h_primneutrons_costheta_fv->Fill(cosTheta,nPrimNeut);
                                	h_totalneutrons_costheta_fv->Fill(cosTheta,nCaptures);
                                	h_pmtvolneutrons_costheta_fv->Fill(cosTheta,nCapturesPMTVol);

				}
				for (int i_part = 0; i_part < (int) true_primary_pdgs->size(); i_part++){
					primaryPdgs->push_back(true_primary_pdgs->at(i_part));
				}
				tout->Fill();

			}
		}
	}

        h_neutrons->SetLineWidth(2);
        h_neutrons->SetStats(0);
        h_neutrons->GetXaxis()->SetTitle("Number of neutrons");
        h_neutrons->GetYaxis()->SetTitle("#");

        h_neutrons_mrdstop->SetLineWidth(2);
        h_neutrons_mrdstop->SetLineColor(kBlack);
        h_neutrons_mrdstop->SetStats(0);
        h_neutrons_mrdstop->GetXaxis()->SetTitle("Number of neutrons");
        h_neutrons_mrdstop->GetYaxis()->SetTitle("#");

        h_neutrons_mrdstop_fv->SetLineWidth(2);
        h_neutrons_mrdstop_fv->SetStats(0);
        h_neutrons_mrdstop_fv->SetLineColor(kRed);
        h_neutrons_mrdstop_fv->GetXaxis()->SetTitle("Number of neutrons");
        h_neutrons_mrdstop_fv->GetYaxis()->SetTitle("#");

        h_neutrons_energy->SetStats(0);
        h_neutrons_energy->GetXaxis()->SetTitle("E_{#mu} [MeV]");
        h_neutrons_energy->GetYaxis()->SetTitle("Number of neutrons");

        h_neutrons_energy_zoom->SetStats(0);
        h_neutrons_energy_zoom->GetXaxis()->SetTitle("E_{#mu} [MeV]");
        h_neutrons_energy_zoom->GetYaxis()->SetTitle("Number of neutrons");

        h_neutrons_energy_fv->SetStats(0);
        h_neutrons_energy_fv->GetXaxis()->SetTitle("E_{#mu} [MeV]");
        h_neutrons_energy_fv->GetYaxis()->SetTitle("Number of neutrons");

        h_neutrons_energy_fv_zoom->SetStats(0);
        h_neutrons_energy_fv_zoom->GetXaxis()->SetTitle("E_{#mu} [MeV]");
        h_neutrons_energy_fv_zoom->GetYaxis()->SetTitle("Number of neutrons");

        h_neutrons_costheta->SetStats(0);
        h_neutrons_costheta->GetXaxis()->SetTitle("cos(#theta)");
        h_neutrons_costheta->GetYaxis()->SetTitle("Number of neutrons");

        h_neutrons_costheta_fv->SetStats(0);
        h_neutrons_costheta_fv->GetXaxis()->SetTitle("cos(#theta)");
        h_neutrons_costheta_fv->GetYaxis()->SetTitle("Number of neutrons");


        h_muon_energy->GetXaxis()->SetTitle("E_{#mu} [MeV]");
        h_muon_energy->GetYaxis()->SetTitle("#");
        h_muon_energy->SetStats(0);
        h_muon_energy->SetLineWidth(2);
        h_muon_energy_fv->SetStats(0);
        h_muon_energy_fv->SetLineColor(kRed);
        h_muon_energy_fv->SetLineWidth(2);


        h_muon_vtx_yz->GetXaxis()->SetTitle("Vertex Z [m]");
        h_muon_vtx_yz->GetYaxis()->SetTitle("Vertex Y [m]");
        h_muon_vtx_yz->SetStats(0);
        h_muon_vtx_xz->GetXaxis()->SetTitle("Vertex Z [m]");
        h_muon_vtx_xz->GetYaxis()->SetTitle("Vertex X [m]");
        h_muon_vtx_xz->SetStats(0);

        h_muon_vtx_x->GetXaxis()->SetTitle("Vertex X [m]");
        h_muon_vtx_x->GetYaxis()->SetTitle("#");
        h_muon_vtx_y->GetXaxis()->SetTitle("Vertex Y [m]");
        h_muon_vtx_y->GetYaxis()->SetTitle("#");
        h_muon_vtx_z->GetXaxis()->SetTitle("Vertex Z [m]");
        h_muon_vtx_z->GetYaxis()->SetTitle("#");

	fout->Write();

	fout->Close();

}
