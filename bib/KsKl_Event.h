#ifndef ClassmuEvent
#define ClassmuEvent

class KsKl_Event
{

public:
        // Particle and Lorentz vectors
        Particle Electron;
        Particle Pi_plus;
        Particle Pi_minus;
        Particle Proton;
        Particle *Photons;
        TLorentzVector vBeam;
        TLorentzVector vMissing;
        TLorentzVector vRestProton;


        // Number of reconstructed particle of each type
        int np = 4; // total number of particles
        int rec_e = 0;
        int rec_pi_p = 0;
        int rec_pi_m = 0;
        int rec_p = 0;

        // Kinematic variables
        float Phi_Mass;
	float Ks_Mass;
        float Kl_Mass;
	float Q2;
        float t;
	float W;

        // run number and trigger bit
        float run;
        int trigger_bit;

        KsKl_Event(int nb_part)
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
                np = nb_part;
                Photons = new Particle[np];
        }

        KsKl_Event()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        void Set_nb_part(int input_np)
        {
                np = input_np;
                Photons = new Particle[np];
        }

        void Set_Particles(hipo::bank PART)
        {

                for (int i = 0; i < np; i++)
                {

                        int pid = PART.getInt("pid", i);
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float beta = PART.getFloat("beta", i);
                        int status = abs(PART.getInt("status", i));
                        int charge = PART.getInt("charge", i);
                        float chi2 = PART.getFloat("chi2pid", i);
                        float vx = PART.getFloat("vx", i);
                        float vy = PART.getFloat("vy", i);
                        float vz = PART.getFloat("vz", i);
                        float vt = PART.getFloat("vt", i);

                        if (pid == -211)
                        {
                                Pi_minus.Vector.SetXYZM(px, py, pz, mPi);
                                Pi_minus.index = i;
                                Pi_minus.pid = -211;
                                Pi_minus.beta = beta;
                                Pi_minus.status = status;
                                Pi_minus.chi2 = chi2;
                                Pi_minus.vertex.x = vx;
                                Pi_minus.vertex.y = vy;
                                Pi_minus.vertex.z = vz;
                                Pi_minus.vt = vt;
                                rec_pi_m++;
                        }

                        if (pid == 211)
                        {
                                Pi_plus.Vector.SetXYZM(px, py, pz, mPi);
                                Pi_plus.index = i;
                                Pi_plus.pid = 211;
                                Pi_plus.beta = beta;
                                Pi_plus.status = status;
                                Pi_plus.chi2 = chi2;
                                Pi_plus.vertex.x = vx;
                                Pi_plus.vertex.y = vy;
                                Pi_plus.vertex.z = vz;
                                Pi_plus.vt = vt;
                                rec_pi_p++;
                        }

                        if (pid == 11)
                        {
                                Electron.Vector.SetXYZM(px, py, pz, me);
                                Electron.index = i;
                                Electron.pid = 11;
                                Electron.beta = beta;
                                Electron.status = status;
                                Electron.chi2 = chi2;
                                Electron.vertex.x = vx;
                                Electron.vertex.y = vy;
                                Electron.vertex.z = vz;
                                Electron.vt = vt;
                                rec_e++;
                        }

                        if (pid == 2212)
                        {
                                Proton.Vector.SetXYZM(px, py, pz, mp);
                                Proton.index = i;
                                Proton.pid = 2212;
                                Proton.beta = beta;
                                Proton.status = status;
                                Proton.chi2 = chi2;
                                Proton.vertex.x = vx;
                                Proton.vertex.y = vy;
                                Proton.vertex.z = vz;
                                Proton.vt = vt;
                                rec_p++;
                        }
                }
        }

        bool pass_topology_cut()
        {
                return (rec_pi_p == 1 && rec_pi_m == 1 && rec_e == 1 && rec_p == 1) && (abs(Ks_Mass-mK)<0.5) && (abs(Kl_Mass-mK)<0.5);
        }

        void Associate_detector_resp(hipo::bank CHE, hipo::bank SCIN, hipo::bank CALO)
        {
                vector<Particle> Particles = {Electron, Pi_plus, Pi_minus, Proton};

                CalorimeterResp Calo;
                CheResp Che;
                ScinResp Scin;
                for (int i = 0; i < 4; i++)
                {

                        for (int c = 0; c < CHE.getRows(); c++)
                        {
                                int Chepindex = CHE.getInt("pindex", c);
                                int Chedetector = CHE.getInt("detector", c);
                                int Chesector = CHE.getInt("sector", c);
                                float Chenphe = CHE.getFloat("nphe", c);
                                float Chetime = CHE.getFloat("time", c);
                                float Chechi2 = CHE.getFloat("chi2", c);
                                float Chex = CHE.getFloat("x", c);
                                float Chey = CHE.getFloat("y", c);
                                float Chez = CHE.getFloat("z", c);

                                if (Chepindex == (Particles[i].index))
                                {
                                        Che.detector = Chedetector;
                                        Che.pindex = Chepindex;
                                        Che.sector = Chesector;
                                        Che.nphe = Chenphe;
                                        Che.time = Chetime;
                                        Che.chi2 = Chechi2;
                                        Che.x = Chex;
                                        Che.y = Chey;
                                        Che.z = Chez;
                                        Particles[i].Cherenkov.push_back(Che);
                                }
                        }

                        for (int c = 0; c < SCIN.getRows(); c++)
                        {
                                int Scindetector = SCIN.getInt("detector", c);
                                int Scinpindex = SCIN.getInt("pindex", c);
                                float Scintime = SCIN.getFloat("time", c);
                                float Scinpath = SCIN.getFloat("path", c);
                                float Scinenergy = SCIN.getFloat("energy", c);
                                int Scinsector = SCIN.getInt("sector", c);

                                if (Scinpindex == (Particles[i].index))
                                {
                                        Scin.detector = Scindetector;
                                        Scin.pindex = Scinpindex;
                                        Scin.t = Scintime;
                                        Scin.path = Scinpath;
                                        Scin.energy = Scinenergy;
                                        Scin.sector = Scinsector;

                                        if (Particles[i].Scintillator.energy < Scinenergy)
                                        {
                                                Particles[i].Scintillator = Scin;
                                        };
                                }
                        }

                        for (int c = 0; c < CALO.getRows(); c++)
                        {
                                int Calopindex = CALO.getInt("pindex", c);
                                int Calosector = CALO.getInt("sector", c);
                                int Calolayer = CALO.getInt("layer", c);
                                int Calodetector = CALO.getInt("detector", c);
                                float Caloenergy = CALO.getFloat("energy", c);
                                float Calox = CALO.getFloat("x", c);
                                float Caloy = CALO.getFloat("y", c);
                                float Caloz = CALO.getFloat("z", c);
                                float Calohx = CALO.getFloat("hx", c);
                                float Calohy = CALO.getFloat("hy", c);
                                float Calohz = CALO.getFloat("hz", c);
                                float Calou = CALO.getFloat("lu", c);
                                float Calov = CALO.getFloat("lv", c);
                                float Calow = CALO.getFloat("lw", c);
                                float Calodu = CALO.getFloat("du", c);
                                float Calodv = CALO.getFloat("dv", c);
                                float Calodw = CALO.getFloat("dw", c);
                                float Calom2u = CALO.getFloat("m2u", c);
                                float Calom2v = CALO.getFloat("m2v", c);
                                float Calom2w = CALO.getFloat("m2w", c);
                                float Calom3u = CALO.getFloat("m3u", c);
                                float Calom3v = CALO.getFloat("m3v", c);
                                float Calom3w = CALO.getFloat("m3w", c);

                                if (Calopindex == (Particles[i].index))
                                {
                                        Calo.detector = Calodetector;
                                        Calo.pindex = Calopindex;
                                        Calo.sector = Calosector;
                                        Calo.layer = Calolayer;
                                        Calo.energy = Caloenergy;
                                        Calo.x = Calox;
                                        Calo.y = Caloy;
                                        Calo.z = Caloz;
                                        Calo.hx = Calohx;
                                        Calo.hy = Calohy;
                                        Calo.hz = Calohz;
                                        Calo.u = Calou;
                                        Calo.v = Calov;
                                        Calo.w = Calow;
                                        Calo.du = Calodu;
                                        Calo.dv = Calodv;
                                        Calo.dw = Calodw;
                                        Calo.m2u = Calom2u;
                                        Calo.m2v = Calom2v;
                                        Calo.m2w = Calom2w;
                                        Calo.m3u = Calom3u;
                                        Calo.m3v = Calom3v;
                                        Calo.m3w = Calom3w;
                                        Particles[i].Calorimeter.push_back(Calo);	
                                }
                        }
                }

                Electron = Particles[0]; 
                Pi_plus = Particles[1]; 
                Pi_minus = Particles[2]; 
                Proton = Particles[3];
        }

        void Get_Kinematics()
        {
                // Kinematic variables
                vMissing = Pi_minus.Vector + Pi_plus.Vector + Electron.Vector + Proton.Vector - vRestProton - vBeam;

                Phi_Mass = (Electron.Vector + Proton.Vector - vRestProton - vBeam).M();
	        Ks_Mass = (Pi_minus.Vector + Pi_plus.Vector).M();
                Kl_Mass = vMissing.M();
	        Q2 = (Electron.Vector - vBeam).M2();
                t = (Proton.Vector - vRestProton).M2();
	        W = (Electron.Vector + Proton.Vector - vRestProton - vBeam).M();
        }


        void Set_Run_Number(int input_run)
        {
                run = (float)(input_run);
        }


};

#endif
