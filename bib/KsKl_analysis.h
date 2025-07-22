#ifndef muCLAS12class
#define muCLAS12class

const float me = 0.000511;
const float mp = 0.938;
const float mMu = 0.105658;
const float mPi = 0.139570;
const float mK = 0.497611;
float ebeam = 10.6;

const float Pi = 3.14159269;

////////////////////Configuration booleans
bool all_Gen_vector = true;
////////////////////

class Track
{
public:
        int detector = 0;
        int pindex = 0;
        int sector = 0;
        int layer = 0;
        float x = 0.;
        float y = 0.;
        float z = 0.;
        float chi2 = 0.0;
};

class Traj
{

public:
        int detector = 0;
        int pindex = 0;
        int layer = 0;
        float x = -10000.;
        float y = -10000.;
        float z = -10000.;
};

class CalorimeterResp
{

public:
        int detector, pindex, sector, layer;
        float energy, x, y, z, hx, hy, hz, u, v, w;
        float du, dv, dw;
        float m2u, m2v, m2w;
        float m3u, m3v, m3w;
};

class CheResp
{

public:
        int detector, pindex, sector;
        float nphe, x, y, z, time, chi2;
};

class ScinResp
{

public:
        int detector, pindex, sector;
        float t, path;
        float energy = 0.0;
};

class Vertex
{
public:
        float x;
        float y;
        float z;

        Vertex()
        {
                x = 0.0;
                y = 0.0;
                z = 0.0;
        }

        Vertex(float x_in, float y_in, float z_in)
        {
                x = x_in;
                y = y_in;
                z = z_in;
        }
};

class Particle
{

public:
        int index = -1;

        int pid;
        float beta;
        float vt;
        int status;
        float chi2;
        bool passEC = true;
        bool passR1 = true;
        bool passR2 = true;
        bool passR3 = true;
        Vertex vertex;
        TLorentzVector Vector;
        ScinResp Scintillator;

        vector<CalorimeterResp> Calorimeter;
        vector<CheResp> Cherenkov;
        vector<Traj> Trajs;

        Vertex cluster_local_PCAL;
        Vertex h_cluster_local_PCAL;
        Vertex cluster_local_ECIN;
        Vertex h_cluster_local_ECIN;

        float nphe(int det)
        {
                float np = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                np += Cherenkov[i].nphe;
                }
                return np;
        }

        float TimeChe(int det)
        {
                float time = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                time = Cherenkov[i].time;
                }
                return time;
        }

        float SectorChe(int det)
        {
                int sector = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                        {
                                sector = Cherenkov[i].sector;
                                break;
                        }
                }
                return sector;
        }

        float SectorCalo(int det, int layer)
        {
                int sector = 0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer)
                        {
                                if (det == Calorimeter[i].detector)
                                {
                                        sector = Calorimeter[i].sector;
                                        break;
                                }
                        }
                }
                return sector;
        }

        float Chi2Che(int det)
        {
                float chi2 = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                chi2 = Cherenkov[i].chi2;
                }
                return chi2;
        }

        float XChe(int det)
        {
                float X = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                X = Cherenkov[i].x;
                }
                return X;
        }

        float YChe(int det)
        {
                float Y = 0;
                for (int i = 0; i < Cherenkov.size(); i++)
                {
                        if (det == Cherenkov[i].detector)
                                Y = Cherenkov[i].y;
                }
                return Y;
        }

        float X_CALO(int layer)
        {
                float X = 0;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && Calorimeter[i].detector == ECAL)
                                X = Calorimeter[i].x;
                        break;
                }
                return X;
        }

        float Y_CALO(int layer)
        {
                float Y = 0;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && Calorimeter[i].detector == ECAL)
                                Y = Calorimeter[i].y;
                        break;
                }
                return Y;
        }

        float Z_CALO(int layer)
        {
                float Z = 0;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Z = Calorimeter[i].z;
                        break;
                }
                return Z;
        }

        Vertex Cluster_CALO(int layer)
        {
                Vertex cluster_calo;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        //cout<<"requiered: "<<layer<<" "<<ECAL<<" size "<<Calorimeter.size()<<endl;
                        //cout<<"list: "<<Calorimeter[i].layer<<" "<<Calorimeter[i].detector<<endl;
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                        {
                                cluster_calo.x = Calorimeter[i].x;
                                cluster_calo.y = Calorimeter[i].y;
                                cluster_calo.z = Calorimeter[i].z;
                                break;
                        }
                        
                }
                return cluster_calo;
        }

        Vertex Matched_Cluster_CALO(int layer)
        {
                Vertex cluster_calo;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                        {
                                cluster_calo.x = Calorimeter[i].hx;
                                cluster_calo.y = Calorimeter[i].hy;
                                cluster_calo.z = Calorimeter[i].hz;
                                break;
                        }
                        
                }
                return cluster_calo;
        }

        float N_strip(int layer)
        {
                float N_strip = 0;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                        {
                                N_strip = Calorimeter[i].du + Calorimeter[i].dv + Calorimeter[i].dw;
                                break;
                        }
                }
                return N_strip;
        }

        float U_CALO(int layer)
        {
                float Y = 0;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].u;
                        break;
                }
                return Y;
        }

        float V_CALO(int layer)
        {
                float Y = 0;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].v;
                        break;
                }
                return Y;
        }

        float W_CALO(int layer)
        {
                float Y = 0;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].w;
                        break;
                }
                return Y;
        }

        float SECTOR_CALO(int layer)
        {
                float Y = 0;
                int ECAL = 7;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer && ECAL == Calorimeter[i].detector)
                                Y = Calorimeter[i].sector;
                        break;
                }
                return Y;
        }

        float Energy(int det, int layer)
        {
                float en = 0.;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer)
                        {
                                if (det == Calorimeter[i].detector)
                                        en = en + Calorimeter[i].energy;
                        }
                }
                return en;
        }

        float M2_ECAL(int layer)
        {
                float m2 = -1.0;
                for (int i = 0; i < Calorimeter.size(); i++)
                {
                        if (layer == Calorimeter[i].layer)
                        {
                                m2 = (Calorimeter[i].m2u + Calorimeter[i].m2v + Calorimeter[i].m2w) / 3.;
                        }
                }
                return m2;
        }

        void Associate_DC_traj_to_Particle(hipo::bank TRAJ)
        {

                for (int t = 0; t < TRAJ.getRows(); t++)
                {
                        Traj new_Traj;
                        int TRAJ_detector = TRAJ.getInt("detector", t);
                        int TRAJ_pindex = TRAJ.getInt("pindex", t);
                        int TRAJ_layer = TRAJ.getInt("layer", t);
                        float TRAJ_X = TRAJ.getFloat("x", t);
                        float TRAJ_Y = TRAJ.getFloat("y", t);
                        float TRAJ_Z = TRAJ.getFloat("z", t);

                        if (TRAJ_pindex == index && TRAJ_detector == DC)
                        {
                                new_Traj.detector = TRAJ_detector;
                                new_Traj.pindex = TRAJ_pindex;
                                new_Traj.layer = TRAJ_layer;
                                new_Traj.x = TRAJ_X;
                                new_Traj.y = TRAJ_Y;
                                new_Traj.z = TRAJ_Z;
                                Trajs.push_back(new_Traj);
                        }
                }
        }

        void Get_local_cluster_CALO(int calo_ID)
        {

                int sector = this->SECTOR_CALO(PCAL);

                Vertex cluster_CALO = Cluster_CALO(calo_ID);
                Vertex h_cluster_CALO = Matched_Cluster_CALO(calo_ID);

                if (calo_ID == PCAL)
                {
                        cluster_local_PCAL = Rotate_to_local(cluster_CALO, sector);
                        h_cluster_local_PCAL = Rotate_to_local(h_cluster_CALO, sector);
                }
                if (calo_ID == ECIN)
                {
                        //cout<<"Here"<<endl;
                        //cout<<"cluster_CALO "<<cluster_CALO.x<<endl;
                        cluster_local_ECIN = Rotate_to_local(cluster_CALO, sector);
                        h_cluster_local_ECIN = Rotate_to_local(h_cluster_CALO, sector);
                }
        }

        Vertex Rotate_to_local(Vertex cluster, int sector)
        {
                TVector3 vector_point(cluster.x, cluster.y, cluster.z);
                vector_point.RotateZ(-TMath::Pi() / 3.0 * (sector - 1));
                vector_point.RotateY(-TMath::Pi() / 180.0 * 25);

                Vertex rotated_cluster(vector_point.x(), vector_point.y(), vector_point.z());
                return rotated_cluster;
        }
};

class ThetaPhi
{
public:
        float theta;
        float phi;
};

class PairValue
{
public:
        float mean;
        float sigma;
};

#endif
