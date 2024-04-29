#include "functions_2.h"


double uni_rand(){
    return (rand()/(double) RAND_MAX);
}

void save_data(std::string outFileName, std::vector<std::vector<double>> data){
    std::fstream outFile;
    outFile.open(outFileName, std::fstream::out);
    if(outFile.is_open()){
        for(int i=0; i < data.size(); i++){
            for(int j=0; j < data[i].size(); j++){
                outFile << data[i][j] << " ";
            }
            outFile << "\n";
        }
    }
    outFile.close();
}

std::vector<std::vector<double>> load_atoms(std::string inFileName){
    
    std::vector<std::vector<double>> xyz;

    std::ifstream inFile;
    inFile.open(inFileName.c_str());

    if(inFile.is_open()){
        while(!inFile.eof()){
            std::vector<double> row;
            row.resize(3,0);
            inFile >> row[0] >> row[1] >> row[2];
            xyz.push_back(row);
            // std::cout<<xyz[i][0]<<" "<<xyz[i][1]<<" "<<xyz[i][2]<<"\n";
        }
        inFile.close();
    }

    return xyz;
}

double cut_function(double r, double R1, double R2){
    double f_cut;
    if(r <= R1){
        f_cut = 1;
    } else if(r > R1 && r <= R2){
        double frac = (double) (r - R1)/(R2 - R1);
        f_cut = 1.0/2.0*(1 + cos( M_PI * frac ));
    } else{
        f_cut = 0;
    }
    return f_cut;
}

double B_val(std::vector<std::vector<double>>& atoms, int i, int j, double r_ij[3], double r_ij_len, double R1, double R2, double a0, double c0, double d0, double delta, bool correction){
    double B_ij = 0;
    double B_ji = 0;
    double ksi_ij = 0;
    double ksi_ji = 0;

    double r_ik[3];

    double atoms_num = atoms.size();

    for(int k=0; k<atoms_num; k++){
        if(k != i && k != j){
            double r_ik[3] = {atoms[k][3] - atoms[i][3], atoms[k][4] - atoms[i][4], atoms[k][5] - atoms[i][5]};
            double r_ik_len = sqrt( pow(r_ik[0],2) + pow(r_ik[1],2) + pow(r_ik[2],2) );
            
            double f_cut_ik = cut_function(r_ik_len, R1, R2);
            if(f_cut_ik != 0){
                double s_prod_ijk = (double) r_ij[0]*r_ik[0] + r_ij[1]*r_ik[1] + r_ij[2]*r_ik[2];
                double len_ijk = (double) r_ij_len * r_ik_len;
                double kosinus = (s_prod_ijk/len_ijk);
                if(correction == true && kosinus > 0){
                    ksi_ji += 10;
                }else{
                    double g_function = a0*( 1 + pow(c0,2)/pow(d0,2) - pow(c0,2)/( pow(d0,2) + pow(1 + kosinus, 2) ) );
                    ksi_ij += f_cut_ik*g_function;
                }
            }

            double r_jk[3] = {atoms[k][3] - atoms[j][3], atoms[k][4] - atoms[j][4], atoms[k][5] - atoms[j][5]};
            double r_jk_len = sqrt( pow(r_jk[0],2) + pow(r_jk[1],2) + pow(r_jk[2],2) );

            double r_ji[3] = {r_ij[0]*(-1.), r_ij[1]*(-1.), r_ij[2]*(-1.)};

            double f_cut_jk = cut_function(r_jk_len, R1, R2);
            if(f_cut_jk != 0){
                double s_prod_ijk = (double) r_ji[0]*r_jk[0] + r_ji[1]*r_jk[1] + r_ji[2]*r_jk[2];
                double len_ijk = (double) r_ij_len * r_jk_len;
                double kosinus = (s_prod_ijk/len_ijk);
                if(correction == true && kosinus > 0){
                    ksi_ji += 10;
                }else{
                    double g_function = a0*( 1 + pow(c0,2)/pow(d0,2) - pow(c0,2)/( pow(d0,2) + pow(1 + kosinus, 2) ) );
                    ksi_ji += f_cut_jk*g_function;
                }
            }

        }
    }

    B_ij = (double) pow((1 + ksi_ij), (-1.*delta));
    B_ji = (double) pow((1 + ksi_ji), (-1.*delta));

    double B = (B_ij + B_ji)/2.0;

    return B;
}

double Brenner_i_energy(std::vector<std::vector<double>>& atoms, int i, bool correction){

    double R0, R1, R2, De, S, lambda, delta, a0, c0, d0;

    R0 = 1.315;
    R1 = 1.70;
    R2 = 2.00;
    De = 6.325;
    S = 1.29;
    lambda = 1.5;
    delta = 0.80469;
    a0 = 0.011304;
    c0 = 19;
    d0 = 2.5;

    double Vi = 0; // potencjał cząstkowy

    double r_i[3] = {atoms[i][3], atoms[i][4], atoms[i][5]};

    int atoms_num = atoms.size();
    for (int j = 0; j < atoms_num; j++){
        if(i == j){
            continue;
        }
        double r_j[3] = {atoms[j][3], atoms[j][4], atoms[j][5]};

        double r_ij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};
        double r_ij_len = sqrt( pow(r_ij[0],2) + pow(r_ij[1],2) + pow(r_ij[2],2) );

        // funkcja odcięcia ij
        double f_cut_ij = cut_function(r_ij_len, R1, R2);
        
        if(f_cut_ij != 0){
            // potencjal odpychania ij
            double VR = 0;
                VR = (double) De/(S-1) * exp(-sqrt(2.*S)*lambda*(r_ij_len - R0));

            // potencjal przyciagania ij
            double VA = 0;
                VA = (double) De*S/(S-1) * exp(-sqrt(2./S)*lambda*(r_ij_len - R0));

            // czynnik skalujacy potencjal przyciagania

            double B = B_val(atoms, i, j, r_ij, r_ij_len, R1, R2, a0, c0, d0, delta, correction);

            Vi += (double) f_cut_ij*( VR - B*VA );
        }
    }
    return Vi;
}

std::vector<double> pair_correlation_function(std::vector<std::vector<double>>& atoms, int M, double& r_mean){

    int atoms_num = atoms.size();

    for(int i = 0; i < atoms_num; i++){
        double r_i = sqrt(pow(atoms[i][3],2) + pow(atoms[i][4],2) + pow(atoms[i][5],2));
        r_mean += (double) r_i/atoms_num;
    }

    double r_max = 2.50 * r_mean;
    double Delta_r = (double) r_max/M;
    std::vector<double> pcf(M);

    for (int i = 0; i < atoms_num; i++){
        for(int j = i+1; j < atoms_num; j++){
            double r_i[3] = {atoms[i][3], atoms[i][4], atoms[i][5]};
            double r_j[3] = {atoms[j][3], atoms[j][4], atoms[j][5]};
            double r_ij = sqrt(pow(r_j[0] - r_i[0],2) + pow(r_j[1] - r_i[1],2) + pow(r_j[2] - r_i[2],2));
            int m = floor(r_ij/Delta_r);
            
            if(m<M){
                pcf[m] += (double) (2*4*M_PI*pow(r_mean,2))/(pow(atoms_num,2)*2*M_PI*r_ij*Delta_r);
            }
        }
    }

    return pcf;
}

void rand_atoms_pos(std::vector<std::vector<double>>& atoms, const int n, double r){
    
    std::vector<double> atom_row(6);
    
    for(int i = 0; i < n; i++){
        double U_phi = uni_rand();
        double U_theta = uni_rand();
        atom_row[0] = r; //promien
        atom_row[1] = 2*M_PI*U_phi; //kat phi
        atom_row[2] = M_PI*U_theta; //kat theta
        atom_row[3] = 0;
        atom_row[4] = 0;
        atom_row[5] = 0;
        atoms.push_back(atom_row);
        // std::cout<<atoms_rad[i][0]<<" "<<atoms_rad[i][1]<<" "<<atoms_rad[i][2]<<"\n";
    }
}

void atom_pos_change(std::vector<std::vector<double>>& atoms, int i, double beta, double w_r, double w_phi, double w_theta, bool correction){
    
    double V_old = Brenner_i_energy(atoms, i, correction);
    
    double U_r = uni_rand();
    double U_phi = uni_rand();
    double U_theta = uni_rand();

    double r = atoms[i][0];
    double phi = atoms[i][1];
    double theta = atoms[i][2];

    double Delta_r = r*(2*U_r-1)*w_r;
    double Delta_phi = phi*(2*U_phi-1)*w_phi;
    double Delta_theta = theta*(2*U_theta-1)*w_theta;

    double r_new = r + Delta_r;
    double phi_new = phi + Delta_phi;
    double theta_new = theta + Delta_theta;

    if(phi_new < 0){
        phi_new = phi_new + 2*M_PI;
    }
    if(phi_new > 2*M_PI){
        phi_new = phi_new - 2*M_PI;
    }

    if(theta_new < 0){
        theta_new = theta;
    }
    if(theta_new > M_PI){
        theta_new = theta;
    }

    atoms[i][0] = r_new;
    atoms[i][1] = phi_new;
    atoms[i][2] = theta_new;

    rad_to_xyz(atoms, i);
    
    double V_new = Brenner_i_energy(atoms, i, correction);

    double p_acc = std::min(1.0, exp((-1.0)*beta*(V_new - V_old)));
    double U_acc = uni_rand();
    if(U_acc > p_acc){
        atoms[i][0] = r;
        atoms[i][1] = phi;
        atoms[i][2] = theta;

        rad_to_xyz(atoms, i);
    }
}

void atoms_global_pos_change(std::vector<std::vector<double>>& atoms , const int n, double beta, double W_all, double& V_tot, bool correction){
    
    double V_old = 0;
    for(int i=0; i<n; i++){
        V_old += Brenner_i_energy(atoms, i, correction)/2.0;
    }

    double U_r = uni_rand();
    double r_rescale = (double) (1+W_all*(2.0*U_r - 1));

    for(int i = 0; i < n; i++){
        atoms[i][0] = (double) atoms[i][0]*r_rescale;
        rad_to_xyz(atoms, i);
    }

    double V_new = 0;
    for(int i=0; i<n; i++){
        V_new += Brenner_i_energy(atoms, i, correction)/2.0;
    }

    double p_acc = std::min(1.0, exp((-1.0)*beta*(V_new - V_old)));
    double U_acc = uni_rand();
    if(U_acc <= p_acc){
        V_tot = V_new;
    } else{
        V_tot = V_old;

        for(int i = 0; i < n; i++){
            atoms[i][0] = (double) atoms[i][0]/r_rescale;
            rad_to_xyz(atoms, i);
        }
    }
}

void rad_to_xyz(std::vector<std::vector<double>>& atoms, int i){
    
    double r = atoms[i][0];
    double phi = atoms[i][1];
    double theta = atoms[i][2];

    atoms[i][3] = r*sin(theta)*cos(phi);
    atoms[i][4] = r*sin(theta)*sin(phi);
    atoms[i][5] = r*cos(theta);
}

double beta_parameter(int it, int it_max, double beta_min, double beta_max, double p){
    
    double beta = (double) beta_min + pow(((double) it / (double) it_max), p)*(beta_max - beta_min);
    return beta;
}

double SA(const int n, const int it_max, int M, double w_r, double w_phi, double w_theta, double W_all, double r_start, double beta_min, double beta_max, double p,
        bool correction, std::ofstream *energyInFile=nullptr, std::ofstream *posInFile=nullptr, std::ofstream *pcfInFile=nullptr, std::ofstream *avogadroInFile=nullptr){

    std::vector<std::vector<double>> atoms;

    rand_atoms_pos(atoms, n, r_start);

    double V_tot;
    int atoms_num = atoms.size();
    double beta = 0;
    for(int it = 1; it <= it_max; it++){
        // std::cout<<"Iteracja: "<<it<<"\n";

        V_tot = 0;

        beta = beta_parameter(it, it_max, beta_min, beta_max, p);

        for(int i = 0; i < n; i++){
            atom_pos_change(atoms, i, beta, w_r, w_phi, w_theta, correction);
        }
        
        atoms_global_pos_change(atoms, n, beta, W_all, V_tot, correction);

        if(it%100 == 0){
            if(energyInFile != nullptr){
                double r_mean = 0;
                std::vector<double> pcf = pair_correlation_function(atoms, M, r_mean);
                std::ostream_iterator<double> output_pcf(*pcfInFile, "\n");
                std::copy(std::begin(pcf), std::end(pcf), output_pcf);

                // std::cout<<"Iteracja: "<<it<<"  Energia: "<<V_tot<<" eV\n";
                *energyInFile << it << " "<< beta << " " << V_tot << " " <<r_mean << "\n";

                std::ostream_iterator<double> output_pos(*posInFile, " ");
                for(const auto& vt : atoms) {
                    std::copy(vt.cbegin(), vt.cend(), output_pos);
                    *posInFile << "\n";
                }
            }
            if(avogadroInFile != nullptr){
                if(it == it_max){
                    std::cout<<"Iteracja: "<<it<<"  Energia: "<<V_tot<<" eV\n";
                    *avogadroInFile << n <<"\n\n";
                    for(int i=0; i<n; i++){
                        *avogadroInFile << "C " << atoms[i][3] << " " << atoms[i][4] << " " << atoms[i][5] << "\n";
                    }
                }
            }
        }
    }
    return V_tot;
}