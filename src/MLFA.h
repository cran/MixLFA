  #ifndef MLFA_HPP
#define MLFA_HPP


#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <random>
#include <Rcpp.h>
#include <RcppEigen.h>


using namespace Rcpp;
using namespace Eigen;
using namespace std;

class MLFA {
public:
    // Attributs
    int nb_facteur ;
    int nb_composantes ;
    int nb_observation;
    MatrixXd Y_resultat;
    MatrixXd X_covariable_fixes;
    MatrixXd Z_covariable_aleatoire;
    vector<MatrixXd> Lam;
    std::vector<Eigen::MatrixXd>  beta;
    vector<MatrixXd> s_b;
    vector<MatrixXd> s_e;
    VectorXd id_patient;
    vector<MatrixXd> tau;
    VectorXd pro;
    int max_it;
    int nb_sujet;
    int J; // var observee
    int p;
    int it = 0 ;
    int t_i ;
    vector<MatrixXd> Xn;
    vector<MatrixXd> Zn;
    MatrixXd BIC = MatrixXd::Zero(max_it,1);
    MatrixXd AIC = MatrixXd::Ones(max_it,1);
    MatrixXd ICL = MatrixXd::Ones(max_it,1);
    Rcpp::IntegerVector fixed_factor;

    MatrixXd riF ;
    int  npar =  (((J * nb_facteur - (nb_facteur * (nb_facteur + 1) / 2)) + J + nb_facteur * p + nb_facteur * (nb_facteur + 1) / 2)) + ((nb_facteur * p * (nb_facteur * p + 1) / 2));
    bool VerifNanPro  = false;



    int before=1;
    int after=2;






    struct Results {
        Eigen::MatrixXd Lam;
        Eigen::VectorXd beta;
        Eigen::VectorXd sig;
        Eigen::MatrixXd S_xi;
        Eigen::MatrixXd S_omega;
        double LogLik;
        Eigen::MatrixXd post_prob;
        double BIC;
        double AIC;
        double ICL;
        bool verif = false ;

    };




    struct factorRes {
        vector<vector<vector<MatrixXd>>> E_eta2_y;
        vector<vector<vector<MatrixXd>>> E_eta_y ;
        MatrixXd ri  ;
        MatrixXd Nc ;
    };

    struct betaRes {
        vector<vector<vector<MatrixXd>>> E_b_y;
        vector<vector<vector<MatrixXd>>> E_eta_y ;
        MatrixXd ri  ;
        MatrixXd Nc ;
    };



    // Constructeur

    MLFA(int nb_facteur, int nb_composantes, MatrixXd Y_resultat, MatrixXd X_covariable_fixes,
            MatrixXd Z_covariable_aleatoire, vector<Eigen::MatrixXd> Lam, vector<MatrixXd> beta,
            vector<MatrixXd> s_b, vector<MatrixXd> s_e, VectorXd id_patient, vector<MatrixXd> tau,
            VectorXd pro, int max_it,int nb_sujet , int J , int p, Rcpp::IntegerVector fixed_factor);

    // Méthodes

    void eStep(int k );
    void mStep(vector<vector<vector<MatrixXd>>> E_eta2_y,vector<vector<vector<MatrixXd>>> E_eta_y,vector<vector<vector<MatrixXd>>> E_b_y,vector<vector<vector<MatrixXd>>>E_bb_y,vector<vector<vector<MatrixXd>>> E_eta_b , int k );
    vector<factorRes> factor_loadings2(vector<vector<vector<MatrixXd>>> E_eta2_y,vector<vector<vector<MatrixXd>>> E_eta_y , MatrixXd ri  , MatrixXd Nc);
    MatrixXd y_seg(int ii);
    vector<MatrixXd>Defilement( int i , bool verif );
    vector<betaRes> beta_estimation2(vector<vector<vector<MatrixXd>>> E_eta_y , vector<vector<vector<MatrixXd>>> E_b_y, MatrixXd ri  , MatrixXd Nc );

    vector<Results>Omega( vector<vector<vector<MatrixXd>>>E_bb_y, vector<vector<vector<MatrixXd>>>E_eta2_y, vector<vector<vector<MatrixXd>>>E_eta_y, vector<vector<vector<MatrixXd>>>E_eta_b,vector<vector<vector<MatrixXd>>> E_b_y,  MatrixXd ri  , MatrixXd Nc);
    VectorXd mean_numba(const MatrixXd &a);
    MatrixXd krone(MatrixXd m1 , MatrixXd m2);
    std::vector<Eigen::MatrixXd> getLam() const { return Lam; }
    std::vector<Eigen::MatrixXd> getBeta() const { return beta; }
    std::vector<Eigen::MatrixXd> getS_b() const { return s_b; }
    std::vector<Eigen::MatrixXd> getS_e() const { return s_e; }
    std::vector<Eigen::MatrixXd> getTau() const { return tau; }
    Eigen::VectorXd getPro() const { return pro; }
    Eigen::MatrixXd getBIC() const { return BIC; }
    Eigen::MatrixXd getAIC() const { return AIC; }
    Eigen::MatrixXd getICL() const { return ICL; }
    Eigen::MatrixXd getRi() const { return riF; }








};

#endif // MLFA_HPP
