#include <Eigen/Dense>
#include "MLFA.h"
#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace Eigen;
using namespace std;
using namespace Rcpp;

// Définition de la struct Retour avec les types appropriés
struct Retour {
    std::vector<Eigen::MatrixXd> Lam;
    std::vector<Eigen::MatrixXd> beta;
    std::vector<Eigen::MatrixXd> S_b;
    std::vector<Eigen::MatrixXd> S_e;
    std::vector<Eigen::MatrixXd> tau;
    Eigen::VectorXd pro;
    Eigen::MatrixXd BIC;
    Eigen::MatrixXd AIC;
    Eigen::MatrixXd ICL;
    Eigen::MatrixXd ri;
    bool Verif;
};

// Fonction Compute
Retour Compute(int nb_composantes, int nb_facteur, Eigen::MatrixXd X_covariable_fixes, Eigen::MatrixXd Y_resultat, Eigen::MatrixXd Z_covariable_aleatoire, Eigen::VectorXd id_patient, std::vector<Eigen::MatrixXd> Lam, std::vector<Eigen::MatrixXd> beta, std::vector<Eigen::MatrixXd> tau, std::vector<Eigen::MatrixXd> s_e, std::vector<Eigen::MatrixXd> s_b, int max_it, Eigen::VectorXd pro2,int nb_sujet , int J , int p, Rcpp::IntegerVector fixed_factor) {
    MLFA model(nb_facteur, nb_composantes, Y_resultat, X_covariable_fixes, Z_covariable_aleatoire, Lam, beta, s_b, s_e, id_patient, tau, pro2, max_it, nb_sujet , J  , p, fixed_factor);



    Retour var;
    model.eStep(max_it);
    var.Lam = model.getLam();
    var.beta = model.getBeta();
    var.S_b = model.getS_b();
    var.S_e = model.getS_e();
    var.tau = model.getTau();
    var.pro = model.getPro();
    var.BIC = model.getBIC();
    var.AIC = model.getAIC();
    var.ICL = model.getICL();
    var.ri = model.getRi();
    var.Verif = model.VerifNanPro;

    return var;

}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List MLFA(int nb_composantes, int nb_facteur, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::MatrixXd& Z, const Eigen::VectorXd& id, const std::vector<Eigen::MatrixXd>& Lam, const std::vector<Eigen::MatrixXd>& beta, const std::vector<Eigen::MatrixXd>& tau, const std::vector<Eigen::MatrixXd>& s_e, const std::vector<Eigen::MatrixXd>& s_b, int max_it, const Eigen::VectorXd& pro2,int nb_sujet , int J , int p, Rcpp::IntegerVector fixed_factor) {


    Retour res = Compute(nb_composantes, nb_facteur, X, Y, Z, id, Lam, beta, tau, s_e, s_b, max_it, pro2,nb_sujet , J , p, fixed_factor);

    return List::create(
        Named("Lam") = res.Lam,
        Named("beta") = res.beta,
        Named("S_b") = res.S_b,
        Named("S_e") = res.S_e,
        Named("tau") = res.tau,
        Named("pro") = res.pro,
        Named("BIC") = res.BIC,
        Named("AIC") = res.AIC,
        Named("ICL") = res.ICL,
        Named("ri") = res.ri,
        Named("VerifNan") = res.Verif
    );
}

