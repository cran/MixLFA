//
//  MLFA.cpp
//  test
//
//  Created by Omar Lahbabi on 15/05/2024.
//
#include "MLFA.h"
#include <iostream>
#include <vector>
using namespace std;


MLFA::MLFA(int nb_facteur, int nb_composantes, MatrixXd Y_resultat, MatrixXd X_covariable_fixes,
           MatrixXd Z_covariable_aleatoire, vector<Eigen::MatrixXd> Lam, vector<MatrixXd> beta,
           vector<MatrixXd> s_b, vector<MatrixXd> s_e, VectorXd id_patient, vector<MatrixXd> tau,
           VectorXd pro, int max_it,int nb_sujet , int J, int p, Rcpp::IntegerVector fixed_factor )
  : nb_facteur(nb_facteur),
    nb_composantes(nb_composantes),
    Y_resultat(Y_resultat),
    X_covariable_fixes(X_covariable_fixes),
    Z_covariable_aleatoire(Z_covariable_aleatoire),
    Lam(Lam),
    beta(beta),
    s_b(s_b),
    s_e(s_e),
    id_patient(id_patient),
    tau(tau),
    pro(pro),
    max_it(max_it),
    nb_sujet(nb_sujet),
    J(J),
    p(p),
    fixed_factor(fixed_factor)
{
  Xn = Defilement(0, true);
  Zn = Defilement(0, false);

}




void MLFA::eStep(int k)
{

  vector<int> indices;
  vector<MatrixXd> S_eta(nb_composantes);
  vector<MatrixXd> S_y(nb_composantes);
  vector<MatrixXd> S_yeta(nb_composantes);
  vector<MatrixXd> Syb(nb_composantes);
  vector<MatrixXd> S_etab(nb_composantes);
  vector<MatrixXd> S_eta_y(nb_composantes);
  vector<MatrixXd> S_b_y(nb_composantes);
  vector<MatrixXd> S_etab_y(nb_composantes);

  vector<vector<vector<MatrixXd>>> E_eta_y(nb_composantes);
  vector<vector<vector<MatrixXd>>> E_eta2_y(nb_composantes);
  vector<vector<vector<MatrixXd>>> E_eta_b(nb_composantes);
  vector<vector<vector<MatrixXd>>> E_b_y(nb_composantes);
  vector<vector<vector<MatrixXd>>> E_bb_y(nb_composantes);
  vector<vector<MatrixXd>> mu_b_y(nb_composantes); // Moyenne des effets aléatoires
  vector<vector<MatrixXd>> mu_eta_y(nb_composantes); // Moyenne facteur latents

  // ------------------ Definitions des matrices covariance -------------------

  for (int c = 0; c < nb_composantes; ++c)
  {

    for (int i = 0; i < nb_sujet; ++i)
    {
      t_i = std::count(id_patient.begin(), id_patient.end(), i + 1);
      MatrixXd x_n = Xn[i] ;
      MatrixXd z_n = Zn[i] ;
      // Indices des observations pour le sujet i
      vector<int> indices;

      for (int idx = 0; idx < id_patient.size(); ++idx)
      {
        if (id_patient[idx] == i + 1)
        {
          indices.push_back(idx);
        }
      }

      VectorXd y(J * t_i);


      for (int j = 0; j < t_i; ++j)
      {

        y.segment(j * J, J) = Y_resultat.row(indices[j]).transpose();
      }
      // Initialiser les matrices E_eta_y, E_eta2_y, etc., pour le sujet i
      E_eta_y[c].resize(nb_sujet);
      E_eta2_y[c].resize(nb_sujet);
      E_eta_b[c].resize(nb_sujet);
      E_b_y[c].resize(nb_sujet);
      E_bb_y[c].resize(nb_sujet);
      mu_b_y[c].resize(nb_sujet);
      mu_eta_y[c].resize(nb_sujet);

      E_eta_y[c][i].resize(t_i);
      E_eta2_y[c][i].resize(t_i);
      E_eta_b[c][i].resize(t_i);
      E_b_y[c][i].resize(t_i);
      E_bb_y[c][i].resize(t_i);
      // ---------------------- Calcul des matrices de covariances ------------------


      MatrixXd It = MatrixXd::Identity(t_i, t_i);
      MatrixXd Ib = MatrixXd::Identity(nb_facteur, nb_facteur); // Matrice identite
      // Matrice de covariance des facteurs latents pour la composante c.

      S_eta[c] = z_n * s_b[c]  * z_n.transpose() + krone(It, s_e[c]).eval();
      VectorXd tau2 =   tau[c].array().square();
      S_y[c] = krone(It, Lam[c]) * S_eta[c] * krone(It, Lam[c]).transpose() + krone(It, tau2.asDiagonal()).eval(); // probb le vr
      // Matrice de covariance des résultats pour la composante c.
      S_yeta[c] = krone(It, Lam[c]) * S_eta[c];
      // Matrice de covariance croisée entre les résultats et les facteurs latents pour la composante c.
      Syb[c] = krone(It, Lam[c]) * z_n * s_b[c];
      // Matrice de covariance croisée entre les résultats et les effets aléatoires pour la composante c.
      S_etab[c] = z_n * s_b[c];
      // Matrice de covariance entre les facteurs latents et les effets aléatoires pour la composante c.

      // ---------------- Calcul des moyennes espérées des facteurs latents et effets aléatoires ----------------

      MatrixXd a = (krone(It, Lam[c]) * x_n * beta[c]);
      // Valeurs attendues des résultats pour le sujet i compte tenu des covariables fixes et des effets fixes.
      MatrixXd b = (x_n * beta[c]);
      // Valeurs attendues des facteurs latents pour le sujet i.


      mu_eta_y[c][i] = b + S_yeta[c].transpose() * S_y[c].inverse() * (y - a);
      S_eta_y[c] = S_eta[c] - S_yeta[c].transpose() * S_y[c].inverse() * S_yeta[c];
      a = (krone(It, Lam[c]) * x_n * beta[c]);
      mu_b_y[c][i] = Syb[c].transpose() * S_y[c].inverse() * (y - a);
      S_b_y[c] = s_b[c] - Syb[c].transpose() * S_y[c].inverse() * Syb[c];
      S_etab_y[c] = S_etab[c] - S_yeta[c].transpose() * S_y[c].inverse() * Syb[c];

      for (int t = 0; t < t_i; ++t)
      {
        E_b_y[c][i][t] = mu_b_y[c][i];
        E_bb_y[c][i][t] = E_b_y[c][i][t] * E_b_y[c][i][t].transpose() + S_b_y[c];
        E_eta_y[c][i][t] = mu_eta_y[c][i].block(t * nb_facteur, 0, nb_facteur, 1);
        E_eta2_y[c][i][t] = E_eta_y[c][i][t] * E_eta_y[c][i][t].transpose() + S_eta_y[c].block(t * nb_facteur, t * nb_facteur, nb_facteur, nb_facteur);
        E_eta_b[c][i][t] = E_eta_y[c][i][t] * E_b_y[c][i][t].transpose() + S_etab_y[c].block(t * nb_facteur, 0, nb_facteur, nb_facteur * Z_covariable_aleatoire.cols());
      }
    }
  }




  k--;

  mStep(E_eta2_y,  E_eta_y , E_b_y , E_bb_y ,E_eta_b , k );
}


void MLFA::mStep(vector<vector<vector<MatrixXd>>> E_eta2_y,vector<vector<vector<MatrixXd>>> E_eta_y,vector<vector<vector<MatrixXd>>> E_b_y,vector<vector<vector<MatrixXd>>> E_bb_y,vector<vector<vector<MatrixXd>>> E_eta_b , int k  )
{

  MatrixXd ri ;
  MatrixXd Nc = MatrixXd::Zero(nb_sujet, nb_composantes);
  vector<Results> results_vector;
  vector<betaRes> beta_res;
  vector<factorRes> factor_res;
  VectorXd LL(nb_composantes);
  MatrixXd x_n;
  MatrixXd z_n;



  if (nb_composantes > 1)
  {
    ri  = MatrixXd::Ones(nb_sujet, nb_composantes) * 0.5;
    for (int i = 0; i < nb_sujet; ++i)
    {

      t_i = std::count(id_patient.begin(), id_patient.end(), i + 1);

      VectorXd LL(nb_composantes);
      for(int c = 0 ; c < nb_composantes ; c++ )
      {
        VectorXd tau2 = tau[c].array().square();
        x_n = Xn[i];
        z_n = Zn[i];
        LL(c) = 1;
        MatrixXd E,V;

        for(int t= 0 ; t < t_i ; t++ )
        {
          E = (Lam[c] * x_n.block(nb_facteur * t , 0 , nb_facteur , x_n.cols() )) * beta[c] ;



          MatrixXd diag_matrix = tau2.asDiagonal();

          V = Lam[c] * (z_n.block(nb_facteur * t , 0 , nb_facteur , z_n.cols()) * s_b[c]* z_n.block(nb_facteur * t , 0 , nb_facteur , z_n.cols()).transpose() + s_e[c]) * Lam[c].transpose() + diag_matrix;



          MatrixXd Y_segment = y_seg(i);


          MatrixXd y_petit = Y_segment.block(t,0,1,Y_segment.cols() ).transpose();

          Eigen::MatrixXd diff = y_petit - E;


          double det_V = V.determinant();
          double exponent_value = -0.5 * (diff.transpose() * V.inverse() * diff).value();

          double exp_term = exp(exponent_value);

          LL(c) = LL(c) *  ((1 / sqrt(det_V)) * exp_term );



        }


      }


      for(int c= 0 ; c < nb_composantes; c++ )
      {
        ri(i,c) = (pro(c)* LL(c))/(pro.array() * LL.array()).sum();
      }

      ri(i, nb_composantes-1) = 1 - (ri.block(i, 0, 1, nb_composantes-1).sum());

    }
  }

  if( nb_composantes == 1 )
  {
    ri = MatrixXd::Ones(nb_sujet, nb_composantes) ;
    VectorXd LL(nb_composantes);
  }



  for(int i = 0 ; i < nb_sujet ; i++)
  {
    t_i = std::count(id_patient.begin(), id_patient.end(), i + 1);

    for(int c= 0 ; c < nb_composantes ; c++ )
    {
      Nc(i,c) = (ri ( i , c) * t_i );
    }
  }



  factor_res = factor_loadings2(E_eta2_y, E_eta_y, ri, Nc);
  beta_res = beta_estimation2( E_eta_y,  E_b_y ,ri, Nc );
  results_vector = Omega( E_bb_y, factor_res[0].E_eta2_y, beta_res[0].E_eta_y, E_eta_b, beta_res[0].E_b_y, beta_res[0].ri, beta_res[0].Nc);

  if(it < max_it - 1)
  {

    it++;
  }




  if( k > 0)
  {
    Rcpp::Rcout << "iteration " << it  ;
    Rcpp::Rcout << "\n" ;
    eStep(k);
  }



}


vector<MLFA::factorRes>MLFA::factor_loadings2(vector<vector<vector<MatrixXd>>> E_eta2_y,vector<vector<vector<MatrixXd>>> E_eta_y , MatrixXd ri  , MatrixXd Nc)
{



  for (int c = 0; c < nb_composantes; c++)
  {
    MatrixXd Lam_i1 = MatrixXd::Zero(E_eta2_y[0][0][0].rows(), E_eta2_y[0][0][0].cols());

    MatrixXd Lam_i2 = MatrixXd::Zero(E_eta_y[0][0][0].rows(), J);

    MatrixXd taui1 = MatrixXd::Zero(J, J);
    MatrixXd taui2 = MatrixXd::Zero(J, J);
    MatrixXd taui3 = MatrixXd::Zero(J, J);

    for (int i = 0; i < nb_sujet; i++) {

      t_i = count(id_patient.begin(), id_patient.end(), i + 1);
      vector<int> indices;
      for (int idx = 0; idx < id_patient.size(); ++idx) {
        if (id_patient[idx] == i + 1) {
          indices.push_back(idx);
        }
      }

      VectorXd y = VectorXd::Zero(J * t_i, 1);
      for (int j = 0; j < t_i; ++j) {
        y.segment(j * J, J) = Y_resultat.row(indices[j]).transpose();
      }

      MatrixXd Lam_1 = MatrixXd::Zero(E_eta2_y[0][0][0].rows(), E_eta2_y[0][0][0].cols());
      MatrixXd Lam_2 = MatrixXd::Zero(E_eta_y[0][0][0].rows(), J);


      MatrixXd tau1 = MatrixXd::Zero(J, J);
      MatrixXd tau2 = MatrixXd::Zero(J, J);
      MatrixXd tau3 = MatrixXd::Zero(J, J);

      for (int t = 0; t < t_i; ++t)
      {

        int start_index = t * J;
        Lam_1 += ri(i, c) * E_eta2_y[c][i][t];


        Lam_2 += ri(i, c) * (E_eta_y[c][i][t] * y.segment(start_index, J).transpose());

        tau1 += ri(i, c) * (y.segment(start_index, J) * y.segment(start_index, J).transpose());
        tau2 -= 2 * ri(i, c) * (y.segment(start_index, J) * E_eta_y[c][i][t].transpose() * Lam[c].transpose());
        tau3 += ri(i, c) * (Lam[c] * E_eta2_y[c][i][t] * Lam[c].transpose());
      }

      Lam_i1 += Lam_1;
      Lam_i2 += Lam_2;

      taui1 += tau1;
      taui2 += tau2;
      taui3 += tau3;
    }

    Lam[c] = (Lam_i1.inverse() * Lam_i2).transpose();
    tau[c] = (1.0 / Nc.col(c).sum() * (taui1 + taui2 + taui3).diagonal()).array().sqrt();



    int pos = 0;

    for (int i = 0; i < fixed_factor.size(); i++) {
      int row = fixed_factor[i];

      Lam[c](row, pos) = 1;

      for (int j = pos + 1; j < nb_facteur; j++) {
        Lam[c](row, j) = 0;
      }

      pos++;
    }


  }






  std::vector<factorRes> results_vector;
  factorRes res;
  res.E_eta2_y = E_eta2_y;
  res.E_eta_y = E_eta_y;
  res.Nc = Nc ;
  res.ri = ri ;


  results_vector.push_back(res);
  return results_vector;

}



MatrixXd MLFA::y_seg(int ii )
{

  vector<int> indices;
  bool test = false ;
  int step = 0  ;
  int step2 = 0 ;

  for (int i = 0; i < id_patient.size(); ++i)
  {



    if(id_patient[i]  != ii && test == false )
    {
      step += 1 ;
    }

    if(id_patient[i]  == ii && test == false )
    {
      if( ii >0 )
      {
        step2+= t_i;
      }
      step2 += step ;
      test = true ;
    }


  }




  MatrixXd y( t_i , J);


  for (int j = 0; j < t_i; ++j)
  {
    int index =   j + step2 ;
    y.row(j) = Y_resultat.row(index);
  }


  return y ;

}




vector<MatrixXd> MLFA::Defilement( int ii , bool verif )
{
  vector<MatrixXd> xn ;
  vector<MatrixXd> zn ;


  for (int i = 0; i < nb_sujet; ++i)
  {
    vector<int> indices;
    t_i = count(id_patient.begin(), id_patient.end(), i + 1);

    for (int idx = 0; idx < id_patient.size(); ++idx)
    {
      if (id_patient[idx] == i + 1)
      {
        indices.push_back(idx);
      }
    }

    MatrixXd x(t_i, X_covariable_fixes.cols());
    MatrixXd z(t_i, Z_covariable_aleatoire.cols());
    //VectorXd y(J * t_i);


    for (int j = 0; j < t_i; ++j)
    {
      x.row(j) = X_covariable_fixes.row(indices[j]);
      z.row(j) = Z_covariable_aleatoire.row(indices[j]);
      // y.segment(j * J, J) = Y_resultat.row(indices[j]).transpose();
    }

    MatrixXd It = MatrixXd::Identity(t_i, t_i);
    MatrixXd Ib = MatrixXd::Identity(nb_facteur, nb_facteur); // Matrice identite

    MatrixXd x_n = krone(Ib , x.row(0)).eval();
    MatrixXd z_n = krone(Ib , z.row(0)).eval();


    if (t_i > 1)
    {
      for (int ob = 1; ob < x.rows(); ++ob)
      {
        // Calculer les nouvelles données
        MatrixXd new_x_data = krone(Ib, x.row(ob)).eval();
        MatrixXd new_z_data = krone(Ib, z.row(ob)).eval();

        // Créer des matrices temporaires pour les nouvelles tailles
        MatrixXd x_n_temp = MatrixXd::Zero(x_n.rows() + new_x_data.rows(), x_n.cols());
        MatrixXd z_n_temp = MatrixXd::Zero(z_n.rows() + new_z_data.rows(), z_n.cols());

        // Copier les anciennes données dans les matrices temporaires
        if (x_n.rows() > 0) {
          x_n_temp.topRows(x_n.rows()) = x_n;
        }
        if (z_n.rows() > 0) {
          z_n_temp.topRows(z_n.rows()) = z_n;
        }

        // Ajouter les nouvelles données dans les matrices temporaires
        x_n_temp.bottomRows(new_x_data.rows()) = new_x_data;
        z_n_temp.bottomRows(new_z_data.rows()) = new_z_data;

        // Remplacer les anciennes matrices par les nouvelles
        x_n = x_n_temp;
        z_n = z_n_temp;
      }
    }



    xn.push_back(x_n);
    zn.push_back(z_n);

  }
  if(verif == true  )
  {

    return xn ;
  }
  else{
    return zn;
  }


}


vector<MLFA::betaRes> MLFA::beta_estimation2(vector<vector<vector<MatrixXd>>> E_eta_y , vector<vector<vector<MatrixXd>>> E_b_y, MatrixXd ri  , MatrixXd Nc )
{


  for(int c = 0 ; c < nb_composantes ; c++ )
  {

    MatrixXd beta_i1 = MatrixXd::Zero(nb_facteur, nb_facteur);
    MatrixXd beta_i2 = MatrixXd::Zero( nb_facteur, 1);

    for (int i = 0; i < nb_sujet; ++i)
    {
      t_i = std::count(id_patient.begin(), id_patient.end(), i + 1);
      MatrixXd x_n = Xn[i];
      MatrixXd z_n = Zn[i];

      if( i == 0)
      {

        beta_i1 = x_n.block(0, 0, nb_facteur, x_n.cols()).transpose() * s_e[c].inverse() * x_n.block(0, 0, nb_facteur, x_n.cols()) * 0.0 ;
        beta_i2 = x_n.block(0, 0, nb_facteur, x_n.cols()).transpose() * s_e[c].inverse() * (E_eta_y[c][0][0] - z_n.block(0, 0, nb_facteur, z_n.cols() ) * E_b_y[c][0][0]) * 0.0;
      }


      MatrixXd beta1 = MatrixXd::Zero(nb_facteur, nb_facteur);
      MatrixXd beta2 = MatrixXd::Zero(nb_facteur, 1);

      beta1 =x_n.block(0, 0, nb_facteur, x_n.cols()).transpose() * s_e[c].inverse() * x_n.block(0, 0, nb_facteur, x_n.cols()) * 0.0;
      beta2 = x_n.block(0, 0, nb_facteur, x_n.cols()).transpose()  * s_e[c].inverse() * (E_eta_y[c][0][0] - z_n.block(0, 0, nb_facteur, z_n.cols()) * E_b_y[c][0][0]) * 0.0;

      for(int t = 0 ; t < t_i ; t++)
      {
        beta1 += ri(i, c) * x_n.block( t * nb_facteur, 0, nb_facteur, x_n.cols()).transpose()* s_e[c].inverse()
        * x_n.block(t * nb_facteur, 0, nb_facteur, x_n.cols());
        VectorXd diff = E_eta_y[c][i][t] - z_n.block(t * nb_facteur, 0, nb_facteur, z_n.cols()) * E_b_y[c][i][t];
        beta2 += ri(i, c)* x_n.block(t * nb_facteur, 0, nb_facteur, x_n.cols()).transpose()* s_e[c].inverse()* diff;
      }


      beta_i1=beta_i1+beta1;
      beta_i2=beta_i2+beta2;
    }

    MatrixXd beta_result = beta_i1.inverse() * beta_i2;
    beta[c] = Map<MatrixXd>(beta_result.data(), nb_facteur * p, 1);



  }



  std::vector<betaRes> results_vector;
  betaRes res;

  res.E_b_y = E_b_y;
  res.E_eta_y = E_eta_y;
  res.Nc = Nc ;
  res.ri = ri ;

  results_vector.push_back(res);

  return results_vector;


}








vector<MLFA::Results> MLFA::Omega( vector<vector<vector<MatrixXd>>>E_bb_y, vector<vector<vector<MatrixXd>>>E_eta2_y, vector<vector<vector<MatrixXd>>>E_eta_y, vector<vector<vector<MatrixXd>>>E_eta_b,vector<vector<vector<MatrixXd>>> E_b_y,  MatrixXd ri  , MatrixXd Nc )
{

  for(int c= 0 ; c < nb_composantes ; c++)
  {
    MatrixXd S_ib = MatrixXd::Zero(E_bb_y[c][0][0].rows(), E_bb_y[c][0][0].cols());
    MatrixXd S_ie = MatrixXd::Zero(s_e[c].rows(), s_e[c].cols());


    for (int i = 0; i < nb_sujet; ++i)
    {


      t_i = std::count(id_patient.begin(), id_patient.end(), i + 1);

      MatrixXd x_n = Xn[i] ;
      MatrixXd z_n = Zn[i] ;
      MatrixXd  Sb1 = MatrixXd::Zero(E_bb_y[c][0][0].rows(), E_bb_y[c][0][0].cols());
      MatrixXd  Se1= MatrixXd::Zero(s_e[c].rows(), s_e[c].cols());
      MatrixXd beta_reshaped = Map<MatrixXd>(beta[c].data(), 1 , nb_facteur * p);
      MatrixXd beta_reshaped2 = Map<MatrixXd>(beta[c].data(), nb_facteur * p , 1 );



      for(int t = 0 ; t < t_i ; t++)
      {


        MatrixXd SS = E_eta2_y[c][i][t]
        - ((E_eta_y[c][i][t] * beta_reshaped) * x_n.block(t*nb_facteur,0,nb_facteur,x_n.cols()).transpose())
        - (E_eta_b[c][i][t] * z_n.block(t*nb_facteur,0,nb_facteur,z_n.cols()).transpose())
        - ((x_n.block(t*nb_facteur,0,nb_facteur,x_n.cols()) * beta_reshaped2) * E_eta_y[c][i][t].transpose())
        - (z_n.block(t*nb_facteur,0,nb_facteur,z_n.cols()) * E_eta_b[c][i][t].transpose())
        + (((x_n.block(t*nb_facteur,0,nb_facteur,x_n.cols()) * beta_reshaped2) * E_b_y[c][i][t].transpose()) * z_n.block(t*nb_facteur,0,nb_facteur,z_n.cols()).transpose())
        + (((z_n.block(t*nb_facteur,0,nb_facteur,z_n.cols()) * E_b_y[c][i][t]) * beta_reshaped) * x_n.block(t*nb_facteur,0,nb_facteur,x_n.cols()).transpose())
        + (((x_n.block(t*nb_facteur,0,nb_facteur,x_n.cols()) * beta_reshaped2) * beta_reshaped) * x_n.block(t*nb_facteur,0,nb_facteur,x_n.cols()).transpose())
        + ((z_n.block(t*nb_facteur,0,nb_facteur,z_n.cols()) * E_bb_y[c][i][t]) * z_n.block(t*nb_facteur,0,nb_facteur,z_n.cols()).transpose());


        Se1 = Se1 + (1 / Nc.col(c).sum()) * ri(i, c) * SS;

      }

      S_ib += ri(i, c) * E_bb_y[c][i][t_i-1];
      S_ie = S_ie + Se1;
    }

    s_b[c] = ( 1 / ri.col(c).sum()) * S_ib;
    s_b[c] = (s_b[c] + s_b[c].transpose()) / 2.0;
    s_e[c] = (S_ie + S_ie.transpose()) / 2.0;
  }

  // -------------------------------------------------------


  pro = mean_numba(ri) ;
  MatrixXd L = MatrixXd::Zero(nb_sujet, nb_composantes);
  int N;



  for(int c =  0 ; c < nb_composantes ; c++)
  {

    N = Y_resultat.rows();


    MatrixXd x_n ;
    MatrixXd z_n ;
    MatrixXd E,V;
    VectorXd tau2 = tau[c].array().square();

    for(int i = 0 ; i < nb_sujet; i++)
    {

      t_i = std::count(id_patient.begin(), id_patient.end(), i + 1);

      x_n = Xn[i];
      z_n = Zn[i];

      L(i,c) = 1 ;




      for( int t = 0 ; t < t_i ; t++ )
      {
        E = (Lam[c] * x_n.block(nb_facteur * t , 0 , nb_facteur , x_n.cols() )) * beta[c] ;
        MatrixXd diag_matrix = tau2.asDiagonal();

        V = Lam[c] * (z_n.block(nb_facteur * t , 0 , nb_facteur , z_n.cols()) * s_b[c]* z_n.block(nb_facteur * t , 0 , nb_facteur , z_n.cols()).transpose() + s_e[c]) * Lam[c].transpose() + diag_matrix;


        MatrixXd Y_segment =y_seg(i);
        MatrixXd y_petit = Y_segment.block(t,0,1,Y_segment.cols() ).transpose();


        Eigen::MatrixXd diff = y_petit - E;

        double det_V = V.determinant();
        double exponent_value = -0.5 * (diff.transpose() * V.inverse() * diff).value();
        double exp_term = exp(exponent_value);

        L(i, c) = L(i,c) *  ((1 / sqrt(2 * M_PI * det_V)) * exp_term );


      }
    }
  }


  VectorXd ll = pro[0] * L.col(0) * 0.0;
  VectorXd cl = pro[0] * L.col(0) * 0.0;

  for(int c = 0 ; c < nb_composantes ; c++ )
  {
    ll = ll + pro[c] * L.col(c);
    VectorXd ent = ri.col(c).array() * ri.col(c).array().log();

    for (int i = 0; i < ent.size(); ++i)
    {
      if (isnan(ent(i)))
      {
        ent(i) = 0;

      }
    }
    cl = cl + ent;
  }

  double Lik = ll.array().log().sum();
  double ICL2 = cl.sum();



  MatrixXd BIC_MLFA = MatrixXd::Zero(max_it,1);
  MatrixXd AIC_MLFA = MatrixXd::Ones(max_it,1);
  MatrixXd ICL_MLFA = MatrixXd::Ones(max_it,1);

  N = Y_resultat.rows();






  BIC(it) = -2 * Lik + (nb_composantes * ((J * nb_facteur - (nb_facteur * (nb_facteur + 1) / 2)) + J + nb_facteur * p + nb_facteur * (nb_facteur + 1) / 2) + nb_composantes - 1) * log(N) + (nb_composantes * (nb_facteur * p * (nb_facteur * p + 1) / 2)) * log(nb_sujet);
  AIC(it) = -2 * Lik + 2 * (nb_composantes * ((J * nb_facteur - (nb_facteur * (nb_facteur + 1) / 2)) + J + nb_facteur * p + nb_facteur * (nb_facteur + 1) / 2) + nb_composantes - 1 + nb_composantes * (nb_facteur * p * (nb_facteur * p + 1) / 2));
  ICL(it) = BIC(it) - 2 * ICL2;



  riF = ri ;
  BIC_MLFA(it) = -2 * Lik + (nb_composantes * ((J * nb_facteur - (nb_facteur * (nb_facteur + 1) / 2)) + J + nb_facteur * p + nb_facteur * (nb_facteur + 1) / 2) + nb_composantes - 1) * log(N) + (nb_composantes * (nb_facteur * p * (nb_facteur * p + 1) / 2)) * log(nb_sujet);
  AIC_MLFA(it) = -2 * Lik + 2 * (nb_composantes * ((J * nb_facteur - (nb_facteur * (nb_facteur + 1) / 2)) + J + nb_facteur * p + nb_facteur * (nb_facteur + 1) / 2) + nb_composantes - 1 + nb_composantes * (nb_facteur * p * (nb_facteur * p + 1) / 2));
  ICL_MLFA(it) = BIC_MLFA(it) - 2 * ICL2;



  std::vector<Results> results_vector;
  Results res;


  // Ajouter la structure au vecteur de résultats
  results_vector.push_back(res);

  return results_vector;

}








VectorXd MLFA:: mean_numba(const MatrixXd &a) {

  int cols = a.cols();
  VectorXd res(cols);
  for (int i = 0; i < cols; ++i)
  {
    res(i) = a.col(i).mean();
  }
  return res;
}




MatrixXd MLFA::krone(MatrixXd m1, MatrixXd m2)
{

  int m , n , pp , q ;
  m = m1.cols();
  n = m1.rows();
  pp = m2.cols();
  q = m2.rows();


  int p2 = m*pp;
  int n2 = n*q;

  MatrixXd resultat = MatrixXd::Zero(n2,p2);

  for(int i = 0 ; i < n ; i++ )
  {
    for(int j = 0 ; j < m ; j++ )
    {
      for(int k = 0 ; k < q ; k++ )
      {
        for(int  t = 0 ; t < pp ; t++ )
        {
          resultat(i * q + k, j * pp + t) =m1(i, j) * m2(k, t);
        }

      }
    }
  }

  return resultat;


}








