#include <RcppArmadillo.h>
#include <Rcpp.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

//' R's implementation of the Moore-Penrose pseudo matrix inverse
//' 
//' @param m matrix
//' @return matrix inverse of m
// [[Rcpp::export]]
arma::mat Rginv(const arma::mat m){
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, m, "dc");
  arma::uvec Positive = arma::find(S > 1E-06 * S(1));
  if(all(Positive)){
    arma::mat D = diagmat(S);
    return V * (1/D * U.t());
  }else if(!any(Positive)){
    return arma::zeros(m.n_rows, m.n_cols);
  }else{
    S.elem(Positive) = 1/S.elem(Positive);
    arma::mat D = diagmat(S);
    return V * D * U.t();
  }
}
// arma::mat Rginv(const arma::mat& m){
//   arma::mat U, V;
//   arma::vec S;
//   arma::svd(U, S, V, m, "dc");
//   arma::uvec Positive = arma::find(S > 0.0);
//   if(Positive.size() == 0){
//     return arma::zeros(m.n_rows, m.n_cols);
//   }else if(all(Positive)){
//     arma::mat D = diagmat(S);
//     return V * (1/D * U.t());
//   }else if(!any(Positive)){
//     return arma::zeros(m.n_rows, m.n_cols);
//   }else{
//     S.elem(Positive) = 1/S.elem(Positive);
//     arma::mat D = diagmat(S);
//     return V * D * U.t();
//   }
// }

//' Generalized matrix inverse
//' 
//' @param m matrix
//' @return matrix inverse of m
// [[Rcpp::export]]
arma::mat gen_inv(arma::mat m){
  arma::mat out(m.n_rows, m.n_cols);
  try{
    out = inv(m);
  }catch(std::exception &ex){
    out = Rginv(m);
  }
  return out;
}

//' Steady State Probabilities
//' 
//' Finds the steady state probabilities from a transition matrix
//' mat = |p_11 p_21 ... p_m1|
//'       |p_12 p_22 ... p_m2|
//'       |...            ...|
//'       |p_1m p_2m ... p_mm|
//' where the columns sum to 1
//' 
//' @param mat square SxS matrix of probabilities with column sums of 1. S
//' represents the number of states
//' @return matrix of dimensions Sx1 with steady state probabilities
//' @examples
//' \dontrun{
//' library(kimfilter)
//' Pm = rbind(c(0.8406, 0.0304), 
//'            c(0.1594, 0.9696))
//' ss_prob(Pm)
//' }
//' @export
// [[Rcpp::export]]
arma::mat ss_prob(arma::mat mat){
  arma::mat Zero(mat.n_rows, 1);
  arma::mat One(1, mat.n_cols);
  One = One.ones();
  arma::mat A1(mat.n_rows, mat.n_rows);
  
  arma::mat A = join_cols(A1.eye() - mat, One);
  arma::mat B = gen_inv(A.t() * A) * A.t();
  return B * join_cols(Zero, One.submat(0, 0, 0, 0));
}

//' Matrix self rowbind
//' 
//' @param mat matrix
//' @param times integer
//' @return matrix
// [[Rcpp::export]]
arma::mat self_rbind(arma::mat mat, int times){
  arma::mat ret(times, mat.n_cols, arma::fill::zeros);
  for(int i = 0; i < times; i++){
    ret.row(i) = mat;
  }
  return ret;
}

//' Check if list contains a name
//' 
//' @param s a string name
//' @param L a list object
//' @return boolean
// [[Rcpp::export]]                                                                                                                                           
bool contains(std::string s, Rcpp::List L){                                                                                                                  
  Rcpp::CharacterVector nv = L.names();                                                                                                                     
  for(int i = 0; i < nv.size(); i++){                                                                                                                         
    if(std::string(nv[i]) == s){                                                                                                                        
      return true;                                                                                                                                      
    }                                                                                                                                                     
  }                                                                                                                                                         
  return false;                                                                                                                                             
} 

//' Kim Filter
//' 
//' @param ssm list describing the state space model, must include names
//' B0 - N_b x 1 x n_state array of matrices, initial guess for the unobserved components 
//' P0 - N_b x N_b x n_state array of matrices, initial guess for the covariance matrix of the unobserved components
//' Dm - N_b x 1 x n_state array of matrices, constant matrix for the state equation
//' Am - N_y x 1 x n_state array of matrices, constant matrix for the observation equation
//' Fm - N_b X p x n_state array of matrices, state transition matrix
//' Hm - N_y x N_b x n_state array of matrices, observation matrix
//' Qm - N_b x N_b x n_state array of matrices, state error covariance matrix
//' Rm - N_y x N_y x n_state array of matrices, state error covariance matrix
//' betaO - N_y x N_o x n_state array of matrices, coefficient matrix for the observation exogenous data
//' betaS - N_b x N_s x n_state array of matrices, coefficient matrix for the state exogenous data
//' Pm - n_state x n_state matrix, state transition probability matrix
//' @param yt N x T matrix of data
//' @param Xo N_o x T matrix of exogenous observation data
//' @param Xs N_s x T matrix of exogenous state 
//' @param weight column matrix of weights, T x 1
//' @param smooth boolean indication whether to run the backwards smoother
//' @return list of cubes and matrices output by the Kim filter
// [[Rcpp::export]]
Rcpp::List kim_filter_cpp(Rcpp::List& ssm, const arma::mat& yt, 
                      Rcpp::Nullable<Rcpp::NumericMatrix> Xo = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericMatrix> Xs = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericMatrix> weight = R_NilValue,
                      bool smooth = false){
  //Initialize matrices
  int n_cols = yt.n_cols;
  int n_rows = yt.n_rows;
  arma::cube B0 = ssm["B0"];
  arma::cube P0 = ssm["P0"];
  arma::cube Dm = ssm["Dm"];
  arma::cube Am = ssm["Am"];
  arma::cube Fm = ssm["Fm"];
  arma::cube Hm = ssm["Hm"];
  arma::cube Qm = ssm["Qm"];
  arma::cube Rm = ssm["Rm"];
  arma::mat Pm = ssm["Pm"];
  arma::cube betaO;
  arma::cube betaS;
  arma::mat X_o;
  arma::mat X_s;
  
  //Set the exogenous matrices
  if(Xo.isNotNull()){
    X_o = Rcpp::as<arma::mat>(Xo);
  }else{
    X_o = arma::zeros(1, n_cols);
  }
  if(Xs.isNotNull()){
    X_s = Rcpp::as<arma::mat>(Xs);
  }else{
    X_s = arma::zeros(1, n_cols);
  }
  
  //Set the exogenous coefficient matrices
  int n_states = Pm.n_cols;
  if(contains("betaO", ssm)){
    betaO = Rcpp::as<arma::cube>(ssm["betaO"]);
  }else{
    betaO = arma::zeros(n_rows, 1, n_states);
  }
  if(contains("betaS", ssm)){
    betaS = Rcpp::as<arma::cube>(ssm["betaS"]);
  }else{
    betaS = arma::zeros(Fm.n_rows, 1, n_states);
  }
  
  //Define the storage lists
  arma::cube B_tts(Fm.n_rows, n_cols, n_states, arma::fill::zeros);
  arma::cube B_tls(Fm.n_rows, n_cols, n_states, arma::fill::zeros);
  arma::cube B_tlss(Fm.n_rows, n_cols, n_states*n_states, arma::fill::zeros);
  arma::cube B_ttss(Fm.n_rows, n_cols, n_states*n_states, arma::fill::zeros);
  arma::cube N_tss(n_rows, n_cols, n_states*n_states, arma::fill::zeros);
  arma::cube N_ts(n_rows, n_cols, n_states, arma::fill::zeros);
  arma::mat N_t(n_rows, n_cols, arma::fill::zeros);
  arma::cube F_tss(n_rows, n_rows, n_states*n_states, arma::fill::zeros);
  arma::cube F_ts(n_rows, n_rows, n_states, arma::fill::zeros);
  arma::cube F_t(n_rows, n_rows, n_cols, arma::fill::zeros);
  arma::cube K_tss(Fm.n_rows, n_rows, n_states*n_states, arma::fill::zeros);
  arma::cube K_ts(Fm.n_rows, n_rows, n_states, arma::fill::zeros);
  arma::cube K_t(Fm.n_rows, n_rows, n_cols, arma::fill::zeros);
  arma::cube PR_VLss(1, 1, n_states*n_states, arma::fill::zeros);
  arma::cube PROss(1, 1, n_states*n_states, arma::fill::zeros);
  
  arma::field<arma::cube> P_tts(n_states);
  arma::field<arma::cube> P_tls(n_states);
  for(int i = 0; i < n_states; i++){
    P_tts(i) = arma::cube(Fm.n_rows, Fm.n_rows, n_cols, arma::fill::zeros);  
    P_tls(i) = arma::cube(Fm.n_rows, Fm.n_rows, n_cols, arma::fill::zeros);  
  }
  arma::cube P_tl(Fm.n_rows, Fm.n_rows, n_cols, arma::fill::zeros);
  arma::cube P_tt(Fm.n_rows, Fm.n_rows, n_cols, arma::fill::zeros);
  
  arma::field<arma::cube> P_tlss(n_states*n_states);
  arma::field<arma::cube> P_ttss(n_states*n_states);
  for(int i = 0; i < n_states*n_states; i++){
    P_tlss(i) = arma::cube(Fm.n_rows, Fm.n_rows, n_cols, arma::fill::zeros);  
    P_ttss(i) = arma::cube(Fm.n_rows, Fm.n_rows, n_cols, arma::fill::zeros);  
  }
  
  arma::mat Pr_tts(n_cols, n_states, arma::fill::zeros);
  arma::mat Pr_tls(n_cols, n_states, arma::fill::zeros);
  arma::mat Pr_tl(n_states, n_states, arma::fill::zeros);
  arma::mat B_tt(Fm.n_rows, n_cols, arma::fill::zeros);
  arma::mat B_tl(Fm.n_rows, n_cols, arma::fill::zeros);
  arma::uvec non_na_idx;
  arma::uvec iv(1);
  
  arma::cube Am_tt(Am.n_rows, Am.n_cols, n_cols, arma::fill::zeros);
  arma::cube Hm_tt(Hm.n_rows, Hm.n_cols, n_cols, arma::fill::zeros);
  arma::cube Rm_tt(Rm.n_rows, Rm.n_cols, n_cols, arma::fill::zeros);
  arma::cube Dm_tt(Dm.n_rows, Dm.n_cols, n_cols, arma::fill::zeros);
  arma::cube Fm_tt(Fm.n_rows, Fm.n_cols, n_cols, arma::fill::zeros);
  arma::cube Qm_tt(Qm.n_rows, Qm.n_cols, n_cols, arma::fill::zeros);
  arma::cube betaO_tt(betaO.n_rows, betaO.n_cols, n_cols, arma::fill::zeros);
  arma::cube betaS_tt(betaS.n_rows, betaS.n_cols, n_cols, arma::fill::zeros);
  arma::mat lnl(1, n_cols, arma::fill::zeros);
  arma::mat temp;
  arma::mat y_tl(n_rows, n_cols);
  arma::mat y_tt(n_rows, n_cols);
  
  //Define some matrix transforms
  arma::mat Fm_t(Fm.n_rows, Fm.n_cols);
  arma::mat Hm_t(Hm.n_rows, Hm.n_cols);
  arma::mat F_tss_inv(F_tss.n_rows, F_tss.n_cols);
  arma::mat F_tss_submat;
  
  //Initialize the filter
  arma::cube B_lls = B0;
  arma::cube P_lls = P0;
  arma::mat Pr = ss_prob(Pm);
  arma::mat w(n_cols, 1);
  int s = 0;
  double Pr_val = 0.0;
  
  //Rescale the weights
  if(weight.isNotNull()){
    w = Rcpp::as<arma::mat>(weight);
  }else{
    w = arma::ones(n_cols, 1);
  }
  w = w * n_cols/arma::as_scalar(sum(w));  
  
  //Kim Filter = Hamilton + Kalman filter routine
  for(int i = 0; i < n_cols; i++){
    Pr_tls.row(i) = (Pm * Pr).t();
    
    // Joint probabilities conditional on t-1
    //Pr[S_t=i,S_{t-1}=j|Y_{t-1}] = Pr[S_t=i|S_{t-1}=j,Y_{t-1}]*Pr[S_{t-1}=j|Y_{t-1}] 
    Pr_tl = Pm % self_rbind(Pr.t(), n_states);
    
    //Reinitialize
    s = 0;
    Pr_val = 0.0;
    
    //Find the non-missing values
    non_na_idx = arma::find_finite(yt.col(i));
    
    for(int stl = 0; stl < n_states; stl++){
      for(int st = 0; st < n_states; st++){
        Fm_t = Fm.slice(st).t();
        Hm_t = Hm.slice(st).t();
        
        //When S_{t-1}=j, S_{t}=i
        //B^{i,j}_{t|t-1} = D_j + Fm %*% B^{j}_{t|t-1}
        B_tlss.slice(s).col(i) = Dm.slice(st) + Fm.slice(st) * B_lls.slice(stl) + betaS.slice(st) * X_s.col(i);
        
        //Initial predictions of the unobserved component
        //P^{i,j}_{t|t-1} = P^{i,j}_{t|t-1} = Fm %*% P^{j}_{t-1|t-1} %*% t(Fm) + Qm
        P_tlss(s).slice(i) = Fm.slice(st) * P_lls.slice(stl) * Fm_t + Qm.slice(st);
        
        //Forecast errors
        //N^{i,j}_{t|t-1} = yti[, j] - Hm %*% B^{i,j}_{t|t-1}
        N_tss.slice(s).col(i) = yt.col(i) - Am.slice(st) - Hm.slice(st) * B_tlss.slice(s).col(i) + betaO.slice(st) * X_o.col(i);
        
        //Variance of forecast errors
        //F^{i,j}_{t|t-1} = Hm %*% P^{i,j}_{t|t-1} %*% t(Hm) + Rm
        F_tss.slice(s) = Hm.slice(st) * P_tlss(s).slice(i) * Hm_t + Rm.slice(st);
        
        //Variance of forecast errors
        //F^{i,j}_{t|t-1} = Hm %*% P^{i,j}_{t|t-1} %*% t(Hm) + Rm
        F_tss_inv = arma::zeros(F_tss.slice(s).n_rows, F_tss.slice(s).n_cols);
        if(!non_na_idx.is_empty()){
          //Variance of the prediction error conditional on t-1
          F_tss.slice(s).submat(non_na_idx, non_na_idx) = Hm.slice(st).rows(non_na_idx) * P_tlss(s).slice(i) * Hm_t.cols(non_na_idx) + Rm.slice(st).submat(non_na_idx, non_na_idx);
          F_tss_submat = F_tss.slice(s).submat(non_na_idx, non_na_idx);
          F_tss_inv.submat(non_na_idx, non_na_idx) = gen_inv(F_tss_submat);
        }
        //Kalman gain conditional on t-1
        K_tss.slice(s) = P_tlss(s).slice(i) * Hm_t * F_tss_inv;

        if(!non_na_idx.is_empty()){
          iv[0] = i;
          //Final estimate of the unobserved values
          B_ttss.slice(s).col(i) = B_tlss.slice(s).col(i) + K_tss.slice(s).cols(non_na_idx) * N_tss.slice(s).submat(non_na_idx, iv);
          //Final estimate of the covariance matrix
          P_ttss(s).slice(i) = P_tlss(s).slice(i) - K_tss.slice(s).cols(non_na_idx) * Hm.slice(st).rows(non_na_idx) * P_tlss(s).slice(i);
          //Update the likelihood
          PR_VLss.slice(s) = ((1/sqrt(det(F_tss.slice(s).submat(non_na_idx, non_na_idx)))) * exp(-0.5*N_tss.slice(s).submat(non_na_idx, iv).t() * F_tss_inv.submat(non_na_idx, non_na_idx) * N_tss.slice(s).submat(non_na_idx, iv))) * Pr_tl(st, stl);
        }else{
          B_ttss.slice(s).col(i) = B_tlss.slice(s).col(i);
          P_ttss(s).slice(i) = P_tlss(s).slice(i);
          PR_VLss.slice(s) = 0.0;
        }
        
        //Joint density conditional on t-1
        Pr_val = Pr_val + arma::as_scalar(PR_VLss.slice(s));
        s++;
      }
    }
    
    for(int s = 0; s < n_states*n_states; s++){
      //Joint probabilities conditional on t
      //Pr[S_t=j,S_{t-1}=i|Y_t] = f[y_t|S_t=j,S_{t-1}=i,Y_{t-1}]*Pr[S_t=j,S{t-1}=i|Yt]/f[y_t|Y_{t-1}]
      //PRO^{i,j} = PR_VL^{i,j}/Pr_val
      if(Pr_val != 0){
        if(arma::as_scalar(Pr_val) > 0){
          PROss.slice(s) = PR_VLss.slice(s)/arma::as_scalar(Pr_val);
        }else{
          PROss.slice(s) = PR_VLss.slice(s)/arma::as_scalar(arma::datum::inf);
        }
      }else{
        PROss.slice(s) = 0.0;
      }
    }
    
    for(int s = 0; s < n_states; s++){
      //Probabilities conditional on t: Pr[S_t=s|Yt]
      //Pr[j] = sum_{i}(PRO^{i,j})
      if(Pr_val != 0.0){
        Pr.row(s) = 0.0;
        for(int j = 0; j < n_states; j++){
          Pr.row(s) += PROss.slice(s + j*n_states);
        }
      }else{
        Pr.row(s) = Pr_tls.row(s).col(s);
      }
    }
    Pr_tts.row(i) = Pr.t();
    
    //Collapsing terms
    K_t.slice(i) = arma::zeros(K_t.slice(i).n_rows, K_t.slice(i).n_cols);
    for(int s = 0; s < n_states; s++){
      B_tts.slice(s).col(i) = arma::zeros(B_tts.slice(s).col(i).n_rows, B_tts.slice(s).col(i).n_cols);
      B_tls.slice(s).col(i) = arma::zeros(B_tls.slice(s).col(i).n_rows, B_tls.slice(s).col(i).n_cols);
      for(int j = 0; j < n_states; j++){
        //B^{j}_{t|t} = (sum_{i}(PRO^{i,j}*B^{i,j}_{t|t}))/Pr[j]
        B_tts.slice(s).col(i) += arma::as_scalar(PROss.slice(s + j*n_states))*B_ttss.slice(s + j*n_states).col(i);
        B_tls.slice(s).col(i) += arma::as_scalar(PROss.slice(s + j*n_states))*B_tlss.slice(s + j*n_states).col(i);
      }
      if(arma::as_scalar(Pr(s, 0)) <= 0){
        B_tts.slice(s).col(i) /= arma::as_scalar(arma::datum::inf);
        B_tls.slice(s).col(i) /= arma::as_scalar(arma::datum::inf);
      }else{
        B_tts.slice(s).col(i) /= arma::as_scalar(Pr.row(s));
        B_tls.slice(s).col(i) /= arma::as_scalar(Pr.row(s));
      }
      
      F_ts.slice(s) = arma::zeros(F_ts.slice(s).n_rows, F_ts.slice(s).n_cols);
      N_ts.slice(s) = arma::zeros(N_ts.slice(s).n_rows, N_ts.slice(s).n_cols);
      K_ts.slice(s) = arma::zeros(K_ts.slice(s).n_rows, K_ts.slice(s).n_cols);
      P_tts(s).slice(i) = arma::zeros(P_tts(s).slice(i).n_rows, P_tts(s).slice(i).n_cols);
      P_tls(s).slice(i) = arma::zeros(P_tls(s).slice(i).n_rows, P_tls(s).slice(i).n_cols);
      for(int j = 0; j < n_states; j++){
        //K^{j} = (sum_{i}(PRO^{i,j}*K^{i,j}))/Pr[j]
        F_ts.slice(s) += arma::as_scalar(PROss.slice(s + j*n_states))*F_tss.slice(s + j*n_states);
        N_ts.slice(s) += arma::as_scalar(PROss.slice(s + j*n_states))*N_tss.slice(s + j*n_states);
        K_ts.slice(s) += arma::as_scalar(PROss.slice(s + j*n_states))*K_tss.slice(s + j*n_states);
        
        //P^{j}_{t|t} = (sum{i}(PRO^{i,j}*(P^{i,j}_{t|t} + (B^{j}_{t|t} - B^{i,j}_{t|t}) %*% t(B^{j}_{t|t} - B^{i,j}_{t|t}))))/Pr[j]
        temp = B_tts.slice(s).col(i) - B_ttss.slice(s + j*n_states).col(i);
        P_tts(s).slice(i) += arma::as_scalar(PROss.slice(s + j*n_states))*(P_ttss(s + j*n_states).slice(i) + temp*temp.t());
        
        temp = B_tls.slice(s).col(i) - B_tlss.slice(s + j*n_states).col(i);
        P_tls(s).slice(i) += arma::as_scalar(PROss.slice(s + j*n_states))*(P_tlss(s + j*n_states).slice(i) + temp*temp.t());
      }
      if(arma::as_scalar(Pr(s, 0)) <= 0){
        F_ts.slice(s) /= arma::as_scalar(arma::datum::inf);
        N_ts.slice(s) /= arma::as_scalar(arma::datum::inf);
        K_ts.slice(s) /= arma::as_scalar(arma::datum::inf);
        P_tts(s).slice(i) /= arma::as_scalar(arma::datum::inf);
        P_tls(s).slice(i) /= arma::as_scalar(arma::datum::inf);
      }else{
        F_ts.slice(s) /= arma::as_scalar(Pr.row(s));
        N_ts.slice(s) /= arma::as_scalar(Pr.row(s));
        K_ts.slice(s) /= arma::as_scalar(Pr.row(s));
        P_tts(s).slice(i) /= arma::as_scalar(Pr.row(s));
        P_tls(s).slice(i) /= arma::as_scalar(Pr.row(s));
      }
     
      K_t.slice(i) += arma::as_scalar(Pr.row(s)) * K_ts.slice(s);
      F_t.slice(i) += arma::as_scalar(Pr.row(s)) * F_ts.slice(s);
      N_t.col(i) += arma::as_scalar(Pr.row(s)) * N_ts.slice(s).col(i);
      Am_tt.slice(i) += arma::as_scalar(Pr.row(s)) * Am.slice(s);
      Hm_tt.slice(i) += arma::as_scalar(Pr.row(s)) * Hm.slice(s);
      Rm_tt.slice(i) += arma::as_scalar(Pr.row(s)) * Rm.slice(s);
      Dm_tt.slice(i) += arma::as_scalar(Pr.row(s)) * Dm.slice(s);
      Fm_tt.slice(i) += arma::as_scalar(Pr.row(s)) * Fm.slice(s);
      Qm_tt.slice(i) += arma::as_scalar(Pr.row(s)) * Qm.slice(s);
      betaO_tt.slice(i) += arma::as_scalar(Pr.row(s)) * betaO.slice(s);
      betaS_tt.slice(i) += arma::as_scalar(Pr.row(s)) * betaS.slice(s);
      B_tt.col(i) += (arma::as_scalar(Pr.row(s)) * B_tts.slice(s).col(i));
      B_tl.col(i) += (arma::as_scalar(Pr.row(s)) * B_tls.slice(s).col(i));
      P_tl.slice(i) += arma::as_scalar(Pr.row(s)) * P_tls(s).slice(i);
      P_tt.slice(i) += arma::as_scalar(Pr.row(s)) * P_tts(s).slice(i);
      
      B_lls.slice(s) = B_tts.slice(s).col(i);
      P_lls.slice(s) = P_tts(s).slice(i);
    }
    
    if(Pr_val != 0){
      lnl.col(i) = -log(Pr_val);
    }else{
      lnl.col(i) = 0.0;
    }
    
    // Get the finalized predictions
    y_tl.col(i) = Am_tt.slice(i) + Hm_tt.slice(i) * B_tl.col(i) + betaO_tt.slice(i) * X_o.col(i);
    y_tt.col(i) = Am_tt.slice(i) + Hm_tt.slice(i) * B_tt.col(i) + betaO_tt.slice(i) * X_o.col(i);
  }
  
  //Rescale the weights to sum to the length of the data
  w = w * n_cols / arma::as_scalar(sum(w));
  
  if(smooth == true){
    //Define variables
    int t = n_cols - 1;
    
    //Define storage lists
    arma::cube B_tTss(Fm.n_rows, 1, n_states*n_states, arma::fill::zeros);
    arma::cube P_tTss(Fm.n_rows, Fm.n_rows, n_states*n_states, arma::fill::zeros);
    arma::cube Pr_tTss(1, 1, n_states*n_states, arma::fill::zeros);
    arma::mat temp;
    
    for(int i = t - 1; i >= 0; i--){
      int s = 0;
      
      for(int st = 0; st < n_states; st++){
        for(int stf = 0; stf < n_states; stf++){
          //Full information inference on unobserved component and its covariance matrix
          //B^{j,k}_{t|T} = B^{j}_{t|t} + P^{j}_{t|t} %*% t(Fm) %*% solve(P^{j,k}_{t+1|t}) %*% (B^{k}_{t+1|T} - B^{j,k}_{t+1|t})
          temp = P_tts(st).slice(i) * Fm_t * Rginv(P_tlss(s).slice(i + 1));
          B_tTss.slice(s) = B_tts.slice(st).col(i) + temp * (B_tts.slice(st).col(i + 1) - B_tlss.slice(s).col(i + 1));
          
          // //P^{j,k} = P^{j}_{t|t} + P^{j}_{t|t} %*% t(Fm) %*% solve(P^{j,k}_{t+1|t}) %*% (P^{k}_{t+t|T} - P^{j,k}_{t+1|t}) %*% t(P^{j}_{t|t} %*% t(Fm) %*% solve(P^{j,k}_{t+1|t}))
          P_tTss.slice(s) = P_tts(st).slice(i) + temp * (P_tts(st).slice(i + 1) - P_tlss(s).slice(i + 1)) * temp.t();
          
          //Full information inference on probabilities
          //Pr[S_t|Y_T]: #Pr[S_{t+1}=k|Y_T]*Pr[S_{t+1}=k|S_t=j]*Pr[S_t=j|Y_t]/Pr[S_{t+1}=k|Y_t]
          Pr_tTss.slice(s) = Pr_tts(i + 1, stf) * Pm(stf, st) * Pr_tts(i, st);
          if(arma::as_scalar(Pr_tls(i + 1, stf)) <= 0){
            Pr_tTss.slice(s) /= arma::as_scalar(arma::datum::inf);
          }else{
            Pr_tTss.slice(s) /= Pr_tls(i + 1, stf);
          }
          s++;
        }
      }
      
      //Collapsing terms
      for(int s = 0; s < n_states; s++){
        //Pr[S_t=d|Y_T]
        //P^{j}_{t|t} = sum_{i}(Pr[S_t|Y_T]^{i,j})
        Pr_tts(i, s) = 0;
        for(int j = n_states - 1; j >= 0; j--){
          Pr_tts(i, s) += arma::as_scalar(Pr_tTss.slice((s + 1)*n_states - j - 1));
        }
        
        //B^{j}_{t|T} = (sum_{i}(Pr[S_t|Y_T]^{i,j}*B^{j,k}_{t|T}))/Pr[j]
        B_tts.slice(s).col(i) = arma::zeros(B_tts.slice(s).col(i).n_rows, B_tts.slice(s).col(i).n_cols);
        for(int j = n_states - 1; j >= 0; j--){
          B_tts.slice(s).col(i) += arma::as_scalar(Pr_tTss.slice((s + 1)*n_states - j - 1)) * B_tTss.slice((s + 1)*n_states - j - 1);
        }
        if(arma::as_scalar(Pr_tts(i, s)) <= 0){
          B_tts.slice(s).col(i) /= arma::as_scalar(arma::datum::inf);
        }else{
          B_tts.slice(s).col(i) /= Pr_tts(i, s);
        }
        
        //P^{j}_{t|T} = (sum_{i}(Pr[S_t|Y_T]^{i,j}*(B^{j}_{t|T} - B^{j,k}_{t|T}) %*% t(B^{j}_{t|T} - B^{j,k}_{t|T})))/Pr[j]
        P_tts(s).slice(i) = arma::zeros(P_tts(s).slice(i).n_rows, P_tts(s).slice(i).n_cols);
        for(int j = n_states - 1; j >= 0; j--){
          arma::mat temp = (B_tts.slice(s).col(i) - B_tTss.slice((s + 1)*n_states - j - 1));
          P_tts(s).slice(i) += arma::as_scalar(Pr_tTss.slice((s + 1)*n_states - j - 1))*(P_tTss.slice((s + 1)*n_states - j - 1) + temp * temp.t());
        }
        if(arma::as_scalar(Pr_tts(i, s)) <= 0){
          P_tts(s).slice(i) /= arma::as_scalar(arma::datum::inf);
        }else{
          P_tts(s).slice(i) /= Pr_tts(i, s);
        }
        
        B_tt.col(i) += (arma::as_scalar(Pr_tts(i, s)) * B_tts.slice(s).col(i));
        Am_tt.slice(i) += arma::as_scalar(Pr_tts(i, s)) * Am.slice(s);
        Hm_tt.slice(i) += arma::as_scalar(Pr_tts(i, s)) * Hm.slice(s);
        Rm_tt.slice(i) += arma::as_scalar(Pr_tts(i, s)) * Rm.slice(s);
        Dm_tt.slice(i) += arma::as_scalar(Pr_tts(i, s)) * Dm.slice(s);
        Fm_tt.slice(i) += arma::as_scalar(Pr_tts(i, s)) * Fm.slice(s);
        Qm_tt.slice(i) += arma::as_scalar(Pr_tts(i, s)) * Qm.slice(s);
      }
      
      // Get the finalized predictions
      y_tt.col(i) = Am_tt.slice(i) + Hm_tt.slice(i) * B_tt.col(i) + betaO_tt.slice(i) * X_o.col(i);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("lnl") = arma::as_scalar(-(lnl * w)),
                            Rcpp::Named("Pr_tl") = Pr_tls,
                            Rcpp::Named("Pr_tt") = Pr_tts,
                            Rcpp::Named("y_tl") = y_tl,
                            Rcpp::Named("y_tt") = y_tt,
                            Rcpp::Named("B_tl") = B_tl,
                            Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("P_tl") = P_tl,
                            Rcpp::Named("P_tt") = P_tt,
                            Rcpp::Named("Am_tt") = Am_tt,
                            Rcpp::Named("Hm_tt") = Hm_tt,
                            Rcpp::Named("Rm_tt") = Rm_tt,
                            Rcpp::Named("Dm_tt") = Dm_tt,
                            Rcpp::Named("Fm_tt") = Fm_tt, 
                            Rcpp::Named("Qm_tt") = Qm_tt,
                            Rcpp::Named("betaO_tt") = betaO_tt,
                            Rcpp::Named("betaS_tt") = betaS_tt,
                            Rcpp::Named("F_t") = F_t,
                            Rcpp::Named("K_t") = K_t);
}
