#include <Rcpp.h>
using namespace Rcpp;

// 1 - Moran statistic function: two stats, one value for NS neighbors and a different for EW.//
 // [[Rcpp::export]]
NumericVector moran(NumericVector y, NumericMatrix NN){
  int n = y.size();
  double ym, dd = 0., res = 0, ss, nn=0, res2=0, nn2=0, ss0;
  NumericVector out(2);
  ym = mean(y);
  // NS neighbors
  for(int i=0; i<n; i++) {
    for(int j=1; j<3; j++) {
      if ( NN(i,j) > 0 ) {
        res += (y[i]-ym)*(y[NN(i,j)-1]-ym);
        nn++;
      }
    }

    // EW neighbors
    for(int k=3; k<5; k++) {
      if ( NN(i,k) > 0 ) {
        res2 += (y[i]-ym)*(y[NN(i,k)-1]-ym);
        nn2++;
      }
    }

    dd = y[i] - ym;
    ss += dd * dd;
  }

  if (ss == 0) {
    out[0] = 0;
    out[1] = 0;
  } else {
  out[0] = (n * res)/(ss*nn);
  out[1] = (n * res2)/(ss*nn2);
  }
  return(out);
}


 // 2.1 - Simulate a WP mrf, in a regular grid with rook structure neightbors
 // [[Rcpp::export]]
NumericMatrix sim_mrfC(int meanpar, double bet0, double bet1, double ens, double eew, int R, int iter, NumericVector a, NumericMatrix NN) {
  int n = NN.nrow(), r1=0, r2=0; // asuming NN is order to have independents chuncks.
  // meanpar define parametrization for the mean: 1 - centered, 2=ID link, and 3-original.

  NumericMatrix out(iter, n); // if we want to return all simulations not only lastone
  NumericVector Ycurr(n), eta(4), vn(5), xx(2); // Yim stores previous iteration, YY stores current iteration (output)
  std::vector< std::string > prs;

  eta[0] = ens;
  eta[1] = ens;
  eta[2] = eew;
  eta[3] = eew;

  if ( meanpar == 1 ) {
  // initial value
  for (int j = 0; j < n; j++) {
    r1  = rpois(1, 1 )[0];
    Ycurr[j] = std::min(r1 , R);
  }

 for (int i = 1; i < iter; i++) {
  // Ycurr = YY(_,i-1);
      for (int j = 0; j < n; j++) {
        vn = NN(j,_); // the first element of vn is an index, not a neighbor
        double etapar1 = 0  ;
        for (int k = 1; k < 5; k++) {
          if (vn[k] > 0 ) {
            etapar1 += eta[k-1]*(Ycurr[vn[k]-1] - (bet0 + bet1*a[vn[k]-1]))/4.0;
          }
        }
        double etapar = log(bet0 + bet1*a[j]) + etapar1 ;
        if ((bet0 + bet1*a[j]) <= 0) {
          Ycurr[j] = 0;
        } else if ( etapar > log(4*R) ) {
            Ycurr[j] = R;
        } else {
        r2  = rpois(1, exp(etapar) )[0];
        Ycurr[j] = std::min(r2,R) ;
        out(i,j) = Ycurr[j];
        }
    }
  }
  } else if ( meanpar == 2 ) {
      // initial value
  for (int j = 0; j < n; j++) {
    // r1  = rpois(1, bet0 + bet1*a[j])[0];
    r1  = rpois(1, 1)[0];
    Ycurr[j] = std::min(r1 , R);
  }

 for (int i = 1; i < iter; i++) {
  // Ycurr = YY(_,i-1);
      for (int j = 0; j < n; j++) {
        vn = NN(j,_); // the first element of vn is an index, not a neighbor
        double etapar = bet0 + bet1*a[j]  ;
        for (int k = 1; k < 5; k++) {
          if (vn[k] > 0 ) {
            etapar += eta[k-1]*(Ycurr[vn[k]-1] - bet0 - bet1*a[vn[k]-1])/4.0;
          }
        }
        if (etapar <0) {
          Ycurr[j] = 0;
//        } else
//        if (etapar > log(R)) {
//          Ycurr[j] = R;
        } else {
          r2  = rpois(1, etapar )[0];
          Ycurr[j] = std::min(r2,R) ;
        out(i,j) = Ycurr[j];
        }
    }
  }
  } else if ( meanpar == 3 ) {
      // initial value
  for (int j = 0; j < n; j++) {
    r1  = rpois(1, 1) [0];
    Ycurr[j] = std::min(r1 , R);
  }

 for (int i = 1; i < iter; i++) {
  // Ycurr = YY(_,i-1);
      for (int j = 0; j < n; j++) {
        vn = NN(j,_); // the first element of vn is an index, not a neighbor
        double etapar = bet0 + bet1*a[j]  ;
        for (int k = 1; k < 5; k++) {
          if (vn[k] > 0 ) {
            etapar += eta[k-1]*(Ycurr[vn[k]-1] - exp( bet0 + bet1*a[vn[k]-1]) )/4.0;
          }
        }
        if ( etapar > log(4*R) ) {
          Ycurr[j] = R;
        } else {
        r2  = rpois(1, exp(etapar) )[0];
        Ycurr[j] = std::min(r2,R) ;
        out(i,j) = Ycurr[j];
        }
    }
  }
  }
  return out;
}

 // 2.2 - Simulate a WP mrf, in a regular grid with rook structure neightbors
 //        with Log parametrization.
 // [[Rcpp::export]]
NumericVector sim_mrfClog(double bet0, double bet1, double ens, double eew, int R, int iter, NumericVector a, NumericMatrix NN) {
  int n = NN.nrow(), r1=0, r2=0; // asuming NN is order to have independents chuncks.
  NumericVector Ycurr(n), eta(4), vn(5), xx(2); // Yim stores previous iteration, YY stores current iteration (output)

  eta[0] = ens;
  eta[1] = ens;
  eta[2] = eew;
  eta[3] = eew;


  return Ycurr;
}





// 3 - Standar ABC basic code

//// [[Rcpp::export]]
//NumericMatrix abc_a2(NumericVector obs, int smp, int R, NumericVector a, NumericMatrix NN) {
//  int n = obs.size(), iter=500, c1, c2, c3, atmp;
//  NumericVector z(n), mrnObs(2),mi(2), theta(4);
//  NumericMatrix out(smp,7);
//  double bet0, bet1, ens, eew, mnObs;
//
//  mnObs = mean(obs);
//  mrnObs = moran(obs, NN);
//
//  for (int i = 0; i < smp; i++) {
//  // simulate from the prior
//  bet0 = rnorm(1, 0, 1)[0];
//  bet1 = runif(1, -1, 1)[0];
//  ens = 2*rbeta(1, 4, 4)[0] - 1;
//  eew = 2*rbeta(1, 4, 4)[0] - 1;
//
//  // simulate data set
//  z = sim_mrfC(bet0, bet1, ens, eew, R, iter, a, NN);
//
//// simulate data set: control if the simulation are all at the truncation value, or that somethings
//// went wrong at cpp and some simulation are negative.
//// c1=0, c2=0, c3=0, atmp=0;
//// while( (c1 == 0 || c2 == 0) & (atmp < 5) ) {
//// z = sim_mrfC(bet, ens, eew, R, iter, a, NN);
////  if( mean(z) == R ) {    c1 = 0; } else {c1 = 1; }
////  if( min(z) < 0 )  { c2 = 0;  } else { c2 = 1;}
////  atmp += 1; }
//
//  // compute statistic and distance value and save results
//  out(i,0)=bet0; out(i,1)=bet1; out(i,2)=ens; out(i,3)=eew;
//  if (mean(z) < 0) {
//    for (int k =4; k< 8; k++) {
//      out(i,k) = -1;
//    }
//  } else {
//  out(i, 4)  =  mean(z);
//  mi  =  moran(z, NN);
//  out(i, 5) = mi[0]; out(i, 6) = mi[1];
//  }
//  }
//  return(out);
//}

