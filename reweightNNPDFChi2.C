double pt[] = {7.02, 7.97, 9.9, 11.6, 13.4, 15.6, 19,
	       22.2, 25.7, 29.7, 34.4, 39.7, 46.3, 53.8};
double ptsyst[] = {0.26, 0.30, 0.36, 0.40, 0.47, 0.5, 0.60,
		   0.63, 0.74, 0.83, 0.95, 1.1, 1.2, 1.5};
double asym[] = {0.00002, -0.0022, 0.0016, 0.0005, 0.0015, 0.0029, 0.0016,
		 0.0044, 0.0050, 0.0036, 0.0169, -0.0049, 0.0122, 0.0018};
double aLL_stat[] = {0.0013, 0.0014, 0.0010, 0.0011, 0.0013, 0.0016, 0.0016,
		     0.0018, 0.0021, 0.0027, 0.0037, 0.0054, 0.0084, 0.0137};
double aLL_syst[] = {0.0001, 0.0005, 0.0001, 0.0001, 0.0001, 0.0002, 0.0001,
		     0.0002, 0.0002, 0.0003, 0.0004, 0.0003, 0.0007, 0.0008};
double ue[] = {0.00036, 0.00033, 0.00031, 0.00028, 0.00027, 0.00026, 0.00025,
		     0.00025, 0.00024, 0.00023, 0.00023, 0.00023, 0.00023, 0.00023};
//di-jet
double minv_a[] = {19.04, 21.62, 26.38, 32.30, 38.42,
		   45.18, 53.42, 63.52, 75.55, 89.12};
double minv_unc_a[] = {0.82, 0.96, 1.08, 1.16, 1.20,
		       1.48, 1.55, 1.80, 2.21, 2.38};
double aLL_a[] = {-0.0128, 0.0090, 0.0079, -0.0012, 0.0101,
		  -0.0013, 0.0048, 0.0052, 0.0363, -0.0218};
double aLL_stat_a[] = {0.0066, 0.0052, 0.0050, 0.0052, 0.0061,
		       0.0064, 0.0081, 0.0108, 0.0167, 0.0264};
double aLL_syst_a[] = {0.0002, 0.0002, 0.0004, 0.0004, 0.0004,
		       0.0006, 0.0009, 0.0010, 0.0020, 0.0048};
double ue_a[] = {0.00025, 0.00058, 0.00062, 0.00056, 0.00051,
		 0.00046, 0.00044, 0.00037, 0.00037, 0.00031};
//
double minv_b[] = {18.80, 21.80, 26.37, 32.24, 38.42,
		   45.83, 54.14, 64.17, 76.06, 89.81, 107.92};
double minv_unc_b[] = {0.81, 0.94, 1.09, 1.07, 1.31,
		       1.41, 1.69, 1.83, 2.13, 2.48, 2.89};
double aLL_b[] = {-0.0023, 0.0041, 0.0016, 0.0029, -0.0063,
		  0.0020, 0.0128, -0.0022, -0.0010, -0.0160, -0.0205};
double aLL_stat_b[] = {0.0053, 0.0036, 0.0033, 0.0034, 0.0040,
		       0.0041, 0.0050, 0.0065, 0.0096, 0.0143, 0.0242};
double aLL_syst_b[] = {0.0002, 0.0003, 0.0003, 0.0004, 0.0005,
		       0.0007, 0.0008, 0.0011, 0.0014, 0.0020, 0.0027};
double ue_b[] = {0.00040, 0.00063, 0.00041, 0.00043, 0.00035,
		 0.00036, 0.00030, 0.00031, 0.00028, 0.00026, 0.00026};
//
double minv_c[] = {19.81, 22.01, 26.16, 32.70, 38.71,
		   46.01, 53.89, 64.75, 77.15, 91.14};
double minv_unc_c[] = {0.81, 1.12, 1.08, 1.19, 1.34,
		       1.36, 1.78, 1.74, 2.27, 2.58};
double aLL_c[] = {0.0058, -0.0006, -0.0043, 0.0049, 0.0046,
		  0.0155, -0.0045, 0.0104, 0.0346, 0.0593};
double aLL_stat_c[] = {0.0085, 0.0066, 0.0062, 0.0065, 0.0077,
		       0.0079, 0.0098, 0.0127, 0.0192, 0.0294};
double aLL_syst_c[] = {0.0002, 0.0006, 0.0004, 0.0006, 0.0007,
		       0.0006, 0.0012, 0.0014, 0.0019, 0.0073};
double ue_c[] = {0.00026, 0.00063, 0.00074, 0.00078, 0.00081,
		 0.00052, 0.00053, 0.00041, 0.00041, 0.00041};
//
double minv_d[] = {20.54, 22.09, 26.31, 31.72, 38.48,
		   45.21, 54.15, 64.30, 76.66, 89.57, 106.86};
double minv_unc_d[] = {2.14, 0.95, 0.98, 1.22, 1.31,
		       1.55, 1.69, 1.92, 2.15, 2.62, 2.97};
double aLL_d[] = {0.0054, 0.0042, 0.0051, -0.0031, -0.0018,
		  -0.0040, 0.0034, 0.0050, 0.0058, 0.0291, -0.0055};
double aLL_stat_d[] = {0.0161, 0.0051, 0.0050, 0.0059, 0.0060,
		       0.0070, 0.0087, 0.0123, 0.0178, 0.0296, 0.0461};
double aLL_syst_d[] = {0.0002, 0.0002, 0.0003, 0.0003, 0.0004,
		       0.0005, 0.0008, 0.0011, 0.0024, 0.0014, 0.0028};
double ue_d[] = {0.00131, 0.00022, 0.00072, 0.00052, 0.00054,
		 0.00044, 0.00044, 0.00036, 0.00036, 0.00036, 0.00029};

//
const char *ver = "v2";
class Grfun{
public:
  Grfun() : mGr(0x0) {}
  Grfun(TGraph *gr) : mGr(gr) {}
  double operator() (double *x, double *p) {
    // function implementation using class data members
    if(mGr) return mGr->Eval(x[0], 0, "S");
    else return 0.;
    //return x[0];
  }
  TGraph *mGr;
};
const int ndi = 14;
const int ntop = 4;
const int top[] = {1, 2, 4, 3};
const int ndj_a = 10, ndj_b =11, ndj_c = 10, ndj_d = 11;
//calculating chi^2 and weights for replica irep
double reweightNNPDFChi2(int irep = 1)
{
  //total dimensions of all data points
  int nd = 56;
  TMatrixD corr = fillMatrix("run12.correlation.v1.txt", nd);
  corr.Print();

  TMatrixD sigm_i = getSigma(14, aLL_stat, aLL_syst, ue); 
  TMatrixD sigm_a = getSigma(10, aLL_stat_a, aLL_syst_a, ue_a); 
  TMatrixD sigm_b = getSigma(11, aLL_stat_b, aLL_syst_b, ue_b); 
  TMatrixD sigm_c = getSigma(10, aLL_stat_c, aLL_syst_c, ue_c); 
  TMatrixD sigm_d = getSigma(11, aLL_stat_d, aLL_syst_d, ue_d); 

  TMatrixD sigm = addMatrixDiag(sigm_i, sigm_a, sigm_b, sigm_c, sigm_d);
  TMatrixD cov = sigm*corr*sigm;
  cov.Print();
  
  TGraph *gr = new TGraph(ndi, pt, asym);
  gr->Print("all");

  TGraph *grpt = new TGraph(ndi, pt, ptsyst);
  grpt->Print("all");

  //nnpdf aLL best fit
  TGraph *grnnpdf = new TGraph("nnpdf1.1.asymData.txt", "%lg %lg %*lg");
  int npg =  grnnpdf->GetN();
  grnnpdf->Print("all");
  //calculating uncertainty due to the uncertainty on x-axis
  TMatrixD ptsm = getPtSyst(grpt, grnnpdf);
  ptsm.Print();
  //di-jet
  TGraph *gr_dj[ntop], *gr_minv[ntop], *gr_nlo[ntop];
  //top a
  gr_dj[0] = new TGraph(ndj_a, minv_a, aLL_a); //gr_dj[0]->Print("all");
  gr_minv[0] = new TGraph(ndj_a, minv_a, minv_unc_a); //gr_minv[0]->Print("all");
  gr_nlo[0] = new TGraph(Form("dijet_nlo/dijet_all0_bin%d.txt", top[0]), "%lg %*lg %*lg %lg %*lg %*lg"); //gr_nlo[0]->Print("all");
  TMatrixD minvsm_a = getPtSyst(gr_minv[0], gr_nlo[0]); minvsm_a.Print();

  //top b
  gr_dj[1] = new TGraph(ndj_b, minv_b, aLL_b); //gr_dj[1]->Print("all");
  gr_minv[1] = new TGraph(ndj_b, minv_b, minv_unc_b); //gr_minv[1]->Print("all");
  gr_nlo[1] = new TGraph(Form("dijet_nlo/dijet_all0_bin%d.txt", top[1]), "%lg %*lg %*lg %lg %*lg %*lg"); //gr_nlo[1]->Print("all");
  TMatrixD minvsm_b = getPtSyst(gr_minv[1], gr_nlo[1]); minvsm_b.Print();

  //top c
  gr_dj[2] = new TGraph(ndj_c, minv_c, aLL_c); //gr_dj[2]->Print("all");
  gr_minv[2] = new TGraph(ndj_c, minv_c, minv_unc_c); //gr_minv[2]->Print("all");
  gr_nlo[2] = new TGraph(Form("dijet_nlo/dijet_all0_bin%d.txt", top[2]), "%lg %*lg %*lg %lg %*lg %*lg"); //gr_nlo[2]->Print("all");
  TMatrixD minvsm_c = getPtSyst(gr_minv[2], gr_nlo[2]); //minvsm_c.Print();

  //top d
  gr_dj[3] = new TGraph(ndj_d, minv_d, aLL_d); //gr_dj[3]->Print("all");
  gr_minv[3] = new TGraph(ndj_d, minv_d, minv_unc_d); //gr_minv[3]->Print("all");
  gr_nlo[3] = new TGraph(Form("dijet_nlo/dijet_all0_bin%d.txt", top[3]), "%lg %*lg %*lg %lg %*lg %*lg"); //gr_nlo[3]->Print("all");
  TMatrixD minvsm_d = getPtSyst(gr_minv[3], gr_nlo[3]); minvsm_d.Print();

  //add x-axis uncertainty to the diagonal entries of the covariance matrix
  TMatrixD covx = addDiag(cov, ptsm, minvsm_a, minvsm_b, minvsm_c, minvsm_d);
  covx.Print();
  //invert covariance matrix
  TMatrixD covi = covx; covi.Invert();
  covi.Print();

  Printf("replica irep = %d", irep);
  TGraph *grep = new TGraph(Form("replica_%dgr_replica_%dData.txt", irep, irep));

  //calculating difference between data and model replcia for inclusive jets
  TMatrixD dim = getDiff(gr, grep);

  TGraph *gdj_rep[ntop];
  for(int it = 0; it < ntop; it++){
    gdj_rep[it] = new TGraph(Form("dijet_nlo/dijet_all%d_bin%d.txt", irep, top[it]), "%lg %*lg %*lg %lg %*lg %*lg");
  }
  //calculating difference between data and model replcia for dijets
  TMatrixD dam = getDiff(gr_dj[0], gdj_rep[0]);
  TMatrixD dbm = getDiff(gr_dj[1], gdj_rep[1]);
  TMatrixD dcm = getDiff(gr_dj[2], gdj_rep[2]);
  TMatrixD ddm = getDiff(gr_dj[3], gdj_rep[3]);
  //chi2 without adding relative luminosity and beam polarization scale uncertainties 
  TMatrixD dm = addMatrix(dim, dam, dbm, dcm, ddm);
  TMatrixD dmt = transpose(dm);
  TMatrixD chi2m = dmt*covi*dm;
  chi2m.Print();
  //chi2 due to relative luminosity
  //vector beta_1
  TMatrixD rlim = getRL(gr); //rlim.Print();
  TMatrixD rlam = getRL(gr_dj[0]); //rlam.Print();
  TMatrixD rlbm = getRL(gr_dj[1]); //rlbm.Print();
  TMatrixD rlcm = getRL(gr_dj[2]); //rlcm.Print();
  TMatrixD rldm = getRL(gr_dj[3]); //rldm.Print();

  TMatrixD rlm = addMatrix(rlim, rlam, rlbm, rlcm, rldm);
  TMatrixD rlmt = transpose(rlm);
  TMatrixD a11m = rlmt*covi*rlm;

  //chi2 due to polarization
  //vector beta_2
  TMatrixD polim = getPol(gr, grep);
  
  TMatrixD polam = getPol(gr_dj[0], gdj_rep[0]);
  TMatrixD polbm = getPol(gr_dj[1], gdj_rep[1]);
  TMatrixD polcm = getPol(gr_dj[2], gdj_rep[2]);
  TMatrixD poldm = getPol(gr_dj[3], gdj_rep[3]);
  
  TMatrixD polm = addMatrix(polim, polam, polbm, polcm, poldm);
  TMatrixD polmt = transpose(polm);
  TMatrixD a22m = polmt*covi*polm;
  
  TMatrixD a12m = rlmt*covi*polm;
  
  TMatrixD a21m = polmt*covi*rlm;
  
  TMatrixD b1m = rlmt*covi*dm;

  TMatrixD b2m = polmt*covi*dm;
  //chi2 correction on the chi2 due to relative luminosity and beam polarization scale uncertainty
  TMatrixD chi2sm = getSystChi2(a11m, a12m, a21m, a22m, b1m, b2m);
  chi2sm.Print();

  double chi2 = chi2m(0,0); double chi2s = chi2sm(0,0);
  double chi2n = chi2-chi2s;
  Printf("chi2=%lg", chi2n);
  if(chi2n < 0) return chi2n;
  double ww = TMath::Power(chi2n, 0.5*nd-1)*TMath::Exp(-0.5*chi2n);
  Printf("ndof=%d", nd);
  Printf("weight=%.3lg", ww);
  return chi2n;
}
TMatrixD getSigma(int n, double stat[], double syst[], double ue[], double rl = 0.00022)
{
  //remove rl contribution from the RL/UE uncertainty presented in the paper
  TMatrixD m(n, n); m.Zero();

  for(int ii = 0; ii < n; ii++){
    double s1 = stat[ii], s2 = syst[ii], s3 = ue[ii];
    double sum = s1*s1+s2*s2+s3*s3-rl*rl;
    sum = TMath::Sqrt(sum);
    m(ii, ii) = sum;
  }
  return m;
}
TMatrixD addMatrixDiag(TMatrixD mi, TMatrixD ma, TMatrixD mb, TMatrixD mc, TMatrixD md)
{
  int ni = mi.GetNrows();
  int na = ma.GetNrows();
  int nb = mb.GetNrows();
  int nc = mc.GetNrows();
  int nd = md.GetNrows();
  int nm = ni+na+nb+nc+nd;
  TMatrixD ms(nm, nm); ms.Zero();
  
  for(int ii = 0; ii < nm; ii++){
    if(ii < ni)
      ms(ii,ii) += mi(ii, ii);
    else if(ii < ni+na)
      ms(ii, ii) += ma(ii-ni, ii-ni);
    else if(ii < ni+na+nb)
      ms(ii, ii) += mb(ii-ni-na, ii-ni-na);
    else if(ii < ni+na+nb+nc)
      ms(ii, ii) += mc(ii-ni-na-nb, ii-ni-na-nb);
    else
      ms(ii, ii) += md(ii-ni-na-nb-nc, ii-ni-na-nb-nc);
  }
  return ms;
}
TMatrixD addMatrix(TMatrixD mi, TMatrixD ma, TMatrixD mb, TMatrixD mc, TMatrixD md)
{
  int ni = mi.GetNrows();
  int na = ma.GetNrows();
  int nb = mb.GetNrows();
  int nc = mc.GetNrows();
  int nd = md.GetNrows();

  int nm = ni+na+nb+nc+nd;
  TMatrixD mm(nm, 1);
  
  for(int ii = 0; ii < nm; ii++){
    if(ii < ni)
      mm(ii, 0) += mi(ii, 0);
    else if(ii < ni+na)
      mm(ii, 0) += ma(ii-ni, 0);
    else if(ii < ni+na+nb)
      mm(ii, 0) += mb(ii-ni-na, 0);
    else if(ii < ni+na+nb+nc)
      mm(ii, 0) += mc(ii-ni-na-nb, 0);
    else
      mm(ii, 0) += md(ii-ni-na-nb-nc, 0);
  }
  return mm;
}
TMatrixD addDiag(TMatrixD mm, TMatrixD mi, TMatrixD ma, TMatrixD mb, TMatrixD mc, TMatrixD md)
{
  int nm = mm.GetNrows();
  int ni = mi.GetNrows();
  int na = ma.GetNrows();
  int nb = mb.GetNrows();
  int nc = mc.GetNrows();
  int nd = md.GetNrows();
  TMatrixD ms(mm);
  if(nm != ni+na+nb+nc+nd) return TMatrix(1,1);
  
  for(int ii = 0; ii < nm; ii++){
    if(ii < ni)
      ms(ii,ii) += mi(ii, ii);
    else if(ii < ni+na)
      ms(ii, ii) += ma(ii-ni, ii-ni);
    else if(ii < ni+na+nb)
      ms(ii, ii) += mb(ii-ni-na, ii-ni-na);
    else if(ii < ni+na+nb+nc)
      ms(ii, ii) += mc(ii-ni-na-nb, ii-ni-na-nb);
    else
      ms(ii, ii) += md(ii-ni-na-nb-nc, ii-ni-na-nb-nc);
  }
  return ms;
}
TMatrixD getPtSyst(TGraph *gr, TGraph *gf)
{
  int np = gr->GetN();
  TMatrixD m(np, np);
  m.Zero();
  double xmin = 5, xmax = 60;

  Grfun *gfun = new Grfun(gf);
  TF1 *fg = new TF1("fg", gfun, xmin, xmax, 0, "Grfun"); 

  for(int ip = 0; ip < np; ip++){
    double xx, yy;
    gr->GetPoint(ip, xx, yy);
    double ep = 1./(xmax-xmin);
    double df = fg->Derivative(xx, 0, ep);

    double e2 = TMath::Power(yy*df, 2.);
    m(ip, ip) = e2; 
  }
  //delete fg; delete gfun;
  return m;
}
TMatrixD getSystChi2(TMatrixD a11, TMatrixD a12, TMatrixD a21, TMatrixD a22, TMatrixD b1, TMatrixD b2)
{
  //calculating inverse of A and B^T*A^-1*B
  TMatrixD aa(2,2);
  double m1 = 1.+a11(0,0); double m2 = a12(0,0); double m3 = a21(0,0); double m4 = 1.+a22(0,0);
  double det = m1*m4-m2*m3;
  aa(0,0) = m4/det; aa(0,1) = -m2/det; aa(1,0) = -m3/det; aa(1,1) = m1/det;
  aa.Print();
  TMatrixD bb(2,1);
  bb(0,0) = b1(0,0); bb(1,0) = b2(0,0);
  bb.Print();
  TMatrixD bbt = transpose(bb);
  
  return bbt*aa*bb;
}
TMatrixD fillMatrix(const char* name, int nn)
{
  TMatrixD m(nn, nn);

  FILE *fp = fopen(name, "r");
  double me;
  for(int ix = 0; ix < nn; ix++){
    for(int iy = 0; iy < nn; iy++){
      fscanf(fp, "%lf", &me);
      m(ix, iy) = me;
    }
  }
  fclose(fp);

  return m;
}
TMatrixD getRL(TGraph *gra, double rl = 0.00022)
{
  int na = gra->GetN();
  TMatrixD m(na, 1);
  for(int ip = 0; ip < na; ip++){
    m(ip, 0) = rl;
  }
  return m;
}
TMatrixD getPol(TGraph *gra, TGraph *grb, double scale = 0.066)
{
  int na = gra->GetN();
  TMatrixD m(na, 1);
  for(int ip = 0; ip < na; ip++){
    double xx, yy;
    gra->GetPoint(ip, xx, yy);
    double yv = grb->Eval(xx, 0, "S");
    //double yv = grb->Eval(xx);
    //Printf("pol yv = %.3lg yy = %.3lg", yv, yy);
    double pol = yv*scale;
    m(ip, 0) = pol;
  }
  return m;
}
TMatrixD getDiff(TGraph *gra, TGraph *grb)
{
  int na = gra->GetN();
  TMatrixD m(na, 1);
  for(int ip = 0; ip < na; ip++){
    double xx, yy;
    gra->GetPoint(ip, xx, yy);
    double yv = grb->Eval(xx, 0, "S");
    //double yv = grb->Eval(xx);
    //Printf("diff yv = %.3lg yy = %.3lg", yv, yy);
    double dd = yy - yv;
    m(ip, 0) = dd;
  }
  return m;
}
TMatrixD transpose(TMatrixD  vv)
{
  int nr = vv.GetNrows();
  TMatrixD m(1, nr);
  m.Zero();
  for(int ir = 0; ir < nr; ir++){
    double vr = vv(ir, 0);
    m(0, ir) = vr;
  }
  return m;
}
