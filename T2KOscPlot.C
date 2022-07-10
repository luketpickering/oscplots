#include "/root/scratch/T2K/oscplots/Prob3-Wrapper/Prob3ppWrapper.hxx"

#include "TColor.h"
#include "colordef.h"
#include "plotutils.h"

#include <array>
#include <limits>

struct OFlux {
  OFlux(TH1 *flux, TH1 *fluxb, double min, double max) {
    TSpline3 spl(flux, 0, min, max);
    int nsteps = 1000;
    gflux = TGraph(nsteps + 1);
    gflux.SetPoint(0, min, 0);
    for (int i = 1; i < nsteps; ++i) {
      double e = min + i * ((max - min) / double(nsteps));
      gflux.SetPoint(i, e, spl.Eval(e));
    }
    gflux.SetPoint(nsteps, max, 0);

    TSpline3 splb(fluxb, 0, min, max);
    gfluxb = TGraph(nsteps + 1);
    gfluxb.SetPoint(0, min, 0);
    for (int i = 1; i < nsteps; ++i) {
      double e = min + i * ((max - min) / double(nsteps));
      gfluxb.SetPoint(i, e, spl.Eval(e));
    }
    gfluxb.SetPoint(nsteps, max, 0);
  }

  void AddOsc(OscillationHelper &oh, bool b = false) {
    int nsteps = 1000;
    gosc.emplace_back(nsteps + 1);
    oscIsB.push_back(b);
    for (int i = 0; i < nsteps + 1; ++i) {
      double e, f;
      (b ? gfluxb : gflux).GetPoint(i, e, f);

      gosc.back().SetPoint(i, e, f * oh.GetWeight(e));
    }
  }

  void AddXSec(TGraph &g, TGraph &gb) {
    int nsteps = 1000;
    for (size_t i = 0; i < gosc.size(); ++i) {
      auto &go = gosc[i];
      gevrate.emplace_back(nsteps + 1);
      bool isb = oscIsB[i];
      for (int i = 0; i < nsteps + 1; ++i) {
        double e, f;
        go.GetPoint(i, e, f);

        double ex, x;
        (isb ? gb : g).GetPoint(i, ex, x);

        gevrate.back().SetPoint(i, e, f * x * ex);
      }
    }
  }
  TGraph gflux;
  TGraph gfluxb;
  std::vector<TGraph> gosc;
  std::vector<bool> oscIsB;
  std::vector<TGraph> gevrate;
};

struct XSecs {
  XSecs(std::initializer_list<TH1 *> xsecs, double min, double max) {
    for (auto h : xsecs) {
      TSpline3 spl(h, 0, min, max);
      int nsteps = 1000;
      gxsecs.emplace_back(nsteps);
      for (int i = nsteps; i >= 0; --i) {
        double e = min + i * ((max - min) / double(nsteps));
        if (spl.Eval(e) < 0) {
          break;
        }
        gxsecs.back().SetPoint(i, e, spl.Eval(e));
      }
    }
  }

  void ForceInc(size_t inc) {
    auto &ginc = gxsecs[inc];
    for (size_t i = 0; i < ginc.GetN(); ++i) {
      double e, incx;
      ginc.GetPoint(i, e, incx);
      double sum = 0;
      for (size_t x = 0; x < gxsecs.size(); ++x) {
        if (x == inc) {
          continue;
        }
        double e, xs;
        gxsecs[x].GetPoint(i, e, xs);
        sum += xs;
      }
      if (sum > incx) {
        std::cout << e << " incx: " << incx << ", sum: " << sum << std::endl;
        ginc.SetPoint(i, e, sum);
      }
    }
  }
  std::vector<TGraph> gxsecs;
};

double GetMaximum(TGraph const &g) {
  double m = -std::numeric_limits<double>::max();
  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    m = std::max(m, y);
  }
  return m;
}

double GetMinimum(TGraph const &g) {
  double m = std::numeric_limits<double>::max();
  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    m = std::min(m, y);
  }
  return m;
}

double GetMaximumX(TGraph const &g) {
  double m = -std::numeric_limits<double>::max();
  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    m = std::max(m, x);
  }
  return m;
}

double GetMinimumX(TGraph const &g) {
  double m = std::numeric_limits<double>::max();
  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    m = std::min(m, x);
  }
  return m;
}

TGraph DummyGraph(double xmin, double xmax, double ymin, double ymax) {
  TGraph dg(4);

  dg.SetPoint(0, xmin, ymin);
  dg.SetPoint(1, xmin, ymax);
  dg.SetPoint(2, xmax, ymax);
  dg.SetPoint(3, xmax, ymin);

  dg.SetLineColorAlpha(0, 0);
  return dg;
}

TGraph Scale(TGraph const &g, double s) {
  TGraph o(g.GetN());

  for (int i = 0; i < g.GetN(); ++i) {
    double x, y;
    g.GetPoint(i, x, y);
    o.SetPoint(i, x, y * s);
  }
  return o;
}

TGraph Clone(TGraph const &g) { return Scale(g, 1); }

double Integrate(TGraph const &g) {
  double ig = 0;
  for (int i = 1; i < g.GetN(); ++i) {
    double x0, y0, x1, y1;
    g.GetPoint(i - 1, x0, y0);
    g.GetPoint(i, x1, y1);

    ig += ((x1 - x0) * (y0)) + ((x1 - x0) * (y1 - y0) / 2.0);
  }
  return ig;
}

TGraph DrawEllipse(std::array<double, 6> osc_params, TGraph gf,
                   TGraph const &xs, TGraph gfb, TGraph const &xsb,
                   double baseline, bool use_xs = true) {

  int nsteps = 100;
  TGraph ell(nsteps);
  for (int dcpi = 0; dcpi < nsteps; ++dcpi) {
    double dcp = dcpi * (2.0 * M_PI / double(nsteps));
    std::cout << dcp << std::endl;

    osc_params[5] = dcp;

    OscillationHelper osc, oscb;
    osc.Setup_baseline(osc_params.data(), baseline);
    osc.SetOscillationChannel(14, 12);

    oscb.Setup_baseline(osc_params.data(), baseline);
    oscb.SetOscillationChannel(-14, -12);

    TGraph gosc = gf;
    TGraph goscb = gfb;

    for (int i = 0; i < gf.GetN(); ++i) {
      double e, x, eb, xb;
      xs.GetPoint(i, e, x);
      xsb.GetPoint(i, eb, xb);
      if (!use_xs) {
        x = 1;
        xb = 1;
      }

      double f, fb;
      gf.GetPoint(i, e, f);
      gfb.GetPoint(i, eb, fb);

      double w = osc.GetWeight(e);
      double wb = oscb.GetWeight(e);

      gosc.SetPoint(i, e, f * x * w);
      goscb.SetPoint(i, eb, fb * xb * wb);
    }
    ell.SetPoint(dcpi, Integrate(gosc), Integrate(goscb));

    if ((dcpi + 1) == nsteps) {
      ell.SetPoint(0, Integrate(gosc), Integrate(goscb));
    }
  }
  return ell;
}

void AxisStyle(TGraph &g) {
  g.SetTitle("");
  g.GetYaxis()->SetLabelSize(0.06);
  g.GetXaxis()->SetLabelSize(0.06);

  g.GetYaxis()->SetTitleSize(0.06);
  g.GetXaxis()->SetTitleSize(0.06);

  g.GetYaxis()->SetNdivisions(505);
  g.GetXaxis()->SetNdivisions(505);

  g.GetYaxis()->SetTitleOffset(0.8);
  g.GetXaxis()->SetTitleOffset(0.9);
}

std::string GetString(std::array<double, 6> osc_params, size_t param) {
  std::stringstream ss("");
  ss.precision(3);
  switch (param) {
  case 0: {
    ss << "\\sin^{2}(\\theta_{12})=" << osc_params[0];
    break;
  }
  case 1: {
    ss << "\\sin^{2}(\\theta_{13})=" << osc_params[1];
    break;
  }
  case 2: {
    ss << "\\sin^{2}(\\theta_{23})=" << osc_params[2];
    break;
  }
  case 3: {
    ss << "#Deltam_{12}^{2}=" << osc_params[3] * 1.0E5 << "#times10^{-5} eV";
    break;
  }
  case 4: {
    ss << "#Deltam_{32}^{2}=" << osc_params[4] * 1.0E3 << "#times10^{-3} eV";
    break;
  }
  case 5: {
    double dcpopi = osc_params[5] / M_PI;
    if (fabs(dcpopi - 0.5) < 1E-5) {
      ss << "#delta_{CP}=#pi/2";
    } else if (fabs(dcpopi - 1) < 1E-5) {
      ss << "#delta_{CP}=#pi";
    } else if (fabs(dcpopi - 1.5) < 1E-5) {
      ss << "#delta_{CP}=3#pi/2";
    } else {
      ss << "#delta_{CP}=" << osc_params[5];
    }
  }
  }
  return ss.str();
}

void DrawOscParameters(std::array<double, 6> osc_params, int baseline = 0,
                       double ymax = 0.75) {

  TLatex ltx;
  ltx.SetTextAlign(12);

  double ydelta = 0.065;

  ltx.DrawLatexNDC(0.5, ymax - 0 * ydelta, GetString(osc_params, 0).c_str());
  ltx.DrawLatexNDC(0.5, ymax - 1 * ydelta, GetString(osc_params, 1).c_str());
  ltx.DrawLatexNDC(0.5, ymax - 2 * ydelta, GetString(osc_params, 2).c_str());
  ltx.DrawLatexNDC(0.2, ymax - 0 * ydelta, GetString(osc_params, 3).c_str());
  ltx.DrawLatexNDC(0.2, ymax - 1 * ydelta, GetString(osc_params, 4).c_str());
  ltx.DrawLatexNDC(0.75, ymax - 2 * ydelta, GetString(osc_params, 5).c_str());

  if (baseline) {
    std::stringstream ss("");
    ss << "L=" << baseline << " km";
    ltx.DrawLatexNDC(0.75, ymax - 0 * ydelta, ss.str().c_str());
  }
}

size_t disp(size_t i) { return i * 4 + 0; }
size_t app(size_t i) { return i * 4 + 1; }
size_t dispb(size_t i) { return i * 4 + 2; }
size_t appb(size_t i) { return i * 4 + 3; }

size_t disp(size_t i, size_t x, size_t n) { return x * n * 4 + i * 4 + 0; }
size_t app(size_t i, size_t x, size_t n) { return x * n * 4 + i * 4 + 1; }
size_t dispb(size_t i, size_t x, size_t n) { return x * n * 4 + i * 4 + 2; }
size_t appb(size_t i, size_t x, size_t n) { return x * n * 4 + i * 4 + 3; }

void T2KOscPlot() {

  int kT2KRed;
  int kT2KGreen;
  int kNOvABlue;

  defc(kNOvABlue, {34, 62, 145}, "kNOvABlue");
  defc(kT2KGreen, {0, 128, 0}, "kT2KGreen");
  defc(kT2KRed, {128, 0, 0}, "kT2KRed");

  int pRed;
  int pOrange;
  int pYellow;
  int pGreen;
  int pBlue;
  int pPurple;
  int pBlack;
  int pBrown;
  int pCyan;
  int pPink;
  int pWhite;
  int ncol = 10;
  int pOffset = 9;

  defc(pRed, 0xff191f, "pRed");
  defc(pOrange, 0xff8b11, "pOrange");
  defc(pYellow, 0xfcff00, "pYellow");
  defc(pGreen, 0x00b600, "pGreen");
  defc(pBlue, 0x120790, "pBlue");
  defc(pPurple, 0x71007d, "pPurple");
  defc(pBlack, 0x000000, "pBlack");
  defc(pBrown, 0x582d0b, "pBrown");
  defc(pCyan, 0x57cfe7, "pCyan");
  defc(pPink, 0xffa6c2, "pPink");
  defc(pWhite, 0xffffff, "pWhite");

  TFile ft2k("fluxes/tuned13av3/run1-9a_v2/"
             "sk_tuned13av3_13anom_run1-9a_numode_fine.root");
  TH1 *ht2k = nullptr;
  ft2k.GetObject("enu_sk_tuned13a_numu", ht2k);
  ht2k->RebinX(2);

  TFile ft2kb("fluxes/tuned13av3/run5c-9d_v2/"
              "sk_tuned13av3_13anom_run5c-9d_antinumode_fine.root");
  TH1 *ht2kb = nullptr;
  ft2kb.GetObject("enu_sk_tuned13a_numu", ht2kb);
  ht2kb->RebinX(2);

  int rebin = 5;

  TFile fxsec("neut/neutsigmaenu.root");
  TH1 *hxsec_ccinc = nullptr;
  fxsec.GetObject("CCInc", hxsec_ccinc);
  hxsec_ccinc->RebinX(rebin);
  hxsec_ccinc->Scale(1.0E39 / double(rebin));

  TH1 *hxsec_0pi = nullptr;
  fxsec.GetObject("CC0Pi", hxsec_0pi);
  hxsec_0pi->RebinX(rebin);
  hxsec_0pi->Scale(1.0E39 / double(rebin));

  TH1 *hxsec_1pi = nullptr;
  fxsec.GetObject("CC1Pi", hxsec_1pi);
  hxsec_1pi->RebinX(rebin);
  hxsec_1pi->Scale(1.0E39 / double(rebin));

  TH1 *hxsec_npi = nullptr;
  fxsec.GetObject("CCNPi", hxsec_npi);
  hxsec_npi->RebinX(rebin);
  hxsec_npi->Scale(1.0E39 / double(rebin));

  TFile fxsecb("neut/neutsigmaenu.nub.root");
  TH1 *hxsecb_ccinc = nullptr;
  fxsecb.GetObject("CCInc", hxsecb_ccinc);
  hxsecb_ccinc->RebinX(rebin);
  hxsecb_ccinc->Scale(1.0E39 / double(rebin));

  TH1 *hxsecb_0pi = nullptr;
  fxsecb.GetObject("CC0Pi", hxsecb_0pi);
  hxsecb_0pi->RebinX(rebin);
  hxsecb_0pi->Scale(1.0E39 / double(rebin));

  TH1 *hxsecb_1pi = nullptr;
  fxsecb.GetObject("CC1Pi", hxsecb_1pi);
  hxsecb_1pi->RebinX(rebin);
  hxsecb_1pi->Scale(1.0E39 / double(rebin));

  TH1 *hxsecb_npi = nullptr;
  fxsecb.GetObject("CCNPi", hxsecb_npi);
  hxsecb_npi->RebinX(rebin);
  hxsecb_npi->Scale(1.0E39 / double(rebin));

  double min = 0.15;
  double max = 4;

  size_t kCCInc = 0;
  size_t kCC0Pi = 1;
  size_t kCC1Pi = 2;
  size_t kCCNPi = 3;

  size_t NXSec = 4;

  XSecs xsgraphs({hxsec_ccinc, hxsec_0pi, hxsec_1pi, hxsec_npi}, 0, 5);
  xsgraphs.ForceInc(0);

  XSecs xsgraphsb({hxsecb_ccinc, hxsecb_0pi, hxsecb_1pi, hxsecb_npi}, 0, 5);
  xsgraphsb.ForceInc(0);

  OFlux T2KFlux(ht2k, ht2kb, 0.15, 4);

  double t2k_baseline = 295;

  // Sin^2(Theta_12)
  double s2th12 = 0.297;
  // Sin^2(Theta_13)
  double s2th13 = 0.0214;
  // Sin^2(Theta_23)
  double s2th23 = 0.56;
  // Dm^2_21
  double dm2_21 = 7.37E-5;
  //|Dm^2_Atm|
  double dm2_atm = 2.49E-3;
  // dcp
  double dcp = 0;

  double dm2_atm_t2kp = 2.59E-3;
  double dm2_atm_t2km = 2.38E-3;
  double s2th23_t2km = 0.44;
  double s2th23_t2kp = 0.6;

  size_t kT2KBF = 0;
  size_t kT2KBF_dm32p = 1;
  size_t kT2KBF_dm32m = 2;

  size_t kT2KBF_th23p = 3;
  size_t kT2KBF_th23m = 4;

  size_t kT2KBF_dcp0 = 5;
  size_t kT2KBF_dcp05 = 6;
  size_t kT2KBF_dcp1 = 7;
  size_t kT2KBF_dcp15 = 8;

  size_t kT2KBF_IH = 9;

  size_t kT2KBF_IH_dcp0 = 10;
  size_t kT2KBF_IH_dcp05 = 11;
  size_t kT2KBF_IH_dcp1 = 12;
  size_t kT2KBF_IH_dcp15 = 13;

  std::vector<std::array<double, 6>> osc_params;
  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, dm2_atm, dcp});

  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, dm2_atm_t2kp, dcp});
  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, dm2_atm_t2km, dcp});

  osc_params.push_back({s2th12, s2th13, s2th23_t2km, dm2_21, dm2_atm, dcp});
  osc_params.push_back({s2th12, s2th13, s2th23_t2kp, dm2_21, dm2_atm, dcp});

  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, dm2_atm, 0});
  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, dm2_atm, M_PI / 2.0});
  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, dm2_atm, M_PI});
  osc_params.push_back(
      {s2th12, s2th13, s2th23, dm2_21, dm2_atm, 3.0 * M_PI / 2.0});

  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, -dm2_atm, dcp});

  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, -dm2_atm, 0});
  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, -dm2_atm, M_PI / 2.0});
  osc_params.push_back({s2th12, s2th13, s2th23, dm2_21, -dm2_atm, M_PI});
  osc_params.push_back(
      {s2th12, s2th13, s2th23, dm2_21, -dm2_atm, 3.0 * M_PI / 2.0});

  std::vector<OscillationHelper> T2KOsc;
  for (auto &op : osc_params) {
    for (auto osc : std::vector<std::pair<int, int>>{
             {14, 14}, {14, 12}, {-14, -14}, {-14, -12}}) {
      T2KOsc.emplace_back();
      T2KOsc.back().Setup_baseline(op.data(), t2k_baseline);
      T2KOsc.back().SetOscillationChannel(osc.first, osc.second);
      T2KFlux.AddOsc(T2KOsc.back(), osc.first < 0);
    }
  }

  size_t NOsc = osc_params.size();

  for (size_t i = 0; i < xsgraphs.gxsecs.size(); ++i) {
    T2KFlux.AddXSec(xsgraphs.gxsecs[i], xsgraphsb.gxsecs[i]);
  }

  TCanvas c1("c1", "", 1400, 1000);

  // T2K CV Disp
  {
    TCanvas c1("c1", "", 1400, 1000);

    double scale = GetMaximum(T2KFlux.gevrate[disp(kT2KBF, kCC0Pi, NOsc)]);
    auto gt2k_dmc =
        Scale(T2KFlux.gevrate[disp(kT2KBF, kCC0Pi, NOsc)], 1.0 / scale);

    gt2k_dmc.GetXaxis()->SetRangeUser(min, 1.25);
    gt2k_dmc.GetYaxis()->SetRangeUser(0, 1.6);
    gt2k_dmc.GetYaxis()->SetTitle("Events (Arbitrary Units)");
    gt2k_dmc.GetXaxis()->SetTitle("E_{#nu} (GeV)");

    AxisStyle(gt2k_dmc);

    c1.SetBottomMargin(0.125);
    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetLeftMargin(0.125);

    gt2k_dmc.SetLineColor(pBrown);
    gt2k_dmc.SetLineWidth(5);

    gt2k_dmc.Draw("AC");

    TLegend leg(0.475, 0.15 + 0.525, 0.95, 0.425 + 0.525);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&gt2k_dmc, GetString(osc_params[kT2KBF], 4).c_str(), "l");

    leg.Draw();

    c1.Print("T2KOscXSecPlot_cvdisp.png");
  }

  // T2K Dm23 Movement
  {
    TCanvas c1("c1", "", 1400, 1000);

    double scale = GetMaximum(T2KFlux.gevrate[disp(kT2KBF, kCC0Pi, NOsc)]);
    auto gt2k_dmc =
        Scale(T2KFlux.gevrate[disp(kT2KBF, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dmp =
        Scale(T2KFlux.gevrate[disp(kT2KBF_dm32p, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dmm =
        Scale(T2KFlux.gevrate[disp(kT2KBF_dm32m, kCC0Pi, NOsc)], 1.0 / scale);

    gt2k_dmc.GetXaxis()->SetRangeUser(min, 1.25);
    gt2k_dmc.GetYaxis()->SetRangeUser(0, 1.6);
    gt2k_dmc.GetYaxis()->SetTitle("Events (Arbitrary Units)");
    gt2k_dmc.GetXaxis()->SetTitle("E_{#nu} (GeV)");

    AxisStyle(gt2k_dmc);

    c1.SetBottomMargin(0.125);
    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetLeftMargin(0.125);

    gt2k_dmc.SetLineColor(pBrown);
    gt2k_dmc.SetLineWidth(5);
    gt2k_dmp.SetLineColor(pCyan);
    gt2k_dmp.SetLineWidth(5);
    gt2k_dmm.SetLineColor(pPink);
    gt2k_dmm.SetLineWidth(5);

    gt2k_dmc.Draw("AC");
    gt2k_dmp.Draw("C");
    gt2k_dmm.Draw("C");

    TLegend leg(0.475, 0.15 + 0.525, 0.95, 0.425 + 0.525);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&gt2k_dmp, GetString(osc_params[kT2KBF_dm32p], 4).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dmc, GetString(osc_params[kT2KBF], 4).c_str(), "l");

    leg.AddEntry(&gt2k_dmm, GetString(osc_params[kT2KBF_dm32m], 4).c_str(),
                 "l");

    leg.Draw();

    c1.Print("T2KOscXSecPlot_dm32.png");
  }

  // T2K th23 Movement
  {
    TCanvas c1("c1", "", 1400, 1000);

    double scale = GetMaximum(T2KFlux.gevrate[disp(kT2KBF, kCC0Pi, NOsc)]);
    auto gt2k_dmc =
        Scale(T2KFlux.gevrate[disp(kT2KBF, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dmp =
        Scale(T2KFlux.gevrate[disp(kT2KBF_th23p, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dmm =
        Scale(T2KFlux.gevrate[disp(kT2KBF_th23m, kCC0Pi, NOsc)], 1.0 / scale);

    gt2k_dmc.GetXaxis()->SetRangeUser(min, 1.25);
    gt2k_dmc.GetYaxis()->SetRangeUser(0, 1.6);
    gt2k_dmc.GetYaxis()->SetTitle("Events (Arbitrary Units)");
    gt2k_dmc.GetXaxis()->SetTitle("E_{#nu} (GeV)");

    AxisStyle(gt2k_dmc);

    c1.SetBottomMargin(0.125);
    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetLeftMargin(0.125);

    gt2k_dmc.SetLineColor(pBrown);
    gt2k_dmc.SetLineWidth(5);
    gt2k_dmp.SetLineColor(pPink);
    gt2k_dmp.SetLineWidth(5);
    gt2k_dmm.SetLineColor(pCyan);
    gt2k_dmm.SetLineWidth(5);

    gt2k_dmc.Draw("AC");
    gt2k_dmp.Draw("C");
    gt2k_dmm.Draw("C");

    TLegend leg(0.475, 0.15 + 0.525, 0.95, 0.425 + 0.525);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&gt2k_dmm, GetString(osc_params[kT2KBF_th23m], 2).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dmc, GetString(osc_params[kT2KBF], 2).c_str(), "l");
    leg.AddEntry(&gt2k_dmp, GetString(osc_params[kT2KBF_th23p], 2).c_str(),
                 "l");

    leg.Draw();

    c1.Print("T2KOscXSecPlot_th23.png");
  }

  // T2K dcp app Movement
  {
    TCanvas c1("c1", "", 1400, 1000);

    double scale = GetMaximum(T2KFlux.gevrate[app(kT2KBF, kCC0Pi, NOsc)]);
    auto gt2k_dcp0 =
        Scale(T2KFlux.gevrate[app(kT2KBF_dcp0, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp05 =
        Scale(T2KFlux.gevrate[app(kT2KBF_dcp05, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp1 =
        Scale(T2KFlux.gevrate[app(kT2KBF_dcp1, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp15 =
        Scale(T2KFlux.gevrate[app(kT2KBF_dcp15, kCC0Pi, NOsc)], 1.0 / scale);

    gt2k_dcp0.GetXaxis()->SetRangeUser(min, 1.25);
    gt2k_dcp0.GetYaxis()->SetRangeUser(0, 1.5 * GetMaximum(gt2k_dcp0));
    gt2k_dcp0.GetYaxis()->SetTitle("Events (Arbitrary Units)");
    gt2k_dcp0.GetXaxis()->SetTitle("E_{#nu} (GeV)");

    AxisStyle(gt2k_dcp0);

    c1.SetBottomMargin(0.125);
    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetLeftMargin(0.125);

    gt2k_dcp0.SetLineColor(pBrown);
    gt2k_dcp0.SetLineWidth(5);
    gt2k_dcp05.SetLineColor(pCyan);
    gt2k_dcp05.SetLineWidth(5);
    gt2k_dcp1.SetLineColor(pPink);
    gt2k_dcp1.SetLineWidth(5);
    gt2k_dcp15.SetLineColor(pGreen);
    gt2k_dcp15.SetLineWidth(5);

    gt2k_dcp0.Draw("AC");
    gt2k_dcp05.Draw("C");
    gt2k_dcp1.Draw("C");
    gt2k_dcp15.Draw("C");

    TLegend leg(0.7, 0.4, 0.95, 0.975);
    leg.SetHeader("#nu_{e}, T2K Best Fit, L=295 km");
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&gt2k_dcp0, GetString(osc_params[kT2KBF_dcp0], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp05, GetString(osc_params[kT2KBF_dcp05], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp1, GetString(osc_params[kT2KBF_dcp1], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp15, GetString(osc_params[kT2KBF_dcp15], 5).c_str(),
                 "l");

    TLegendEntry *header =
        static_cast<TLegendEntry *>(leg.GetListOfPrimitives()->First());
    header->SetTextAlign(32);

    leg.Draw();

    c1.Print("T2KOscXSecPlot_app.png");
  }

  // T2K dcp appb Movement
  {
    TCanvas c1("c1", "", 1400, 1000);

    double scale = GetMaximum(T2KFlux.gevrate[app(kT2KBF, kCC0Pi, NOsc)]);
    auto gt2k_dcp0 =
        Scale(T2KFlux.gevrate[appb(kT2KBF_dcp0, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp05 =
        Scale(T2KFlux.gevrate[appb(kT2KBF_dcp05, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp1 =
        Scale(T2KFlux.gevrate[appb(kT2KBF_dcp1, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp15 =
        Scale(T2KFlux.gevrate[appb(kT2KBF_dcp15, kCC0Pi, NOsc)], 1.0 / scale);

    gt2k_dcp0.GetXaxis()->SetRangeUser(min, 1.25);
    gt2k_dcp0.GetYaxis()->SetRangeUser(0, 1.5 * GetMaximum(gt2k_dcp0));
    gt2k_dcp0.GetYaxis()->SetTitle("Events (Arbitrary Units)");
    gt2k_dcp0.GetXaxis()->SetTitle("E_{#nu} (GeV)");

    AxisStyle(gt2k_dcp0);

    c1.SetBottomMargin(0.125);
    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetLeftMargin(0.125);

    gt2k_dcp0.SetLineColor(pBrown);
    gt2k_dcp0.SetLineWidth(5);
    gt2k_dcp05.SetLineColor(pCyan);
    gt2k_dcp05.SetLineWidth(5);
    gt2k_dcp1.SetLineColor(pPink);
    gt2k_dcp1.SetLineWidth(5);
    gt2k_dcp15.SetLineColor(pGreen);
    gt2k_dcp15.SetLineWidth(5);

    gt2k_dcp0.Draw("AC");
    gt2k_dcp05.Draw("C");
    gt2k_dcp1.Draw("C");
    gt2k_dcp15.Draw("C");

    TLegend leg(0.7, 0.4, 0.95, 0.975);
    leg.SetHeader("#bar{#nu}_{e}, T2K Best Fit, L=295 km");
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&gt2k_dcp0, GetString(osc_params[kT2KBF_dcp0], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp05, GetString(osc_params[kT2KBF_dcp05], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp1, GetString(osc_params[kT2KBF_dcp1], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp15, GetString(osc_params[kT2KBF_dcp15], 5).c_str(),
                 "l");

    TLegendEntry *header =
        static_cast<TLegendEntry *>(leg.GetListOfPrimitives()->First());
    header->SetTextAlign(32);

    leg.Draw();

    c1.Print("T2KOscXSecPlot_appb.png");
  }

  // T2K dcp app Movement IH
  {
    TCanvas c1("c1", "", 1400, 1000);

    double scale = GetMaximum(T2KFlux.gevrate[app(kT2KBF_IH, kCC0Pi, NOsc)]);
    auto gt2k_dcp0 =
        Scale(T2KFlux.gevrate[app(kT2KBF_IH_dcp0, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp05 =
        Scale(T2KFlux.gevrate[app(kT2KBF_IH_dcp05, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp1 =
        Scale(T2KFlux.gevrate[app(kT2KBF_IH_dcp1, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp15 =
        Scale(T2KFlux.gevrate[app(kT2KBF_IH_dcp15, kCC0Pi, NOsc)], 1.0 / scale);

    gt2k_dcp0.GetXaxis()->SetRangeUser(min, 1.25);
    gt2k_dcp0.GetYaxis()->SetRangeUser(0, 1.5 * GetMaximum(gt2k_dcp0));
    gt2k_dcp0.GetYaxis()->SetTitle("Events (Arbitrary Units)");
    gt2k_dcp0.GetXaxis()->SetTitle("E_{#nu} (GeV)");

    AxisStyle(gt2k_dcp0);

    c1.SetBottomMargin(0.125);
    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetLeftMargin(0.125);

    gt2k_dcp0.SetLineColor(pBrown);
    gt2k_dcp0.SetLineWidth(5);
    gt2k_dcp0.SetLineStyle(2);
    gt2k_dcp05.SetLineColor(pCyan);
    gt2k_dcp05.SetLineWidth(5);
    gt2k_dcp05.SetLineStyle(2);
    gt2k_dcp1.SetLineColor(pPink);
    gt2k_dcp1.SetLineWidth(5);
    gt2k_dcp1.SetLineStyle(2);
    gt2k_dcp15.SetLineColor(pGreen);
    gt2k_dcp15.SetLineWidth(5);
    gt2k_dcp15.SetLineStyle(2);

    gt2k_dcp0.Draw("AC");
    gt2k_dcp05.Draw("C");
    gt2k_dcp1.Draw("C");
    gt2k_dcp15.Draw("C");

    TLegend leg(0.7, 0.4, 0.95, 0.975);
    leg.SetHeader("#nu_{e}, T2K Best Fit, L=295 km");
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&gt2k_dcp0, GetString(osc_params[kT2KBF_IH_dcp0], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp05, GetString(osc_params[kT2KBF_IH_dcp05], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp1, GetString(osc_params[kT2KBF_IH_dcp1], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp15, GetString(osc_params[kT2KBF_IH_dcp15], 5).c_str(),
                 "l");

    TLegendEntry *header =
        static_cast<TLegendEntry *>(leg.GetListOfPrimitives()->First());
    header->SetTextAlign(32);

    leg.Draw();

    c1.Print("T2KOscXSecPlot_appIH.png");
  }

  // T2K dcp appb Movement IH
  {
    TCanvas c1("c1", "", 1400, 1000);

    double scale = GetMaximum(T2KFlux.gevrate[app(kT2KBF_IH, kCC0Pi, NOsc)]);
    auto gt2k_dcp0 =
        Scale(T2KFlux.gevrate[appb(kT2KBF_IH_dcp0, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp05 = Scale(
        T2KFlux.gevrate[appb(kT2KBF_IH_dcp05, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp1 =
        Scale(T2KFlux.gevrate[appb(kT2KBF_IH_dcp1, kCC0Pi, NOsc)], 1.0 / scale);
    auto gt2k_dcp15 = Scale(
        T2KFlux.gevrate[appb(kT2KBF_IH_dcp15, kCC0Pi, NOsc)], 1.0 / scale);

    gt2k_dcp0.GetXaxis()->SetRangeUser(min, 1.25);
    gt2k_dcp0.GetYaxis()->SetRangeUser(0, 1.5 * GetMaximum(gt2k_dcp0));
    gt2k_dcp0.GetYaxis()->SetTitle("Events (Arbitrary Units)");
    gt2k_dcp0.GetXaxis()->SetTitle("E_{#nu} (GeV)");

    AxisStyle(gt2k_dcp0);

    c1.SetBottomMargin(0.125);
    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetLeftMargin(0.125);

    gt2k_dcp0.SetLineColor(pBrown);
    gt2k_dcp0.SetLineWidth(5);
    gt2k_dcp0.SetLineStyle(2);
    gt2k_dcp05.SetLineColor(pCyan);
    gt2k_dcp05.SetLineWidth(5);
    gt2k_dcp05.SetLineStyle(2);
    gt2k_dcp1.SetLineColor(pPink);
    gt2k_dcp1.SetLineWidth(5);
    gt2k_dcp1.SetLineStyle(2);
    gt2k_dcp15.SetLineColor(pGreen);
    gt2k_dcp15.SetLineWidth(5);
    gt2k_dcp15.SetLineStyle(2);

    gt2k_dcp0.Draw("AC");
    gt2k_dcp05.Draw("C");
    gt2k_dcp1.Draw("C");
    gt2k_dcp15.Draw("C");

    TLegend leg(0.7, 0.4, 0.95, 0.975);
    leg.SetHeader("#bar{#nu}_{e}, T2K Best Fit, L=295 km");
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&gt2k_dcp0, GetString(osc_params[kT2KBF_IH_dcp0], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp05, GetString(osc_params[kT2KBF_IH_dcp05], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp1, GetString(osc_params[kT2KBF_IH_dcp1], 5).c_str(),
                 "l");
    leg.AddEntry(&gt2k_dcp15, GetString(osc_params[kT2KBF_IH_dcp15], 5).c_str(),
                 "l");

    TLegendEntry *header =
        static_cast<TLegendEntry *>(leg.GetListOfPrimitives()->First());
    header->SetTextAlign(32);

    leg.Draw();

    c1.Print("T2KOscXSecPlot_appbIH.png");
  }

  // ellipse event rate
  {

    TCanvas c1("c1", "", 1400, 1000);

    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetBottomMargin(0.125);
    c1.SetLeftMargin(0.15);

    double scale = GetMaximum(T2KFlux.gevrate[app(kT2KBF, kCC0Pi, NOsc)]);

    auto ell =
        DrawEllipse(osc_params[kT2KBF], Scale(T2KFlux.gflux, 1.0 / scale),
                    xsgraphs.gxsecs[kCC0Pi], Scale(T2KFlux.gfluxb, 1.0 / scale),
                    xsgraphsb.gxsecs[kCC0Pi], t2k_baseline);

    auto ell_IH =
        DrawEllipse(osc_params[kT2KBF_IH], Scale(T2KFlux.gflux, 1.0 / scale),
                    xsgraphs.gxsecs[kCC0Pi], Scale(T2KFlux.gfluxb, 1.0 / scale),
                    xsgraphsb.gxsecs[kCC0Pi], t2k_baseline);

    double xmin = 0.8 * std::min(GetMinimumX(ell), GetMinimumX(ell_IH));
    double xmax = 1.1 * std::max(GetMaximumX(ell), GetMaximumX(ell_IH));

    double ymin = 0.8 * std::min(GetMinimum(ell), GetMinimum(ell_IH));
    double ymax = 1.1 * std::max(GetMaximum(ell), GetMaximum(ell_IH));

    auto dg = DummyGraph(xmin, xmax, ymin, ymax);

    AxisStyle(dg);

    dg.GetYaxis()->SetRangeUser(ymin, ymax);
    dg.GetYaxis()->SetTitle("#bar{#nu}_{e} Events (Arbitrary Units)");
    dg.GetYaxis()->SetTitleOffset(1.2);

    dg.GetXaxis()->SetRangeUser(xmin, xmax);
    dg.GetXaxis()->SetTitle("#nu_{e} Events (Arbitrary Units)");

    double n, nb;
    ell.GetPoint(0, n, nb);
    TMarker mdcp0(n, nb, 20);
    mdcp0.SetMarkerColor(pBrown);
    mdcp0.SetMarkerSize(3);

    ell.GetPoint(25, n, nb);
    TMarker mdcp05(n, nb, 21);
    mdcp05.SetMarkerColor(pCyan);
    mdcp05.SetMarkerSize(3);

    ell.GetPoint(50, n, nb);
    TMarker mdcp1(n, nb, 33);
    mdcp1.SetMarkerColor(pPink);
    mdcp1.SetMarkerSize(4);

    ell.GetPoint(75, n, nb);
    TMarker mdcp15(n, nb, 34);
    mdcp15.SetMarkerColor(pGreen);
    mdcp15.SetMarkerSize(3);


    dg.Draw("AC");
    ell.SetLineWidth(5);
    ell.Draw("C");

    mdcp0.Draw();
    mdcp05.Draw();
    mdcp1.Draw();
    mdcp15.Draw();

    TLegend leg(0.7, 0.4, 0.95, 0.975);
    leg.SetHeader("#nu_{e}, T2K Best Fit, L=295 km");
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&mdcp0, GetString(osc_params[kT2KBF_dcp0], 5).c_str(), "p");
    leg.AddEntry(&mdcp05, GetString(osc_params[kT2KBF_dcp05], 5).c_str(), "p");
    leg.AddEntry(&mdcp1, GetString(osc_params[kT2KBF_dcp1], 5).c_str(), "p");
    leg.AddEntry(&mdcp15, GetString(osc_params[kT2KBF_dcp15], 5).c_str(), "p");

    TLegendEntry *header =
        static_cast<TLegendEntry *>(leg.GetListOfPrimitives()->First());
    header->SetTextAlign(32);

    leg.Draw();

    c1.Print("T2KOscXSecPlot_ellipse_NH.png");
  }

  // ellipse event rate
  {

    TCanvas c1("c1", "", 1400, 1000);

    c1.SetTopMargin(0.03);
    c1.SetRightMargin(0.03);
    c1.SetBottomMargin(0.125);
    c1.SetLeftMargin(0.15);

    double scale = GetMaximum(T2KFlux.gevrate[app(kT2KBF, kCC0Pi, NOsc)]);

    auto ell =
        DrawEllipse(osc_params[kT2KBF], Scale(T2KFlux.gflux, 1.0 / scale),
                    xsgraphs.gxsecs[kCC0Pi], Scale(T2KFlux.gfluxb, 1.0 / scale),
                    xsgraphsb.gxsecs[kCC0Pi], t2k_baseline);

    auto ell_IH =
        DrawEllipse(osc_params[kT2KBF_IH], Scale(T2KFlux.gflux, 1.0 / scale),
                    xsgraphs.gxsecs[kCC0Pi], Scale(T2KFlux.gfluxb, 1.0 / scale),
                    xsgraphsb.gxsecs[kCC0Pi], t2k_baseline);

    double xmin = 0.8 * std::min(GetMinimumX(ell), GetMinimumX(ell_IH));
    double xmax = 1.1 * std::max(GetMaximumX(ell), GetMaximumX(ell_IH));

    double ymin = 0.8 * std::min(GetMinimum(ell), GetMinimum(ell_IH));
    double ymax = 1.1 * std::max(GetMaximum(ell), GetMaximum(ell_IH));

    auto dg = DummyGraph(xmin, xmax, ymin, ymax);

    AxisStyle(dg);

    dg.GetYaxis()->SetRangeUser(ymin, ymax);
    dg.GetYaxis()->SetTitle("#bar{#nu}_{e} Events (Arbitrary Units)");
    dg.GetYaxis()->SetTitleOffset(1.2);

    dg.GetXaxis()->SetRangeUser(xmin, xmax);
    dg.GetXaxis()->SetTitle("#nu_{e} Events (Arbitrary Units)");

    double n, nb;
    ell.GetPoint(0, n, nb);
    TMarker mdcp0(n, nb, 20);
    mdcp0.SetMarkerColor(pBrown);
    mdcp0.SetMarkerSize(3);

    ell.GetPoint(25, n, nb);
    TMarker mdcp05(n, nb, 21);
    mdcp05.SetMarkerColor(pCyan);
    mdcp05.SetMarkerSize(3);

    ell.GetPoint(50, n, nb);
    TMarker mdcp1(n, nb, 33);
    mdcp1.SetMarkerColor(pPink);
    mdcp1.SetMarkerSize(4);

    ell.GetPoint(75, n, nb);
    TMarker mdcp15(n, nb, 34);
    mdcp15.SetMarkerColor(pGreen);
    mdcp15.SetMarkerSize(3);

    ell_IH.GetPoint(0, n, nb);
    TMarker mdcp0_IH(n, nb, 20);
    mdcp0_IH.SetMarkerColor(pBrown);
    mdcp0_IH.SetMarkerSize(3);

    ell_IH.GetPoint(25, n, nb);
    TMarker mdcp05_IH(n, nb, 21);
    mdcp05_IH.SetMarkerColor(pCyan);
    mdcp05_IH.SetMarkerSize(3);

    ell_IH.GetPoint(50, n, nb);
    TMarker mdcp1_IH(n, nb, 33);
    mdcp1_IH.SetMarkerColor(pPink);
    mdcp1_IH.SetMarkerSize(4);

    ell_IH.GetPoint(75, n, nb);
    TMarker mdcp15_IH(n, nb, 34);
    mdcp15_IH.SetMarkerColor(pGreen);
    mdcp15_IH.SetMarkerSize(3);

    dg.Draw("AC");

    ell.SetLineWidth(5);
    ell_IH.SetLineWidth(5);

    ell.Draw("C");
    ell_IH.SetLineStyle(2);
    ell_IH.Draw("C");

    mdcp0.Draw();
    mdcp05.Draw();
    mdcp1.Draw();
    mdcp15.Draw();

    mdcp0_IH.Draw();
    mdcp05_IH.Draw();
    mdcp1_IH.Draw();
    mdcp15_IH.Draw();

    TLegend leg(0.7, 0.4, 0.95, 0.975);
    leg.SetHeader("#nu_{e}, T2K Best Fit, L=295 km");
    leg.SetBorderSize(0);
    leg.SetTextSize(0.065);
    leg.SetNColumns(1);
    leg.SetFillColorAlpha(0, 0);
    leg.AddEntry(&mdcp0, GetString(osc_params[kT2KBF_dcp0], 5).c_str(), "p");
    leg.AddEntry(&mdcp05, GetString(osc_params[kT2KBF_dcp05], 5).c_str(), "p");
    leg.AddEntry(&mdcp1, GetString(osc_params[kT2KBF_dcp1], 5).c_str(), "p");
    leg.AddEntry(&mdcp15, GetString(osc_params[kT2KBF_dcp15], 5).c_str(), "p");

    TLegendEntry *header =
        static_cast<TLegendEntry *>(leg.GetListOfPrimitives()->First());
    header->SetTextAlign(32);

    leg.Draw();

    TLegend leg2(0.15, 0.2, 0.6, 0.4);
    leg2.SetBorderSize(0);
    leg2.SetTextSize(0.065);
    leg2.SetNColumns(1);
    leg2.SetFillColorAlpha(0, 0);
    leg2.AddEntry(&ell, "Normal Hierarchy", "l");
    leg2.AddEntry(&ell_IH, "Inverted Hierarchy", "l");
    leg2.Draw();

    c1.Print("T2KOscXSecPlot_ellipse_both.png");
  }
}
