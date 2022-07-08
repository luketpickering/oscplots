#pragma once

#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <array>

void defc(int &c, std::array<int, 3> cs, std::string const &n) {
  static int currC = 51;
  TColor *col = gROOT->GetColor(currC);
  if (!col) {
    col = new TColor(currC, Float_t(cs[0]) / 255.0, Float_t(cs[1]) / 255.0,
                     Float_t(cs[2]) / 255.0, n.c_str());
  } else {
    col->SetRGB(Float_t(cs[0]) / 255.0, Float_t(cs[1]) / 255.0,
                Float_t(cs[2]) / 255.0);
  }
  c = currC;
  currC++;
}

void GetRGB(long in, int &R, int &G, int &B) {
  B = (in >> (8 * 0)) & 0xff;
  G = (in >> (8 * 1)) & 0xff;
  R = (in >> (8 * 2)) & 0xff;
}

void defc(int &c, long in, std::string const &n) {
  std::array<int, 3> cs;
  GetRGB(in, cs[0], cs[1], cs[2]);
  defc(c, cs, n);
}

static int kDUNEOrange;
static int kDUNEBlue;
static int kMSUGreen;
static int kMSUPurple;

static int kT2KRed;
static int kT2KGreen;
static int kNOvABlue;

static size_t const npastel = 11;
static int pastelWheel[npastel];

static size_t const nSORNBright = 7;
static int SORNBrightWheel[nSORNBright];

static size_t const nSORNVibrant = 7;
static int SORNVibrantWheel[nSORNVibrant];

static size_t const nSORNMuted = 10;
static int SORNMutedWheel[nSORNMuted];

inline void DeclareColors() {
  defc(kDUNEOrange, {242, 159, 84}, "kDUNEOrange");
  defc(kDUNEBlue, {125, 172, 213}, "kDUNEBlue");
  defc(kMSUGreen, {13, 177, 75}, "kMSUGreen");
  defc(kMSUPurple, {110, 0, 95}, "kMSUPurple");

  defc(kNOvABlue, {34, 62, 145}, "kNOvABlue");
  defc(kT2KGreen, {0, 128, 0}, "kT2KGreen");
  defc(kT2KRed, {128, 0, 0}, "kT2KRed");

  defc(pastelWheel[0], 0xff191f, "pRed");
  defc(pastelWheel[1], 0xff8b11, "pOrange");
  defc(pastelWheel[2], 0xfcff00, "pYellow");
  defc(pastelWheel[3], 0x00b600, "pGreen");
  defc(pastelWheel[4], 0x120790, "pBlue");
  defc(pastelWheel[5], 0x71007d, "pPurple");
  defc(pastelWheel[6], 0x000000, "pBlack");
  defc(pastelWheel[7], 0x582d0b, "pBrown");
  defc(pastelWheel[8], 0x57cfe7, "pCyan");
  defc(pastelWheel[9], 0xffa6c2, "pPink");
  defc(pastelWheel[10], 0xffffff, "pWhite");

  defc(SORNBrightWheel[0], 0x4477AA, "SORNBrightWheel_0");
  defc(SORNBrightWheel[1], 0x66CCEE, "SORNBrightWheel_1");
  defc(SORNBrightWheel[2], 0x228833, "SORNBrightWheel_2");
  defc(SORNBrightWheel[3], 0xCCBB44, "SORNBrightWheel_3");
  defc(SORNBrightWheel[4], 0xEE6677, "SORNBrightWheel_4");
  defc(SORNBrightWheel[5], 0xAA3377, "SORNBrightWheel_5");
  defc(SORNBrightWheel[6], 0xBBBBBB, "SORNBrightWheel_6");

  defc(SORNVibrantWheel[0], 0x0077BB, "SORNVibrantWheel_0");
  defc(SORNVibrantWheel[1], 0x33BBEE, "SORNVibrantWheel_1");
  defc(SORNVibrantWheel[2], 0x009988, "SORNVibrantWheel_2");
  defc(SORNVibrantWheel[3], 0xEE7733, "SORNVibrantWheel_3");
  defc(SORNVibrantWheel[4], 0xCC3311, "SORNVibrantWheel_4");
  defc(SORNVibrantWheel[5], 0xEE3377, "SORNVibrantWheel_5");
  defc(SORNVibrantWheel[6], 0xBBBBBB, "SORNVibrantWheel_6");

  defc(SORNMutedWheel[0], 0x332288, "SORNMutedWheel_0");
  defc(SORNMutedWheel[1], 0x88CCEE, "SORNMutedWheel_1");
  defc(SORNMutedWheel[2], 0x44AA99, "SORNMutedWheel_2");
  defc(SORNMutedWheel[3], 0x117733, "SORNMutedWheel_3");
  defc(SORNMutedWheel[4], 0x999933, "SORNMutedWheel_4");
  defc(SORNMutedWheel[5], 0xDDCC77, "SORNMutedWheel_5");
  defc(SORNMutedWheel[6], 0xCC6677, "SORNMutedWheel_6");
  defc(SORNMutedWheel[7], 0x882255, "SORNMutedWheel_7");
  defc(SORNMutedWheel[8], 0xAA4499, "SORNMutedWheel_8");
  defc(SORNMutedWheel[9], 0xDDDDDD, "SORNMutedWheel_9");
}

Int_t BirdPalette[100];
Int_t W2RPalette[100];
Int_t BWRPalette[100];

int *GetBirdPalette() {
  static bool first = true;
  if (first) {
    first = false;
    Double_t Red[] = {0., 1.0};
    Double_t Green[] = {0., 1.0};
    Double_t Blue[] = {1., 0.0};
    Double_t Length[] = {0., 1.0};
    Int_t FI =
        TColor::CreateGradientColorTable(2, Length, Red, Green, Blue, 100);
    for (int i = 0; i < 100; i++) {
      BirdPalette[i] = FI + i;
    }
  }
  return BirdPalette;
}

void SetBirdPalette() {
  gStyle->SetPalette(100, GetBirdPalette());
  gStyle->SetNumberContours(100);
}

int *GetW2RPalette(double maxpoint = 1) {
  static bool first = true;
  if (first) {
    first = false;
    Double_t Red[] = {1., 1.0, 1.0};
    Double_t Green[] = {1., 0.0, 0.0};
    Double_t Blue[] = {1., 0.0, 0.0};
    Double_t Length[] = {0., maxpoint, 1.0};
    Int_t FI =
        TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 100);
    for (int i = 0; i < 100; i++) {
      W2RPalette[i] = FI + i;
    }
  }
  return W2RPalette;
}

void SetW2RPalette(double maxpoint = 1) {
  gStyle->SetPalette(100, GetW2RPalette(maxpoint));
  gStyle->SetNumberContours(100);
}

int *GetBWRPalette() {
  static bool first = true;
  if (first) {
    first = false;
    Double_t Red[] = {0., 1.0, 1.0};
    Double_t Green[] = {0., 1.0, 0.0};
    Double_t Blue[] = {1., 1.0, 0.0};
    Double_t Length[] = {0., 0.5, 1.0};
    Int_t FI =
        TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 100);
    for (int i = 0; i < 100; i++) {
      BWRPalette[i] = FI + i;
    }
  }
  return BWRPalette;
}

void SetBWRPalette() {
  gStyle->SetPalette(100, GetBWRPalette());
  gStyle->SetNumberContours(100);
}