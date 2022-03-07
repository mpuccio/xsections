#include "G4AntiProton.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentBarNucleonNucleusXsc.hh"
#include "G4ChipsProtonInelasticXS.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"

#include <TFile.h>
#include <TObjString.h>
#include <TString.h>
#include <TGraph.h>

void g4xsection() {
  auto proton = G4Proton::Definition();
  auto antiproton = G4AntiProton::Definition();

  G4ComponentBarNucleonNucleusXsc matterXS;
  G4ChipsProtonInelasticXS matterXS1;
  G4ComponentAntiNuclNuclearXS antimatterXS;
  matterXS.BuildPhysicsTable(*proton);

  TFile protonData("xsections.root");
  TFile output("g4xsection.root", "recreate");
  auto antidir = output.mkdir("antiproton");
  auto dir = output.mkdir("proton");
  dir->cd();
  for (auto key : *(protonData.GetListOfKeys())) {
    TString name{key->GetName()};
    auto zS = (TObjString*)((name.Tokenize("-"))->At(0));
    auto aS = (TObjString*)((name.Tokenize("-"))->At(1));
    int z = zS->String().Atoi();
    int a = aS->String().Atoi();
    std::vector<double> energy, mom, xs, xs1, antixs;
    for (int iE{0}; iE <= 1000; ++iE) {
      energy.push_back(100 + iE * 5);
      mom.push_back(std::sqrt(energy.back() * (energy.back() + 2. * 938.272)));
      xs.push_back(matterXS.GetInelasticElementCrossSection(proton, energy.back() * megaelectronvolt, z, a) / barn); ///A is dummy here
      if (a) {
        xs1.push_back(matterXS1.GetChipsCrossSection(mom.back(), z, a - z, 2212) / barn); ///A is dummy here
        antixs.push_back(antimatterXS.GetInelasticIsotopeCrossSection(antiproton, energy.back() * megaelectronvolt, z, a) /barn);
      }
    }
    TGraph gh(energy.size(), energy.data(), xs.data());
    TGraph ghP(energy.size(), mom.data(), xs.data());
    dir->cd();
    gh.Write(zS->String());
    ghP.Write(zS->String() + "_mom");
    if (a) {
      TGraph ghPalt(energy.size(), mom.data(), xs1.data());
      ghPalt.Write(zS->String() + "-" + aS->String() + "_mom1");
      antidir->cd();
      TGraph antigh(energy.size(), energy.data(), antixs.data());
      antigh.Write(zS->String() + "-" + aS->String());
      TGraph antighmom(energy.size(), mom.data(), antixs.data());
      antighmom.Write(zS->String() + "-" + aS->String() + "_mom");
    }
  }
}