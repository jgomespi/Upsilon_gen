#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include <math.h>
using namespace Pythia8;
using namespace std;

int main() { 
  gStyle->SetOptStat(1111111); //Opção da caixa de estatística

  Pythia pythia;
  // Colisão pp numa energia de centro de massa de 7 TeV
  pythia.readString("Beams:eCM = 7000."); // energia do CM
  pythia.readString("Beams:idA = 2212");  // um feixe de prótons
  pythia.readString("Beams:idB = 2212");  // outro feixe de prótons

  pythia.readString("Bottomonium:all=on"); // esta é a única flag ligada, que ativa a produção de todas as ressonâncias de bottomonium 
  pythia.init();
  Event *event = &pythia.event;

  // Declarando os histrogramas

  // Variáveis cinemáticas do Upsilon (eta, pT e phi) sem e com corte no pT e eta dos múons (pT>3,5 GeV e |eta|<2,1).
  TH1D* UpsilonPt = new TH1D("UpsilonPt","Upsilon pT",30,0,30);
  TH1D* UpsilonCutPt = new TH1D("UpsilonCutPt","Upsilon withCuts pT",30,0,30);
  TH1D* UpsilonEta = new TH1D("UpsilonEta","Upsilon Eta",20,-10.,10.);
  TH1D* UpsilonCutEta = new TH1D("UpsilonCutEta","Upsilon withCut Eta",20,-10.,10.);
  TH1D* UpsilonPhi = new TH1D("UpsilonPhi","Upsilon ",64,-3.2,3.2);
  TH1D* UpsilonCutPhi = new TH1D("UpsilonCutPhi","UpsilonCutPt",64,-3.2,3.2);

  // pT do Upsilon com os pesos referentes à polarização, para alpha=1 (pol1) e alpha=-1 (pol2), também sem e com corte nas variáveis dos múons
  TH1D* UpsilonPt_polar1 = new TH1D("UpsilonPt_polar1","Upsilon pT I = 1+cos^2#theta",30,0,30);
  TH1D* UpsilonCutPt_polar1 = new TH1D("UpsilonCutPt_polar1","Upsilon withCut pT I = 1+cos^2#theta",30,0,30);
  TH1D* UpsilonPt_polar2 = new TH1D("UpsilonPt_polar2","Upsilon pT I = 1-cos^2#theta",30,0,30);
  TH1D* UpsilonCutPt_polar2 = new TH1D("UpsilonCutPt_polar2","Upsilon withCut pT I = 1-cos^2#theta",30,0,30);

  // Variáveis cinemáticas (pT, eta e phi) dos múons positivos e negativos       
  TH1D* munPt = new TH1D("munPt","mu- pT",100,0,30);
  TH1D* mupPt = new TH1D("mupPt","mu+ pT",100,0,30);
  TH1D* munEta = new TH1D("munEta","mu- Eta",100,-10.,10.);
  TH1D* mupEta = new TH1D("mupEta","mu+ Eta",100,-10.,10.);
  TH1D* munPhi = new TH1D("munPhi","mu- Phi",100,-3.2,3.2);
  TH1D* mupPhi = new TH1D("mupPhi","mu+ Phi",100,-3.2,3.2);

  double acept(0.),aceptpol1(0),aceptpol2(0); // Cálculo da aceptância para alpha =0 (acept), =1 (aceptpol1) e =-1 (aceptpol2).
  int ncut(0),ncutpol1(0),ncutpol2(0); // Contador de eventos com corte na variável dos múons
  int ntotal=0; // Contador do número total de eventos
  int nev = 30000000; // Número total de eventos gerados
  int UpsilonDaughter1(0), UpsilonDaughter2(0); // Variáveis que vão receber as partículas filhas do Upsilon.
  int indexUpsilon(0); // Índice do Upsilon no evento gerado
  int munIndex(0), mupIndex(0); // Índice do par de múons no evento gerado
  double thetastar(0), Ipol1(0), Ipol2(0); // theta* e os pesos correspondentes a cada polarização (alpha=1 (pol1) e =-1 (pol2))

  // Começa o loop de eventos
  for (int iEvent = 0; iEvent < nev; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < 1) {pythia.info.list(); pythia.event.list();} // Imprime o primeiro evento

    indexUpsilon = -1; // Índice do Upsilon recebe uma flag -1
    // Começa o loop de partículas
    for (int i = 0; i < pythia.event.size(); ++i){

      Particle& theParticle = pythia.event[i];

      // Procura por um Upsilon (id=553)
      if (abs(theParticle.id()) == 553) {
        indexUpsilon = i; // Pega o índice do Upsilon
        break; // Se acha o Upsilon, sai do loop de partículas
      }     
    } // Fim do loop de partículas.
    if (indexUpsilon == -1) continue; // Se não houve identificação de Upsilon no evento, segue para o próximo evento.

    //Encontra as filhas do Upsilon
    UpsilonDaughter1 = pythia.event[indexUpsilon].daughter1();
    UpsilonDaughter2 = pythia.event[indexUpsilon].daughter2();

    // Zerando o índice dos múons
    munIndex=0;
    mupIndex=0;

    // Procurando por um par de múon-antimúon (id=+/-13) entre as filhas do Upsilon
    if (UpsilonDaughter1<UpsilonDaughter2) {
      for (int i=UpsilonDaughter1; i<=UpsilonDaughter2; ++i) {
        if (pythia.event[i].id()==13)  munIndex=i;
        if (pythia.event[i].id()==-13) mupIndex=i;
      }
    }
    // Checando se 2 filhas do Upsilon são um par de múon-antimúon
    if (munIndex!=0 && mupIndex!=0) {
      // Associando cada partícula: muon+, muon- e Upsilon a um TLorentzVector
      TLorentzVector MuonP(pythia.event[mupIndex].px(),pythia.event[mupIndex].py(),pythia.event[mupIndex].pz(),pythia.event[mupIndex].e());
      TLorentzVector MuonN(pythia.event[munIndex].px(),pythia.event[munIndex].py(),pythia.event[munIndex].pz(),pythia.event[munIndex].e());
      TLorentzVector MuonP_CM(pythia.event[mupIndex].px(),pythia.event[mupIndex].py(),pythia.event[mupIndex].pz(),pythia.event[mupIndex].e()); // A seguir faremos um boost
      TLorentzVector MuonN_CM(pythia.event[munIndex].px(),pythia.event[munIndex].py(),pythia.event[munIndex].pz(),pythia.event[munIndex].e()); // A seguir faremos um boost
      TLorentzVector Upsilon(pythia.event[munIndex].px(),pythia.event[munIndex].py(),pythia.event[munIndex].pz(),pythia.event[munIndex].e());

      TVector3 Upsilon_CM = -(MuonP+MuonN).BoostVector(); // Obter o centro de massa do par de múons.

      // Realizando um boost, ou seja, obtendo as variáveis cinemáticas dos múons no sistema de referência do seu CM (onde o  Upsilon está em repouso)
      MuonP_CM.Boost(Upsilon_CM);
      MuonN_CM.Boost(Upsilon_CM);

      // Imprime na tela informações sobre o evento:
      cout << "Event number " << iEvent << endl;
      cout << "Found an event " << pythia.event[indexUpsilon].name() << " -> " << pythia.event[munIndex].name() << " " << pythia.event[mupIndex].name() << endl;
      cout << "Mu+ 4-mom = " << pythia.event[munIndex].p() << endl;
      cout << "Mu- 4-mom = " << pythia.event[mupIndex].p() << endl;
      
      // Preenchendo os histogramas com as variáveis cinemáticas dos múons
      mupPt->Fill(pythia.event[mupIndex].pT());
      munPt->Fill(pythia.event[munIndex].pT());   
      mupEta->Fill(pythia.event[mupIndex].phi());
      munEta->Fill(pythia.event[munIndex].phi()); 
      mupPhi->Fill(pythia.event[mupIndex].eta());
      munPhi->Fill(pythia.event[munIndex].eta());  

      // ângulo theta*, que é o ângulo entre o momentum do múon no sistema de repouso do Upsilon e o momentum do Upsilon no sistema de laboratório.
      thetastar = abs(MuonP_CM.Theta() - pythia.event[indexUpsilon].theta());

      // cálculo do peso devido a polarização
      Ipol1 = 1 + cos(thetastar)*cos(thetastar); // alpha = 1
      Ipol2 = 1 - cos(thetastar)*cos(thetastar); // alpha = -1

      // Preenchendo os histogramas com as variáveis cinemáticas do Upsilon
      UpsilonPt->Fill(pythia.event[indexUpsilon].pT());
      UpsilonPt_polar1->Fill(Ipol1*pythia.event[indexUpsilon].pT());
      UpsilonPt_polar2->Fill(Ipol2*pythia.event[indexUpsilon].pT());            
      UpsilonEta->Fill(pythia.event[indexUpsilon].eta());
      UpsilonPhi->Fill(pythia.event[indexUpsilon].phi());

      // Verificar se o Upsilon não polarizado está na região de aceptância
      if (pythia.event[munIndex].pT() > 3.5      && pythia.event[mupIndex].pT() > 3.5 &&
        abs(pythia.event[munIndex].eta()) < 2.1 && abs(pythia.event[mupIndex].eta()) < 2.1) {
        ncut++;
        UpsilonCutPt->Fill(pythia.event[indexUpsilon].pT());                             
        UpsilonCutEta->Fill(pythia.event[indexUpsilon].eta());
        UpsilonCutPhi->Fill(pythia.event[indexUpsilon].phi());
      }
      // Verificar se o Upsilon polarizado com alpha=1 está na região de aceptância
      if (Ipol1*pythia.event[munIndex].pT() > 3.5      && Ipol1*pythia.event[mupIndex].pT() > 3.5 &&
        abs(Ipol1*pythia.event[munIndex].eta()) < 2.1 && abs(Ipol1*pythia.event[mupIndex].eta()) < 2.1) {
        ncutpol1++;
        UpsilonCutPt_polar1->Fill(Ipol1*pythia.event[indexUpsilon].pT());          
      }      
      // Verificar se o Upsilon polarizado com alpha=-1 está na região de aceptância
      if (Ipol2*pythia.event[munIndex].pT() > 3.5      && Ipol2*pythia.event[mupIndex].pT() > 3.5 &&
        abs(Ipol2*pythia.event[munIndex].eta()) < 2.1 && abs(Ipol2*pythia.event[mupIndex].eta()) < 2.1) {
        ncutpol2++;
        UpsilonCutPt_polar2->Fill(Ipol2*pythia.event[indexUpsilon].pT());          
      }       
      ntotal++;
    }
  } // Fim do loop de eventos

  // Cálculo da aceptância para as diferentes polarizações assumidas
  acept = double(ncut)/double(ntotal);
  aceptpol1 = double(ncutpol1)/double(ntotal);
  aceptpol2 = double(ncutpol2)/double(ntotal);

  // Imprimir o resultado na tela
  if (ntotal!=0) cout << "Acceptance = #cut/#total = " << ncut << "/" << ntotal << " = " << acept << endl;
  if (ntotal!=0) cout << "Acceptancepol1 = #cut/#total = " << ncutpol1 << "/" << ntotal << " = " << acept << endl;
  if (ntotal!=0) cout << "Acceptancepol2 = #cut/#total = " << ncutpol2 << "/" << ntotal << " = " << acept << endl;

  // Informação sobre a estatística da geração dos eventos
  pythia.stat();
  // Criar o arquivo de output onde os histogramas serão armazenados
  TFile* outFile = new TFile("gen.root","RECREATE");
  UpsilonPhi->Write(); UpsilonCutPhi->Write();
  UpsilonEta->Write(); UpsilonCutEta->Write();
  UpsilonPt->Write();  UpsilonCutPt->Write();
  UpsilonPt_polar1->Write(); UpsilonCutPt_polar1->Write();    
  UpsilonPt_polar2->Write(); UpsilonCutPt_polar2->Write();     
  munPt->Write();  mupPt->Write();
  munEta->Write(); mupEta->Write();
  munPhi->Write(); mupPhi->Write();    

  // Calculando a aceptância em função do pT para diferentes polarizações
  TH1 *Accept = (TH1*)UpsilonCutPt->Clone("Accept"); 
  TH1 *Acceptpol1 = (TH1*)UpsilonCutPt_polar1->Clone("Acceptpol1"); 
  TH1 *Acceptpol2 = (TH1*)UpsilonCutPt_polar2->Clone("Acceptpol2"); 

  Accept->Divide(UpsilonCutPt,UpsilonPt);
  Accept->Write();
  Acceptpol1->Divide(UpsilonCutPt_polar1,UpsilonPt_polar1);
  Acceptpol1->Write();
  Acceptpol2->Divide(UpsilonCutPt_polar2,UpsilonPt_polar2);
  Acceptpol2->Write();

  TCanvas* c1 = new TCanvas("c1","",800,800);
  Acceptpol1->SetLineColor(2);
  Acceptpol2->SetLineColor(3);

  Acceptpol1->SetMarkerStyle(2);
  Acceptpol1->Draw("P");
  Accept->SetMarkerStyle(3);  
  Accept->Draw("P""SAME");
  Acceptpol2->SetMarkerStyle(4);
  Acceptpol2->Draw("P""SAME");  
  c1->SaveAs("CMS_acceptance.pdf","pdf");

  outFile->Close();
  return 0;
}
