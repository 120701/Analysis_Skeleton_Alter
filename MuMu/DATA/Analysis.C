//look at eta,pt,invmass,phi of truth level Zee

// Header guard to ensure file is imported properly
#ifndef Analysis
#define Analysis

// Include the file that lets the program know about the data
#include "backend/CLoop.h"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <fstream>
#include <cmath>

double acop_angle(const double del_phi){
    if(del_phi < M_PI){
        return M_PI - del_phi;
    }else{
        return del_phi - M_PI;
    }
}

double phi_star(double phi_1, double phi_2, double eta_1, double eta_2)
{
    double del_phi, phi_acop, phi_star;
    del_phi = std::abs((phi_1) - (phi_2));
    phi_acop = acop_angle(del_phi);
    double sin_theta = sqrt(1.0 - pow(tanh(((eta_1) - (eta_2)) / 2.0), 2));
    phi_star = tan(phi_acop / 2.0) * sin_theta;
    return phi_star;
}

bool fiducial(double phi_1, double phi_2, double eta_1, double eta_2, double pt_1, double pt_2)
{
    double inv_mass;
    inv_mass=sqrt(2 * (pt_1) * (pt_2) * (cosh((eta_1) - (eta_2)) - cos((phi_1) - (phi_2))));

    if (((pt_1) > 27) && ((pt_2) > 27) && (abs(eta_1) < 2.5) && (abs(eta_2) < 2.5) && (inv_mass < 116) && (inv_mass > 66))
        return true;
    else
        return false;
}

//declare the histograms
void CLoop::Book(double lumFactor) {

    // construct the variable bins
    double binEdges[155] = {0};
    int j = 0;
    for(j = 0; j <= 50; j++) binEdges[j] = j*0.001;
    for(j = 1; j <= 25; j++) binEdges[j+50] = 0.05 + j*0.002;
    for(j = 1; j <= 25; j++) binEdges[j+75] = 0.1 + j*0.008;
    for(j = 1; j <= 10; j++) binEdges[j+100] = 0.3 + j*0.01;
    for(j = 1; j <= 20; j++) binEdges[j+110] = 0.4 + j*0.02;
    for(j = 1; j <= 4; j++) binEdges[j+130] = 0.8 + j*0.05;
    for(j = 1; j <= 5; j++) binEdges[j+134] = 1.0 + j*0.2;
    for(j = 1; j <= 4; j++) binEdges[j+139] = 2.0 + j*0.5;
    for(j = 1; j <= 6; j++) binEdges[j+143] = 4.0 + j*1.0;
    for(j = 1; j <= 5; j++) binEdges[j+149] = 10.0 + j*2.0;

    h_A = new TH1F("A","phi* distribution", 154, binEdges);
    h_B = new TH1F("B","phi* distribution", 154, binEdges);
    h_C = new TH1F("C","phi* distribution", 154, binEdges);
    h_D = new TH1F("D","phi* distribution", 154, binEdges);

}

//filling histograms
void CLoop::Fill(double weight, int z_sample) {
  
    //MATH CONSTANTS
    double pi = TMath::Pi();

    //ESSENTIAL CALCULATIONS
    double inv_mass=sqrt(2*muon_0_p4->Pt()*muon_1_p4->Pt()*(cosh(muon_0_p4->Eta()-muon_1_p4->Eta())-cos(muon_0_p4->Phi()-muon_1_p4->Phi())));
    
    //IDENTIFICATION flags
    bool muon_medium_id = muon_0_id_medium && muon_1_id_medium;

    //ISOLATION flags
    bool muon_loose_iso = muon_0_iso_Loose_FixedRad * muon_1_iso_Loose_FixedRad;
    bool muon_tight_iso = muon_0_iso_TightTrackOnly_FixedRad * muon_1_iso_TightTrackOnly_FixedRad;

    //EVENT SELECTION flags
    //oppositely charged
    bool muon_opposite_charge = (muon_0_q == -1*muon_1_q);
    //dimuon
    bool dimuon = (n_muons == 2);
    //matched
    bool matched = (muon_1_matched && muon_0_matched);

    //SELECTION AND FILLING
    if(muon_medium_id && dimuon && fiducial(muon_0_p4->Phi(), muon_1_p4->Phi(), muon_0_p4->Eta(), muon_1_p4->Eta(), muon_0_p4->Pt(), muon_1_p4->Pt())){
        //Filling histos
        if(muon_opposite_charge){
            if(muon_loose_iso){
                h_A->Fill(phi_star(muon_0_p4->Phi(), muon_1_p4->Phi(), muon_0_p4->Eta(), muon_1_p4->Eta()), weight);
            }
            else
                h_B->Fill(phi_star(muon_0_p4->Phi(), muon_1_p4->Phi(), muon_0_p4->Eta(), muon_1_p4->Eta()), weight);
        }else{
            if(muon_loose_iso){
                h_C->Fill(phi_star(muon_0_p4->Phi(), muon_1_p4->Phi(), muon_0_p4->Eta(), muon_1_p4->Eta()), weight);
            }else{
                h_D->Fill(phi_star(muon_0_p4->Phi(), muon_1_p4->Phi(), muon_0_p4->Eta(), muon_1_p4->Eta()), weight);
            }
        }
    }
}

//Write histograms
void CLoop::Style(double lumFactor) {
    h_A->Write();
    h_B->Write();
    h_C->Write();
    h_D->Write();
}
#endif // End header guard

