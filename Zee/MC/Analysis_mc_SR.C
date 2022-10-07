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
    // double binEdges[155] = {0};
    // int j = 0;
    // for(j = 0; j <= 50; j++) binEdges[j] = j*0.001;
    // for(j = 1; j <= 25; j++) binEdges[j+50] = 0.05 + j*0.002;
    // for(j = 1; j <= 25; j++) binEdges[j+75] = 0.1 + j*0.008;
    // for(j = 1; j <= 10; j++) binEdges[j+100] = 0.3 + j*0.01;
    // for(j = 1; j <= 20; j++) binEdges[j+110] = 0.4 + j*0.02;
    // for(j = 1; j <= 4; j++) binEdges[j+130] = 0.8 + j*0.05;
    // for(j = 1; j <= 5; j++) binEdges[j+134] = 1.0 + j*0.2;
    // for(j = 1; j <= 4; j++) binEdges[j+139] = 2.0 + j*0.5;
    // for(j = 1; j <= 6; j++) binEdges[j+143] = 4.0 + j*1.0;
    // for(j = 1; j <= 5; j++) binEdges[j+149] = 10.0 + j*2.0;

    // double binEdges[139] = {0};
    // int j = 0;
    // for(j = 0; j <= 50; j++) binEdges[j] = j*0.001;
    // for(j = 1; j <= 25; j++) binEdges[j+50] = 0.05 + j*0.002;
    // for(j = 1; j <= 25; j++) binEdges[j+75] = 0.1 + j*0.008;
    // for(j = 1; j <= 10; j++) binEdges[j+100] = 0.3 + j*0.01;
    // for(j = 1; j <= 5; j++) binEdges[j+110] = 0.4 + j*0.02;
    // for(j = 1; j <= 6; j++) binEdges[j+115] = 0.5 + j*0.05;
    // for(j = 1; j <= 2; j++) binEdges[j+121] = 0.8 + j*0.1;
    // for(j = 1; j <= 5; j++) binEdges[j+123] = 1.0 + j*0.2;
    // for(j = 1; j <= 4; j++) binEdges[j+128] = 2.0 + j*0.5;
    // for(j = 1; j <= 2; j++) binEdges[j+132] = 4.0 + j*1.0;
    // for(j = 1; j <= 2; j++) binEdges[j+134] = 6.0 + j*2.0;
    // for(j = 1; j <= 2; j++) binEdges[j+136] = 10.0 + j*5.0;

    double binEdges[77] = {0};
    int j = 0;
    for(j = 0; j <= 25; j++) binEdges[j] = j*0.002;
    for(j = 1; j <= 10; j++) binEdges[j+25] = 0.05 + j*0.005;
    for(j = 1; j <= 15; j++) binEdges[j+35] = 0.1 + j*0.02;
    for(j = 1; j <= 5; j++) binEdges[j+50] = 0.4 + j*0.04;
    for(j = 1; j <= 4; j++) binEdges[j+55] = 0.6 + j*0.05;
    for(j = 1; j <= 2; j++) binEdges[j+59] = 0.8 + j*0.1;
    for(j = 1; j <= 5; j++) binEdges[j+61] = 1.0 + j*0.2;
    for(j = 1; j <= 4; j++) binEdges[j+66] = 2.0 + j*0.5;
    for(j = 1; j <= 2; j++) binEdges[j+70] = 4.0 + j*1.0;
    for(j = 1; j <= 2; j++) binEdges[j+72] = 6.0 + j*2.0;
    for(j = 1; j <= 2; j++) binEdges[j+74] = 10.0 + j*5.0;

    h_A = new TH1F("A","phi* distribution", 76, binEdges);
    h_B = new TH1F("B","phi* distribution", 76, binEdges);
    h_C = new TH1F("C","phi* distribution", 76, binEdges);
    h_D = new TH1F("D","phi* distribution", 76, binEdges);
    h_A_reside = new TH1F("A_reside","phi* distribution reside", 76, binEdges);
    h_A_flow = new TH1F("A_flow","phi* distribution flow", 76, binEdges);
    h_A_misid = new TH1F("A_misid","phi* distribution misid", 76, binEdges);
    h_B_reside = new TH1F("B_reside","phi* distribution reside", 76, binEdges);
    h_B_flow = new TH1F("B_flow","phi* distribution flow", 76, binEdges);
    h_B_misid = new TH1F("B_misid","phi* distribution misid", 76, binEdges);
    h_C_reside = new TH1F("C_reside","phi* distribution reside", 76, binEdges);
    h_C_flow = new TH1F("C_flow","phi* distribution flow", 76, binEdges);
    h_C_misid = new TH1F("C_misid","phi* distribution misid", 76, binEdges);
    h_D_reside = new TH1F("D_reside","phi* distribution reside", 76, binEdges);
    h_D_flow = new TH1F("D_flow","phi* distribution flow", 76, binEdges);
    h_D_misid = new TH1F("D_misid","phi* distribution misid", 76, binEdges);

}

//filling histograms
void CLoop::Fill(double weight, int z_sample) {

    //MATH CONSTANTS
    double pi = TMath::Pi();

    //ESSENTIAL CALCULATIONS
    double inv_mass=sqrt(2*elec_0_p4->Pt()*elec_1_p4->Pt()*(cosh(elec_0_p4->Eta()-elec_1_p4->Eta())-cos(elec_0_p4->Phi()-elec_1_p4->Phi())));
    
    //IDENTIFICATION flags
    bool elec_medium_id = elec_0_id_medium && elec_1_id_medium;

    //ISOLATION flags
    bool elec_loose_iso = elec_0_iso_FCLoose*elec_1_iso_FCLoose;
    bool elec_tight_iso = elec_0_iso_FCTight*elec_1_iso_FCTight;

    //EVENT SELECTION flags
    //oppositely charged
    bool elec_opposite_charge = (elec_0_q == -1*elec_1_q);
    //dielec
    bool dielec = (n_electrons == 2);
    //matched
    bool matched = (elec_1_matched && elec_0_matched);

    //SELECTION AND FILLING
    if(dielec && elec_medium_id && fiducial(elec_0_p4->Phi(), elec_1_p4->Phi(), elec_0_p4->Eta(), elec_1_p4->Eta(), elec_0_p4->Pt(), elec_1_p4->Pt())){
        double phi_star_detector = phi_star(elec_0_p4->Phi(), elec_1_p4->Phi(), elec_0_p4->Eta(), elec_1_p4->Eta());
        if(matched){
            //decide event type
            double phi_star_truth = phi_star(elec_0_matched_p4->Phi(), elec_1_matched_p4->Phi(), elec_0_matched_p4->Eta(), elec_1_matched_p4->Eta());
            bool flag_reside = false;
            if (phi_star_truth < 0.05) {if(floor(phi_star_detector*2000) == floor(phi_star_truth*2000)) flag_reside = true;}
            else if (phi_star_truth < 0.1) {if(floor((phi_star_detector-0.05)*1000) == floor((phi_star_truth-0.05)*1000)) flag_reside = true;}
            else if (phi_star_truth < 0.3) {if(floor((phi_star_detector-0.1)*250) == floor((phi_star_truth-0.1)*250)) flag_reside = true;}
            else if (phi_star_truth < 0.4) {if(floor((phi_star_detector-0.3)*100) == floor((phi_star_truth-0.3)*100)) flag_reside = true;}
            else if (phi_star_truth < 0.8) {if(floor((phi_star_detector-0.4)*50) == floor((phi_star_truth-0.4)*50)) flag_reside = true;}
            else {if(floor((phi_star_detector-0.8)*20) == floor((phi_star_truth-0.8)*20)) flag_reside = true;}
            //filling histos
            if(flag_reside){
                if(elec_opposite_charge){
                    if(elec_loose_iso){
                        h_A_reside->Fill(phi_star_detector, weight);
                        h_A->Fill(phi_star_detector, weight);
                    }else{
                        h_B_reside->Fill(phi_star_detector, weight);
                        h_B->Fill(phi_star_detector, weight);
                    }
                }else{
                    if(elec_loose_iso){
                        h_C_reside->Fill(phi_star_detector, weight);
                        h_C->Fill(phi_star_detector, weight);
                    }else{
                        h_D_reside->Fill(phi_star_detector, weight);
                        h_D->Fill(phi_star_detector, weight);
                    }
                }
            }else{
                if(elec_opposite_charge){
                    if(elec_loose_iso){
                        h_A_flow->Fill(phi_star_detector, weight);
                        h_A->Fill(phi_star_detector, weight);
                    }else{
                        h_B_flow->Fill(phi_star_detector, weight);
                        h_B->Fill(phi_star_detector, weight);
                    }
                }else{
                    if(elec_loose_iso){
                        h_C_flow->Fill(phi_star_detector, weight);
                        h_C->Fill(phi_star_detector, weight);
                    }else{
                        h_D_flow->Fill(phi_star_detector, weight);
                        h_D->Fill(phi_star_detector, weight);
                    }
                }
            }
        }else{
            if(elec_opposite_charge){
                if(elec_loose_iso){
                    h_A_misid->Fill(phi_star_detector, weight);
                    h_A->Fill(phi_star_detector, weight);
                }else{
                    h_B_misid->Fill(phi_star_detector, weight);
                    h_B->Fill(phi_star_detector, weight);
                }
            }else{
                if(elec_loose_iso){
                    h_C_misid->Fill(phi_star_detector, weight);
                    h_C->Fill(phi_star_detector, weight);
                }else{
                    h_D_misid->Fill(phi_star_detector, weight);
                    h_D->Fill(phi_star_detector, weight);
                }
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
    h_A_reside->Write();
    h_A_flow->Write();
    h_A_misid->Write();
    h_B_reside->Write();
    h_B_flow->Write();
    h_B_misid->Write();
    h_C_reside->Write();
    h_C_flow->Write();
    h_C_misid->Write();
    h_D_reside->Write();
    h_D_flow->Write();
    h_D_misid->Write();
}
#endif // End header guard