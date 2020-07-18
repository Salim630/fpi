//
// Created by salim on 19/10/19.
//

#ifndef DNS_TOOL_H
#define DNS_TOOL_H


#include <DoubleVect.h>
#include <DoubleTab.h>

#include <string.h>
#include <string>
#include <IntTab.h>



class Tool {
public:
    static  std::string myCode;
    static IntTab isFirstCollision;
    static int formule_mu; // viscosite du fluide appliquer sur les maille diphasique et les aretes
    static DoubleTab memorisedElongation;
    static DoubleVect myOrigine;
    static DoubleVect myLongueurs;
    static IntVect myNb_Noeuds;
    static double myRayon;
    static double mySigma;
    static double d_desactivation_lubrification;
    static int compteur_;

    static DoubleTab F_old;
    static DoubleTab F_now;
    static DoubleTab raideur;
    static DoubleTab e_eff;
    static double vitessRelImp;

    static DoubleTab vitesses_compo;
    static DoubleTab positions_compo;
    static DoubleVect num_compo_;
    //static Champ_Inc num_compo;
private:
    static double myMuPhase1;
    static double myMuPhase0;
    static DoubleVect myIndic;
    static DoubleTab myNormaleInterfaceElem;


public:

    //Utilitaire
    static double calcMyViscLam(int elem1,int elem2, int elem3, int elem4);  //calcule de visclam a l'arrette avec ponderation harmonique

    //Setters and Getters

    static const DoubleVect &getMyIndic();
    static void setMyIndic(const DoubleVect &myIndic);

    static double getMyMuPhase1();
    static void setMyMuPhase1(double myMuPhase1);

    static double getMyMuPhase0();
    static void setMyMuPhase0(double myMuPhase0);

    static const DoubleTab &getMyNormaleInterfaceElem();
    static void setMyNormaleInterfaceElem(const DoubleTab &myNormaleInterfaceElem);

    static double calc_positions_bords(ArrOfDouble &positions_bords);


    //--
    static double calc_positions_bords2(ArrOfDouble &positions_bords);

    static double module_vecteur(DoubleTab &vecteur);

    static double prod_scal(DoubleTab &A, DoubleTab &B);

    static int checkForDuplicates(ArrOfInt &vector);


    static void backup_myVariables();
    static void load_myVariables();


    static int modele_collision;
    static int decalage_bords;
    static DoubleVect valeurs_decalage;
    static int transport_vitesse_cg;
    static int force_sur_elem_diphasiques;
};



#endif //DNS_TOOL_H
