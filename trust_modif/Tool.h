//
// Created by salim on 19/10/19.
//

#ifndef DNS_TOOL_H
#define DNS_TOOL_H


#include <DoubleVect.h>
#include <DoubleTab.h>

#include <string.h>
#include <string>

#include <Ref_Zone_VDF.h>
#include <Zone_VDF.h>

#include <Zone_VF.h>
#include <Ref_Zone_VF.h>

class Tool {
public:
    static  std::string myCode;
    static REF(Zone_VDF)    ma_zone_VDF_;
    static REF(Zone_VF)    myZone_vf_;

    static DoubleTab myVitesse; //vitesse discritise
    static DoubleTab myVitesseFaces;
    static  DoubleTab& myCalculer_vitesse_faces(DoubleTab& v_faces_stockage); //vitesse aux faces
    static void print_doubletab(const DoubleTab& v_faces_stockage, Nom nom_fichier );
    static DoubleTab myVitesseSommets; //vitesse aux sommets

    static double myTime;
    static double myOldTime;

    static double myMuPhase1;
    static double myMuPhase0;
    static DoubleVect myIndic;
    static DoubleTab myNormaleInterfaceElem;
    static int myNiter;

public:

    //Utilitaire
    static double calcMyViscLam(int elem1,int elem2, int elem3, int elem4);  //calcule de visclam a l'arrette avec ponderation harmonique
    static void   calcMyVisc_areteInterne(int elem1,int elem2, int elem3, int elem4, double& mu_a, double& mu_h); // renvois les deux ponderation
    static void   calcMyVisc_fa7Elem(int elem, double& mu_a, double& mu_h); // renvois les deux ponderation pour
    //Setters and Getters

    static const DoubleVect &getMyIndic();
    static void setMyIndic(const DoubleVect &myIndic);

    static double getMyMuPhase1();
    static void setMyMuPhase1(double myMuPhase1);

    static double getMyMuPhase0();
    static void setMyMuPhase0(double myMuPhase0);

    static const DoubleTab &getMyNormaleInterfaceElem();
    static void setMyNormaleInterfaceElem(const DoubleTab &myNormaleInterfaceElem);

    //--
};



#endif //DNS_TOOL_H
