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
    static DoubleTab memorisedElongation;
    static DoubleVect myOrigine;
    static DoubleVect myLongueurs;
    static IntVect myNb_Noeuds;
    static double myRayon;
    static double d_desactivation_lubrification;

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
};



#endif //DNS_TOOL_H
