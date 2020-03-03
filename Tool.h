//
// Created by salim on 19/10/19.
//

#ifndef DNS_TOOL_H
#define DNS_TOOL_H


#include <DoubleVect.h>
#include <DoubleTab.h>

#include <string.h>
#include <string>



class Tool {
public:
    static  std::string myCode;
    static int isFirstCollision;
    static double memorisedElongation;
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

    //--
};



#endif //DNS_TOOL_H
