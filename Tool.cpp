//
// Created by salim on 19/10/19.
//

#include <Tool.h>

//declaration des membre donnees
std::string Tool::myCode="commicode 000std";

double Tool::myMuPhase1=-1;
double Tool::myMuPhase0=-1;
DoubleVect Tool::myIndic;
DoubleTab Tool::myNormaleInterfaceElem;



//implementation des membre fonctions

//utilitaire
double Tool::calcMyViscLam(int elem1,int elem2, int elem3, int elem4){

	  double indicArete = 0.25*(myIndic[elem1] + myIndic[elem2]
	                            +myIndic[elem3] + myIndic[elem4]);

	  double myViscLam  = (myMuPhase0 * myMuPhase1)/(myMuPhase1 - indicArete * ( myMuPhase1 - myMuPhase0));

    //printf("elem{1,2,3,4} ={%d|%d|%d|%d}\n",elem1,elem2,elem3,elem4);
    //printf("indic{1,2,3,4={%f|%f|%f|%f}\n",myIndic[elem1],myIndic[elem2],myIndic[elem3],myIndic[elem4]);
    //printf("indic arete=%f; mu0=%f, mu1=%f, visclam=%f\n", indicArete,myMuPhase0,myMuPhase1,myViscLam);

	return myViscLam;
}
// setters and getters
double Tool::getMyMuPhase1() {
    return myMuPhase1;
}

void Tool::setMyMuPhase1(double myMuPhase1) {
    Tool::myMuPhase1 = myMuPhase1;
}



double Tool::getMyMuPhase0() {
    return myMuPhase0;
}

void Tool::setMyMuPhase0(double myMuPhase0) {
    Tool::myMuPhase0 = myMuPhase0;
}



const DoubleVect &Tool::getMyIndic() {
    return myIndic;
}

void Tool::setMyIndic(const DoubleVect &myIndic) {
    Tool::myIndic = myIndic;
}



const DoubleTab &Tool::getMyNormaleInterfaceElem() {
    return myNormaleInterfaceElem;
}

void Tool::setMyNormaleInterfaceElem(const DoubleTab &myNormaleInterfaceElem) {
    Tool::myNormaleInterfaceElem = myNormaleInterfaceElem;
}



//--

//--