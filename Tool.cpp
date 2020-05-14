#include <Tool.h>
#include <IntTab.h>
#include <Pave.h>

//declaration des membre donnees
std::string Tool::myCode="commicode 015.3rmp";

double Tool::myMuPhase1=-1;
double Tool::myMuPhase0=-1;
double Tool::myRayon=-1;
double Tool::d_desactivation_lubrification =0;

DoubleVect Tool::myIndic;
DoubleTab Tool::myNormaleInterfaceElem;

IntTab Tool::isFirstCollision;
DoubleTab Tool::memorisedElongation;

DoubleVect Tool::myOrigine(3);
DoubleVect Tool::myLongueurs(3);
IntVect Tool::myNb_Noeuds(3);
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

double Tool::calc_positions_bords(ArrOfDouble &positions_bords)
{

    double dx =myLongueurs(0)/(myNb_Noeuds(0)-1);
    double epsi =dx/4;
    positions_bords[0] = myOrigine(0) + epsi;
    positions_bords[1] = myOrigine(1) + epsi;
    positions_bords[2] = myOrigine(2) + epsi;
    positions_bords[3] = myOrigine(0)+myLongueurs(0) - epsi;
    positions_bords[4] = myOrigine(1)+myLongueurs(1) - epsi;
    positions_bords[5] = myOrigine(2)+myLongueurs(2) - epsi;
    //for (int i = 0; i < 6; i++)
    //{
    //    printf("  (a) positions_bords[%d]=%f\n",i,positions_bords[i]);
    //}
return dx;
}



//--

//--
