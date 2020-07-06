#include <Tool.h>
#include <IntTab.h>
#include <Pave.h>


#include <communications.h>
#include <SFichier.h>
#include <EcrFicCollecte.h>
#include <EcrFicPartage.h>
#include <EFichier.h>
#include <LecFicDiffuse.h>
#include <LecFicDistribue.h>
#include <LecFicDistribueBin.h>
#include <EcrFicCollecteBin.h>
#include <Champ_Inc.h>


//declaration des membre donnees
std::string Tool::myCode="commicode 015.6rpo";
int dimension = 3;
double Tool::myMuPhase1=-1;
double Tool::myMuPhase0=-1;
DoubleVect Tool::myIndic;
DoubleTab Tool::myNormaleInterfaceElem;



DoubleTab Tool::memorisedElongation;


double Tool::myRayon=-1;
double Tool::mySigma=-1;
double Tool::d_desactivation_lubrification =0;

DoubleVect Tool::myOrigine(3);
DoubleVect Tool::myLongueurs(3);
IntVect Tool::myNb_Noeuds(3);

int Tool::compteur_=0;
DoubleTab Tool::F_old;
DoubleTab Tool::F_now;
DoubleTab Tool::raideur;
DoubleTab Tool::e_eff;

double Tool::vitessRelImp;
IntTab Tool::isFirstCollision;
int Tool::isMuPhaseFluide=0;
DoubleTab Tool::vitesses_compo;
DoubleTab Tool::positions_compo;
DoubleVect Tool::num_compo_;
//Champ_Inc Tool::num_compo;
//implementation des membre fonctions

//utilitaire
double Tool::calcMyViscLam(int elem1,int elem2, int elem3, int elem4){

	  double indicArete = 0.25*(myIndic[elem1] + myIndic[elem2]
	                            +myIndic[elem3] + myIndic[elem4]);

	  double myViscLam  =isMuPhaseFluide? myMuPhase1 : (myMuPhase0 * myMuPhase1)/(myMuPhase1 - indicArete * ( myMuPhase1 - myMuPhase0));

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

double Tool::calc_positions_bords2(ArrOfDouble &positions_bords)
{

    double dx =myLongueurs(0)/(myNb_Noeuds(0)-1);
    double epsi =dx/4;
    positions_bords[0] = myOrigine(0)- myRayon + epsi;
    positions_bords[1] = myOrigine(1)- myRayon + epsi;
    positions_bords[2] = myOrigine(2)- myRayon + epsi;
    positions_bords[3] = myOrigine(0)+myLongueurs(0) +myRayon - epsi;
    positions_bords[4] = myOrigine(1)+myLongueurs(1) +myRayon - epsi;
    positions_bords[5] = myOrigine(2)+myLongueurs(2) +myRayon - epsi;
    //for (int i = 0; i < 6; i++)
    //{
    //    printf("  (a) positions_bords[%d]=%f\n",i,positions_bords[i]);
    //}
    return dx;
}

double Tool::module_vecteur(DoubleTab &vecteur)
{
    int dim =3;
    double module = 0;
    for (int d = 0; d < dim; d++)
    {
        double tmp = vecteur(d);
        tmp *= tmp;
        module += tmp;
    }
    module = sqrt(module);
    return module;
}

double Tool::prod_scal(DoubleTab &A, DoubleTab &B)
{
    int dim = 3 ;
    double prod = 0 ;
    for (int d = 0; d < dim; d++) prod += A(d) * B(d) ;
    return prod;
}

int Tool::checkForDuplicates(ArrOfInt &vector)
{
    int flag =0;
    ArrOfInt copy_vector(vector);
    const int size = copy_vector.size_array();
    copy_vector.ordonne_array();
    for (int i = 0; i < size-1; i++)
    {
        if (copy_vector(i)==copy_vector(i+1))
        {
            flag = 1;
            Cerr << copy_vector(i) << " is duplicate !!" << finl ;
        }

       // Cerr << vector[i]<<" after sort : " << copy_vector[i]<<finl;
    }

    return flag;
}


void Tool::backup_myVariables()
{
    if (Process::je_suis_maitre())
    {
    SFichier os("backup.dat");
    os <<compteur_<<finl;
    os<<myRayon<<finl;
    os<<mySigma<<finl;
    os<<myOrigine<<myLongueurs<<myNb_Noeuds;
    os<<positions_compo<<vitesses_compo;
    os <<F_old<<raideur<< e_eff<<finl;
    }

    //EcrFicCollecteBin pos ("backup.lat");
    //pos <<num_compo_;

    /*
     SFichier os("backup.dat");
    EcrFicCollecte os1("file.data");
    EcrFicPartage os2("file.data");

    */

}

void Tool::load_myVariables()
{
    LecFicDiffuse is("backup.dat");
    is >> Tool::compteur_ >>myRayon>>mySigma;
    is >> Tool::myOrigine>>Tool::myLongueurs>>myNb_Noeuds;
    is >> Tool::positions_compo>>Tool::vitesses_compo;
    is >> Tool::F_old;
    is >> Tool::raideur;
    is >> Tool::e_eff;

    //LecFicDistribueBin pis ("backup.lat");
    //pis >> num_compo_;
    /*
    EFichier is ("file.data");
    LecFicDiffuse is1 ("file.data");
    LecFicDistribue is2 ("file.data");
    */

}

//--

//--
