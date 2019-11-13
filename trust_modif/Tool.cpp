//
// Created by salim on 19/10/19.
//

#include <Tool.h>
#include <Sortie_Fichier_base.h>

//declaration des membre donnees
std::string Tool::myCode="commicode 010dpv";

double Tool::myMuPhase1=-1;
double Tool::myMuPhase0=-1;
DoubleVect Tool::myIndic;
DoubleTab Tool::myNormaleInterfaceElem;
DoubleTab Tool::myVitesse;
DoubleTab Tool::myVitesseSommets;
DoubleTab Tool::myVitesseFaces;

double Tool::myTime;
double Tool::myOldTime=0.0;
int Tool::myNiter=-1;

REF(Zone_VDF)    Tool::ma_zone_VDF_;
REF(Zone_VF)     Tool::myZone_vf_;

//implementation des membre fonctions

//utilitaire
// penderation harmonique de la viscosié a l'arete
double Tool::calcMyViscLam(int elem1,int elem2, int elem3, int elem4){

	  double indicArete = 0.25*(myIndic[elem1] + myIndic[elem2] +myIndic[elem3] + myIndic[elem4]);

	  double myViscLam  = (myMuPhase0 * myMuPhase1)/(myMuPhase1 - indicArete * ( myMuPhase1 - myMuPhase0));

    //printf("elem{1,2,3,4} ={%d|%d|%d|%d}\n",elem1,elem2,elem3,elem4);
    //printf("indic{1,2,3,4={%f|%f|%f|%f}\n",myIndic[elem1],myIndic[elem2],myIndic[elem3],myIndic[elem4]);
    //printf("indic arete=%f; mu0=%f, mu1=%f, visclam=%f\n", indicArete,myMuPhase0,myMuPhase1,myViscLam);

	return myViscLam;
}

void Tool::calcMyVisc_areteInterne(int elem1,int elem2, int elem3, int elem4, double& mu_a, double& mu_h) {

    double indicArete = 0.25 * (myIndic[elem1] + myIndic[elem2]
                                + myIndic[elem3] + myIndic[elem4]);

    // ponderation aritmetique
     mu_a = indicArete * (myMuPhase1 - myMuPhase0) + myMuPhase0;
    // ponderation harmonique
     mu_h = (myMuPhase0 * myMuPhase1) / (myMuPhase1 - indicArete * (myMuPhase1 - myMuPhase0));

    //printf("elem{1,2,3,4} ={%d|%d|%d|%d}\n",elem1,elem2,elem3,elem4);
    //printf("indic{1,2,3,4={%f|%f|%f|%f}\n",myIndic[elem1],myIndic[elem2],myIndic[elem3],myIndic[elem4]);
    //printf("indic arete=%f; mu0=%f, mu1=%f, visclam=%f\n", indicArete,myMuPhase0,myMuPhase1,myViscLam);

}

void Tool::calcMyVisc_fa7Elem(int elem, double& mu_a, double& mu_h) {

    double indicArete = myIndic[elem];

    // ponderation aritmetique
    mu_a = indicArete * (myMuPhase1 - myMuPhase0) + myMuPhase0;
    // ponderation harmonique
    mu_h = (myMuPhase0 * myMuPhase1) / (myMuPhase1 - indicArete * (myMuPhase1 - myMuPhase0));

    //printf("elem{1,2,3,4} ={%d|%d|%d|%d}\n",elem1,elem2,elem3,elem4);
    //printf("indic{1,2,3,4={%f|%f|%f|%f}\n",myIndic[elem1],myIndic[elem2],myIndic[elem3],myIndic[elem4]);
    //printf("indic arete=%f; mu0=%f, mu1=%f, visclam=%f\n", indicArete,myMuPhase0,myMuPhase1,myViscLam);

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

 DoubleTab& Tool::myCalculer_vitesse_faces(DoubleTab& v_faces_stockage)
{
    const Zone_VDF&    zone_VDF      = ma_zone_VDF_.valeur();
    const IntVect&     orientation   = zone_VDF.orientation();
    const IntTab&      faces_voisins = zone_VDF.face_voisins();
    const DoubleVect& volumes       = zone_VDF.volumes();  // volumes des elements
    const IntTab&      elem_faces    = zone_VDF.elem_faces();
    const DoubleTab&   v_faces  = myVitesse;
    const int       dim      = Objet_U::dimension;
    const int       nb_faces = v_faces.dimension(0);
    v_faces_stockage.resize(nb_faces, dim);
    int i_face;
    ArrOfDouble composante_vitesse(3);
    for (i_face = 0; i_face < nb_faces; i_face++)
    {
        const int orientation_face = orientation(i_face);
        composante_vitesse=0;
        int composante;

        // Numeros des deux elements voisins de la face (-1 si face de bord)
        int elem[2];
        elem[0] = faces_voisins(i_face, 0);
        elem[1] = faces_voisins(i_face, 1);

        // Volumes de ces deux elements (0. si pas de voisin)
        double volume_elem[2] = {0., 0.};
        if (elem[0] >= 0)
            volume_elem[0] = volumes(elem[0]);
        if (elem[1] >= 0)
            volume_elem[1] = volumes(elem[1]);

        const double i_volume_total = 1. / (volume_elem[0] + volume_elem[1]);

        for (composante = 0; composante < dim; composante++)
        {
            if (composante == orientation_face)
            {
                composante_vitesse[composante] = v_faces(i_face);
            }
            else
            {
                // Calcul de la moyenne des vitesses sur les faces voisines
                // qui ont la bonne orientation:
                composante_vitesse[composante] = 0.;
                int i_elem;
                for (i_elem = 0; i_elem < 2; i_elem++)
                {
                    if (elem[i_elem] >= 0)
                    {
                        const int element = elem[i_elem];
                        const int face1 = elem_faces(element, composante);
                        const int face2 = elem_faces(element, composante + dim);
                        const double v1 = v_faces(face1);
                        const double v2 = v_faces(face2);
                        const double p = volume_elem[i_elem]; // Ponderation par le volume
                        composante_vitesse[composante] += (v1 + v2) * 0.5 * p;
                    }
                }
                composante_vitesse[composante] *= i_volume_total;
            }
        }
        for (composante = 0; composante < dim; composante++)
        {
            v_faces_stockage(i_face, composante) = composante_vitesse[composante];
        }
    }
    myVitesseFaces=v_faces_stockage;
    return v_faces_stockage;
}

 void Tool::print_doubletab( const DoubleTab& tab, Nom nom_fichier ){
	 //std::string nom_fichier="mavariable.txt";
     //Nom nom_fic = nom_fichier;
	 Cerr<<"\t\t!! writing in file: "<< nom_fichier <<finl;
	 ofstream fout;
 	 fout.open(nom_fichier,ios::app);
     fout << "Time: " << Tool::myTime<< std::endl;
     fout.setf(ios::scientific);
    int nb_dim=tab.nb_dim();

    switch (nb_dim){
    case 1:
    {
        int n=tab.dimension(0);
        for(int i=0; i<n; i++){
    	    fout <<i<<"\t\t "<< tab(i) << std::endl;
        }
    }
    	break;
    case 2:
    {
        int n=tab.dimension(0);
        int m=tab.dimension(1);
        for(int i=0; i<n; i++){
    	    fout <<i<<"\t\t ";
    	    for(int j=0; j<m; j++){
    	      fout << tab(i,j) <<"\t\t ";
    	    }
    	    fout<< std::endl;
        }
    }
    	break;
    default :
    	Cerr<<"not yet programmed for this size"<<std::endl;
    }





      //fout << "Done!" << std::endl;
      fout.close();
 }
