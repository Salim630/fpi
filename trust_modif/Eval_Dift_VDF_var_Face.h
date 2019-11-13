/****************************************************************************
* Copyright (c) 2018, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//
// File:        Eval_Dift_VDF_var_Face.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Eval_Dift_VDF_var_Face_included
#define Eval_Dift_VDF_var_Face_included

#include <Eval_Dift_VDF_var.h>
#include <Eval_VDF_Face.h>
#include <Ref_Turbulence_paroi_base.h>
#include <Ref_Mod_turb_hyd_base.h>

#include <Tool.h>
//
// .DESCRIPTION class Eval_Dift_VDF_var_Face
//
// Evaluateur VDF pour la diffusion totale (laminaire et turbulente)
// Le champ diffuse est un Champ_Face
// Le champ de diffusivite n'est pas constant.

//
// .SECTION voir aussi Eval_Dift_VDF_var

class Eval_Dift_VDF_var_Face : public Eval_Dift_VDF_var, public Eval_VDF_Face
{

public:

  inline Eval_Dift_VDF_var_Face();
  void associer_modele_turbulence(const Mod_turb_hyd_base& );
  void mettre_a_jour( );

  inline int calculer_arete_fluide() const ;
  inline int calculer_arete_paroi() const ;
  inline int calculer_arete_paroi_fluide() const ;
  inline int calculer_arete_symetrie() const ;
  inline int calculer_arete_interne() const ;
  inline int calculer_arete_mixte() const ;
  inline int calculer_fa7_sortie_libre() const ;
  inline int calculer_arete_periodicite() const;
  inline int calculer_arete_symetrie_paroi() const;
  inline int calculer_arete_symetrie_fluide() const;

  // Fonctions qui servent a calculer le flux de grandeurs scalaires
  // Elles renvoient le flux calcule

  inline double flux_fa7_sortie_libre(const DoubleTab&, int , const Neumann_sortie_libre&, int ) const;
  inline double flux_fa7_elem(const DoubleTab&, int, int, int) const ;
  inline double flux_arete_interne(const DoubleTab&, int, int, int, int) const ;
  inline double flux_arete_mixte(const DoubleTab&, int, int, int, int) const ;
  inline double flux_arete_symetrie(const DoubleTab&, int, int, int, int) const ;
  inline double flux_arete_paroi(const DoubleTab&, int, int, int, int) const ;
  inline void flux_arete_fluide(const DoubleTab&, int, int,
                                int, int, double& , double& ) const ;
  inline void flux_arete_paroi_fluide(const DoubleTab&, int, int,
                                      int, int, double& , double& ) const ;
  inline void flux_arete_periodicite(const DoubleTab&, int, int, int, int,
                                     double&, double&) const ;
  inline void flux_arete_symetrie_fluide(const DoubleTab&, int, int, int, int,
                                         double&, double&) const ;
  inline double flux_arete_symetrie_paroi(const DoubleTab&, int, int, int, int) const ;

  // Fonctions qui servent a calculer le flux de grandeurs vectorielles
  // Elles sont de type void et remplissent le tableau flux

  inline void flux_fa7_elem(const DoubleTab&, int, int, int, DoubleVect& flux) const;
  inline void flux_fa7_sortie_libre(const DoubleTab&, int , const Neumann_sortie_libre&,
                                    int, DoubleVect& flux) const;
  inline void flux_arete_interne(const DoubleTab&, int, int, int,
                                 int, DoubleVect& flux) const ;
  inline void flux_arete_mixte(const DoubleTab&, int, int, int,
                               int, DoubleVect& flux) const ;
  inline void flux_arete_symetrie(const DoubleTab&, int, int, int,
                                  int, DoubleVect& flux) const ;
  inline void flux_arete_paroi(const DoubleTab&, int, int, int,
                               int, DoubleVect& flux) const ;
  inline void flux_arete_fluide(const DoubleTab&, int, int, int,
                                int, DoubleVect& , DoubleVect&) const ;
  inline void flux_arete_paroi_fluide(const DoubleTab&, int, int,
                                      int, int, DoubleVect& , DoubleVect&) const;
  inline void flux_arete_symetrie_fluide(const DoubleTab&, int, int,
                                         int, int, DoubleVect&, DoubleVect&) const;
  inline void flux_arete_symetrie_paroi(const DoubleTab&, int, int, int,
                                        int, DoubleVect& flux) const ;

  inline double tau_tan(int face,int k) const;
  inline void flux_arete_periodicite(const DoubleTab&, int, int,
                                     int, int, DoubleVect&, DoubleVect& ) const ;

  // Fonctions qui servent a calculer les coefficients de la matrice pour des grandeurs
  // scalaires.

  inline void coeffs_fa7_elem(int, int, int, double& aii, double& ajj) const;
  inline void coeffs_fa7_sortie_libre(int, const Neumann_sortie_libre&, double& aii, double& ajj ) const;
  inline void coeffs_arete_interne(int, int, int, int, double& aii, double& ajj) const;
  inline void coeffs_arete_mixte(int, int, int, int, double& aii, double& ajj) const;
  inline void coeffs_arete_symetrie(int, int, int, int, double& aii1_2, double& aii3_4, double& ajj1_2) const;
  inline void coeffs_arete_paroi(int, int, int, int, double& aii1_2, double& aii3_4, double& ajj1_2) const;
  inline void coeffs_arete_fluide(int, int, int, int, double& aii1_2, double& aii3_4, double& ajj1_2) const;
  inline void coeffs_arete_paroi_fluide(int, int, int, int, double& aii1_2, double& aii3_4, double& ajj1_2) const;
  inline void coeffs_arete_periodicite(int, int, int, int, double& aii, double& ajj) const;
  inline void coeffs_arete_symetrie_fluide(int, int, int, int, double& aii1_2, double& aii3_4, double& ajj1_2) const;
  inline void coeffs_arete_symetrie_paroi(int, int, int, int, double& aii1_2,
                                          double& aii3_4, double& ajj1_2) const;

  // Fonctions qui servent a calculer la contribution des conditions limites
  // au second membre pour l'implicite dans le cas scalaire.

  inline double secmem_fa7_elem( int, int, int) const;
  inline double secmem_fa7_sortie_libre(int, const Neumann_sortie_libre&, int ) const;
  inline double secmem_arete_interne(int, int, int, int) const;
  inline double secmem_arete_mixte(int, int, int, int) const;
  inline double secmem_arete_symetrie(int, int, int, int) const;
  inline double secmem_arete_paroi(int, int, int, int ) const;
  inline void secmem_arete_fluide(int, int, int, int, double&, double&) const;
  inline void secmem_arete_paroi_fluide(int, int, int, int, double&, double&) const;
  inline void secmem_arete_periodicite(int, int, int, int, double&, double&) const;
  inline void secmem_arete_symetrie_fluide(int, int, int, int, double&, double&) const;
  inline double secmem_arete_symetrie_paroi(int, int, int, int ) const;

  // Fonctions qui servent a calculer les coefficients de la matrice pour des grandeurs
  // vectorielles.

  inline  void coeffs_fa7_elem(int, int, int, DoubleVect& aii, DoubleVect& ajj) const;
  inline void coeffs_fa7_sortie_libre(int , const Neumann_sortie_libre&, DoubleVect& aii, DoubleVect& ajj) const;
  inline void coeffs_arete_interne(int, int, int, int, DoubleVect& aii, DoubleVect& ajj) const;
  inline void coeffs_arete_mixte(int, int, int, int, DoubleVect& aii, DoubleVect& ajj) const;
  inline void coeffs_arete_symetrie(int, int, int, int, DoubleVect& aii1_2, DoubleVect& aii3_4, DoubleVect& ajj1_2) const;
  inline void coeffs_arete_paroi(int, int, int, int, DoubleVect& aii1_2, DoubleVect& aii3_4, DoubleVect& ajj1_2) const;
  inline void coeffs_arete_fluide(int, int, int, int, DoubleVect& aii1_2, DoubleVect& aii3_4, DoubleVect& ajj1_2) const;
  inline void coeffs_arete_paroi_fluide(int, int, int, int, DoubleVect& aii1_2, DoubleVect& aii3_4, DoubleVect& ajj1_2) const;
  inline void coeffs_arete_periodicite(int, int, int, int, DoubleVect& aii, DoubleVect& ajj) const;
  inline void coeffs_arete_symetrie_fluide(int, int, int, int,
                                           DoubleVect& aii1_2, DoubleVect& aii3_4, DoubleVect& ajj1_2) const;
  inline void coeffs_arete_symetrie_paroi(int, int, int, int, DoubleVect& aii1_2, DoubleVect& aii3_4,     DoubleVect& ajj1_2) const;

  // Fonctions qui servent a calculer la contribution des conditions limites
  // au second membre pour l'implicite dans le cas vectoriel.

  inline void secmem_fa7_elem(int, int, int, DoubleVect& flux) const;
  inline void secmem_fa7_sortie_libre(int , const Neumann_sortie_libre&, int, DoubleVect& flux) const;
  inline void secmem_arete_interne(int, int, int, int, DoubleVect& flux) const;
  inline void secmem_arete_mixte(int, int, int, int, DoubleVect& flux) const;
  inline void secmem_arete_symetrie(int, int, int, int, DoubleVect& ) const;
  inline void secmem_arete_paroi(int, int, int, int, DoubleVect& ) const;
  inline void secmem_arete_fluide(int, int, int, int, DoubleVect&, DoubleVect&) const;
  inline void secmem_arete_paroi_fluide(int, int, int, int, DoubleVect&, DoubleVect&) const;
  inline void secmem_arete_periodicite(int, int, int, int, DoubleVect&, DoubleVect&) const;
  inline void secmem_arete_symetrie_fluide(int, int, int, int, DoubleVect&, DoubleVect&) const;
  inline void secmem_arete_symetrie_paroi(int, int, int, int, DoubleVect& ) const;

private:

  REF(Mod_turb_hyd_base) le_modele_turbulence;
  REF(Turbulence_paroi_base) loipar;
  DoubleTab tau_tan_;
  DoubleTab k_;
  int indic_bas_Re, indic_lp_neg;
};

//
// Fonctions inline de la classe Eval_Dift_VDF_var_Face
//
inline double Eval_Dift_VDF_var_Face::tau_tan(int face, int k) const
{
  int nb_faces = la_zone.valeur().nb_faces();
  const ArrOfInt& ind_faces_virt_bord = la_zone.valeur().ind_faces_virt_bord();
  int f;
  if(face>=tau_tan_.dimension(0))
    f = ind_faces_virt_bord[face-nb_faces];
  else
    f=face;
  if(f>=tau_tan_.dimension_tot(0))
    {
      Cerr << "Erreur dans tau_tan " << finl;
      Cerr << "dimension : " << tau_tan_.dimension(0) << finl;
      Cerr << "dimension_tot : " << tau_tan_.dimension_tot(0) << finl;
      Cerr << "face : " << face << finl;
      Process::exit();
    }
  return tau_tan_(f,k);
}

inline Eval_Dift_VDF_var_Face::Eval_Dift_VDF_var_Face()
{}


//// calculer_arete_fluide
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_fluide() const
{
  return 1;
}


//// calculer_arete_paroi
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_paroi() const
{
  return 1;
}


//// calculer_arete_paroi_fluide
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_paroi_fluide() const
{
  return 0;
}

//// calculer_arete_periodicite
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_periodicite() const
{
  return 1;
}

//// calculer_arete_symetrie
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_symetrie() const
{
  return 0;
}

//// calculer_arete_interne
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_interne() const
{
  return 1;
}

////  calculer_arete_mixte
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_mixte() const
{
  return 1;
}

////   calculer_fa7_sortie_libre
//

inline int Eval_Dift_VDF_var_Face::calculer_fa7_sortie_libre() const
{
  return 1;
}

//// calculer_arete_symetrie_paroi
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_symetrie_paroi() const
{
  return 1;
}

//// calculer_arete_symetrie_fluide
//

inline int Eval_Dift_VDF_var_Face::calculer_arete_symetrie_fluide() const
{
  return 1;
}

// Fonctions de calcul des flux pour une inconnue scalaire

//// flux_fa7_sortie_libre
//

inline double Eval_Dift_VDF_var_Face::flux_fa7_sortie_libre(const DoubleTab&, int face,
                                                            const Neumann_sortie_libre&, int ) const
{
  double flux=0;
  /*
    double k_elem;
    int element = elem(face,0);
    if (elem(face,0) == -1) element = elem(face,1);
    if (k.nb_dim() == 1)
    k_elem = k(element);
    else if (k.nb_dim() == 2)
    k_elem = k(element,0);
    flux = - 2./3.*k_elem*surface(face);
  */
  return flux;
}

//// coeffs_fa7_sortie_libre
//

inline void Eval_Dift_VDF_var_Face::coeffs_fa7_sortie_libre(int, const Neumann_sortie_libre&,
                                                            double& aii, double& ajj ) const
{
  aii = ajj = 0;
}

//// secmem_fa7_sortie_libre
//

inline double Eval_Dift_VDF_var_Face::secmem_fa7_sortie_libre(int, const Neumann_sortie_libre&,
                                                              int ) const
{
  return 0;
}


//// flux_arete_interne
//

inline double Eval_Dift_VDF_var_Face::flux_arete_interne(const DoubleTab& inco, int fac1,
                                                         int fac2, int fac3, int fac4) const
{
  double flux;
  int ori1 = orientation(fac1);
  int ori3 = orientation(fac3);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  int elem3 = elem_(fac4,0);
  int elem4 = elem_(fac4,1);

  const int dim = 3;
  // la somme des trois orientaion = 3, pour trouver la deuxieme on fait 3 - les deux autre
  int oria = 3-ori1-ori3;

  const Zone_VF& zvf = Tool::myZone_vf_;
  const DoubleTab& v_fac=Tool::myVitesseFaces;
  const DoubleTab& v_som=Tool::myVitesseSommets;

  // algorythme pour retrouver les sommet de l'arete commune aux 2 faces,
  // Attention!! fonction seulement si les tableaux sont ordonnee

  const int nb_s_pf=4;          //nombre de sommet par face : 4 en 3D
  int som_fac1[nb_s_pf];        //les sommets de la face 1
  int som_fac2[nb_s_pf];        //les sommets de la face 2
  int som_art[2];               // vecteur pour stocker les 2 sommets par arete en 3D

// remplissage des vecteurs avec les sommets des face 1 et 2
  for(int s=0; s<nb_s_pf; s++)
    {
      som_fac1[s]=zvf.face_sommets(fac1,s);
      som_fac2[s]=zvf.face_sommets(fac2,s);
    }

  // recherche des sommets commnun aux deux faces
  //
  // TODO: cree une conectivite  face num arete commune et arete/sommets en VDF.
  // il faudra penser a optimiser ce bout de code calculer a chaque fois les sommet commun n'est franchement pas efficiant,
  // il vaut mieu cree une conectivite 2 face -> arete -> sommet une seule fois au debut du calcule

  int ii = 0, jj = 0, kk=0;
  while (ii < nb_s_pf && jj < nb_s_pf)
    {
      if (som_fac1[ii] > som_fac2[jj])
        {
          jj++;
        }
      else if (som_fac2[jj] > som_fac1[ii])
        {
          ii++;
        }
      else
        {
          // quand som_fac2[j] == som_fac1[i]
          som_art[kk]=som_fac1[ii];
          //Cerr << som_art[kk] << " ";
          ii++;
          jj++;
          kk++;
        }
    }

  // calcule de la distance entre deux sommets
  const DoubleTab& coord_som = zvf.zone().domaine().les_sommets();
  double som1=som_art[0];
  double som2=som_art[1];

  //remplissage du vecteur distance (delta x)
  double dx[dim];
  dx[ori1]=dist_face(fac3,fac4,ori1);
  dx[ori3]=dist_face(fac1,fac2,ori3);
  dx[oria]=coord_som(som2,oria)-coord_som(som1,oria);


  //remplissage du tableau des diffrences de vitesse delta U
  // compo pour vitesse u,v,w et
  DoubleTab du(dim,dim);
  for(int compo =0; compo<dim; compo++)
    {
      du(compo,ori1)=v_fac(fac4,compo)-v_fac(fac3,compo);  // coresspendance avec l'orientation des distances
      du(compo,ori3)=v_fac(fac2,compo)-v_fac(fac1,compo);  //
      du(compo,oria)=v_som(som2,compo)-v_som(som1,compo);  //
    }

//remplisage du tableau aux 9 derivee en 3D
  DoubleTab dudx(dim,dim);
  for(int compo =0; compo<dim; compo++)
    {
      for(int o=0; o<dim; o++)
        {
          dudx(compo, o) = du(compo, o) / dx[o]; //
        }
    }


  // remplissage du vecteur de normale
  double norme=0, i_norme=0, n[dim];
  DoubleTab n_elem = Tool::myNormaleInterfaceElem;

  for (int compo = 0; compo < dim; compo++)
    {
      // interpolation des composante de la normal a l'arete
      n[compo] =  0.25*(n_elem(elem1,compo) + n_elem(elem2,compo) +n_elem(elem3,compo) + n_elem(elem4,compo));
      norme+= n[compo] * n[compo];
    }
  //normalisation de la normale a l'arete
  i_norme=1./sqrt(norme);
  for (int compo = 0; compo < dim ; compo++)
    {
      n[compo]*=i_norme;
    }



// calcule de tau_c en trois partie  tau_c=t1+t2+t3 avec
//t1=(dui/dxk+duk/dxi)*nk*nj
//t2=(duj/dxk+duk/dxj)*nk*ni
//t3=-2*(duk/dxm+dum/dxk)*nx*nm*ni*nj

  double t1=0, t2=0, t3=0, tau_c=0;
  for (int k = 0; k < dim ; k++)
    {
      t1+=( dudx(ori3,k) + dudx(k,ori3) )*n[k]*n[ori1];
      t2+=( dudx(ori1,k) + dudx(k,ori1) )*n[k]*n[ori3];

      for (int l = 0; l < dim ; l++)
        {
          t3+=-2*( dudx(k,l) + dudx(l,k)  )*n[k]*n[l]*n[ori3]*n[ori1];
        }
    }
  tau_c=t1+t2+t3;

  //Cerr<<tau_c<<finl;
  // old expression
  //double visc_lam = 0.25*(dv_diffusivite(elem1) + dv_diffusivite(elem2)
  //                        +dv_diffusivite(elem3) + dv_diffusivite(elem4));
  double mu_a, mu_h;
  Tool::calcMyVisc_areteInterne(elem1,elem2,elem3,elem4,mu_a,mu_h);

  double visc_lam=mu_a; // ponderation aritmetique via l'indicatrice interpolee a l'arete
  double visc_turb = 0.25*(dv_diffusivite_turbulente(elem1) + dv_diffusivite_turbulente(elem2)
                           +dv_diffusivite_turbulente(elem3) + dv_diffusivite_turbulente(elem4));

  double tau = (inco[fac4]-inco[fac3])/dist_face(fac3,fac4,ori1); //dui/dxj=du/dz (dudx(ori3,ori1)
  double tau_tr = (inco[fac2]-inco[fac1])/dist_face(fac1,fac2,ori3);  //duj/dxi=dw/dx
  double reyn = (tau + tau_tr)*visc_turb;

  flux = 0.25*(reyn + visc_lam*(tau+tau_tr)+(mu_h-mu_a)*tau_c)*(surface(fac1)+surface(fac2))
         *(porosite(fac1)+porosite(fac2));

  // TODO: coder proprement la double ponderation de viscosite, et pourquoi pas la rondre compatible 2D
  //  pour la teste en sequenciel

  //  flux = 0.25*(reyn + visc_lam*(tau+tau_tr)+(mu_h-mu_a)*tau_c)*(surface(fac1)+surface(fac2))
  //         *(porosite(fac1)+porosite(fac2));

  //validation :
// if(fac1 == 93 && fac2==94 && fac3== 56 && fac4== 62 && Tool::myNiter==1)
  if(fac1 == 93 && fac2==94 && fac3== 56 && fac4== 62)
    {
      //combine pour ne pas imprimer le temps 0, et faire une seul impression par pas de temps

      Cerr<< " \n je suis pass ici \n" <<finl;
      Tool::print_doubletab(du,"du_arete.txt");
      Tool::print_doubletab(dudx,"dudx_arete.txt");

      //donne relatif aux somets
      int fac[]= {fac1, fac2, fac3, fac4};
      DoubleTab cf = zvf.xv();
// faces
      for (int j = 0; j < 4 ; j++)
        {
          Cerr <<j << "\t" <<fac[j] << "\t" << orientation(fac[j]) <<"\t" ;
          for (int compo = 0; compo < dim; compo++)
            {
              Cerr << cf(fac[j],compo) << "\t"  ;
            }

          for (int compo = 0; compo < dim; compo++)
            {
              Cerr  << v_fac(fac[j],compo) <<"\t" ;
            }
          Cerr << finl;
        }
      //sommets
      Cerr << finl;

      for (int j = 0; j < 2 ; j++)
        {
          Cerr <<j << "\t" <<som_art[j] << "\t" << oria <<"\t" ;
          for (int compo = 0; compo < dim; compo++)
            {
              Cerr << coord_som(som_art[j],compo) << "\t"  ;
            }

          for (int compo = 0; compo < dim; compo++)
            {
              Cerr  << v_som(som_art[j],compo) <<"\t" ;
            }
          Cerr << finl;
        }
      //--
      Cerr << finl;

      // donnes des normales
      Cerr << n_elem(elem1,0)<<"\t"<< n_elem(elem2,0)<<"\t"<< n_elem(elem3,0)<<"\t"<< n_elem(elem4,0)<<"\t"<<n[0]<< finl;
      Cerr << n_elem(elem1,1)<<"\t"<< n_elem(elem2,1)<<"\t"<< n_elem(elem3,1)<<"\t"<< n_elem(elem4,1)<<"\t"<<n[1]<< finl;
      Cerr << n_elem(elem1,2)<<"\t"<< n_elem(elem2,2)<<"\t"<< n_elem(elem3,2)<<"\t"<< n_elem(elem4,2)<<"\t"<<n[2]<< finl;
      double moy_indic = 0.25*(Tool::myIndic(elem1)+Tool::myIndic(elem2)+Tool::myIndic(elem3)+Tool::myIndic(elem4));
      Cerr << Tool::myIndic(elem1)<<"\t"<<Tool::myIndic(elem2)<<"\t"<<Tool::myIndic(elem3)<<"\t"<<Tool::myIndic(elem4)<<"\t"<<moy_indic<< finl;
      Cerr << finl;

      Cerr <<"norme:\t"<< norme<<finl;
      Cerr <<"T1:\t"<< t1<<finl;
      Cerr <<"T2:\t"<< t2<<finl;
      Cerr <<"T3:\t"<< t3<<finl;
      Cerr <<"indic:\t"<<moy_indic<<finl;
      Cerr <<"mu_h\t"<<mu_h <<finl;
      Cerr <<"mu_a\t"<<mu_a <<finl;
      Cerr <<"tau:\t"<<tau <<finl;
      Cerr <<"tau_tr:\t"<<tau_tr <<finl;
      Cerr <<"tau_c:\t"<<tau_c <<finl;
      Cerr <<"flux:\t"<<flux <<finl;
    }
  return flux;
}

//// coeffs_arete_interne
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_interne(int fac1, int fac2, int fac3, int fac4,
                                                         double& aii, double& ajj) const
{

  int ori1 = orientation(fac1);
  int ori3 = orientation(fac3);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  int elem3 = elem_(fac4,0);
  int elem4 = elem_(fac4,1);

  double visc_lam = 0.25*(dv_diffusivite(elem1) + dv_diffusivite(elem2)
                          +dv_diffusivite(elem3) + dv_diffusivite(elem4));
  double visc_turb = 0.25*(dv_diffusivite_turbulente(elem1) + dv_diffusivite_turbulente(elem2)
                           +dv_diffusivite_turbulente(elem3) + dv_diffusivite_turbulente(elem4));

  double tau = 1/dist_face(fac3,fac4,ori1);
  double tau_tr = 1/dist_face(fac1,fac2,ori3);
  double reyn = (tau + tau_tr)*visc_turb;

  aii = ajj = 0.25*(reyn + visc_lam*(tau+tau_tr))*(surface(fac1)+surface(fac2))
              *(porosite(fac1)+porosite(fac2));
}

//// secmem_arete_interne
//

inline double Eval_Dift_VDF_var_Face::secmem_arete_interne(int fac1, int fac2, int fac3, int fac4) const
{
  return 0;
}


//// flux_arete_mixte
//

// Sur les aretes mixtes les termes croises du tenseur de Reynolds
// sont nuls: il ne reste donc que la diffusion laminaire

inline double Eval_Dift_VDF_var_Face::flux_arete_mixte(const DoubleTab& inco, int fac1,
                                                       int fac2, int fac3, int fac4) const
{
  double flux=0;
  if (inco[fac4]*inco[fac3] != 0)
    {
      double visc_lam=0;
      int element;

      if ((element=elem_(fac3,0)) != -1)
        visc_lam += dv_diffusivite(element);
      if ((element=elem_(fac3,1)) != -1)
        visc_lam += dv_diffusivite(element);
      if ((element=elem_(fac4,0)) != -1)
        visc_lam += dv_diffusivite(element);
      if ((element=elem_(fac4,1)) != -1)
        visc_lam += dv_diffusivite(element);

      visc_lam/=3.0;

      int ori=orientation(fac1);
      double tau = (inco[fac4]-inco[fac3])/dist_face(fac3,fac4,ori);
      flux = 0.25*tau*(surface(fac1)+surface(fac2))*
             visc_lam*(porosite(fac1)+porosite(fac2));
    }
  return flux;
}

//// coeffs_arete_mixte
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_mixte(int fac1, int fac2, int fac3, int fac4,
                                                       double& aii, double& ajj) const
{

  if (inconnue->valeurs()[fac4]*inconnue->valeurs()[fac3] != 0)
    {
      double visc_lam=0;
      int element;

      if ((element=elem_(fac3,0)) != -1)
        visc_lam += dv_diffusivite(element);
      if ((element=elem_(fac3,1)) != -1)
        visc_lam += dv_diffusivite(element);
      if ((element=elem_(fac4,0)) != -1)
        visc_lam += dv_diffusivite(element);
      if ((element=elem_(fac4,1)) != -1)
        visc_lam += dv_diffusivite(element);

      visc_lam/=3.0;

      int ori=orientation(fac1);
      aii = ajj= 0.25*(surface(fac1)+surface(fac2))*visc_lam*
                 (porosite(fac1)+porosite(fac2))/dist_face(fac3,fac4,ori);
    }
  else
    {
      aii=ajj=0;
    }
}

//// secmem_arete_mixte
//

inline double Eval_Dift_VDF_var_Face::secmem_arete_mixte(int fac1, int fac2, int fac3, int fac4) const
{
  return 0;
}


//// flux_arete_fluide
//

inline void Eval_Dift_VDF_var_Face::flux_arete_fluide(const DoubleTab& inco, int fac1,
                                                      int fac2, int fac3, int signe,
                                                      double& flux3,double& flux1_2) const
{
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  int ori= orientation(fac3);
  double vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang1,ori,la_zcl)+
                        Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang2,ori,la_zcl));

  double dist = dist_norm_bord(fac1);
  double tau = signe * (vit_imp - inco[fac3])/dist;
  double tau_tr = (inco[fac2] - inco[fac1])/dist_face(fac1,fac2,ori);
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = (tau + tau_tr)*visc_turb;
  double coef = ((tau + tau_tr)*visc_lam + reyn);
  flux3 = coef*surf*poros;
  flux1_2 = ((tau + tau_tr)*visc_lam + reyn)*surface(fac3)*porosite(fac3);

}

//// coeffs_arete_fluide
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_fluide(int fac1, int fac2, int fac3, int signe,
                                                        double& aii1_2, double& aii3_4,
                                                        double& ajj1_2) const
{
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  int ori= orientation(fac3);

  double dist = dist_norm_bord(fac1);
  double tau = signe/dist;
  double tau_tr = 1/dist_face(fac1,fac2,ori);
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = (tau + tau_tr)*visc_turb;
  double coef = ((tau + tau_tr)*visc_lam + reyn);

  // Calcul de aii3_4
  aii3_4 = coef*surf*poros;

  // Calcul de aii1_2 et ajj1_2
  aii1_2 = ajj1_2  = ((tau + tau_tr)*visc_lam + reyn)*surface(fac3)*porosite(fac3);
}

//// secmem_arete_fluide
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_fluide(int fac1, int fac2, int fac3, int signe,
                                                        double& flux3, double& flux1_2) const
{
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  int ori= orientation(fac3);
  double vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang1,ori,la_zcl)+
                        Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang2,ori,la_zcl));

  double dist = dist_norm_bord(fac1);
  double tau = signe*vit_imp/dist;
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = tau*visc_turb;
  double coef = (tau*visc_lam + reyn);

  flux3 = coef*surf*poros;
  flux1_2 = 0;
}


//// flux_arete_paroi
//

inline double Eval_Dift_VDF_var_Face::flux_arete_paroi(const DoubleTab& inco, int fac1,
                                                       int fac2, int fac3, int signe ) const
{
  double flux;
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int ori = orientation(fac3);
  double vit = inco(fac3);
  double vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang1,ori,la_zcl)+
                        Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang2,ori,la_zcl));

  if ( (indic_bas_Re==1) || (indic_lp_neg==1) )
    {
      int elem1 = elem_(fac3,0);
      int elem2 = elem_(fac3,1);
      if (elem1==-1)
        elem1 = elem2;
      else if (elem2==-1)
        elem2 = elem1;
      double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
      double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                              + dv_diffusivite_turbulente(elem2));

      double dist = dist_norm_bord(fac1);
      double tau  = signe*(vit_imp - inco[fac3])/dist;
      double surf = 0.5*(surface(fac1)+surface(fac2));
      flux = tau*surf*(visc_lam+visc_turb);
    }
  else
    {
      // On calcule u_star*u_star*surf sur chaque partie de la facette de Qdm
      int signe_terme;
      if ( vit < vit_imp )
        signe_terme = -1;
      else
        signe_terme = 1;

      //30/09/2003  YB : influence de signe terme eliminee, signe pris en compte dans la loi de paroi
      signe_terme = 1;

      double tau1 = tau_tan(rang1,ori)*0.5*surface(fac1);
      double tau2 = tau_tan(rang2,ori)*0.5*surface(fac2);
      double coef = tau1+tau2;
      flux = signe_terme*coef;
    }
  return flux;

}

//// coeffs_arete_paroi
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_paroi(int fac1, int fac2, int fac3, int signe, double& aii1_2, double& aii3_4, double& ajj1_2) const
{
  if ( (indic_bas_Re==1) || (indic_lp_neg==1) )
    {
      int elem1 = elem_(fac3,0);
      int elem2 = elem_(fac3,1);
      if (elem1==-1)
        elem1 = elem2;
      else if (elem2==-1)
        elem2 = elem1;
      double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
      double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                              + dv_diffusivite_turbulente(elem2));

      double dist = dist_norm_bord(fac1);
      double surf = 0.5*(surface(fac1)+surface(fac2));
      aii3_4 = signe*surf*(visc_lam+visc_turb)/dist;
      aii1_2 = 0;
      ajj1_2 = 0;
    }
  else
    {
      aii3_4 = 0;
      aii1_2 = 0;
      ajj1_2 = 0;
    }
}


//// secmem_arete_paroi
//

inline double Eval_Dift_VDF_var_Face::secmem_arete_paroi(int fac1, int fac2, int fac3, int signe) const
{
  double flux;
  int ori = orientation(fac3);
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int k= orientation(fac3);
  const DoubleTab& inco = inconnue->valeurs();
  double vit = inco(fac3);
  double vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang1,k,la_zcl)+
                        Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang2,k,la_zcl));
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));

  if ( (indic_bas_Re==1) || (indic_lp_neg==1) )
    {
      double dist = dist_norm_bord(fac1);
      double tau  = signe*(vit_imp - inco[fac3])/dist;
      double surf = 0.5*(surface(fac1)+surface(fac2));
      flux = tau*surf*(visc_lam+visc_turb);
    }
  else
    {
      // On calcule u_star*u_star*surf sur chaque partie de la facette de Qdm
      int signe_terme;
      if ( vit < vit_imp )
        signe_terme = -1;
      else
        signe_terme = 1;

      //30/09/2003  YB : influence de signe terme eliminee, signe pris en compte dans la loi de paroi
      signe_terme = 1;

      double tau1 = tau_tan(rang1,ori)*0.5*surface(fac1);
      double tau2 = tau_tan(rang2,ori)*0.5*surface(fac2);
      double coef = tau1+tau2;
      flux = signe_terme*coef;
    }
  return flux;
}

//// flux_arete_paroi_fluide
//

inline void Eval_Dift_VDF_var_Face::flux_arete_paroi_fluide(const DoubleTab& inco, int fac1,
                                                            int fac2 , int fac3, int signe,
                                                            double& flux3, double& flux1_2) const
{
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  int ori= orientation(fac3);

  // On ne sait pas qui de fac1 ou de fac2 est la face de paroi
  double vit_imp;
  if (est_egal(inco[fac1],0)) // fac1 est la face de paroi
    vit_imp = Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang2,ori,la_zcl);
  else  // fac2 est la face de paroi
    vit_imp = Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang1,ori,la_zcl);

  double dist = dist_norm_bord(fac1);
  double tau = signe * (vit_imp - inco[fac3])/dist;
  double tau_tr = (inco[fac2] - inco[fac1])/dist_face(fac1,fac2,ori);
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = (tau + tau_tr)*visc_turb;
  double coef = ((tau + tau_tr)*visc_lam + reyn);
  flux3 = coef*surf*poros;
  flux1_2 = ((tau + tau_tr)*visc_lam + reyn)*surface(fac3)*porosite(fac3);

}

//// coeffs_arete_paroi_fluide
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_paroi_fluide(int fac1, int fac2, int fac3, int signe,
                                                              double& aii1_2, double& aii3_4,
                                                              double& ajj1_2) const
{
  double dist;
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  int ori= orientation(fac3);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));

  //Calcul des aii et ajj 3_4
  // On ne sait pas qui de fac1 ou de fac2 est la face de paroi

  dist = dist_norm_bord(fac1);
  double tau = signe/dist;
  double tau_tr = 1/dist_face(fac1,fac2,ori);
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = (tau + tau_tr)*visc_turb;
  double coef = ((tau + tau_tr)*visc_lam + reyn);

  aii3_4 = coef*surf*poros;
  aii1_2 = ajj1_2 =((tau + tau_tr)*visc_lam + reyn)*surface(fac3)*porosite(fac3);
}


//// secmem_arete_paroi_fluide
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_paroi_fluide(int fac1, int fac2, int fac3, int signe,
                                                              double& flux3, double& flux1_2) const
{
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam =  0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  int ori= orientation(fac3);

  // On ne sait pas qui de fac1 ou de fac2 est la face de paroi

  double vit_imp;
  if (est_egal(inconnue->valeurs()[fac1],0)) // fac1 est la face de paroi
    vit_imp = Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang2,ori,la_zcl);
  else  // fac2 est la face de paroi
    vit_imp = Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang1,ori,la_zcl);

  double dist = dist_norm_bord(fac1);
  double tau = signe*vit_imp/dist;
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = tau*visc_turb;
  double coef = (tau*visc_lam + reyn);

  flux3 = coef*surf*poros;
  flux1_2 = 0;
}

//// flux_arete_periodicite
//

inline void Eval_Dift_VDF_var_Face::flux_arete_periodicite(const DoubleTab& inco,
                                                           int fac1, int fac2, int fac3, int fac4,
                                                           double& flux3_4, double& flux1_2) const
{
  double flux;
  int ori1 = orientation(fac1);
  int ori3 = orientation(fac3);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  int elem3 = elem_(fac4,0);
  int elem4 = elem_(fac4,1);
  double dist3_4 = dist_face_period(fac3,fac4,ori1);
  double dist1_2 = dist_face_period(fac1,fac2,ori3);

  double visc_lam = 0.25*(dv_diffusivite(elem1) + dv_diffusivite(elem2)
                          +dv_diffusivite(elem3) + dv_diffusivite(elem4));
  double visc_turb = 0.25*(dv_diffusivite_turbulente(elem1) + dv_diffusivite_turbulente(elem2)
                           +dv_diffusivite_turbulente(elem3) + dv_diffusivite_turbulente(elem4));

  double tau = (inco[fac4]-inco[fac3])/dist3_4;
  double tau_tr = (inco[fac2]-inco[fac1])/dist1_2;
  double reyn = (tau + tau_tr)*visc_turb;

  flux = 0.25*(reyn + visc_lam*(tau + tau_tr))*(surface(fac1)+surface(fac2))
         *(porosite(fac1)+porosite(fac2));

  flux3_4 = flux;

  flux = 0.25*(reyn + visc_lam*(tau + tau_tr))*(surface(fac3)+surface(fac4))
         *(porosite(fac3)+porosite(fac4));

  flux1_2 = flux;
}

//// coeffs_arete_periodicite
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_periodicite(int fac1, int fac2, int fac3, int fac4,
                                                             double& aii, double& ajj) const
{
  int ori1 = orientation(fac1);
  int ori3 = orientation(fac3);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  int elem3 = elem_(fac4,0);
  int elem4 = elem_(fac4,1);
  double dist3_4 = dist_face_period(fac3,fac4,ori1);
  double dist1_2 = dist_face_period(fac1,fac2,ori3);

  double visc_lam = 0.25*(dv_diffusivite(elem1) + dv_diffusivite(elem2)
                          +dv_diffusivite(elem3) + dv_diffusivite(elem4));
  double visc_turb = 0.25*(dv_diffusivite_turbulente(elem1) + dv_diffusivite_turbulente(elem2)
                           +dv_diffusivite_turbulente(elem3) + dv_diffusivite_turbulente(elem4));

  double tau = 1/dist3_4;
  double tau_tr = 1/dist1_2;
  double reyn = (tau + tau_tr)*visc_turb;

  aii = ajj =0.25*(reyn + visc_lam*(tau + tau_tr))*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
}


//// secmem_arete_periodicite
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_periodicite(int fac1, int fac2, int fac3, int fac4,
                                                             double& flux3_4, double& flux1_2) const
{
  ;
}

//// flux_arete_symetrie
//

inline double Eval_Dift_VDF_var_Face::flux_arete_symetrie(const DoubleTab&, int,
                                                          int, int, int) const
{
  return 0;
}

//// coeffs_arete_symetrie
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_symetrie(int, int, int, int, double&, double&, double&) const
{
  ;
}

//// secmem_arete_symetrie
//

inline double Eval_Dift_VDF_var_Face::secmem_arete_symetrie(int, int, int, int) const
{
  return 0;
}


//// flux_fa7_elem
//

inline double Eval_Dift_VDF_var_Face::flux_fa7_elem(const DoubleTab& inco, int elem,
                                                    int fac1, int fac2) const
{
  double flux,k_elem;
  int ori=orientation(fac1);
  double tau = (inco[fac2]-inco[fac1])/dist_face(fac1,fac2,ori);
  double surf = 0.5*(surface(fac1)*porosite(fac1)+surface(fac2)*porosite(fac2));
  if (k_.nb_dim() == 1)
    k_elem = k_(elem);
  else if (k_.nb_dim() == 2)
    k_elem = k_(elem,0);
  else
    {
      assert(0);
      Process::exit();
      k_elem=-1;
    }
  double reyn = 2./3.*k_elem - 2.*dv_diffusivite_turbulente(elem)*tau ;
  // On verifie que les termes diagonaux du tenseur de reynolds sont bien positifs
  // Sinon on annulle :
  if (reyn < 0) reyn=0. ;


  //HMS ----
  const Zone_VF& zvf = Tool::myZone_vf_;


  int const nb_fac_pElem =	6; // a  generaliser pour la 2D
  const int dim =3;
  double fac[nb_fac_pElem];
  for (int i = 0; i < nb_fac_pElem; i++)
    {
      fac[i]=zvf.elem_faces(elem,i);
    }

  double dx[dim];
  for (int o = 0; o < dim ; o++)
    {
      dx[o]= dist_face(fac[o],fac[o+dim],o) ;
    }

  const DoubleTab& v_fac=Tool::myVitesseFaces;

  DoubleTab du(dim,dim); //du( compostante, orientation )

  for (int o = 0; o < dim ; o++)
    {
      for(int compo =0; compo<dim; compo++)
        {
          du(compo,o)=v_fac(fac[o+dim],compo)-v_fac(fac[o],compo);  // coresspendance avec l'orientation des distances
//
        }
    }
//  Tool::print_doubletab(du,"du.txt");
  //remplisage du tableau aux 9 derivee en 3D
  DoubleTab dudx(dim,dim);
  for(int compo =0; compo<dim; compo++)
    {
      for(int o=0; o<dim; o++)
        {
          dudx(compo, o) = du(compo, o) / dx[o]; //
        }
    }
//  Tool::print_doubletab(dudx,"dudx.txt");
  double n[dim];
  DoubleTab n_elem = Tool::myNormaleInterfaceElem;

  for (int compo = 0; compo < dim; compo++)
    {
      // interpolation des composante de la normal a l'arete
      n[compo] =  n_elem(elem,compo);
    }

  double t1=0, t2=0, tau_c=0;
  for (int k = 0; k < dim ; k++)
    {
      t1+=2*( dudx(ori,k) + dudx(k,ori) )*n[k]*n[ori];

      for (int m = 0; m < dim ; m++)
        {
          t2+=-2*( dudx(k,m) + dudx(m,k)  )*n[k]*n[m]*n[ori]*n[ori];
        }
    }
  tau_c=t1+t2;
////
///



  //---
  // Cerr <<tau_c <<finl;
  double mu_a, mu_h;
  Tool::calcMyVisc_fa7Elem(elem, mu_a, mu_h);
  //---
  flux = (-reyn + mu_a*2*tau+(mu_h-mu_a)*tau_c) * surf ;

  if(elem==3 && fac1==39 &&fac2==93 )
    {
      //  //impression pour verification
//// donnees relatif aux faces
      Tool::myOldTime=Tool::myTime;
      Tool::print_doubletab(du,"du_fa7.txt");
      Tool::print_doubletab(dudx,"dudx_fa7.txt");
//
      DoubleTab cf = zvf.xv();

      for (int j = 0; j < 6 ; j++)
        {
          Cerr <<j<< "\t"<<fac[j] << "\t" << orientation(fac[j]) <<"\t" ;
          for (int compo = 0; compo < 3; compo++)
            {
              Cerr << cf(fac[j],compo) << "\t" ;
            }

          for (int compo = 0; compo < 3; compo++)
            {
              Cerr << v_fac(fac[j],compo) <<"\t" ;
            }

          Cerr << finl;
        }


      // donnes des normales
      Cerr <<"nx:\t" <<n[0]<< finl;
      Cerr <<"ny:\t" <<n[1]<< finl;
      Cerr <<"nz:\t" <<n[2]<< finl;
      Cerr <<"T1:\t"<< t1<<finl;
      Cerr <<"T2:\t"<< t2<<finl;
      Cerr <<"indic:\t"<<Tool::myIndic(elem)<<finl;
      Cerr <<"mu_h\t"<<mu_h <<finl;
      Cerr <<"mu_a\t"<<mu_a <<finl;
      Cerr <<"tau:\t"<<tau <<finl;
      Cerr <<"tau_c:\t"<<tau_c <<finl;
      Cerr <<"flux:\t"<<flux <<finl;

    }
  return flux;
}


//// coeffs_fa7_elem
//

inline void Eval_Dift_VDF_var_Face::coeffs_fa7_elem(int elem ,int fac1, int fac2, double& aii, double& ajj) const
{
  double k_elem;
  int ori=orientation(fac1);
  double tau = 1/dist_face(fac1,fac2,ori);
  double surf = 0.5*(surface(fac1)*porosite(fac1)+surface(fac2)*porosite(fac2));
  if (k_.nb_dim() == 1)
    k_elem = k_(elem);
  else if (k_.nb_dim() == 2)
    k_elem = k_(elem,0);
  else
    {
      assert(0);
      Process::exit();
      k_elem=-1;
    }
  double reyn = 2./3.*k_elem - 2.*dv_diffusivite_turbulente(elem)*tau ;
  // On verifie que les termes diagonaux du tenseur de reynolds sont bien positifs
  // Sinon on annulle :
  if (reyn < 0) reyn=0. ;

  aii = ajj =(-reyn +dv_diffusivite(elem)*tau) * surf;
}


//// secmem_fa7_elem
//

inline double Eval_Dift_VDF_var_Face::secmem_fa7_elem(int elem ,int fac1, int fac2) const
{
  return 0;
}

//// flux_arete_symetrie_fluide
//

inline void Eval_Dift_VDF_var_Face::flux_arete_symetrie_fluide(const DoubleTab& inco, int fac1,
                                                               int fac2, int fac3, int signe,
                                                               double& flux3,double& flux1_2) const
{
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  int ori= orientation(fac3);
  double vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord_sym(inco,inconnue->temps(),rang1,ori,la_zcl)+
                        Champ_Face_get_val_imp_face_bord_sym(inco,inconnue->temps(),rang2,ori,la_zcl));

  double dist = dist_norm_bord(fac1);
  double tau = signe * (vit_imp - inco[fac3])/dist;
  double tau_tr = (inco[fac2] - inco[fac1])/dist_face(fac1,fac2,ori);
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = (tau + tau_tr)*visc_turb;
  double coef = ((tau + tau_tr)*visc_lam + reyn);
  flux3 = coef*surf*poros;
  flux1_2 = ((tau + tau_tr)*visc_lam + reyn)*surface(fac3)*porosite(fac3);
}

//// coeffs_arete_symetrie_fluide
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_symetrie_fluide(int fac1, int fac2, int fac3, int signe,
                                                                 double& aii1_2, double& aii3_4,
                                                                 double& ajj1_2) const
{
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  int ori= orientation(fac3);

  double dist = dist_norm_bord(fac1);
  double tau = signe/dist;
  double tau_tr = 1/dist_face(fac1,fac2,ori);
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = (tau + tau_tr)*visc_turb;
  double coef = ((tau + tau_tr)*visc_lam + reyn);

  // Calcul de aii3_4
  aii3_4 = coef*surf*poros;

  // Calcul de aii1_2 et ajj1_2
  aii1_2 = ajj1_2  = ((tau + tau_tr)*visc_lam + reyn)*surface(fac3)*porosite(fac3);
}

//// secmem_arete_symetrie_fluide
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_symetrie_fluide(int fac1, int fac2, int fac3, int signe,
                                                                 double& flux3, double& flux1_2) const
{
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam = 0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  int ori= orientation(fac3);
  double vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang1,ori,la_zcl)+
                        Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang2,ori,la_zcl));

  double dist = dist_norm_bord(fac1);
  double tau = signe*vit_imp/dist;
  double surf = 0.5*(surface(fac1)+surface(fac2));
  double poros = 0.5*(porosite(fac1)+porosite(fac2));
  double reyn = tau*visc_turb;
  double coef = (tau*visc_lam + reyn);

  flux3 = coef*surf*poros;
  flux1_2 = 0;
}

//// flux_arete_symetrie_paroi
//

inline double Eval_Dift_VDF_var_Face::flux_arete_symetrie_paroi(const DoubleTab& inco, int fac1,
                                                                int fac2, int fac3, int signe ) const
{
  double flux;
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam =  0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));

  int ori = orientation(fac3);
  //  double vit = inco(fac3);
  double vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord_sym(inco,inconnue->temps(),rang1,ori,la_zcl)+
                        Champ_Face_get_val_imp_face_bord_sym(inco,inconnue->temps(),rang2,ori,la_zcl));
  if ((indic_bas_Re==1) || (indic_lp_neg==1))
    {
      double dist = dist_norm_bord(fac1);
      double tau  = signe*(vit_imp - inco[fac3])/dist;
      double surf = 0.5*(surface(fac1)+surface(fac2));
      flux = tau*surf*(visc_lam+visc_turb);
    }
  else
    {
      // solution retenue pour le calcul du frottement sachant que sur l'une des faces
      // tau_tan vaut 0 car c'est une face qui porte une condition de symetrie
      // On calcule u_star*u_star*surf sur chaque partie de la facette de Qdm
      double tau1 = tau_tan(rang1,ori)*0.5*surface(fac1);
      double tau2 = tau_tan(rang2,ori)*0.5*surface(fac2);
      double coef = tau1+tau2;
      flux = coef;
    }
  return flux;
}

//// coeffs_arete_symetrie_paroi
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_symetrie_paroi(int fac1, int fac2,
                                                                int fac3, int signe,
                                                                double& aii1_2, double& aii3_4, double& ajj1_2) const
{
  if ( (indic_bas_Re==1) || (indic_lp_neg==1) )
    {
      int elem1 = elem_(fac3,0);
      int elem2 = elem_(fac3,1);
      double visc_lam =  0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
      double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                              + dv_diffusivite_turbulente(elem2));

      double dist = dist_norm_bord(fac1);
      double surf = 0.5*(surface(fac1)+surface(fac2));
      aii3_4 = signe*surf*(visc_lam+visc_turb)/dist;
      aii1_2 = 0;
      ajj1_2 = 0;
    }
  else
    {
      aii3_4 = 0;
      aii1_2 = 0;
      ajj1_2 = 0;
    }
}


//// secmem_arete_symetrie_paroi
//

inline double Eval_Dift_VDF_var_Face::secmem_arete_symetrie_paroi(int fac1, int fac2, int fac3, int signe) const
{
  double flux;
  int rang1 = (fac1-premiere_face_bord);
  int rang2 = (fac2-premiere_face_bord);
  int ori = orientation(fac3);
  int elem1 = elem_(fac3,0);
  int elem2 = elem_(fac3,1);
  double visc_lam =  0.5*(dv_diffusivite(elem1)+dv_diffusivite(elem2));
  double visc_turb = 0.5*(dv_diffusivite_turbulente(elem1)
                          + dv_diffusivite_turbulente(elem2));
  double vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang1,ori,la_zcl)+
                        Champ_Face_get_val_imp_face_bord(inconnue->temps(),rang2,ori,la_zcl));
  if ( (indic_bas_Re==1) || (indic_lp_neg==1) )
    {
      const DoubleTab& inco = inconnue->valeurs();
      double dist = dist_norm_bord(fac1);
      double tau  = signe*(vit_imp - inco[fac3])/dist;
      double surf = 0.5*(surface(fac1)+surface(fac2));
      flux = tau*surf*(visc_lam+visc_turb);
    }
  else
    {
      // solution retenue pour le calcul du frottement sachant que sur l'une des faces
      // tau_tan vaut 0 car c'est une face qui porte une condition de symetrie
      // On calcule u_star*u_star*surf sur chaque partie de la facette de Qdm
      double tau1 = tau_tan(rang1,ori)*0.5*surface(fac1);
      double tau2 = tau_tan(rang2,ori)*0.5*surface(fac2);
      double coef = tau1+tau2;
      flux = coef;
    }
  return flux;
}

// Fonctions de calcul des flux pour une inconnue vectorielle
// Elles ne sont pas codees pour l'instant

//// flux_fa7_sortie_libre
//

inline void Eval_Dift_VDF_var_Face::flux_fa7_sortie_libre(const DoubleTab&, int ,
                                                          const Neumann_sortie_libre&,
                                                          int, DoubleVect& ) const
{
  // A coder !
}

////coeffs_fa7_sortie_libre
//

inline void Eval_Dift_VDF_var_Face::coeffs_fa7_sortie_libre(int , const Neumann_sortie_libre&,
                                                            DoubleVect&, DoubleVect& ) const
{
  // A Coder!
}

////secmem_fa7_sortie_libre
//

inline void Eval_Dift_VDF_var_Face::secmem_fa7_sortie_libre(int , const Neumann_sortie_libre&,
                                                            int, DoubleVect&  ) const
{
  // A Coder!
}


//// flux_arete_interne
//

inline void Eval_Dift_VDF_var_Face::flux_arete_interne(const DoubleTab&, int, int,
                                                       int, int, DoubleVect& ) const
{
  // A coder !
}

//// coeffs_arete_interne
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_interne(int , int , int , int ,
                                                         DoubleVect&, DoubleVect& ) const
{
  // A coder!
}

//// secmem_arete_interne
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_interne(int ,
                                                         int , int , int ,
                                                         DoubleVect& ) const
{
  // A coder!
}

//// flux_arete_mixte
//

inline void Eval_Dift_VDF_var_Face::flux_arete_mixte(const DoubleTab&, int, int,
                                                     int, int, DoubleVect& ) const
{
  // A coder !
}

//// coeffs_arete_mixte
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_mixte(int , int , int , int ,
                                                       DoubleVect& , DoubleVect& ) const
{
  // A coder!
}

//// secmem_arete_mixte
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_mixte(int, int , int , int ,
                                                       DoubleVect& ) const
{
  // A coder!
}

//// flux_arete_fluide
//

inline void Eval_Dift_VDF_var_Face::flux_arete_fluide(const DoubleTab&, int, int,
                                                      int, int, DoubleVect&, DoubleVect& ) const
{
  // A coder !
}

//// coeffs_arete_fluide
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_fluide(int ,int ,
                                                        int , int , DoubleVect&, DoubleVect&, DoubleVect&) const
{
  // A Coder!
}

//// secmem_arete_fluide
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_fluide(int ,int ,int , int ,
                                                        DoubleVect& , DoubleVect& ) const
{
  // A Coder!
}

//// flux_arete_paroi
//

inline void Eval_Dift_VDF_var_Face::flux_arete_paroi(const DoubleTab&, int, int,
                                                     int, int, DoubleVect& ) const
{
  // A coder !
}

//// coeffs_arete_paroi
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_paroi(int , int , int , int ,
                                                       DoubleVect&, DoubleVect&, DoubleVect&) const
{
  // A coder!
}

//// secmem_arete_paroi
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_paroi(int, int , int , int ,
                                                       DoubleVect& ) const
{
  // A coder!
}

//// flux_arete_paroi_fluide
//

inline void Eval_Dift_VDF_var_Face::flux_arete_paroi_fluide(const DoubleTab&, int, int,
                                                            int, int, DoubleVect&, DoubleVect&) const
{
  // A coder !
}

//// coeffs_arete_paroi_fluide
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_paroi_fluide(int, int , int , int ,
                                                              DoubleVect&, DoubleVect&, DoubleVect& ) const
{
  // A Coder!
}


//// secmem_arete_paroi_fluide
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_paroi_fluide(int, int , int , int ,
                                                              DoubleVect& , DoubleVect& ) const
{
  // A Coder!
}

//// flux_arete_periodicite
//

inline void Eval_Dift_VDF_var_Face::flux_arete_periodicite(const DoubleTab& , int ,
                                                           int , int , int ,
                                                           DoubleVect& , DoubleVect& ) const
{
  ;
}

//// coeffs_arete_periodicite
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_periodicite(int , int , int , int ,
                                                             DoubleVect& , DoubleVect& ) const
{
  // A Coder!
}

//// secmem_arete_periodicite
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_periodicite(int , int , int , int ,
                                                             DoubleVect& , DoubleVect& ) const
{
  // A Coder!
}

//// flux_arete_symetrie
//

inline void Eval_Dift_VDF_var_Face::flux_arete_symetrie(const DoubleTab&, int, int,
                                                        int, int, DoubleVect& ) const
{
  // A coder !
}

//// coeffs_arete_symetrie
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_symetrie(int, int, int, int, DoubleVect&, DoubleVect&, DoubleVect&) const
{
  // A Coder!
}

//// secmem_arete_symetrie
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_symetrie(int, int,
                                                          int, int, DoubleVect& ) const
{
  // A Coder!
}


//// flux_fa7_elem
//

inline void Eval_Dift_VDF_var_Face::flux_fa7_elem(const DoubleTab&, int, int,
                                                  int, DoubleVect& ) const
{
  // A coder !
}

//// coeffs_fa7_elem
//

inline void Eval_Dift_VDF_var_Face::coeffs_fa7_elem(int , int ,int , DoubleVect&, DoubleVect& ) const
{
  // A Coder!
}

//// secmem_fa7_elem
//

inline void Eval_Dift_VDF_var_Face::secmem_fa7_elem(int , int ,
                                                    int , DoubleVect& ) const
{
  // A Coder!
}

//// flux_arete_symetrie_fluide
//

inline void Eval_Dift_VDF_var_Face::flux_arete_symetrie_fluide(const DoubleTab&, int, int,
                                                               int, int, DoubleVect&, DoubleVect& ) const
{
  // A coder !
}

inline void Eval_Dift_VDF_var_Face::coeffs_arete_symetrie_fluide(int ,int ,
                                                                 int , int , DoubleVect&, DoubleVect&, DoubleVect&) const
{
  // A Coder!
}

//// secmem_arete_symetrie_fluide
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_symetrie_fluide(int ,int ,int , int ,
                                                                 DoubleVect& , DoubleVect& ) const
{
  // A Coder!
}

//// flux_arete_symetrie_paroi
//

inline void Eval_Dift_VDF_var_Face::flux_arete_symetrie_paroi(const DoubleTab&, int, int,
                                                              int, int, DoubleVect& ) const
{
  // EMPTY //
  ;
}

//// coeffs_arete_symetrie_paroi
//

inline void Eval_Dift_VDF_var_Face::coeffs_arete_symetrie_paroi(int , int , int , int ,
                                                                DoubleVect&, DoubleVect&, DoubleVect&) const
{
  // A coder!
}

//// secmem_arete_symetrie_paroi
//

inline void Eval_Dift_VDF_var_Face::secmem_arete_symetrie_paroi(int, int , int , int ,
                                                                DoubleVect& ) const
{
  // A coder!
}

#endif
