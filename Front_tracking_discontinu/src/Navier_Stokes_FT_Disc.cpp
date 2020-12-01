/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Navier_Stokes_FT_Disc.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/55
//
//////////////////////////////////////////////////////////////////////////////
#include <Navier_Stokes_FT_Disc.h>
#include <Probleme_FT_Disc_gen.h>
#include <Modele_turbulence_hyd_nul.h>
#include <Discret_Thyd.h>
#include <Ref_Navier_Stokes_FT_Disc.h>
#include <Operateur_Diff_base.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>
#include <Assembleur_base.h>
#include <Ref_Champ_base.h>
#include <Schema_Temps_base.h>
#include <Zone_VDF.h>
#include <Debog.h>
// Ces includes seront a retirer quand on aura clairement separe les operations
// specifiques au VDF et VEF
#include <Zone_VF.h>
#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Ref_Convection_Diffusion_Temperature_FT_Disc.h>
#include <Terme_Source_Constituant_Vortex_VEF_Face.h>
#include <DoubleTrav.h>
#include <Matrice_Morse_Sym.h>
#include <Matrice_Bloc.h>
#include <Param.h>
#include <Tool.h>
#include <Statistiques.h>
#include <Connex_components_FT.h>
#include <communications.h>
#include <Connex_components.h>
#include <fstream>
#include <iomanip>
#include <SFichier.h>
#include <EFichier.h>

// #include <chrono>
// using namespace std::chrono;

Implemente_instanciable_sans_constructeur_ni_destructeur(Navier_Stokes_FT_Disc,"Navier_Stokes_FT_Disc",Navier_Stokes_Turbulent);

Implemente_ref(Navier_Stokes_FT_Disc);
Declare_vect(REF(Transport_Interfaces_FT_Disc));
Implemente_vect(REF(Transport_Interfaces_FT_Disc));

class Navier_Stokes_FT_Disc_interne
{
public:
  Navier_Stokes_FT_Disc_interne() :
    correction_courbure_ordre_(0), // Par defaut, pas de prise en compte de la courbure pour corriger le champ etendu delta_vitesse
    mpoint_inactif(0),   // Par defaut, mpoint cree un saut de vitesse
    mpointv_inactif(0),   // Par defaut, mpointv  cree un saut de vitesse
    matrice_pression_invariante(0),   // Par defaut, recalculer la matrice pression
    clipping_courbure_interface(1e40),// Par defaut, pas de clipping
    terme_gravite_(GRAVITE_GRAD_I),   // Par defaut terme gravite ft sans courants parasites
    is_explicite(1),                  // Par defaut, calcul explicite de vpoint etape predicition
    is_penalized(0),                  // Par defaut, pas de penalisation L2 du forcage
    eta(1.0),                         // Par defaut, coefficient de penalisation L2 = 1.
    p_ref_pena(-1.e40),               // Par defaut, pas de penalisation L2 de la pression sinon valeur reference
    is_pfl_flottant(0),               // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
    x_pfl_imp(-1.e40),                // Par defaut, x, y, z du point de modification de la pression fluide
    y_pfl_imp(-1.e40),
    z_pfl_imp(-1.e40)
  {};
  int correction_courbure_ordre_;
  int mpoint_inactif;
  int mpointv_inactif;

  Champ_Fonc second_membre_projection;
  Champ_Fonc derivee_u_etoile;
  Champ_Fonc gradient_pression;
  Champ_Fonc terme_diffusion;
  Champ_Fonc terme_convection;
  Champ_Fonc terme_source;
  Champ_Fonc terme_source_interfaces;

  Champ_Fonc indicatrice_p1b;
  Champ_Fonc gradient_indicatrice;
  Champ_Fonc potentiel_faces;
  Champ_Fonc potentiel_elements;
  // delta_u_interface = la partie "saut de vitesse" du champ de vitesse a l'interface
  Champ_Inc delta_u_interface;
  Champ_Fonc laplacien_d;
  Champ_Fonc mpoint;
  Champ_Fonc mpoint_vap;
  // Variation temporelle indicatrice de phase
  Champ_Fonc derivee_temporelle_indicatrice;

  //-----------------
  Champ_Fonc terme_source_collisions;
  Champ_Fonc  num_compo;
  Motcle modele_collisions_;
  //-----------------

  LIST(REF(Champ_base)) liste_champs_compris;

  // Si matrice_pression_invariante != 0,
  //   on ne recalcule pas la matrice de pression a chaque pas de temps.
  int matrice_pression_invariante;
  // Si on veut ajouter une interface a vitesse imposee :
  //  reference a l'equation de transport correspondante :
  VECT(REF(Transport_Interfaces_FT_Disc)) ref_eq_interf_vitesse_imposee;
  // Si le fluide est diphasique, c'est l'indicatrice de l'equation suivante
  // qui est utilisee pour determiner les proprietes du fluide:
  // (masse volumique, viscosite, tension superficielle, ...)
  REF(Transport_Interfaces_FT_Disc) ref_eq_interf_proprietes_fluide;
  // Si le fluide est diphasique, la reference au fluide:
  REF(Fluide_Diphasique) ref_fluide_diphasique;

  REF(Convection_Diffusion_Temperature_FT_Disc) ref_equation_mpoint_;
  REF(Convection_Diffusion_Temperature_FT_Disc) ref_equation_mpoint_vap_;

  // Valeur maximale de courbure autorisee pour calculer le
  // terme source de tension de surface (clipping si valeur superieur)
  double clipping_courbure_interface;

  enum Terme_Gravite { GRAVITE_RHO_G, GRAVITE_GRAD_I };
  Terme_Gravite terme_gravite_;
  Noms equations_concentration_source_fluide_;
  // Si is_explicite != 0,
  //   on calcul vpoint de facon explicite dans l etape de prediction des vitesses.
  int is_explicite;
  // Si is_penalized != 0,
  //   on penalise L2 le terme de forcage.
  int is_penalized;
  // Valeur de l'inverse du coefficient de penalisation L2 du terme de forcage.
  double eta;
  // Valeur pour la penalisation L2 de la pression.
  double p_ref_pena;
  // Point de penalisation L2 de la pression du fluide
  int is_pfl_flottant; // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
  double x_pfl_imp;
  double y_pfl_imp;
  double z_pfl_imp;
};

// Description:
//  Calcul de champ_rho_faces_ et champ_rho_elem_ en fonction de l'indicatrice:
//   rho_elem_ = indicatrice * ( rho(phase_1) - rho(phase_0) ) + rho(phase_0)
//   rho_faces_ = 0.5 * (rho_elem(voisin_0) + rho_elem(voisin_1))
//  Le calcul des viscosites cinematique et dynamique est pour l'instant le suivant:
//   nu_elem = indicatrice * ( nu(phase_1) - nu(phase_0) ) + nu(phase_0)
//   mu_elem = nu_elem * rho_elem
static void FT_disc_calculer_champs_rho_mu_nu_dipha(const Zone_dis_base&      zone_dis_base,
                                                    const Fluide_Diphasique& fluide,
                                                    const DoubleVect&         indicatrice_elem,
                                                    DoubleVect& rho_elem,
                                                    DoubleVect& nu_elem,
                                                    DoubleVect& mu_elem,
                                                    DoubleVect& rho_faces)
{
  const Fluide_Incompressible& phase_0 = fluide.fluide_phase(0);
  const Fluide_Incompressible& phase_1 = fluide.fluide_phase(1);
  const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique().valeurs();
  const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique().valeurs();
  const double rho_phase_0 = tab_rho_phase_0(0,0);
  const double rho_phase_1 = tab_rho_phase_1(0,0);
  const double delta_rho = rho_phase_1 - rho_phase_0;
  const DoubleTab& tab_nu_phase_0 = phase_0.viscosite_cinematique().valeur().valeurs();
  const DoubleTab& tab_nu_phase_1 = phase_1.viscosite_cinematique().valeur().valeurs();
  const double nu_phase_0 = tab_nu_phase_0(0,0);
  const double nu_phase_1 = tab_nu_phase_1(0,0);
  const double delta_nu = nu_phase_1 - nu_phase_0;
  const double mu_phase_0 = nu_phase_0 * rho_phase_0;
  const double mu_phase_1 = nu_phase_1 * rho_phase_1;
  const double delta_mu = mu_phase_1 - mu_phase_0;
  double mu = 0.;
  const int formule_mu = fluide.formule_mu();

  Tool::setMyMuPhase0(mu_phase_0);
  Tool::setMyMuPhase1(mu_phase_1);
  Tool::setMyIndic(indicatrice_elem);
  // Calcul de rho, nu, mu aux elements
  {
    const int n = indicatrice_elem.size();
    assert(n == rho_elem.size());
    assert(n == nu_elem.size());
    assert(n == mu_elem.size());
    for (int i = 0; i < n; i++)
      {
        const double indic = indicatrice_elem[i];
        const double rho = indic * delta_rho + rho_phase_0;
        rho_elem[i] = rho;
        const double nu  = indic * delta_nu  + nu_phase_0; // HMS: faux, car nu = mu/rho
        nu_elem[i]  = nu;
        Tool::formule_mu=formule_mu;
        switch(formule_mu)
          {
          case 0: // standard default method (will be removed)
            {
              //Cerr << "!! Attention !! mu phase_1 ( fluide ) sur elem diphasiques !!" << finl;
              mu  = indic == 0 ? mu_phase_0 : mu_phase_1 ;

            }
            break;
          case 1: // Arithmetic average
            {
              //Cerr << "arithmetic" << finl;
              mu  = indic * delta_mu  + mu_phase_0;
            }
            break;
          case 2: // Harmonic average
            {
              //  Cerr << "harmonic" << finl;
              mu  = (mu_phase_0 * mu_phase_1)/(mu_phase_1 - indic * delta_mu);
            }
            break;
          default:
            {
              Cerr << "The method specified for formule_mu in not recognized. \n" << finl;
              Cerr << "you can choose : standard, arithmetic or harmonic. \n" << finl;
              //Cerr << "We should not be here Navier_Stokes_FT_Disc::FT_disc_calculer_champs_rho_mu_nu_dipha" << finl;
              Process::exit();
            }
          }
        mu_elem[i]  = mu;
        nu_elem[i]  = mu / rho; // correction de nu, fin HMS
      }

    rho_elem.echange_espace_virtuel();
    mu_elem.echange_espace_virtuel();
    nu_elem.echange_espace_virtuel();
  }

// Calcul de rho aux faces (on suppose que la vitesse est aux faces)
  {
    assert(rho_elem.size() == zone_dis_base.nb_elem());
    const IntTab& face_voisins = zone_dis_base.face_voisins();
    const int nfaces = face_voisins.dimension(0);
    assert(rho_faces.size() == nfaces);
    for (int i = 0; i < nfaces; i++)
      {
        const int elem0 = face_voisins(i, 0);
        const int elem1 = face_voisins(i, 1);
        double rho = 0.;
        if (elem0 >= 0)
          rho = rho_elem[elem0];
        if (elem1 >= 0)
          rho += rho_elem[elem1];
        if (elem0 >= 0 && elem1 >= 0)
          rho *= 0.5;
        rho_faces[i] = rho;
      }
    rho_faces.echange_espace_virtuel();
  }
}

static void FT_disc_calculer_champs_rho_mu_nu_mono(const Zone_dis_base& zdis,
                                                   const Fluide_Incompressible& fluide,
                                                   Champ_Fonc& champ_rho_elem_,
                                                   Champ_Don& champ_mu_,
                                                   Champ_Don& champ_nu_,
                                                   Champ_Fonc& champ_rho_faces_)
{

  if (sub_type(Champ_Uniforme,champ_rho_elem_.valeur()) && (sub_type(Champ_Uniforme,champ_nu_.valeur())))
    {
      const DoubleTab& tab_rho_phase_0 = fluide.masse_volumique().valeurs();
      const double rho = tab_rho_phase_0(0,0);
      const DoubleTab& tab_nu_phase_0 = fluide.viscosite_cinematique().valeur().valeurs();
      const double nu = tab_nu_phase_0(0,0);
      const double mu = nu * rho;
      champ_rho_elem_. valeur().valeurs() = rho;
      champ_nu_.       valeur().valeurs() = nu;
      champ_mu_.       valeur().valeurs() = mu;
      champ_rho_faces_.valeur().valeurs() = rho;
    }
  else
    {
      const int nb_el = zdis.zone().nb_elem();
      const Zone_VF& zvf = ref_cast(Zone_VF,zdis);
      const IntTab& face_vois = zvf.face_voisins();
      const DoubleVect& vol = zvf.volumes();
      const int nb_faces = zvf.nb_faces();
      int elem1,elem2;
      double volume;
      //


      IntVect les_polys(nb_el);
      for (int i=0; i<nb_el; i++)
        les_polys(i) = i;

      const DoubleTab& cg = zvf.xp();
      DoubleTab& val_rho = champ_rho_elem_.valeur().valeurs();
      DoubleTab& val_nu = champ_nu_.valeur().valeurs();
      DoubleTab& val_mu = champ_mu_.valeur().valeurs();
      DoubleTab& val_rho_faces = champ_rho_faces_.valeur().valeurs();

      fluide.masse_volumique()->valeur_aux_elems(cg,les_polys,val_rho);
      fluide.viscosite_dynamique()->valeur_aux_elems(cg,les_polys,val_mu);
      val_rho.echange_espace_virtuel();
      val_mu.echange_espace_virtuel();

      for (int el=0; el<nb_el; el++)
        val_nu(el) = val_mu(el)/val_rho(el);
      val_nu.echange_espace_virtuel();

      for (int fac=0; fac<nb_faces; fac++)
        {
          elem1 = face_vois(fac,0);
          elem2 = face_vois(fac,1);
          val_rho_faces(fac) = 0.;
          volume = 0.;
          if (elem1!=-1)
            {
              val_rho_faces(fac) += val_rho(elem1)*vol(elem1);
              volume += vol(elem1);
            }
          if (elem2!=-1)
            {
              val_rho_faces(fac) += val_rho(elem2)*vol(elem2);
              volume += vol(elem2);
            }
          val_rho_faces(fac) /=volume;
        }

      val_rho_faces.echange_espace_virtuel();
    }
}
// Description: constructeur par defaut
Navier_Stokes_FT_Disc::Navier_Stokes_FT_Disc()
{
  variables_internes_ = new Navier_Stokes_FT_Disc_interne;
  // terme_diffusif.associer_diffusivite(champ_mu_);
  is_repulsion=0;
}

// Description: le destructeur qui va avec
Navier_Stokes_FT_Disc::~Navier_Stokes_FT_Disc()
{
  delete variables_internes_;
}

Sortie& Navier_Stokes_FT_Disc::printOn(Sortie& os) const
{
  return Navier_Stokes_std::printOn(os);
}

Entree& Navier_Stokes_FT_Disc::readOn(Entree& is)
{
  terme_diffusif.associer_diffusivite(champ_mu_);
  return Navier_Stokes_Turbulent::readOn(is);
}

void Navier_Stokes_FT_Disc::set_param(Param& param)
{
  Navier_Stokes_Turbulent::set_param(param);
  param.ajouter_flag("matrice_pression_invariante",&variables_internes().matrice_pression_invariante);
  param.ajouter_non_std("equation_interfaces_vitesse_imposee",(this) );
  param.ajouter_non_std("equations_interfaces_vitesse_imposee",(this) );
  param.ajouter_non_std("equation_interfaces_proprietes_fluide",(this));
  param.ajouter("clipping_courbure_interface",&variables_internes().clipping_courbure_interface);
  param.ajouter_non_std("equation_temperature_mpoint",(this));
  param.ajouter_non_std("equation_temperature_mpoint_vapeur",(this));
  param.ajouter_non_std("terme_gravite",(this));
  param.ajouter("equations_concentration_source_vortex",&variables_internes().equations_concentration_source_fluide_);
  param.ajouter_non_std("repulsion_aux_bords",(this));
  param.ajouter_non_std("penalisation_forcage",(this));
  param.ajouter_flag("mpoint_inactif_sur_qdm",&variables_internes().mpoint_inactif);
  param.ajouter_flag("mpoint_vapeur_inactif_sur_qdm",&variables_internes().mpointv_inactif);
  param.ajouter("correction_courbure_ordre", &variables_internes().correction_courbure_ordre_);

  param.ajouter("rayon", &Tool::myRayon,Param::REQUIRED);
  param.ajouter("dist_act", &Tool::mySigma);
  param.ajouter("modele_lubrification", &Tool::modele_lubrification);
  param.ajouter("d_desactivation_lubrification", &Tool::d_desactivation_lubrification);
  param.ajouter("modele_collisions", &variables_internes().modele_collisions_,Param::REQUIRED);
  param.ajouter_arr_size_predefinie("param_oscillateur", &Tool::param_osillateur);

  param.ajouter("decalage_bords", &Tool::decalage_bords);
  param.ajouter("force_sur_elem_diphasiques", &Tool::force_sur_elem_diphasiques);
  param.ajouter_arr_size_predefinie("valeurs_decalage", &Tool::valeurs_decalage);

  param.ajouter_arr_size_predefinie("Origine", &Tool::myOrigine,Param::REQUIRED);
  param.ajouter_arr_size_predefinie("Nombre_de_Noeuds", &Tool::myNb_Noeuds,Param::REQUIRED);
  param.ajouter_arr_size_predefinie("Longueurs", &Tool::myLongueurs,Param::REQUIRED);
}
int Navier_Stokes_FT_Disc::modele_collisions() const
{
  const Motcle& le_mot_cle =variables_internes().modele_collisions_;
  if (le_mot_cle == "ressort_amorti")
    return 0;
  else if (le_mot_cle == "mohagheg")
    return 1;
  else if (le_mot_cle == "hybrid" )
    return 2;
  else if (le_mot_cle == "breugem" )
    return 3;
  else
    return -1;
}
int Navier_Stokes_FT_Disc::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="diffusion")
    {
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      Cerr << "Reading and typing of the diffusion operator : " << finl;
      // Si on a lu le modele de turbulence et qu'il est nul,
      // alors on utilise l'operateur de diffusion standard.
      if (le_modele_turbulence.non_nul() // L'operateur a ete type (donc lu)
          && sub_type(Modele_turbulence_hyd_nul, le_modele_turbulence.valeur()))
        {
          is >> terme_diffusif; // Operateur de diffusion standard (non turbulent)
          if (Process::je_suis_maitre())
            {
              Cerr << " WARNING : standard diffusion operator used for front_tracking\n";
              Cerr << "           the transposed term grad(v) is missing !!!" << finl;
            }
        }
      else
        lire_op_diff_turbulent(is);

      // Le coefficient de diffusion est une viscosite dynamique.
      // Il faut le diviser par rho pour calculer le pas de temps de stabilite.
      terme_diffusif.valeur().
      associer_champ_masse_volumique(champ_rho_elem_.valeur());

      // Pareil avec le modele de turbulence
      if (le_modele_turbulence.non_nul())
        le_modele_turbulence.valeur().
        associer_champ_masse_volumique(champ_rho_elem_.valeur());
      return 1;
    }
  else if ((mot=="equation_interfaces_vitesse_imposee") || (mot=="equation_interfaces_proprietes_fluide"))
    {
      Motcle nom_equation;
      is >> nom_equation;
      Cerr << mot << " equation : " << nom_equation << finl;
      const Transport_Interfaces_FT_Disc& eq =
        probleme_ft().equation_interfaces(nom_equation);

      if (mot=="equation_interfaces_vitesse_imposee")
        {
          if (variables_internes().ref_eq_interf_vitesse_imposee.size()>0)
            {
              Cerr<<"Error: You have already a equation_interfaces_vitesse_imposee defined."<<finl;
              Cerr<<"Since 1.6.4 TRUST version, you need to use the following syntax when"<<finl;
              Cerr<<"using several equation_interfaces_vitesse_imposee equations:"<<finl;
              Cerr<<"equations_interfaces_vitesse_imposee number_of_equations equation_name_one equation_name_two ..."<<finl;
              exit();
            }
          else
            {
              Cerr << "===============================================================================" << finl;
              Cerr << "Warning: You are using a future obsolete syntax for defining a solid interface:" << finl;
              Cerr << "equation_interfaces_vitesse_imposee " << nom_equation << finl;
              Cerr << "Should be written:" << finl;
              Cerr << "equations_interfaces_vitesse_imposee 1 " << nom_equation << finl;
              Cerr << "===============================================================================" << finl;
            }
          variables_internes().ref_eq_interf_vitesse_imposee.add(eq);
        }
      else if (mot=="equation_interfaces_proprietes_fluide")
        variables_internes().ref_eq_interf_proprietes_fluide = eq;
      return 1;
    }
  else if (mot=="equations_interfaces_vitesse_imposee")
    {
      Noms na;
      is >> na;
      variables_internes().ref_eq_interf_vitesse_imposee.reset();
      for (int i=0; i<na.size(); i++)
        {
          const Transport_Interfaces_FT_Disc& eq =
            probleme_ft().equation_interfaces(na[i]);
          variables_internes().ref_eq_interf_vitesse_imposee.add(eq);
        }
      return 1;
    }
  else if (mot=="equation_temperature_mpoint")
    {
      Nom nom_equation;
      is >> nom_equation;
      Cerr << " equation : " << nom_equation << finl;
      const Equation_base& eq =
        probleme_ft().get_equation_by_name(nom_equation);
      if (!sub_type(Convection_Diffusion_Temperature_FT_Disc, eq))
        {
          Cerr << " Error : equation is not of type Convection_Diffusion_Temperature_FT_Disc"
               << finl;
          exit();
        }
      variables_internes().ref_equation_mpoint_ =
        ref_cast(Convection_Diffusion_Temperature_FT_Disc, eq);
      return 1;
    }
  else if (mot=="equation_temperature_mpoint_vapeur")
    {
      Nom nom_equation;
      is >> nom_equation;
      Cerr << " equation : " << nom_equation << finl;
      const Equation_base& eq =
        probleme_ft().get_equation_by_name(nom_equation);
      if (!sub_type(Convection_Diffusion_Temperature_FT_Disc, eq))
        {
          Cerr << " Error : equation is not of type Convection_Diffusion_Temperature_FT_Disc"
               << finl;
          Process::exit();
        }
      variables_internes().ref_equation_mpoint_vap_ =
        ref_cast(Convection_Diffusion_Temperature_FT_Disc, eq);
      return 1;
    }
  else if (mot=="terme_gravite")
    {
      Motcles mots;
      mots.add("rho_g"); // Terme source volumique traditionnel avec courants parasites
      mots.add("grad_i"); // Terme source "front-tracking" sans courants parasites mais pbs en VEF.
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading terme_gravite : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          variables_internes().terme_gravite_ = Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G;
          break;
        case 1:
          variables_internes().terme_gravite_ = Navier_Stokes_FT_Disc_interne::GRAVITE_GRAD_I;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis <<" has been found."<< finl;
          barrier();
          exit();
        }
      return 1;
    }
  else if (mot=="repulsion_aux_bords")
    {
      is_repulsion=1;
      is >> minx;
      is >> maxx;
      is >> pente;
      Cerr <<"Interfaces repulsion on the boundaries for : "<<minx<<" "<<maxx<<finl;
      return 1;
    }
  else if (mot=="penalisation_forcage")
    {
      variables_internes().is_penalized = 1;
      variables_internes().is_explicite = 0;
      variables_internes().eta = 1.e-12;
      Cerr << "Navier_Stokes_FT_Disc : option penalisation_forcage" << finl;
      Motcles mots;
      mots.add("pression_reference"); // Valeur reference pour penalisation L2 pression
      mots.add("domaine_flottant_fluide"); // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
      Motcle motbis;
      is >> motbis;
      Motcle accouverte = "{" , accfermee = "}" ;
      if (motbis == accouverte)
        {
          is >> motbis;
          while (motbis != accfermee)
            {
              int rang = mots.search(motbis);
              switch(rang)
                {
                case 0:
                  {
                    is >> variables_internes().p_ref_pena;
                    Cerr << "Lecture de la pression de reference : " << variables_internes().p_ref_pena << finl;
                    break;
                  }
                case 1:
                  {
                    variables_internes().is_pfl_flottant = 1;
                    is >> variables_internes().x_pfl_imp ;
                    is >> variables_internes().y_pfl_imp ;
                    if( Objet_U::dimension == 3 )
                      is >> variables_internes().z_pfl_imp ;
                    Cerr <<"Domaine flottant fluide"<< finl;
                    Cerr <<"Lecture de la position du point de reference pression fluide : "<<variables_internes().x_pfl_imp <<" "<<variables_internes().y_pfl_imp<<" "<<variables_internes().z_pfl_imp<< finl;
                    break;
                  }
                default:
                  Cerr << "Erreur, on attendait " << mots << "On a trouve : " << motbis << finl;
                  barrier();
                  exit();
                }
              is >> motbis;
            }
        }
      else
        {
          Cerr << "Erreur a la lecture des parametres de la penalisation du forcage " << finl;
          Cerr << "On attendait : " << accouverte << finl;
          exit();
        }
    }
  else
    return Navier_Stokes_Turbulent::lire_motcle_non_standard(mot,is);
  return 1;
}

const Champ_Don& Navier_Stokes_FT_Disc::diffusivite_pour_transport()
{
  return champ_mu_;
}

void Navier_Stokes_FT_Disc::associer_pb_base(const Probleme_base& un_probleme)
{
  if (! sub_type(Probleme_FT_Disc_gen, un_probleme))
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::associer_pb_base\n";
      Cerr << " Navier_Stokes_FT_Disc equation must be associated to\n";
      Cerr << " a Probleme_FT_Disc_gen problem type\n";
      exit();
    }
  probleme_ft_ = ref_cast(Probleme_FT_Disc_gen, un_probleme);
  Navier_Stokes_std::associer_pb_base(un_probleme);
}

const Fluide_Diphasique& Navier_Stokes_FT_Disc::fluide_diphasique() const
{
  const REF(Fluide_Diphasique) & fluide_dipha = variables_internes().ref_fluide_diphasique;
  if (! fluide_dipha.non_nul())
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::fluide_diphasique()\n";
      Cerr << " The fluid has not been associated\n";
      Cerr << " (check that the fluid is of Fluide_Diphasique type)" << finl;
      assert(0);
      exit();
    }
  return fluide_dipha.valeur();
}

void Navier_Stokes_FT_Disc::associer_milieu_base(const Milieu_base& un_fluide)
{
  if (sub_type(Fluide_Diphasique, un_fluide))
    variables_internes().ref_fluide_diphasique = ref_cast(Fluide_Diphasique, un_fluide);
  else
    Navier_Stokes_Turbulent::associer_milieu_base(un_fluide);
}

const Milieu_base& Navier_Stokes_FT_Disc::milieu() const
{
  if (variables_internes().ref_fluide_diphasique.non_nul())
    return variables_internes().ref_fluide_diphasique.valeur();
  else
    return Navier_Stokes_Turbulent::milieu();
}

Milieu_base& Navier_Stokes_FT_Disc::milieu()
{
  if (variables_internes().ref_fluide_diphasique.non_nul())
    return variables_internes().ref_fluide_diphasique.valeur();
  else
    return Navier_Stokes_Turbulent::milieu();
}

// Description:
//  Discretisation des champs utilises dans cette equation.
//  Fonction appelee par Probleme_base::discretiser.
//  B. Mathieu : a titre experimental, au lieu de dupliquer les noms
//               des champs ici et dans "a_pour_champ_Fonc", on stocke
//               les champs dans une liste. (voir a_pour_champ_fonc).
void Navier_Stokes_FT_Disc::discretiser()
{
  Navier_Stokes_Turbulent::discretiser();
  const Discretisation_base& dis = discretisation();
  const double temps = schema_temps().temps_courant();
  const Zone_dis_base& ma_zone_dis = zone_dis().valeur();
  LIST(REF(Champ_base)) & champs_compris = variables_internes().liste_champs_compris;

  dis.discretiser_champ("champ_elem", ma_zone_dis,
                        "diffusivite", "m^2/s",
                        1 /* composantes */, temps,
                        champ_nu_);
  champs_compris.add(champ_nu_.valeur());
  //Nouvelle formulation
  champs_compris_.ajoute_champ(champ_nu_);

  dis.discretiser_champ("champ_elem", ma_zone_dis,
                        "viscosite_dynamique", "kg/m.s",
                        1 /* composantes */, temps,
                        champ_mu_);
  champs_compris.add(champ_mu_.valeur());
  champs_compris_.ajoute_champ(champ_mu_);

  dis.discretiser_champ("champ_elem", ma_zone_dis,
                        "masse_volumique", "kg/m3",
                        1 /* composantes */, temps,
                        champ_rho_elem_);
  champs_compris.add(champ_rho_elem_.valeur());
  champs_compris_.ajoute_champ(champ_rho_elem_);

  // La masse volumique discretisee sur les volumes de controle de la vitesse
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "masse_volumique_vnodes", "kg/m3",
                        1 /* composantes */, temps,
                        champ_rho_faces_);
  champs_compris.add(champ_rho_faces_.valeur());
  champs_compris_.ajoute_champ(champ_rho_faces_);

  // Variables internes
  dis.discretiser_champ("pression", ma_zone_dis,
                        "second_membre_projection", "",
                        1 /* composantes */, temps,
                        variables_internes().second_membre_projection);
  champs_compris.add(variables_internes().second_membre_projection.valeur());
  champs_compris_.ajoute_champ(variables_internes().second_membre_projection);
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "gradient_pression", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().gradient_pression);
  champs_compris.add(variables_internes().gradient_pression.valeur());
  champs_compris_.ajoute_champ(variables_internes().gradient_pression);
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "derivee_u_etoile", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().derivee_u_etoile);
  champs_compris.add(variables_internes().derivee_u_etoile.valeur());
  champs_compris_.ajoute_champ(variables_internes().derivee_u_etoile);
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "terme_diffusion_vitesse", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_diffusion);
  champs_compris.add(variables_internes().terme_diffusion.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_diffusion);
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "terme_convection_vitesse", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_convection);
  champs_compris.add(variables_internes().terme_convection.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_convection);
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "terme_source_vitesse", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_source);
  champs_compris.add(variables_internes().terme_source.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_source);
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "terme_source_interfaces", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_source_interfaces);
  champs_compris.add(variables_internes().terme_source_interfaces.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_source_interfaces);
  if (dis.que_suis_je() == "VEFPreP1B")
    {
      dis.discretiser_champ("pression", ma_zone_dis,
                            "indicatrice_p1b", "",
                            1 /* composantes */, temps,
                            variables_internes().indicatrice_p1b);
      champs_compris.add(variables_internes().indicatrice_p1b.valeur());
      champs_compris_.ajoute_champ(variables_internes().indicatrice_p1b);
    }

  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "gradient_indicatrice", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().gradient_indicatrice);
  champs_compris.add(variables_internes().gradient_indicatrice.valeur());
  champs_compris_.ajoute_champ(variables_internes().gradient_indicatrice);
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "potentiel_faces", "",
                        1 /* composante */, temps,
                        variables_internes().potentiel_faces);
  champs_compris.add(variables_internes().potentiel_faces.valeur());
  champs_compris_.ajoute_champ(variables_internes().potentiel_faces);
  dis.discretiser_champ("champ_elem", ma_zone_dis,
                        "potentiel_elements", "",
                        1 /* composante */, temps,
                        variables_internes().potentiel_elements);
  champs_compris.add(variables_internes().potentiel_elements.valeur());
  champs_compris_.ajoute_champ(variables_internes().potentiel_elements);
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "vitesse_delta_interface", "m/s",
                        -1 /* nb composantes par defaut */, 1, temps,
                        variables_internes().delta_u_interface);
  // Pour pouvoir faire filtrer_L2:
  variables_internes().delta_u_interface.associer_eqn(*this);
  champs_compris.add(variables_internes().delta_u_interface.valeur());
  champs_compris_.ajoute_champ(variables_internes().delta_u_interface);
  dis.discretiser_champ("pression", ma_zone_dis,
                        "pression_laplacien_d", "",
                        1 /* composante */, temps,
                        variables_internes().laplacien_d);
  champs_compris.add(variables_internes().laplacien_d.valeur());
  champs_compris_.ajoute_champ(variables_internes().laplacien_d);
  dis.discretiser_champ("temperature", ma_zone_dis,
                        "temperature_mpoint", "",
                        1 /* composante */, temps,
                        variables_internes().mpoint);
  champs_compris.add(variables_internes().mpoint.valeur());
  champs_compris_.ajoute_champ(variables_internes().mpoint);
  dis.discretiser_champ("temperature", ma_zone_dis,
                        "temperature_mpointv", "",
                        1 /* composante */, temps,
                        variables_internes().mpoint_vap);
  champs_compris.add(variables_internes().mpoint_vap.valeur());
  champs_compris_.ajoute_champ(variables_internes().mpoint_vap);
  // Pour variation temporelle dI_dt
  dis.discretiser_champ("pression", ma_zone_dis,
                        "derivee_temporelle_indicatrice", "",
                        1 /* composante */, temps,
                        variables_internes().derivee_temporelle_indicatrice);
  champs_compris.add(variables_internes().derivee_temporelle_indicatrice.valeur());
  champs_compris_.ajoute_champ(variables_internes().derivee_temporelle_indicatrice);
  //HMS
  dis.discretiser_champ("vitesse", ma_zone_dis,
                        "terme_source_collisions", "",
                        -1 /* nb composantes par defaut */, temps,
                        variables_internes().terme_source_collisions);
  champs_compris.add(variables_internes().terme_source_collisions.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_source_collisions);

  dis.discretiser_champ("pression", ma_zone_dis,
                        "num_compo", "",
                        1 /* composante */ , temps,
                        variables_internes().num_compo);
  champs_compris.add(variables_internes().num_compo.valeur());
  champs_compris_.ajoute_champ(variables_internes().num_compo);

  //dis.discretiser_champ("pression", ma_zone_dis,
  //                      "num_compo", "",
  //                      1 /* composante */,1, temps,
  //                      Tool::num_compo);
  //champs_compris.add(Tool::num_compo.valeur());
  //champs_compris_.ajoute_champ(Tool::num_compo);
}

// Description:
//  Methode surchargee de Navier_Stokes_std, appelee par
//   Navier_Stokes_std::discretiser().
//  L'assembleur pression est particulier pour le front-tracking
//  en VEF (en attendant qu'on factorise tous ces assembleurs pression)
void Navier_Stokes_FT_Disc::discretiser_assembleur_pression()
{
  const Discretisation_base& dis = discretisation();
  if (dis.que_suis_je() == "VDF")
    {
      // Assembleur pression generique (prevu pour rho_variable)
      assembleur_pression_.typer("Assembleur_P_VDF");
    }
  else if (dis.que_suis_je() == "VEFPreP1B")
    {
      // Assembleur pression generique (prevu pour rho_variable)
      assembleur_pression_.typer("Assembleur_P_VEFPreP1B");
    }

  Assembleur_base& assembleur = assembleur_pression_.valeur();
  assembleur.associer_zone_dis_base(zone_dis().valeur());
  // B. Mathieu, 07_2004 : premier jet de la methode, on resout en pression :
  assembleur.set_resoudre_increment_pression(0);
  assembleur.set_resoudre_en_u(1);
}

// Description: methode appelee par Navier_Stokes_std::preparer_calcul...
void Navier_Stokes_FT_Disc::projeter()
{
  if (Process::je_suis_maitre() && limpr())
    Cerr << "Navier_Stokes_FT_Disc::projeter does nothing" << finl;
}

// Description: methode appelee par Probleme_base::preparer_calcul()
int Navier_Stokes_FT_Disc::preparer_calcul()
{
  Cerr << "Navier_Stokes_FT_Disc::preparer_calcul()" << finl;
  Equation_base::preparer_calcul();

  le_modele_turbulence.preparer_calcul();


  {
    // Si le mot cle equation_interfaces_proprietes_fluide a ete specifie,
    // l'indicatrice de phase correspondant a l'equation sert a calculer
    // les proprietes physiques du milieu.
    // Sinon, le milieu doit etre de type Fluide_Incompressible.

    REF(Transport_Interfaces_FT_Disc) & ref_equation =
      variables_internes().ref_eq_interf_proprietes_fluide;

    if (ref_equation.non_nul())
      {

        // Couplage Navier-Stokes / Fluide : les interfaces determinent la position des phases.
        // On utilise les proprietes du fluide diphasique

        if (Process::je_suis_maitre())
          {
            Cerr << "Initialisation of the physical properties (rho, mu, ...)\n"
                 << " based on the indicatrice field of the equation " << ref_equation.valeur().le_nom() << finl;

          }
        FT_disc_calculer_champs_rho_mu_nu_dipha(zone_dis().valeur(),
                                                fluide_diphasique(),
                                                ref_equation.valeur().
                                                get_update_indicatrice().valeurs(), // indicatrice
                                                champ_rho_elem_.valeur().valeurs(),
                                                champ_nu_.valeur().valeurs(),
                                                champ_mu_.valeur().valeurs(),
                                                champ_rho_faces_.valeur().valeurs());
      }
    else
      {

        // Pas de couplage Navier-Stokes / Fluide : le fluide est monophasique.

        const Fluide_Incompressible& phase_0 = ref_cast(Fluide_Incompressible,milieu());
        const Zone_dis_base& zdis = zone_dis().valeur();
        FT_disc_calculer_champs_rho_mu_nu_mono(zdis,
                                               phase_0,
                                               champ_rho_elem_,
                                               champ_mu_,
                                               champ_nu_,
                                               champ_rho_faces_);
      }

    // Premiere version :les termes sources sont homogenes a rho*v
    //   sources().associer_champ_rho(champ_rho_elem_.valeur());
    // Nouvelle version: les termes "sources()" sont homogenes a v.
    //   on ne fait rien.
  }

  DoubleTab& tab_vitesse = inconnue().valeurs();

  // On assemble la matrice de pression pour la premiere fois.
  assembleur_pression().valeur().assembler_rho_variable(matrice_pression_,
                                                        champ_rho_faces_.valeur());
  // Informe le solveur que la matrice a change :
  solveur_pression().valeur().reinit();
  // Calcul du second membre :
  //  div(vpoint) a l'interieur du domaine,
  //  prise en compte des conditions aux limites de pression/vitesse
  DoubleTab& secmem = variables_internes().second_membre_projection.valeurs();
  divergence.calculer(tab_vitesse, secmem);
  secmem *= -1;
  // Il faut faire ceci car on ne resout pas en "increment de pression":
  assembleur_pression_.valeur().modifier_secmem(secmem);



  // Ajout pour la sauvegarde au premier pas de temps si reprise
  la_pression.changer_temps(schema_temps().temps_courant());
  // HMS Tool:: moi aussi je vais faire pareil tien !
  //variables_internes().num_compo.changer_temps(schema_temps().temps_courant());



  // Resolution du systeme en pression : calcul de la_pression
  solveur_pression_.resoudre_systeme(matrice_pression_.valeur(),
                                     secmem,
                                     la_pression.valeurs()
                                    );
  assembleur_pression_.modifier_solution(la_pression);
  // Calcul d(u)/dt = vpoint + 1/rho*grad(P)
  DoubleTab& gradP = variables_internes().gradient_pression.valeurs();
  gradient.calculer(la_pression.valeur().valeurs(), gradP);
  solveur_masse.appliquer(gradP);
  // Correction de la vitesse :
  if (projection_a_faire()) // Temporaire pour permettre de ne pas resoudre NS avec mettant operateurs nuls et projection_initiale 0
    {
      int i, j;
      const DoubleTab& rho_faces = champ_rho_faces_.valeur().valeurs();
      const int n = tab_vitesse.dimension(0);
      if (tab_vitesse.nb_dim() == 1)
        {
          for (i = 0; i < n; i++)
            tab_vitesse(i) -= gradP(i) / rho_faces(i);
        }
      else
        {
          const int m = tab_vitesse.dimension(1);
          for (i = 0; i < n; i++)
            {
              double x = rho_faces(i);
              for (j = 0; j < m; j++)
                {
                  tab_vitesse(i,j) -= gradP(i,j) / x;
                }
            }
        }
      tab_vitesse.echange_espace_virtuel();
    }

  if (le_traitement_particulier.non_nul())
    {
      if (le_traitement_particulier.valeur().support_ok())
        le_traitement_particulier.valeur().
        associer_champ_masse_volumique(champ_rho_faces_.valeur());
      le_traitement_particulier.preparer_calcul_particulier();
    }

  return 1;
}

const int& Navier_Stokes_FT_Disc::get_is_penalized() const
{
  return variables_internes().is_penalized;
}

void Navier_Stokes_FT_Disc::preparer_pas_de_temps()
{
}

void Navier_Stokes_FT_Disc::mettre_a_jour(double temps)
{
  Navier_Stokes_Turbulent::mettre_a_jour(temps);

  champ_rho_elem_.mettre_a_jour(temps);
  champ_rho_faces_.mettre_a_jour(temps);
  champ_mu_.mettre_a_jour(temps);
  champ_nu_.mettre_a_jour(temps);
}



// Description:
//  Calcul des forces de tension superficielles (F_sigma) et de la partie
//  a rotationnel non nul de la gravite (G_rot) (si GRAVITE_GRAD_I) :
//  F_sigma = INTEGRALE sur le volume de controle (
//            sigma_aux_faces * courbure_aux_faces * gradient(indicatrice)
//            + gradient_sigma )
//  G_rot   = INTEGRALE sur le volume de controle (
//               phi * gradient(rho) )   (avec phi = potentiel de pesanteur)
// Parametre: gradient_indicatrice
// Signification: le gradient de l'indicatrice issu de l'operateur "gradient",
//  donc homogene a l'integrale du gradient sur les volumes de controle de la vitesse.
// Parametre: potentiel_faces
// Signification: un champ aux faces a une composante, ou on stocke le "potentiel aux faces"
// Parametre: champ
// Signification: le champ aux faces (meme discretisation que la vitesse)
//  ou on stocke le terme source des forces superficielles.
void Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles(const Maillage_FT_Disc& maillage,
                                                                 const Champ_base& gradient_indicatrice,
                                                                 Champ_base& potentiel_elements,
                                                                 Champ_base& potentiel_faces,
                                                                 Champ_base& champ) const
{
  const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
  // Nombre de faces
  const int nb_faces = zone_vf.nb_faces();
  {
    // Le champ et le gradient de l'indicatrice doivent etre aux faces
    assert(champ.valeurs().dimension(0) == nb_faces);
    assert(potentiel_faces.valeurs().dimension(0) == nb_faces);
    assert(gradient_indicatrice.valeurs().dimension(0) == nb_faces);
  }
  const int dim = Objet_U::dimension;

  // Calcul du "potentiel aux sommets du maillage"
  ArrOfDouble potentiel_sommets;
  potentiel_sommets.resize_array(maillage.nb_sommets(), Array_base::NOCOPY_NOINIT);

  // (ce potentiel est constant sur chaque portion connexe d'interface si
  //  l'interface est a l'equilibre).
  //  courbure * sigma + potentiel_gravite_phase_1 - potentiel_gravite_phase_0
  {
    // Pour l'instant, la tension de surface est constante
    const Fluide_Diphasique& fluide_dipha = fluide_diphasique();
    const double sigma = fluide_dipha.sigma();

    const int n = maillage.nb_sommets();
    const ArrOfDouble& courbure_sommets = maillage.get_update_courbure_sommets();
    //Calcul du "potentiel aux sommets"
    int i;
    {
      const double clipping_courbure_max = variables_internes().clipping_courbure_interface;
      int clip_counter = 0;
      for (i = 0; i < n; i++)
        {
          double c = courbure_sommets[i];
          // Clipping de la courbure: si la courbure est superieure a la
          // valeur maxi autorisee, on limite (permet de ne pas plomber le
          // pas de temps s'il y a une singularite geometrique dans le maillage)
          if (fabs(c) > clipping_courbure_max)
            {
              clip_counter++;
              c = ((c > 0) ? 1. : -1.) * clipping_courbure_max;
            }
          potentiel_sommets[i] = c * sigma;
        }
      clip_counter = mp_sum(clip_counter);
      if (clip_counter > 0 && je_suis_maitre())
        {
          Cerr << "Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles : clip_count "
               << clip_counter << finl;
        }
    }

    //Ajout des termes de gravite
    if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_GRAD_I)
      {
        if (milieu().a_gravite())
          {
            const double rho_0 = fluide_dipha.fluide_phase(0).masse_volumique().valeurs()(0,0);
            const double rho_1 = fluide_dipha.fluide_phase(1).masse_volumique().valeurs()(0,0);
            const double delta_rho = rho_1 - rho_0;

            // Pour l'instant : gravite uniforme g => phi(s) = - x scalaire g
            const DoubleTab& gravite = milieu().gravite().valeurs();
            if (gravite.nb_dim() != 2 || gravite.dimension(1) != dim)
              {
                Cerr << "Error for the method calculer_champ_forces_superficielles\n";
                Cerr << "gravite.dimension(1) != Objet_U::dimension" << finl;
                Process::exit();
              }
            const DoubleTab& sommets = maillage.sommets();

            for (i = 0; i < n; i++)
              {
                double p = 0.;
                for (int j = 0; j < dim; j++)
                  p += sommets(i,j) * gravite(0,j);

                if(is_repulsion)
                  {
                    double dx=0.;
                    double px=sommets(i,0);
                    if(px>maxx)
                      dx=px-maxx;
                    else if (px<minx)
                      dx=minx-px;

                    double dy=0.;
                    double py=sommets(i,1);
                    if(py>maxx)
                      dy=py-maxx;
                    else if (py<minx)
                      dy=minx-py;

                    p+=sqrt(dx*dx+dy*dy)*pente;
                  }

                potentiel_sommets[i] -= p * delta_rho;
              }
          }
      }
  }
#if 0
  // Calcul du "potentiel aux faces" : interpolation sur les faces des elements
  // traverses par l'interface.
  {
    // Tableau des facettes du maillage interfaces
    const IntTab& facettes = maillage.facettes();
    // Tableau des faces des elements :
    const IntTab& elem_faces = zone_vf.elem_faces();
    const int nb_faces_par_element = elem_faces.dimension(1);
    // Le tableau "potentiel aux faces" a remplir :
    DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    valeurs_potentiel_faces = 0.;
    // Somme des poids des contributions ajoutees aux faces
    ArrOfDouble poids(nb_faces);
    poids = 0.;

    // Boucle sur les faces de l'interface
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_facette_elem = intersections.index_facette();
    double pot[3] = {0., 0., 0.};

    int fa7;
    const int nb_facettes = maillage.nb_facettes();
    for (fa7 = 0; fa7 < nb_facettes; fa7++)
      {
        // Potentiel aux sommets de la facette:
        pot[0] = potentiel_sommets[facettes(fa7,0)];
        pot[1] = potentiel_sommets[facettes(fa7,1)];
        if (dim == 3)
          pot[2] = potentiel_sommets[facettes(fa7,2)];
        // Boucle sur les elements traverses par la facette
        int index = index_facette_elem[fa7];
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
            // Calcul du potentiel au centre de l'intersection
            double p = 0.;
            int i;
            for (i = 0; i < 3; i++)
              p += data.barycentre_[i] * pot[i];
            // La contribution aux faces est ponderee par la surface d'intersection
            p *= data.surface_intersection_;

            // Ajout de la contribution aux faces de l'element
            const int element = data.numero_element_;
            for (i = 0; i < nb_faces_par_element; i++)
              {
                const int face = elem_faces(element, i);
                poids[face] += data.surface_intersection_;
                valeurs_potentiel_faces(face) += p;
              }
            index = data.index_element_suivant_;
          }
      }

    // Il reste a diviser les valeurs aux face par la somme des poids
    for (int face = 0; face < nb_faces; face++)
      {
        double p = poids[face];
        if (p > 0.)
          valeurs_potentiel_faces(face) /= p;
      }
    valeurs_potentiel_faces.echange_espace_virtuel();
  }
#else
  // Calcul du "potentiel aux elements" :
  DoubleTab poids(potentiel_elements.valeurs());
  {
    // Tableau des facettes du maillage interfaces
    const IntTab& facettes = maillage.facettes();
    // Le tableau "potentiel aux faces" a remplir :
    DoubleTab& valeurs_potentiel_elements = potentiel_elements.valeurs();
    const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
    valeurs_potentiel_elements = 0.;
    // Somme des poids des contributions ajoutees aux faces
    poids = 0.;

    // Boucle sur les faces de l'interface
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_facette_elem = intersections.index_facette();
    double pot[3] = {0., 0., 0.};

    int fa7;
    const int nb_facettes = maillage.nb_facettes();
    for (fa7 = 0; fa7 < nb_facettes; fa7++)
      {
        // Potentiel aux sommets de la facette:
        pot[0] = potentiel_sommets[facettes(fa7,0)];
        pot[1] = potentiel_sommets[facettes(fa7,1)];
        if (dim == 3)
          pot[2] = potentiel_sommets[facettes(fa7,2)];
        // Boucle sur les elements traverses par la facette
        int index = index_facette_elem[fa7];
        const double surface_facette = surface_facettes[fa7];
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
            // Calcul du potentiel au centre de l'intersection
            double p = 0.;
            int i;
            for (i = 0; i < 3; i++)
              p += data.barycentre_[i] * pot[i];
            // La contribution aux faces est ponderee par la surface d'intersection
#ifdef AVEC_BUG_SURFACES
            const double surf = data.surface_intersection_;
#else
            const double surf = data.fraction_surface_intersection_ * surface_facette;
#endif
            p *= surf;

            // Ajout de la contribution a l'element
            const int element = data.numero_element_;
            valeurs_potentiel_elements(element) += p;
            poids(element) += surf;

            index = data.index_element_suivant_;
          }
      }
    valeurs_potentiel_elements.echange_espace_virtuel();
    poids.echange_espace_virtuel();
    if (champ.valeurs().nb_dim() > 1)   // VEF ?
      {
        // Compute a potential on elements that are neighbours of
        // elements containing an interface :
        //  For VEF, the gradient of the indicator function can be
        //  non zero on the faces of these elements.
        int element;
        const int nb_elements = poids.dimension(0);
        const IntTab& face_voisins = la_zone_dis.valeur().valeur().face_voisins();
        const Zone_VF& zonevf = ref_cast(Zone_VF, la_zone_dis.valeur().valeur());
        const IntTab& elem_faces = zonevf.elem_faces();
        const int nb_faces_par_element = elem_faces.dimension(1);
        DoubleVect copie_valeurs_potentiel_elements(valeurs_potentiel_elements);
        DoubleVect copie_poids(poids);
        for (element = 0; element < nb_elements; element++)
          {
            double potential = 0.; // sum of potentials of neighbouring elements
            double p = 0.;   // sum of weights of neighbouring elements
            int i_face;
            for (i_face = 0; i_face < nb_faces_par_element; i_face++)
              {
                const int face = elem_faces(element, i_face);
                const int elem_voisin_0 = face_voisins(face, 0);
                const int elem_voisin_1 = face_voisins(face, 1);
                const int elem_voisin = elem_voisin_0 + elem_voisin_1 - element;
                if (elem_voisin >= 0)   // Not a boundary of the domain ?
                  {
                    potential += copie_valeurs_potentiel_elements(elem_voisin);
                    p += copie_poids(elem_voisin);
                  }
              }
            const double old_weight = copie_poids(element);
            // Do not change values for elements that contain an interface:
            if (p > 0. && old_weight == 0.)
              {
                // Decrease weight so that values on faces adjacent to elements
                // containing an interface are almost untouched.
                static double reduction_factor = 0.1;
                potential = potential * reduction_factor;
                p = p * reduction_factor;
                // Assign value
                valeurs_potentiel_elements(element) = potential;
                poids(element) = p;
              }
          }
        valeurs_potentiel_elements.echange_espace_virtuel();
        poids.echange_espace_virtuel();
        Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles poids:",poids);
      }
  }
  {
    // Boucle sur les faces
    int face;
    DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    valeurs_potentiel_faces = 0.;
    const DoubleTab& valeurs_potentiel_elements = potentiel_elements.valeurs();
    const int nb_faces_pot = valeurs_potentiel_faces.dimension(0);
    const IntTab& face_voisins = la_zone_dis.valeur().valeur().face_voisins();
    for (face = 0; face < nb_faces_pot; face++)
      {
        double p = 0.; // Somme des poids des deux elements voisins
        double pot = 0.; // Somme des potentiels
        // Boucle sur les deux elements voisins de la face
        int i;
        for (i = 0; i < 2; i++)
          {
            int element = face_voisins(face, i);
            if (element >= 0)
              {
                // If neighbour exists (not a boundary face)
                p += poids(element);
                pot += valeurs_potentiel_elements(element);
              }
          }
        if (p > 0.)
          {
            valeurs_potentiel_faces(face) = pot / p;
          }
      }
    valeurs_potentiel_faces.echange_espace_virtuel();
    Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles valeurs_potentiel_faces:",valeurs_potentiel_faces);
  }
#endif
  // Derniere operation : calcul du champ
  //   champ = potentiel_faces * gradient_indicatrice
  {
    const DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    const DoubleTab& valeurs_gradient_i = gradient_indicatrice.valeurs();
    DoubleTab& valeurs_champ = champ.valeurs();
    if (valeurs_champ.nb_dim() == 1)   // VDF (1 composante par face)
      {
        for (int face = 0; face < nb_faces; face++)
          {
            const double p = valeurs_potentiel_faces(face);
            valeurs_champ(face) = valeurs_gradient_i(face) * p;

          }
      }
    else     // VEF
      {
        const int nb_compo = valeurs_champ.dimension(1);
        for (int face = 0; face < nb_faces; face++)
          {
            const double p = valeurs_potentiel_faces(face);
            for (int i = 0; i < nb_compo; i++)
              {
                valeurs_champ(face, i) = valeurs_gradient_i(face, i) * p;
              }
          }
      }
    valeurs_champ.echange_espace_virtuel();
    Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles valeurs_champ:",valeurs_champ);
  }
}

//<editor-fold desc="Ancien modele"> Todo : a supprimer

void Navier_Stokes_FT_Disc::calculer_champ_forces_collisions(const DoubleTab& indicatrice, DoubleTab& valeurs_champ,
                                                             double& FC)
{
  Tool::compteur_+=1;
  Cerr << "Navier_Stokes_FT_Disc::calculer_champ_forces_collisions" << finl;
  //auto start = high_resolution_clock::now();

  REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
  const Zone_VF& zone_vdf = ref_cast(Zone_VDF, zone_dis().valeur());
  // recuperation des vitesse compo et des centres de gravites
  DoubleTab& positions = eq_transport.getPositionsCompo();
  DoubleTab& vitesses = eq_transport.getVitessesCompo();


  const int nb_facettes = maillage.nb_facettes();
  ArrOfInt compo_connexes_facettes(nb_facettes); // Init a zero
  int n = search_connex_components_local_FT(maillage, compo_connexes_facettes); //
  int nb_compo_tot = compute_global_connex_components_FT(maillage, compo_connexes_facettes, n); //

  assert(nb_compo_tot !=0);

  // si vitesses na pas les bonnes dimension (premiere pas de temps, en calcule la position des centre de gravite et
  // on initialise les vitesse des inclusions a zero
  if ((vitesses.dimension(0) != nb_compo_tot) || (vitesses.dimension(1) != dimension))
    if (1)
      {
        positions.resize(nb_compo_tot, dimension);
        vitesses.resize(nb_compo_tot, dimension);
        vitesses = 0;

        // calcule des centre de gravite pour les composantes
        assert(nb_compo_tot == positions.dimension(0));

        const int dim = positions.dimension(1); //
        const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
        const IntTab& facettes = maillage.facettes();
        const DoubleTab& sommets = maillage.sommets();
        assert(facettes.dimension(1) == dim);

        // Surface totale de chaque composante connexe, initialise a zero
        ArrOfDouble surfaces_compo(nb_compo_tot);
        positions = 0.;

        // Calcul du centre de gravite de la composante connexe
        //  (centre de gravite de la surface, pas du volume)
        const int nb_facettes_tot = facettes.dimension_tot(0);
        {
          for (int i = 0; i < nb_facettes_tot; i++)
            {
              if (maillage.facette_virtuelle(i)) continue;

              const int compo = compo_connexes_facettes[i];
              const double surface = surface_facettes[i];
              surfaces_compo[compo] += surface;
              // Centre de gravite de la facette, pondere par la surface
              for (int j = 0; j < dim; j++)
                {
                  // Indice du sommet
                  const int s = facettes(i, j);
                  for (int k = 0; k < dim; k++)
                    {
                      // On divisera par dim a la fin:
                      positions(compo, k) += surface * sommets(s, k);
                    }
                }
            }
          mp_sum_for_each_item(surfaces_compo);
          mp_sum_for_each_item(positions);

          positions *= (1. / dim);
          DoubleVect s; // tab_divide prend DoubleVect, pas ArrOfDouble...
          s.ref_array(surfaces_compo);
          tab_divide_any_shape(positions, s);
          //---fin calcul surface_compo et position_centre_gravite_compo ( a rassembler dans une focntion)---//
        }

        Tool::isFirstCollision.resize(nb_compo_tot);
        Tool::isFirstCollision = 1;

        Tool::memorisedElongation.resize(nb_compo_tot);
        Tool::memorisedElongation = -100;

      }


  ArrOfInt elem_cg ;
  zone_vf.zone().chercher_elements(positions, elem_cg);
  Tool::F_now.resize(nb_compo_tot,nb_compo_tot);
  Tool::F_now=0;
  if (Tool::compteur_ == 1)
    {
      Tool::F_old.resize(nb_compo_tot,nb_compo_tot);

      Tool::raideur.resize(nb_compo_tot,nb_compo_tot);
      Tool::e_eff.resize(nb_compo_tot,nb_compo_tot);

      Tool::F_old=0;

      Tool::raideur=-1;
      Tool::e_eff=-1;
    }
  else
    {
      // on rearange les numero de particules pour avoir la meme numerotation qu'au pas de temps precedant
      DoubleVect old_num_compo = variables_internes().num_compo.valeur().valeurs();
      // DoubleVect old_num_compo = Tool::num_compo.valeur().valeurs();
      DoubleTab copy_vitesses = vitesses;
      DoubleTab copy_positions = positions;
      for (int bad_compo = 0; bad_compo < nb_compo_tot; bad_compo++)
        {
          if (elem_cg[bad_compo] != -1)
            {
              int good_compo = old_num_compo[elem_cg[bad_compo]];
              for (int d = 0; d < dimension; d++)
                {
                  vitesses(bad_compo, d) = copy_vitesses(good_compo, d);
                  positions(bad_compo, d) = copy_positions(good_compo, d);
                }
            }

        }

    }
  //auto stop = high_resolution_clock::now();
  //auto duration = duration_cast<milliseconds>(stop - start);
  //Cerr << "Calcul position vitesse des centres : " << duration.count() << " ms" << finl;
  //start = high_resolution_clock::now();


// ----------



  //const DoubleTab& cgf = zone_vf.xv();
  const IntVect& orientation = zone_vdf.orientation();






  //  ----- proprietes particules  -------------
  ArrOfDouble rayons_compo(nb_compo_tot);
  ArrOfDouble volumes_compo(nb_compo_tot);
  ArrOfDouble masse_compo(nb_compo_tot);

//--------------------------------------------------
  const Fluide_Diphasique& mon_fluide =fluide_diphasique();
  int nb_bord = 2 * dimension;
  ArrOfDouble positions_bords(nb_bord);
  const double dt = schema_temps().pas_de_temps();
  const double temps = schema_temps().temps_courant();

  const double rho_solide =  mon_fluide.fluide_phase(0).masse_volumique().valeurs()(0, 0);
  const double rho_fluide =  mon_fluide.fluide_phase(1).masse_volumique().valeurs()(0, 0);
  const double mu_fluide = rho_fluide * mon_fluide.fluide_phase(1).viscosite_cinematique().valeur().valeurs()(0, 0);


  //double mu_fluide = 0.09;
  int indic_phase_fluide = 1; //, indic_phase_solide=0; // TODO a generaliser


  //printf("mu fluide = %f \n" , mu_fluide);
  //calcule du rayon et du volume pour chaque composante

  //double dx = rayon/16;
//double epsi = dx / 4;

//positions_bords[0] = 0 + epsi;
//positions_bords[1] = 0 + epsi;
//positions_bords[2] = 0 + epsi;
//positions_bords[3] = 38.1e-3 - epsi;
//positions_bords[4] = 76.2e-3 - epsi;
//positions_bords[5] = 38.1e-3 - epsi;

  double dx= Tool::calc_positions_bords(positions_bords);
// -------------------------------------------------


// parametre du modele -----------------------------
  double rayon = Tool::myRayon, rugosite=rayon*1e-4;
  //printf("(r) rayon=%f\n",rayon);
  double ed = 0.97;
  int Nc = 8;
  double d_act = 2*dx/rayon;
  double d_sat = rugosite/rayon;
  double d_des = Tool::d_desactivation_lubrification;
  double sig_max=0.015*rayon;
  // FLAGS
  double FS = 1, FD = 0, FL = 0, Fb = 1, Fp = 1, TR = 0, print = 2, CV =1, myPI=3.14159265358979;
  int isPurementSolide=0;
//----------------------------------------------------
//pour debeug jdd rebond
  if (0)
    {
      d_des = -1e-2;
      // mobile : particule p
      positions(0,0)=5e-4;
      positions(0,1)=6.09e-3;
      positions(0,2)=5e-4;

      vitesses(0,0)=-0.1;
      vitesses(0,1)=-0.1;
      vitesses(0,2)=0.2;

      // fixe bord ou particule q
      positions(1,0)=0;
      positions(1,1)=9e-3;
      positions(1,2)=0;

      vitesses(1,0)=0.2;
      vitesses(1,1)=-0.1;
      vitesses(1,2)=-0.2;
    }

  for (int compo = 0; compo < nb_compo_tot; compo++)
    {
      rayons_compo[compo] = rayon;
      volumes_compo[compo] = 4 * myPI * pow(rayons_compo[compo], 3) / 3;
      masse_compo[compo] = volumes_compo[compo] * rho_solide;
    }

  DoubleTab impression_d_int(nb_compo_tot);
  impression_d_int = 6.6666666;
  DoubleTab impression_dist_int_tr(nb_compo_tot);
  impression_dist_int_tr = 7.777777;
  DoubleTab impression_next_dist_int_tr(nb_compo_tot);
  impression_next_dist_int_tr = 8.888888;
  ArrOfInt FP(nb_compo_tot), FB(nb_compo_tot);
  DoubleTab force_collision_bord(nb_compo_tot, dimension);
  DoubleTab force_lubrification_bord(nb_compo_tot, dimension);
  DoubleTab forces_bord(nb_compo_tot, dimension);
  DoubleTab force_lubrification_particule(nb_compo_tot, dimension);
  DoubleTab force_collision_particule(nb_compo_tot, dimension);
  DoubleTab forces_particule(nb_compo_tot, dimension);
  DoubleTab forces_solide(nb_compo_tot, dimension);

  //{
  //  stop = high_resolution_clock::now();
  //  duration = duration_cast<milliseconds>(stop - start);
  //  Cerr << "parametrage : " << duration.count() << " ms"<< finl;
  //  start = high_resolution_clock::now();
  //}


  for (int compo = 0; compo < nb_compo_tot; compo++)
    {

      double raideur_b = (masse_compo[compo] * (myPI * myPI + pow(log(ed), 2))) / pow(Nc * dt, 2);
      double amortis_b = (2 * masse_compo[compo] * log(ed)) / (Nc * dt);

      for (int d = 0; d < dimension; d++)
        {
          //--------------- contribution des bords --------------------------
          for (int bord = 0; bord < 2 * dimension; bord++)
            {

              // force de collision
              // les bord sont perpondiculaire a la direction considere ?
              int ori = bord < dimension ? bord : bord - dimension;
              if (ori == d)
                {
                  double dist_cg = positions(compo, d) - positions_bords[bord];
                  double dist_int = fabs(dist_cg) - rayons_compo[compo];
                  double d_int = dist_int / rayons_compo[compo];

                  double normBord = dist_cg / fabs(dist_cg);
                  double vitRelNorm = vitesses(compo, d) - 0.0;

                  if (d_int <= d_act)
                    {

                      //if (Process::je_suis_maitre()) printf("(c%d-b%d) d:%d, n:%f d_int:%f Vr:%f \n", compo, bord,d, normBord,  d_int,vitRelNorm);


                      if (d_sat < d_int && d_int <= d_act)
                        {
                          FB(compo) += 1;
                          // force de lubrification
                          double lambda = 1 / d_int - log(d_int) / 5 - d_int * log(d_int) / 21;
                          double lambda_act = 1 / d_act - log(d_act) / 5 - d_act * log(d_act) / 21;
                          //Force en N
                          force_lubrification_bord(compo, d) += (-6 * myPI * mu_fluide * rayons_compo[compo] *
                                                                 vitRelNorm * (lambda - lambda_act));

                        }

                      if (d_des < d_int && d_int <= d_sat)
                        {
                          FB(compo) += 2;
                          // force de saturation sature
                          double lambda_sat = 1 / d_sat - log(d_sat) / 5 - d_sat * log(d_sat) / 21;
                          double lambda_act = 1 / d_act - log(d_act) / 5 - d_act * log(d_act) / 21;
                          // Force en N
                          force_lubrification_bord(compo, d) += (-6 * myPI * mu_fluide * rayons_compo[compo] *
                                                                 vitRelNorm * (lambda_sat - lambda_act));

                        }

                      if (d_int <= 0)
                        {
                          FB(compo) += 4;
                          FC = 1;
                          //euler semi implicite
                          double next_positionCg = positions(compo, d) + vitesses(compo, d) * dt;
                          double next_dist_int =
                            fabs(next_positionCg - positions_bords[bord]) - rayons_compo[compo]; // delta^(n+1)

                          if (Tool::isFirstCollision(compo) == 1)
                            {
                              Tool::memorisedElongation(compo) = TR * dist_int;
                              //Tool::memorisedElongation(compo) = 0;
                              Tool::isFirstCollision(compo) = -1;
                            }

                          double next_dist_int_translate = next_dist_int - Tool::memorisedElongation(compo);
                          // collision
                          double F_spring = raideur_b * fabs(next_dist_int_translate) * normBord ;
                          double F_dashpo = -1 * fabs(amortis_b * normBord) * vitRelNorm;
                          //force en N/m^3
                          force_collision_bord(compo, d) += FS * F_spring + 1 * F_dashpo;

                        }
                      else
                        {
                          Tool::isFirstCollision(compo) = 1;
                        }
                    }

                }
              // sommations des forces
            } // end bord
          // force en N/m^3
          forces_bord(compo, d) = force_collision_bord(compo, d) / volumes_compo[compo] +
                                  1 * force_lubrification_bord(compo, d) / volumes_compo[compo];



          //-------------- contribution des particules -----------------------
          for (int parti = compo + 1; parti < nb_compo_tot; parti++)
            {
              double m_e = 1 / (1 / masse_compo[compo] + 1 / masse_compo(parti));
              double raideur_p = (m_e * (myPI * myPI + pow(log(ed), 2))) / pow(Nc * dt, 2);
              double amortis_p = (2*m_e * log(ed)) / (Nc * dt);
              //distance entre les centre de gravites des particules
              double dist_cg = 0;
              for (int j = 0; j < dimension; j++)
                {
                  double tmp = positions(compo, j) - positions(parti, j);
                  tmp *= tmp;
                  dist_cg += tmp;
                }
              dist_cg = sqrt(dist_cg);
              double dist_int = fabs(dist_cg) - (rayons_compo[compo] + rayons_compo[parti]);
              double d_int = dist_int / rayons_compo[compo];
              double Norm_d = (positions(compo, d) - positions(parti, d)) / dist_cg;

              double vitRelNorm_d, prod = 0;
              for (int j = 0; j < dimension; j++)
                {
                  prod += 1 * (vitesses(compo, j) - vitesses(parti, j)) *  // TODO verifier d'ou sort le 0.5 ...
                          (positions(compo, j) - positions(parti, j));
                }
              // composante "d" de la vitesse relative apres projection sur la normal
              vitRelNorm_d = (prod / dist_cg) * Norm_d;



              if (d_int <= d_act)
                {

                  // if (Process::je_suis_maitre()) printf("(c%d-p%d) d:%d, n:%f d_int:%f Vr:%f \n",compo,parti, d, Norm_d,  d_int,vitRelNorm);

                  if (d_sat < d_int && d_int <= d_act)
                    {
                      FP(compo) += 1;
                      FP(parti) += 1;
                      // force de lubrification
                      double lambda = 0.5 / d_int - 9 * log(d_int) / 20 - 3 * d_int * log(d_int) / 56;
                      double lambda_act = 0.5 / d_act - 9 * log(d_act) / 20 - 3 * d_act * log(d_act) / 56;
                      force_lubrification_particule(compo, d) += (-6 * myPI * mu_fluide * rayons_compo[compo] *
                                                                  vitRelNorm_d * (lambda - lambda_act));
                      force_lubrification_particule(parti, d) += -force_lubrification_particule(compo, d);

                    }

                  if (d_des < d_int && d_int <= d_sat)
                    {
                      FP(compo) += 2;
                      FP(parti) += 2;
                      // force de saturation sature
                      double lambda_sat = 0.5 / d_sat - 9 * log(d_sat) / 20 - 3 * d_sat * log(d_sat) / 56;
                      double lambda_act = 0.5 / d_act - 9 * log(d_act) / 20 - 3 * d_act * log(d_act) / 56;

                      force_lubrification_particule(compo, d) += (-6 * myPI * mu_fluide * rayons_compo[compo] *
                                                                  vitRelNorm_d * (lambda_sat - lambda_act));
                      force_lubrification_particule(parti, d) += -force_lubrification_particule(compo, d);

                    }

                  if (d_int <= 0)
                    {
                      Tool::F_now(compo,parti) = 1;
                      FP(compo) += 4;
                      FP(parti) += 4;
                      double next_dist_cg = 0;
                      for (int j = 0; j < dimension; j++)
                        {
                          double tmp = positions(compo, j) - positions(parti, j) +
                                       (vitesses(compo, j) - vitesses(parti, j)) * dt;
                          tmp *= tmp;
                          next_dist_cg += tmp;
                        }
                      next_dist_cg = sqrt(next_dist_cg);
                      double next_dist_int =(fabs(next_dist_cg) - (rayons_compo[compo] + rayons_compo[parti]))*Norm_d;

                      if(Tool::F_now(compo,parti)>Tool::F_old(compo,parti))
                        {
                          DoubleTab Norm(dimension), vitRel(dimension), vitRelNorm(dimension);
                          for (int j = 0; j < dimension; j++)
                            {
                              Norm(j)=(positions(compo, j) - positions(parti, j)) / dist_cg;
                              //vitRel(j)=(vitesses(compo, j) - vitesses(parti, j));
                              vitRelNorm(j) = (prod / dist_cg) * Norm(j);
                            }

                          double module_vitRelNorm = 0;
                          for (int j = 0; j < dimension; j++)
                            {
                              double tmp = vitRelNorm(j);
                              tmp *= tmp;
                              module_vitRelNorm += tmp;
                            }
                          module_vitRelNorm = sqrt(module_vitRelNorm);

                          double Stb=rho_solide*2*rayon*module_vitRelNorm/(9*mu_fluide);

                          Tool::raideur(compo,parti)=m_e*(module_vitRelNorm/sig_max)*(module_vitRelNorm/sig_max);
                          Tool::e_eff(compo,parti)=ed*exp(-35/Stb);


                        }

                      //double F_spring =  -1 * raideur_p * next_dist_int ;
                      //raideur_p =Tool::raideur(compo,parti);
                      double e_eff = prod <=0 ? 1 : Tool::e_eff(compo,parti) ;


                      double F_spring =  -1 * e_eff*e_eff * raideur_p * next_dist_int ;
                      double F_dashpo = amortis_p * vitRelNorm_d;
                      double F = FS * F_spring + FD * F_dashpo;

                      force_collision_particule(compo, d) += +F;
                      force_collision_particule(parti, d) += -F;


                      // test sur le sens des forces

                      //Cerr << "    >> next_dist_int: " << next_dist_int << finl;
                      //Cerr << "    >> F_spring: " << F_spring << finl;
                      //Cerr << "    >> F_dashpo: " << F_dashpo << finl;


                    }

                  forces_particule(compo, d) +=
                    (force_collision_particule(compo, d) + FL * force_lubrification_particule(compo, d)) /
                    volumes_compo[compo];
                  forces_particule(parti, d) +=
                    (force_collision_particule(parti, d) + FL * force_lubrification_particule(parti, d)) /
                    volumes_compo[parti];

                }
              Tool::F_old(compo,parti)=Tool::F_now(compo,parti);
            } // fin boucle parti
          forces_solide(compo, d) += Fb * forces_bord(compo, d) + Fp * forces_particule(compo, d); //N/m^3
        } // fin boucle dimension
    } // fin boucle compo
// ---fin calcule force de collision pour chaque composante-------------------------------------------------------------

  {


    //Cerr << "\n\n force_collision_bord: ";
    //for (int i = 0; i < dimension; i++) Cerr << " | " << force_collision_bord(0, i);
    //Cerr << "\n force_lubrification_bord: ";
    //for (int i = 0; i < dimension; i++) Cerr << " | " << force_lubrification_bord(0, i);
//
    //Cerr << "\n\n force_collision_particule: ";
    //for (int i = 0; i < dimension; i++) Cerr << " | " << force_collision_particule(0, i);
    //Cerr << "\n force_lubrification_particule: ";
    //for (int i = 0; i < dimension; i++) Cerr << " | " << force_lubrification_particule(0, i);
//
    //Cerr << "\n\n forces_bord: ";
    //for (int i = 0; i < dimension; i++) Cerr << " | " << forces_bord(0, i);
    //Cerr << "\n forces_particule: ";
    //for (int i = 0; i < dimension; i++) Cerr << " | " << forces_particule(0, i);
    //Cerr << "\n forces_solide: ";
    //for (int i = 0; i < dimension; i++) Cerr << " | " << forces_solide(0, i);
  }


// -- Application de la force de collision aux faces euleriennes correspandantes a chaque composante

  const DoubleVect& volumes_entrelaces = zone_vf.volumes_entrelaces();
  const int nb_elem = zone_vf.zone().nb_elem();

  IntVect num_compo;
  zone_vf.zone().creer_tableau_elements(num_compo);

  DoubleVect num_compo_copie;
  zone_vf.zone().creer_tableau_elements(num_compo_copie);

  const DoubleTab& xp = zone_vf.xp();

  ArrOfInt  num_lagrange(nb_compo_tot);



  for (int elem = 0; elem < nb_elem; elem++)
    {
      if (isPurementSolide)
        {
          //les elements diphasiques ne sont pas compris dans compo
          num_compo[elem] = (indicatrice[elem] !=  1-indic_phase_fluide) ? -1 : 1;
        }
      else
        {
          //les elements diphasiques sont compris dans compo
          num_compo[elem] = (indicatrice[elem] == indic_phase_fluide) ? -1 : 1;
        }
    }
  num_compo.echange_espace_virtuel();


  const IntTab& elem_faces = zone_vf.elem_faces();
  const IntTab& faces_elem = zone_vf.face_voisins();
  const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();

  //const DoubleTab& cgf = zone_vf.xv();
  const int nb_local_connex_components = search_connex_components_local(elem_faces, faces_elem, num_compo);
  const int nb_connex_components = compute_global_connex_components(num_compo, nb_local_connex_components);

  // on s'assure que le numero eulerien corresepond au bon numero lagrangien
  //remplissage du tableau de corespandance indice interface lagrange
  for (int compo = 0; compo < nb_compo_tot; compo++)
    {
      if (elem_cg[compo] != -1)
        {
          int num_euler = num_compo[elem_cg[compo]];
          num_lagrange[num_euler] = compo;
        }
    }
  mp_max_for_each_item(num_lagrange);


  //syncronisation entre les indice lagrangien et eulerien dans num_compo
  //on parcours toutes les cellule du maillage, on identifie son indice eulerien
  // on utilise l'indice eulerien  pour trouver l'indice lagrangien correspandant
  // on remplace l'indice eulerien par l'indice lagrengien

  const DoubleVect& volumes_maille = zone_vf.volumes();
  DoubleTab volumes_euler(nb_compo_tot);

  for (int elem = 0; elem < nb_elem; elem++)
    {
      if(num_compo[elem]!=-1)
        {
          num_compo[elem]=num_lagrange[num_compo[elem]];
        }
      num_compo_copie[elem]=num_compo[elem]; //conversion d'un int en double pour le postraitement

      if (indicatrice(elem)!=1)
        {
          //volume au elements
          volumes_euler[0]+=(1-indicatrice[elem])*volumes_maille[elem];
          //Cerr << finl << elem << " " << indicatrice[elem] << " " << volumes_maille[elem] << " ";
          //volume au face
          //?
        };
    }





// valeurs_champ.resize(nb_faces);
  // boucle sur les elements du domaine
  for (int elem = 0; elem < nb_elem; elem++)
    {
      int compo = num_compo[elem];
      if (compo != -1)
        {

          for (int ifac = 0; ifac < 2 * dimension; ifac++)
            {
              int face = elem_faces(elem, ifac);
              int ori = orientation(face);
              valeurs_champ(face) =
                (1 - indicatrice_faces(face)) * volumes_entrelaces(face) * forces_solide(compo, ori)*CV*
                //volume discret
                //volumes_compo(compo) / volumes_euler(compo);
                //volume theorique
                1;

            }
        }
    }
  valeurs_champ.echange_espace_virtuel();
  variables_internes().num_compo.valeur().valeurs()=num_compo_copie;  // champ pret pour postraitement
  variables_internes().num_compo.changer_temps(schema_temps().temps_courant());// champ pret pour postraitement
  //Tool::num_compo.valeur().valeurs()=num_compo_copie;  // champ pret pour postraitement




  // =============================================================================================================
  //  impression
  // =============================================================================================================
  //printf("temps: %f pc:%d size compo_connexes_facettes: %d\n",temps,Process::me(),compo_connexes_facettes.size_array() );

  if (0)
    {
      Cout << nb_connex_components <<finl ;
      std::ofstream fout3;
      fout3.open("num_compo.txt", ios::app);
      if (Process::je_suis_maitre())
        {
          fout3 << std::scientific << std::showpos;
          fout3 << "TIME: " << temps << "\t" << std::endl;
          for (int elem = 0; elem < nb_elem; elem++)
            {
              if (num_compo[elem] != -1)
                {
                  fout3 << elem << "\t\t";
                  fout3 << xp(elem, 1) << ",\t" << xp(elem, 0) << "," << xp(elem, 2) << ",";
                  fout3 << num_compo[elem] << ",\t";
                  for (int fac = 0; fac < 6; fac++)
                    {
                      fout3 << valeurs_champ(elem_faces(elem, fac)) << ",\t";
                    }
                  fout3 << std::endl;
                }
            }
          fout3 << std::endl;
        }
      fout3.close();
    }

  if (0)
    {
      Cerr << "2end pass \t" << finl;

      Cerr << "compo: 0 \t";
      for (int d = 0; d < 3; d++)
        {
          Cerr << positions(0, d) << "\t";
        }
      Cerr << finl;
      Cerr << "compo: 1 \t";
      for (int d = 0; d < 3; d++)
        {
          Cerr << positions(1, d) << "\t";
        }
      Cerr << finl;
    }
  if (0)
    {

      for (int compo = 0; compo < nb_compo_tot; compo++) //Todo
        {
          int elem = elem_cg(compo);
          if (compo != num_compo[elem])
            {
              Cerr << "ERROR ";
            }
          Cerr << "(N): " << compo << "\t" << num_compo[elem] << "\t" << elem_cg[compo] << "-" << elem << "\t"
               << temps << "\t\t";
          for (int d = 0; d < 3; d++)
            Cerr << positions(compo, d) << "\t";

          for (int d = 0; d < 3; d++)
            Cerr << xp(elem, d) << "\t";

          Cerr << finl;
        }





    }
// =============================================================================================================
  if(0)
    {
      Cerr << "(V):T " << volumes_compo(0) << "\t(V):SM" << Nom(volumes_euler[0], "%20.14g")  << "\t" ;
      Cerr << finl;

      Cerr << "(F):T ";
      for (int j = 0; j < dimension; j++) Cerr << Nom(forces_solide(0,j)*volumes_compo(0), "%20.14g") << "\t";
      Cerr << finl;

    }
  if (1)
    {
      std::ofstream fout;
      fout.open("forces_composantes_connexes.txt", ios::app);
      //fout << "TEMPS: " << temps << std::endl;
      for (int compo = 0; compo < nb_compo_tot; compo++)
        {
          //if (FB(compo) != 0 || FP(compo) != 0)
          if (print == 1)
            {
              fout << std::scientific << std::showpos;
              fout << "compo: " << compo << "\tFB: " << FB(compo) << "\tFP: " << FP(compo);
              fout << "\t\tP: ";
              for (int i = 0; i < dimension; i++) fout << " " << positions(compo, i);
              fout << "\t\tV: ";
              for (int i = 0; i < dimension; i++) fout << " " << vitesses(compo, i);
              fout << "\t\tFLB: ";
              for (int i = 0; i < dimension; i++) fout << " " << force_lubrification_bord(compo, i);
              fout << "\t\tFCB: ";
              for (int i = 0; i < dimension; i++) fout << " " << force_collision_bord(compo, i);
              fout << "\t\tFB: ";
              for (int i = 0; i < dimension; i++) fout << " " << forces_bord(compo, i);
              fout << "\t\tFP: ";
              for (int i = 0; i < dimension; i++) fout << " " << forces_particule(compo, i);
              fout << std::endl;
            }

          if (print == 2 && Process::je_suis_maitre())
            {
              //int bord =1, d=1;
              //double posCompo_np1 = positions(compo,d) + vitesses(compo,d)*dt ;
              //double dist_int_np1 = fabs(posCompo_np1 - positions_bords[bord]) - rayons_compo[compo];

              fout << std::scientific << std::showpos;
              fout << compo << " [" << FP(compo)/3 << ",\t" << temps << ",\t" << positions(compo, 1) << ",\t"
                   << vitesses(compo, 1) << ",\t";
              fout << forces_solide(compo, 1) << ",\t" << impression_d_int(compo) << ",\t"
                   << impression_dist_int_tr(compo) << ",\t" << impression_next_dist_int_tr(compo) << "], ";
              fout << Tool::isFirstCollision(compo) << ",\t" << Tool::memorisedElongation(compo);
              fout << std::endl;
            }
        }
      fout.close();

      if (1)
        {
          std::ofstream fout2;
          fout2.open("profil_compo.txt", ios::app);
          if (Process::je_suis_maitre())
            {
              fout2 << std::scientific << std::showpos;
              fout2 << temps;
              for (int compo = 0; compo < nb_compo_tot; compo++)
                {
                  fout2 << '\t';
                  for (int i = 0; i < dimension; i++) fout2 << " " << positions(compo, i);
                  fout2 << "  ";
                  for (int i = 0; i < dimension; i++) fout2 << " " << vitesses(compo, i);

                  //fout2 << "compo:" << compo << '\t' << positions(compo, 1) << '\t' << vitesses(compo, 1) << '\t';
                }

              fout2 << std::endl;

            }
          fout2.close();
        }
    }

  // =============================================================================================================
// fin
}

//</editor-fold>

void Navier_Stokes_FT_Disc::calculer_champ_forces_collisions2(const DoubleTab& indicatrice, DoubleTab& valeurs_champ,
                                                              double& FC)
{
  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_forces_collisions", 0);
  statistiques().begin_count(count);

//

//<editor-fold desc="Determiniations des vitesses et positions de chaques particule">
  Tool::compteur_ += 1;
  Cerr << "Navier_Stokes_FT_Disc::calculer_champ_forces_collisions" << finl;
  //auto start = high_resolution_clock::now();

  REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
  const Zone_VF& zone_vdf = ref_cast(Zone_VDF, zone_dis().valeur());
  // recuperation des vitesse compo et des centres de gravites

  DoubleTab& positions = Tool::positions_compo;
  DoubleTab& vitesses = Tool::vitesses_compo;

  const int nb_facettes = maillage.nb_facettes();
  ArrOfInt compo_connexes_facettes(nb_facettes); // Init a zero
  int n = search_connex_components_local_FT(maillage, compo_connexes_facettes); //
  int nb_compo_tot = compute_global_connex_components_FT(maillage, compo_connexes_facettes, n); //

  int nb_bord = 2 * dimension;
  ArrOfDouble positions_bords(nb_bord);

  assert(nb_compo_tot != 0);
  // Cerr << " Zone 1" << finl ;
  // si vitesses na pas les bonnes dimension (premiere pas de temps, en calcule la position des centre de gravite et
  // on initialise les vitesse des inclusions a zero
  if ((vitesses.dimension(0) != nb_compo_tot) || (vitesses.dimension(1) != dimension))
    if (Tool::compteur_ == 1)
      {
        // on initialise les vitesse avec 0
        positions.resize(nb_compo_tot, dimension);
        vitesses.resize(nb_compo_tot, dimension);
        vitesses = 0;

        // calcule des centre de gravite pour les composantes
        assert(nb_compo_tot == positions.dimension(0));

        const int dim = positions.dimension(1); //
        const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
        const IntTab& facettes = maillage.facettes();
        const DoubleTab& sommets = maillage.sommets();
        assert(facettes.dimension(1) == dim);

        // Surface totale de chaque composante connexe, initialise a zero
        ArrOfDouble surfaces_compo(nb_compo_tot);
        positions = 0.;

        // Calcul du centre de gravite de la composante connexe
        //  (centre de gravite de la surface, pas du volume)
        const int nb_facettes_tot = facettes.dimension_tot(0);
        {
          for (int i = 0; i < nb_facettes_tot; i++)
            {
              if (maillage.facette_virtuelle(i)) continue;

              const int compo = compo_connexes_facettes[i];
              const double surface = surface_facettes[i];
              surfaces_compo[compo] += surface;
              // Centre de gravite de la facette, pondere par la surface
              for (int j = 0; j < dim; j++)
                {
                  // Indice du sommet
                  const int s = facettes(i, j);
                  for (int k = 0; k < dim; k++)
                    {
                      // On divisera par dim a la fin:
                      positions(compo, k) += surface * sommets(s, k);
                    }
                }
            }
          mp_sum_for_each_item(surfaces_compo);
          mp_sum_for_each_item(positions);

          positions *= (1. / dim);
          DoubleVect s; // tab_divide prend DoubleVect, pas ArrOfDouble...
          s.ref_array(surfaces_compo);
          tab_divide_any_shape(positions, s);
          //---fin calcul surface_compo et position_centre_gravite_compo ( a rassembler dans une focntion)---//
        }


      }


  ArrOfInt elem_cg;
  zone_vf.zone().chercher_elements(positions, elem_cg);


  //initialisation
  DoubleTab correct_vitesses(nb_compo_tot, dimension), correct_positions(nb_compo_tot, dimension);
  ArrOfInt copy_elem_cg, correctNum(nb_compo_tot);
  Tool::F_now.resize(nb_compo_tot, nb_compo_tot + nb_bord);
  Tool::F_now = 0;
  if (Tool::compteur_ == 1)
    {
      Tool::F_old.resize(nb_compo_tot, nb_compo_tot + nb_bord);
      Tool::raideur.resize(nb_compo_tot, nb_compo_tot + nb_bord);
      Tool::e_eff.resize(nb_compo_tot, nb_compo_tot + nb_bord);

      Tool::F_old = 0;
      Tool::raideur = -1;
      Tool::e_eff = -1;
    }
  else
    {
      // correspandance entre le bon numero eulerian au temps precedant et le mauvais numero lagrangien actuel
      DoubleVect old_num_compo = variables_internes().num_compo.valeur().valeurs();

      {

        zone_vf.zone().chercher_elements(positions, copy_elem_cg);

        for (int wrong_num = 0; wrong_num < nb_compo_tot; wrong_num++)
          {
            int elem = elem_cg[wrong_num];
            int correct_num = elem == -1 ? -1 : old_num_compo[elem]; // TODO  verefier conversion int <-> double
            correctNum(wrong_num) = correct_num;

          }
        mp_max_for_each_item(correctNum);
        int isduplicateValue = Tool::checkForDuplicates(correctNum);
        if (isduplicateValue == 1) exit("ERROR: duplicate value");
      }

      for (int bad_compo = 0; bad_compo < nb_compo_tot; bad_compo++)
        {

          int good_compo = correctNum[bad_compo];

          for (int d = 0; d < dimension; d++)
            {


              correct_vitesses(good_compo, d) = vitesses(bad_compo, d);
              correct_positions(good_compo, d) = positions(bad_compo, d);
            }
        }

      positions = correct_positions;
      vitesses = correct_vitesses;
      zone_vf.zone().chercher_elements(positions, elem_cg);
    }
  //</editor-fold>


//<editor-fold desc="Recuperation des proprietes diphasiques">
  const IntVect& orientation = zone_vdf.orientation();
  //  ----- proprietes particules  -------------
  ArrOfDouble rayons_compo(nb_compo_tot);
  ArrOfDouble volumes_compo(nb_compo_tot);
  ArrOfDouble masse_compo(nb_compo_tot);

//--------------------------------------------------
  const Fluide_Diphasique& mon_fluide = fluide_diphasique();

  const double dt = schema_temps().pas_de_temps();
  const double temps = schema_temps().temps_courant();

  const double rho_solide = mon_fluide.fluide_phase(0).masse_volumique().valeurs()(0, 0);
  const double rho_fluide = mon_fluide.fluide_phase(1).masse_volumique().valeurs()(0, 0);
  const double mu_fluide =
    rho_fluide * mon_fluide.fluide_phase(1).viscosite_cinematique().valeur().valeurs()(0, 0);

  int indic_phase_fluide = 1; //, indic_phase_solide=0; // TODO a generaliser
  double longeur_maille = Tool::calc_positions_bords2(positions_bords);
// -------------------------------------------------


// parametre du modele -----------------------------
  double rayon = Tool::myRayon;
  double sigma = Tool::mySigma; //distance de  d'enfoncement max (Ec=0, Epe=max)
  double ed = 0.97;
  //int Nc = 8;
  // FLAGS
  double CV = 1, myPI = 3.14159265358979;
  double print_next_dist_int=-1;
//----------------------------------------------------

  for (int compo = 0; compo < nb_compo_tot; compo++)
    {
      rayons_compo[compo] = rayon;
      volumes_compo[compo] = 4 * myPI * pow(rayons_compo[compo], 3) / 3;
      masse_compo[compo] = volumes_compo[compo] * rho_solide;
    }


  std::ofstream fout;
  fout.open("details_colisions.txt", ios::app);
  fout << std::scientific;
  //</editor-fold>


//<editor-fold desc="Calcule des forces de collisions">
  DoubleTab forces_solide(nb_compo_tot, dimension);
  for (int compo = 0; compo < nb_compo_tot; compo++)
    {

      for (int voisin = compo + 1; voisin < nb_compo_tot + nb_bord; voisin++)
        {

          //<editor-fold desc="Calcule de la distance entre deux interfaces">
          DoubleTab dX(dimension), dU(dimension);
          dX = 0;
          dU = 0;
          int CollisionParticuleParticule = voisin < nb_compo_tot;
          if (CollisionParticuleParticule)
            {
              for (int d = 0; d < dimension; d++)
                {
                  dX(d) = positions(compo, d) - positions(voisin, d);
                  dU(d) = vitesses(compo, d) - vitesses(voisin, d);
                }

            }
          else
            {
              int bord = voisin - nb_compo_tot;
              int ori = bord < dimension ? bord : bord - dimension;
              dX(ori) = positions(compo, ori) - positions_bords(bord);
              for (int d = 0; d < dimension; d++) dU(d) = vitesses(compo, d);
            }

          double dist_cg = Tool::module_vecteur(dX);
          if (dist_cg == 0)
            {
              Cerr << "ERROR : dist_cg = 0 entre " << compo << " et " << voisin << finl;
              exit();
            }

          double dist_int = CollisionParticuleParticule ? dist_cg - (rayons_compo(compo) + rayons_compo(voisin)) :
                            dist_cg - 2 * rayons_compo(compo);
          //</editor-fold>


          //<editor-fold desc="Calcule de la norme et de la vitesse relative normale">
          DoubleTab norm(dimension);
          for (int d = 0; d < dimension; d++) norm(d) = dX(d) / dist_cg;
          double prod_sacl = Tool::prod_scal(dX, dU);
          DoubleTab dUn(dimension);
          for (int d = 0; d < dimension; d++) dUn(d) = (prod_sacl / dist_cg) * norm(d);
          double vitesseRelNorm = Tool::module_vecteur(dUn);
          //</editor-fold>

          //----------------------------

          //<editor-fold desc="Modele de lubrification">
          int isModeleLubrification =Tool::modele_lubrification ;
          double d_int = fabs(dist_int) / rayons_compo[compo];
          double d_act = 2.0 * longeur_maille / rayons_compo[compo];
          double d_sat = 0.1 * longeur_maille / rayons_compo[compo];

          if (isModeleLubrification && d_int <= d_act )
            {


              double lambda     = CollisionParticuleParticule ? 0.5 / d_int - 9 * log(d_int) / 20 - 3 * d_int * log(d_int) / 56 : 1 / d_int - log(d_int) / 5 - d_int * log(d_int) / 21;
              double lambda_act = CollisionParticuleParticule ? 0.5 / d_act - 9 * log(d_act) / 20 - 3 * d_act * log(d_act) / 56 : 1 / d_act - log(d_act) / 5 - d_act * log(d_act) / 21;
              double lambda_sat = CollisionParticuleParticule ? 0.5 / d_sat - 9 * log(d_sat) / 20 - 3 * d_sat * log(d_sat) / 56 : 1 / d_sat - log(d_sat) / 5 - d_sat * log(d_sat) / 21;

              double delta_lambda=0;

              if (d_sat < d_int && d_int <= d_act)
                delta_lambda= (lambda - lambda_act);


              if (0 < d_int && d_int <= d_sat)
                delta_lambda= (lambda_sat - lambda_act);


              for (int d = 0; d < dimension; d++)
                {
                  double force_lubrification = -6 * myPI * mu_fluide * rayons_compo[compo] * dUn(d) * delta_lambda;
                  Cerr << "(*) d: " << d << "  force_lubrification:  " <<force_lubrification <<finl;
                  continue;
                  forces_solide(compo, d) += +force_lubrification / volumes_compo(compo);
                  if (!CollisionParticuleParticule) continue; //collision avec un bord -> pas de force sur le bord
                  forces_solide(voisin, d) += -force_lubrification / volumes_compo(voisin);
                }


            }
          //</editor-fold>


          Tool::F_now(compo, voisin) = 0;
          if (dist_int <= 0) //contact
            {

              Tool::F_now(compo, voisin) = 1;
              int isFirstStepOfCollision = Tool::F_now(compo, voisin) > Tool::F_old(compo, voisin);

              //<editor-fold desc="schema semi implicte">
              DoubleTab next_dX(dimension);
              for (int d = 0; d < dimension; d++) next_dX(d) = dX(d) + dt * dU(d);
              double next_dist_cg = Tool::module_vecteur(next_dX);
              double next_dist_int = CollisionParticuleParticule ? next_dist_cg -
                                     (rayons_compo(compo) + rayons_compo(voisin)) :
                                     next_dist_cg - 2 * rayons_compo(compo);
              print_next_dist_int= next_dist_int;
              //</editor-fold>

              //new writing
              if(1)
                {
                  double rayon_eff = CollisionParticuleParticule ? 0.5 *
                                     (rayons_compo(compo) + rayons_compo(voisin))
                                     : rayons_compo(compo); //bord

                  double masse_eff = CollisionParticuleParticule ? 1 / (1 / masse_compo(compo) +
                                                                        1 / masse_compo(voisin))
                                     : masse_compo(compo); //bord

                  double Stb = rho_solide * 2 * rayon_eff * vitesseRelNorm / (9 * mu_fluide);

                  DoubleTab force_contact(dimension);
                  const int modele=modele_collisions();
                  switch (modele)
                    {

                    case 0:
                      {
                        //Cerr << "Modele ressort_amorti. \n" << finl;

                        double raideur = Tool::param_osillateur(0) ;
                        double amortisseur = Tool::param_osillateur(1) ;

                        for (int d = 0; d < dimension; d++)
                          force_contact(d)=-1 * raideur * next_dist_int * norm(d) -1*amortisseur*dUn(d);
                      }
                      break;

                    case 1:
                      {
                        //Cerr << "Modele Mohagheg. \n" << finl;

                        //int isFirstStepOfCollision = Tool::F_now(compo, voisin) > Tool::F_old(compo, voisin);
                        if (isFirstStepOfCollision)
                          {
                            Tool::e_eff(compo, voisin) = Stb > 18 ? ed * (1 - 8.65 * pow(Stb, -0.75)) : 0;
                            Tool::raideur(compo, voisin) = masse_eff * pow(vitesseRelNorm / sigma, 2);
                          }
                        int isPhasePenetration = prod_sacl <= 0;
                        double e_eff = isPhasePenetration ? 1 : Tool::e_eff(compo, voisin);
                        double raideur = Tool::raideur(compo, voisin);

                        for (int d = 0; d < dimension; d++)
                          force_contact(d)=-1 * e_eff * e_eff * raideur * next_dist_int * norm(d);
                      }
                      break;

                    case 2:
                      {
                        //Cerr << "Modele hybrid. \n" << finl;

                        //int isFirstStepOfCollision = Tool::F_now(compo, voisin) > Tool::F_old(compo, voisin);
                        if (isFirstStepOfCollision)
                          {
                            double  tau_c = 8 * dt ;
                            Tool::e_eff(compo, voisin) = ed * exp(-35 / (Stb + 1e-6));
                            Tool::raideur(compo, voisin) = (masse_eff * (myPI * myPI + pow(log(ed), 2))) / pow(tau_c, 2);
                          }
                        int isPhasePenetration = prod_sacl <= 0;
                        double e_eff = isPhasePenetration ? 1 : Tool::e_eff(compo, voisin);
                        double raideur = Tool::raideur(compo, voisin);

                        for (int d = 0; d < dimension; d++)
                          force_contact(d)=-1 * e_eff * e_eff * raideur * next_dist_int * norm(d);

                      }
                      break;
                    case 3:
                      {
                        //Cerr << "Modele Breugem. \n" << finl;
                        int Nc=8;
                        double raideur = (masse_eff * (myPI * myPI + pow(log(ed), 2))) / pow(Nc * dt, 2);
                        double amortisseur = -1*(masse_eff * log(ed)) / (Nc * dt);

                        for (int d = 0; d < dimension; d++)
                          force_contact(d)=-1 * raideur * next_dist_int * norm(d) -1*amortisseur*dUn(d);
                      }
                      break;


                    default:
                      Cerr << "The method specified for modele_collision in not recognized. \n" << finl;
                      Cerr << "pls use one of : ressort_amorti,  double_raideurs,  hybrid ! your input is  :"
                           << modele << " \n" << finl;
                      Process::exit();
                    }

                  //<editor-fold desc="Impression en debut de contact solide">
                  if (isFirstStepOfCollision && Process::je_suis_maitre())
                    {

                      fout << " # START!!! P" <<compo <<"V" <<voisin<<"  \n# t_imp= "<<temps<<" \n# Vrn_imp= "<<vitesseRelNorm<<"  \n# St_imp= "<< Stb;

                      fout<<std::endl;
                      fout << "'IMP',"<<compo<<","<<voisin<<","<<temps<<","<<vitesseRelNorm << ","<<Stb <<","<<Tool::e_eff(compo, voisin) <<","<<Tool::raideur(compo, voisin) <<std::endl;
                    }
                  //</editor-fold>

                  for (int d = 0; d < dimension; d++)
                    {

                      //impression sans application a supprimer une fois l'implementation fonctionelle
                      Cerr << "(*) d: " << d << "  force_contact:  " <<force_contact(d) <<finl;
                      // continue;

                      forces_solide(compo, d) += +force_contact(d) / volumes_compo(compo);
                      if (!CollisionParticuleParticule) continue; // collision avec un bord -> pas de force sur le bord
                      forces_solide(voisin, d) += -force_contact(d) / volumes_compo(voisin);
                    }
                }

              //<editor-fold desc="Ipression collision encours ...">
              if(Process::je_suis_maitre() && compo==0 && !isFirstStepOfCollision )
                {
                  fout << "#  COL!!! P" <<compo <<"V" <<voisin<<"   t= "<<temps;
                  DoubleTab fc(dimension); // variable pour le calcule du module : la fonction prend un vecteur a entree unique
                  for (int d = 0; d < dimension; d++) fc(d) = forces_solide(compo,d);
                  fout <<"     DistInt= "<<dist_int  <<"     nDistInt= "<<next_dist_int  <<"     collision_force= "<<Tool::module_vecteur(fc)<< std::endl;
                }
              //</editor-fold>

            }


          //<editor-fold desc="Fin du contact solide">
          int isLastCollision = Tool::F_now(compo, voisin) < Tool::F_old(compo, voisin);
          if (isLastCollision)
            {
              //double  e_calc = vitesseRelNorm/Tool::vitessRelImp ;
              if (Process::je_suis_maitre())
                {
                  fout << "#  END!!! P" <<compo <<"V" <<voisin<<"  \n# t_reb= "<<temps<<" \n# Vrn_reb= "<<vitesseRelNorm  <<std::endl;
                  fout << "'REB',"<<compo<<","<<voisin<<","<<temps<<","<<vitesseRelNorm << std::endl;
                  //fout << "e_app= "<<e_calc << std::endl;
                  //fout << "error= " << 100*fabs(Tool::e_eff(compo,voisin)-e_calc)/Tool::e_eff(compo,voisin)<<std::endl;
                }

            }
          //</editor-fold>

          Tool::F_old(compo, voisin) = Tool::F_now(compo, voisin);
        } // fin boucle parti

    } // fin boucle compo
// ---fin calcule force de collision pour chaque composante-------------------------------------------------------------
  fout.close();
  //</editor-fold>


//<editor-fold desc="Application des forces de collisions sur les faces du domaine">
  const DoubleVect& volumes_entrelaces = zone_vf.volumes_entrelaces();
  const int nb_elem = zone_vf.zone().nb_elem();

  IntVect num_compo;
  zone_vf.zone().creer_tableau_elements(num_compo);

  DoubleVect num_compo_copie;
  zone_vf.zone().creer_tableau_elements(num_compo_copie);

  // const DoubleTab& xp = zone_vf.xp();

  ArrOfInt num_lagrange(nb_compo_tot);


  for (int elem = 0; elem < nb_elem; elem++)
    {
      if (Tool::force_sur_elem_diphasiques)
        {
          //les elements diphasiques sont compris dans num compo
          num_compo[elem] = (indicatrice[elem] == indic_phase_fluide) ? -1 : 1;
        }
      else
        {
          //les elements diphasiques ne sont pas compris dans num compo
          num_compo[elem] = (indicatrice[elem] != 1 - indic_phase_fluide) ? -1 : 1;
        }
    }
  num_compo.echange_espace_virtuel();

  const IntTab& elem_faces = zone_vf.elem_faces();
  const IntTab& faces_elem = zone_vf.face_voisins();
  const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();
  //const DoubleTab& cgf = zone_vf.xv();
  const int nb_local_connex_components = search_connex_components_local(elem_faces, faces_elem, num_compo);
  const int nb_connex_components = compute_global_connex_components(num_compo, nb_local_connex_components);

  // on s'assure que le numero eulerien corresepond au bon numero lagrangien
  //remplissage du tableau de corespandance indice interface lagrange
  for (int compo = 0; compo < nb_compo_tot; compo++)
    {
      int elem = elem_cg[compo];
      if (elem != -1)
        {
          int num_euler = num_compo[elem];
          num_lagrange[num_euler] = compo;
        }
    }
  mp_max_for_each_item(num_lagrange);


  //syncronisation entre les indice lagrangien et eulerien dans num_compo
  //on parcours toutes les cellule du maillage, on identifie son indice eulerien
  // on utilise l'indice eulerien  pour trouver l'indice lagrangien correspandant
  // on remplace l'indice eulerien par l'indice lagrengien

  const DoubleVect& volumes_maille = zone_vf.volumes();
  DoubleTab volumes_euler(nb_compo_tot);

  for (int elem = 0; elem < nb_elem; elem++)
    {
      if (num_compo[elem] != -1)
        {
          num_compo[elem] = num_lagrange[num_compo[elem]];
        }
      num_compo_copie[elem] = num_compo[elem]; //conversion d'un int en double pour le postraitement

      if (indicatrice(elem) != 1)
        {
          //volume au elements
          volumes_euler[0] += (1 - indicatrice[elem]) * volumes_maille[elem];
          //Cerr << finl << elem << " " << indicatrice[elem] << " " << volumes_maille[elem] << " ";
          //volume au face
          //?
        }
    }




// valeurs_champ.resize(nb_faces);
  // boucle sur les elements du domaine
  for (int elem = 0; elem < nb_elem; elem++)
    {
      int compo = num_compo[elem];
      if (compo != -1)
        {

          for (int ifac = 0; ifac < 2 * dimension; ifac++)
            {
              int face = elem_faces(elem, ifac);
              int ori = orientation(face);
              valeurs_champ(face) =
                (1 - indicatrice_faces(face)) * volumes_entrelaces(face) * forces_solide(compo, ori) * CV *
                //volume discret
                //volumes_compo(compo) / volumes_euler(compo);
                //volume theorique
                1;

            }
        }
    }
  valeurs_champ.echange_espace_virtuel();
  //Tool::num_compo_=num_compo_copie;
  variables_internes().num_compo.valeur().valeurs() = num_compo_copie;  // champ pret pour postraitement
  //variables_internes().num_compo.changer_temps(temps);
  //Cerr << "temps de numm compo: " <<variables_internes().num_compo.temps() <<finl;
  //</editor-fold>


//<editor-fold desc="Impression">
  if (0)
    {
      if (Process::je_suis_maitre() && 1)
        fout << std::scientific;
      Cerr << nb_connex_components << longeur_maille << print_next_dist_int << finl;
    }

  if (Process::je_suis_maitre())
    {
      {
        std::ofstream fout2;
        fout2.open("profil_compo.txt", ios::app);
        fout2 << std::scientific << std::showpos;
        fout2 << temps;
        for (int compo = 0; compo < nb_compo_tot; compo++)
          {
            fout2 << '\t';
            for (int i = 0; i < dimension; i++) fout2 << " " << positions(compo, i);
            fout2 << "  ";
            for (int i = 0; i < dimension; i++) fout2 << " " << vitesses(compo, i);
            //fout2 << "compo:" << compo << '\t' << positions(compo, 1) << '\t' << vitesses(compo, 1) << '\t';
          }
        fout2 << std::endl;
        fout2.close();
      }

    }
  //</editor-fold>

// fin
  statistiques().end_count(count);
}

// Description:
//  Calcul du gradient de l'indicatrice.
//  Ce gradient est utilise pour calculer le second membre de l'equation de qdm,
//  contenant les termes de tension de surface.
//  En VEF, on commence par creer un champ P1B a partir du champ P0
//  et on calcule le gradient.
//  Design de classe a revoir pour separer VDF et VEF...
//
// Parametre : indicatrice
// Signification : un champ aux elements (l'espace virtuel doit etre a jour)
// Parametre : gradient_i
// Signification : un champ discretise comme la vitesse dans lequel
//  on met gradient(indicatrice).

void Navier_Stokes_FT_Disc::calculer_gradient_indicatrice(
  const Champ_base& indicatrice,
  const DoubleTab& distance_interface_sommets,
  Champ_base& gradient_i)
{
  if (gradient_i.que_suis_je() == "Champ_Fonc_Face")
    {

      gradient.calculer(indicatrice.valeurs(), gradient_i.valeurs());

    }
  else
    {
      int i;
      const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
      const Zone& zone = zone_vf.zone();
      const int nb_elem_tot = zone.nb_elem_tot();
      //const int nb_sommets = zone.nb_som();
      //const IntTab & les_elems = zone.les_elems();
      //const int nb_sommets_par_element = les_elems.dimension(1);

      // Calcul d'une indicatrice p1bulle
      DoubleTab& indic_p1b = variables_internes().indicatrice_p1b;
      // Verification du support du champ indicatrice_p1b
      if (dimension==2 && indic_p1b.size_totale()!=nb_elem_tot+zone.nb_som_tot())
        {
          Cerr << "The method Navier_Stokes_FT_Disc::calculer_gradient_indicatrice is developped" << finl;
          Cerr << "only for an indicatrice field discretized like P0+P1 for the 2D dimension case."<<finl;
          Cerr << "Please change the discretization." << finl;
          Process::exit();
        }
      if (dimension==3 && indic_p1b.size_totale()!=nb_elem_tot+zone.nb_som_tot()+zone.nb_aretes_tot() && indic_p1b.size_totale()!=nb_elem_tot+zone.nb_som_tot())
        {
          Cerr << "The method Navier_Stokes_FT_Disc::calculer_gradient_indicatrice is developped" << finl;
          Cerr << "only for an indicatrice field discretized like P0+P1 or P0+P1+Pa for the 3D dimension case."<<finl;
          Cerr << "Please change the discretization." << finl;
          Process::exit();
        }

      // Le champ P1bulle contient
      //   nelements valeurs au centre des elements
      // suivies de
      //   nsommets valeurs aux sommets
      indic_p1b = 0.;
      // Valeurs aux centres des elements = indicatrice a l'element
      // On recopie des valeurs sur les elements virtuels car on en
      // a besoin pour le calcul des valeurs aux sommets.
      for (i = 0; i < nb_elem_tot; i++)
        {
          indic_p1b(i) = indicatrice(i);
        }
#if 0
      // Premiere version : pas terrible parce que le stencil est tres large
      //  et il faut calculer la "courbure de l'interface" sur une grande largeur
      //  autour de l'interface.

      // Valeurs aux sommets = moyenne des valeurs des elements adjacents
      //  ponderee par le volume des elements.
      // Il y a la un choix a faire qui peut etre important...
      // (autre idee : projection L2 du champ discontinu sur l'espace
      // P1bulle)
      const DoubleVect& volume = zone_vf.volumes();

      // Somme des poids aux sommets:
      ArrOfDouble poids(nb_sommets);
      poids = 0.;
      // Boucle sur les elements reels et virtuels :
      for (i = 0; i < nb_elem_tot; i++)
        {
          // Ajout d'une contribution a chaque sommet reel de l'element
          double p = volume(i);
          double x = indicatrice(i);
          int j;
          for (j = 0; j < nb_sommets_par_element; j++)
            {
              const int sommet = les_elems(i, j);
              if (sommet < nb_sommets)   // sommet reel ?
                {
                  // Le tableau des ddl contient nb_elem_tot valeurs aux elements
                  // suivies de nb_sommets_tot valeurs aux sommets.
                  // L'inconnue associee au "sommet" se trouve a l'indice
                  // nb_elem_tot + sommet.
                  indic_p1b(nb_elem_tot + sommet) += p * x;
                  poids[sommet] += p;
                }
            }
        }
      // Division par le poids
      for (i = 0; i < nb_sommets; i++)
        {
          const double p = poids[i];
          indic_p1b(nb_elem_tot + i) /= p;
        }
#else
#if 0
      // Calcul de l'indicatrice aux sommets. Deuxieme version, pour
      // limiter le stencil de gradient_indicatrice :
      // l'indicatrice vaut 1 si le sommet est adjacent a un element
      // de phase 1, 0 si le sommet est adjacent a un element de phase 0

      // 1) on met une valeur invalide dans les inconnues:
      for (i = 0; i < nb_sommets; i++)
        indic_p1b(nb_elem_tot + i) = -1;
      // 2) On met 1. ou 0. pour les sommets adjacents a un element
      //    dont l'indicatrice vaut 1. ou 0.
      for (i = 0; i < nb_elem_tot; i++)
        {
          const double indic_element = indic_p1b(i);
          if (indic_element == 0. || indic_element == 1.)
            {
              int j;
              for (j = 0; j < nb_sommets_par_element; j++)
                {
                  const int sommet = les_elems(i,j);
                  indic_p1b(nb_elem_tot + sommet) = indic_element;
                }
            }
        }
      // 3) Pour les sommets indecis, on prend 0 ou 1 selon le signe
      //    de la fonction distance pour ce sommet.
      static const double valeur_distance_invalide = -1e30;
      int error_count = 0;
      for (i = 0; i < nb_sommets; i++)
        {
          const double valeur = indic_p1b(nb_elem_tot + i);
          if (valeur < 0.)
            {
              double newval;
              const double d = distance_interface_sommets(i);
              if (d < valeur_distance_invalide)
                {
                  error_count++;
                  newval = 0.;
                }
              else
                {
                  newval = (d >= 0.) ? 1. : 0.;
                }
              indic_p1b(nb_elem_tot + i) = newval;
            }
        }
      if (error_count > 0)
        {
          Process::Journal()
              << "Navier_Stokes_FT_Disc::calculer_gradient_indicatrice\n"
              << error_count << " sommets ont une valeur indeterminee" << finl;
        }
#else
      // On met l'indicatrice sur P0 et 0 sur P1
#endif
#endif
      indic_p1b.echange_espace_virtuel();

      gradient.calculer(indic_p1b, gradient_i.valeurs());
    }
}

// Description:
//  Calcul du saut de vitesse a l'interface du au changement de phase
//  phase_pilote = -1: u+u0 = champ de vitesse de deplacement de l'interface
//  phase_pilote = 0 : u+u0 = champ de vitesse de la phase 0
//  phase_pilote = 1 : u+u0 = champ de vitesse de la phase 1
//  ordre = 0 : pas de prise en compte de la correction en courbure
//  ordre = 1 : prise en compte de la correction en courbure a l'ordre 1
//  ordre = 2 : prise en compte de la correction en courbure a l'ordre 2
void Navier_Stokes_FT_Disc::calculer_delta_u_interface(Champ_base& champ_u0,
                                                       int phase_pilote,
                                                       int ordre)
{
  //  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_delta_u_interface", 0);

  //  statistiques().begin_count(count);

  // GB : REMPLACEMENT TEMPORAIRE DE :
  //   variables_internes().ref_equation_mpoint_.valeur().calculer_mpoint(variables_internes().mpoint.valeur());
  //   const DoubleTab & mpoint = variables_internes().mpoint.valeur().valeurs();
  // PAR  (voir un peu plus bas)

  REF(Transport_Interfaces_FT_Disc) & refeq_transport =
    variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  // Distance a l'interface discretisee aux elements:
  const DoubleTab& dist = eq_transport.get_update_distance_interface().valeurs();
  DoubleTab phi = calculer_div_normale_interface().valeurs();

  // GB : debut de mon rajout pour permettre de calculer mpoint, mais sans l'utiliser pour le deplacement de l'interface.
  // GB : Pour ce faire, il suffit de ne pas mettre le mot cle " equation_temperature_mpoint temp" dans le jdd.
  const int nn = phi.dimension(0);
  DoubleTab mpoint;
  DoubleTab mpoint_vap;
  mpoint.resize(nn);
  mpoint_vap.resize(nn);
  mpoint=0.; //pour initialiser
  mpoint_vap=0.; //pour initialiser
  if (variables_internes().ref_equation_mpoint_.non_nul())
    {
      variables_internes().ref_equation_mpoint_.valeur().calculer_mpoint(variables_internes().mpoint.valeur());
      // Si inactif, on ne prend pas en compte sa contribution dans le calcul de delta_u:
      if (!variables_internes().mpoint_inactif)
        mpoint = variables_internes().mpoint.valeur().valeurs();
    }
  if (variables_internes().ref_equation_mpoint_vap_.non_nul())
    {
      variables_internes().ref_equation_mpoint_vap_.valeur().calculer_mpoint(variables_internes().mpoint_vap.valeur());
      // Si inactif, on ne prend pas en compte sa contribution dans le calcul de delta_u:
      if (!variables_internes().mpointv_inactif)
        mpoint_vap = variables_internes().mpoint_vap.valeur().valeurs();
    }
  // GB : fin de mon rajout
  if (variables_internes().ref_equation_mpoint_vap_.non_nul())
    mpoint +=mpoint_vap;

  DoubleTab& u0 = champ_u0.valeurs();

  {
    const Fluide_Diphasique& fluide_dipha = fluide_diphasique();
    const Fluide_Incompressible& phase_0 = fluide_dipha.fluide_phase(0);
    const Fluide_Incompressible& phase_1 = fluide_dipha.fluide_phase(1);
    const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique().valeurs();
    const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique().valeurs();
    const double rho_0 = tab_rho_phase_0(0,0);
    const double rho_1 = tab_rho_phase_1(0,0);
    //const double delta_un_sur_rho = 1. / rho_1 - 1. / rho_0;
    const double un_sur_rho_0 = 1. / rho_0;
    const double un_sur_rho_1 = 1. / rho_1;
    //const Zone_VF & zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
    //const DoubleTab & xp = zone_vf.xp();

    const int n = phi.dimension(0);
    for (int i = 0; i < n; i++)
      {
        double d = dist(i);
        double p = 0.;
        if (d >= -1e20)
          {
            const double div_n = phi(i);
            // Distance calculee pour cet element ?
            // Calcul de la fonction phi pour cet element :
            const double mp = mpoint[i];
            switch (ordre)
              {
              case 0:
                // Pas de prise en compte de la correction en courbure
                p = d * mp;
                break;
              case 1:
                //  Prise en compte de la correction en courbure a l'ordre 1
                p = d * (1. - 0.5 * div_n * d) * mp;
                break;
              case 2:
                //  Prise en compte de la correction en courbure a l'ordre 2
                p = d * (1. - 0.5 * div_n * d + div_n * div_n * d * d / 6.) * mp;
                break;
              default:
                Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface ordre" << finl;
                exit();
              }
            switch (phase_pilote)
              {
              case -1:
                // Champ de vitesse tel que u + u0 soit continu et egal a la
                // vitesse de deplacement de l'interface
                if (d < 0)
                  p *= un_sur_rho_0;
                else
                  p *= un_sur_rho_1;
                break;
              case 0:
                // Champ de vitesse tel que u + u0 soit continu et
                // u+u0 = la vitesse de la phase 0 dans la phase 0
                if (d < 0)
                  p = 0.;
                else
                  p *= (un_sur_rho_1 - un_sur_rho_0);
                break;
              case 1:
                // Champ de vitesse tel que u + u0 soit continu et
                // u+u0 = la vitesse de la phase 1 dans la phase 1
                if (d < 0)
                  p *= (un_sur_rho_0 - un_sur_rho_1);
                else
                  p = 0.;
                break;
              default:
                Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface phase_pilote" << finl;
                exit();
              }
          }
        phi(i) = p;
      }
    phi.echange_espace_virtuel();
  }

  // Gradient de phi:
  if (champ_u0.que_suis_je() == "Champ_Face")
    {

      gradient.calculer(phi, u0);

    }
  else
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface\n"
           << " Non code pour " << champ_u0.que_suis_je() << finl;
      exit();
    }

  // On annule la vitesse calculee pour les faces adjacentes a un element
  // invalide.
  {
    const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
    const IntTab&   face_voisins = zone_vf.face_voisins();
    const int nb_faces = zone_vf.nb_faces();
    for (int i = 0; i < nb_faces; i++)
      {
        for (int j = 0; j < 2; j++)
          {
            const int elem = face_voisins(i, j);
            if (elem >= 0 && dist(elem) < -1e20)
              {
                u0(i) = 0.;
                break;
              }
          }
      }
  }
  u0.echange_espace_virtuel();
  solveur_masse.appliquer(u0);
}

// Calcul de l'integrale de dI_dt sur chaque element du maillage.
// Le tableau dI_dt doit avoir la bonne structure. L'espace virtuel n'est pas
// mis a jour.
void Navier_Stokes_FT_Disc::calculer_dI_dt(DoubleVect& dI_dt) const
{
  const double rho_0 = fluide_diphasique().fluide_phase(0).masse_volumique().valeurs()(0,0);
  const double rho_1 = fluide_diphasique().fluide_phase(1).masse_volumique().valeurs()(0,0);
  const double delta_rho = rho_0-rho_1;

  double rho_0_delta_rho = 0.;
  if (delta_rho != 0)
    rho_0_delta_rho = rho_0 / delta_rho;

  const DoubleTab& tab_vitesse = inconnue().valeurs();
  const int vef = (tab_vitesse.nb_dim() == 2);
  const IntTab& face_voisins = zone_dis().valeur().face_voisins();
  const DoubleTab& indic = variables_internes().ref_eq_interf_proprietes_fluide.valeur().inconnue().valeur().valeurs();
  DoubleTab tmp(tab_vitesse); // copie du tableau des vitesses
  const int dim = tab_vitesse.line_size();
  const int n = tab_vitesse.dimension(0);
  for (int i = 0; i < n; i++)
    {
      const int elem0 = face_voisins(i, 0);
      const int elem1 = face_voisins(i, 1);
      double indic_0 = (elem0 >= 0) ? indic[elem0] : indic[elem1];
      double indic_1 = (elem1 >= 0) ? indic[elem1] : indic[elem0];
      // Decentrage : si la face est adjacente a une face monophasique,
      //  prendre la phase pure pour toute la face:
      if (indic_0 == 0. || indic_0 == 1.)
        indic_1 = indic_0;
      if (indic_1 == 0. || indic_1 == 1.)
        indic_0 = indic_1;
      const double indic_face = (indic_0 + indic_1) * 0.5;
      const double x = rho_0_delta_rho - indic_face;
      if (vef)
        for (int j = 0; j < dim; j++)
          tmp(i,j) *= x;
      else
        tmp(i) *= x;
    }

  // Question: il y a un assert_espace_virtuel_vect dans divergence.calculer,
  //  mais l'operateur n'a normalement pas besoin de l'espace virtuel !
  //  La ligne suivante devrait pouvoir etre retiree:
  tmp.echange_espace_virtuel();

  // On cree un tableau avec la meme structure que la pression
  DoubleTab resu;
  resu.copy(variables_internes().second_membre_projection.valeurs(), Array_base::NOCOPY_NOINIT);

  //On utilise un operateur de divergence temporaire et pas celui porte par l equation
  //pour ne pas modifier les flux_bords_ rempli au cours de Navier_Stokes_std::mettre_a_jour

  Operateur_Div div_tmp;
  div_tmp.associer_eqn(*this);
  div_tmp.typer();
  div_tmp.l_op_base().associer_eqn(*this);
  div_tmp->completer();
  div_tmp->calculer(tmp,resu);

  // Extraction des valeurs
  const int nb_elem = zone_dis().valeur().nb_elem();
  assert(nb_elem == dI_dt.size());
  if (!vef)
    {
      // Simple copie
      dI_dt.inject_array(resu, nb_elem);
    }
  else
    {
      // L'integrale de div sur l'element est dimension * la valeur aux elements
      //  renvoyee par l'operateur.
      // Les valeurs aux elements sont au debut du tableau resu.
      dI_dt.inject_array(resu, nb_elem);
      dI_dt *= dimension;
    }
}

// Description:
//  Calcul de la derivee en temps de la vitesse.
DoubleTab& Navier_Stokes_FT_Disc::derivee_en_temps_inco(DoubleTab& vpoint)
{
  // Preparation des champs utilises pour le calcul des derivees en temps
  // S'il n'y a pas d'equation de transport des interfaces entre phases fluides,
  // on ne recalcule pas les proprietes.
  {
    REF(Transport_Interfaces_FT_Disc) & refeq_transport =
      variables_internes().ref_eq_interf_proprietes_fluide;
    if (refeq_transport.non_nul())
      {
        FT_disc_calculer_champs_rho_mu_nu_dipha(zone_dis().valeur(),
                                                fluide_diphasique(),
                                                refeq_transport.valeur().inconnue().valeur().valeurs(),
                                                // (indicatrice)
                                                champ_rho_elem_.valeur().valeurs(),
                                                champ_nu_.valeur().valeurs(),
                                                champ_mu_.valeur().valeurs(),
                                                champ_rho_faces_.valeur().valeurs());
      }
    else
      {
        const Fluide_Incompressible& phase_0 = ref_cast(Fluide_Incompressible,milieu());
        const Zone_dis_base& zdis = zone_dis().valeur();
        FT_disc_calculer_champs_rho_mu_nu_mono(zdis,
                                               phase_0,
                                               champ_rho_elem_,
                                               champ_mu_,
                                               champ_nu_,
                                               champ_rho_faces_);
      }
  }

  vpoint = 0.;

  // =====================================================================
  // Methode de projection :
  // Premiere etape : calcul de u_etoile (tous les termes de N.S. sauf la pression)

  // Contribution des operateurs diffusion et convection :
  // Operateur de diffusion : valeurs discretes homogenes a
  //                                             / d             \    //
  //                INTEGRALE                    | -- (rho * v)  |    //
  //                (sur le volume de controle)  \ dt            /    //
  // B.M. 08/2004 : on envoie la vitesse "v" a l'operateur qui calcule
  //                div (mu * (grad(v)+tr(grad(v))))
  //                (on a associe "mu" a la "diffusivite" de l'operateur,
  //                 voir Navier_Stokes_FT_Disc::lire)
  terme_diffusif.calculer(la_vitesse.valeurs(),
                          variables_internes().terme_diffusion.valeur().valeurs());
  solveur_masse.appliquer(variables_internes().terme_diffusion.valeur().valeurs());

  // Termes sources : gravite et tension de surface,
  // valeurs discretes homogenes a
  //                                             / d             \    //
  //                INTEGRALE                    | -- (rho * v)  |    //
  //                (sur le volume de controle)  \ dt            /    //
  //DoubleTab terme_source_collisions(la_vitesse.valeurs());
  DoubleTab& terme_source_collisions=variables_internes().terme_source_collisions.valeur().valeurs() ;
  terme_source_collisions=0;
  double isCollision =0 ;// -----------------
  {
    // Si une equation de transport est associee aux proprietes du fluide,
    // on ajoute le terme de tension de surface.
    REF(Transport_Interfaces_FT_Disc) & refeq_transport =
      variables_internes().ref_eq_interf_proprietes_fluide;

    if (refeq_transport.non_nul())
      {
        const Champ_base& indicatrice = refeq_transport.valeur().get_update_indicatrice();
        Champ_base& gradient_i = variables_internes().gradient_indicatrice.valeur();
        // Note:
        // On appelle la version const de maillage_interface() (qui est publique) car
        // on passe par const Transport_Interfaces_FT_Disc :
        const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
        const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
        const DoubleTab& distance_interface_sommets =
          eq_transport.get_update_distance_interface_sommets();

        calculer_gradient_indicatrice(indicatrice,
                                      distance_interface_sommets,
                                      gradient_i);
        calculer_champ_forces_superficielles(maillage,
                                             gradient_i,
                                             variables_internes().potentiel_elements,
                                             variables_internes().potentiel_faces,
                                             variables_internes().terme_source_interfaces);

        // calcule du terme source du au collisions  //--------------------------------
        //        variables_internes().terme_source_collisions.valeurs() = 0;
        calculer_champ_forces_collisions2(indicatrice.valeurs(), terme_source_collisions, isCollision);
      }
    else
      {
        variables_internes().terme_source_interfaces.valeurs() = 0;
      }
  }
  solveur_masse.appliquer(variables_internes().terme_source_interfaces.valeur().valeurs());
  solveur_masse.appliquer(terme_source_collisions);




  // Autres termes sources (acceleration / repere mobile)
  //  Valeurs homogenes a
  //                                             / d       \          //
  //                INTEGRALE                    | -- (v)  |          //
  //                (sur le volume de controle)  \ dt      /          //
  // (voir "preparer_calcul", commentaire sources().associer_rho...)
  variables_internes().terme_source.valeur().valeurs() = 0.;
  les_sources.ajouter(variables_internes().terme_source.valeur().valeurs());
  solveur_masse.appliquer(variables_internes().terme_source.valeur().valeurs());

  // Operateur de convection : valeurs discretes homogenes a
  //                                             / d       \          //
  //                INTEGRALE                    | -- (v)  |          //
  //                (sur le volume de controle)  \ dt      /          //
  // B.M. 08/2004 : on transporte "v" et non "rho_v"...

  DoubleTab& terme_convection_valeurs = variables_internes().terme_convection.valeur().valeurs();
  if (schema_temps().diffusion_implicite())
    {
      terme_convection_valeurs=0;
      derivee_en_temps_conv(terme_convection_valeurs,la_vitesse.valeurs());
    }
  else
    {
      terme_convectif.calculer(la_vitesse.valeurs(),terme_convection_valeurs);
    }
  solveur_masse.appliquer(variables_internes().terme_convection.valeur().valeurs());

  // Ajout des differentes contributions a vpoint :
  const DoubleTab& tab_rho_faces = champ_rho_faces_.valeur().valeurs();
  const DoubleVect& volumes_entrelaces = ref_cast(Zone_VF, zone_dis().valeur()).volumes_entrelaces();
  const DoubleTab& tab_diffusion = variables_internes().terme_diffusion.valeur().valeurs();
  const DoubleTab& termes_sources_interf = variables_internes().terme_source_interfaces.valeur().valeurs();
  //const DoubleTab& termes_sources_collisions = variables_internes().terme_source_collisions.valeur().valeurs();
  const DoubleTab& termes_sources = variables_internes().terme_source.valeur().valeurs();
  const DoubleTab& tab_convection = variables_internes().terme_convection.valeur().valeurs();
  const int n = vpoint.dimension(0);
  const int nbdim1 = (vpoint.nb_dim() == 1);
  const int m = (nbdim1 ? 0 : vpoint.dimension(1));

  DoubleTab gravite_face(inconnue().valeurs());
  if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G
      && milieu().a_gravite())
    {
      ArrOfDouble g(dimension);
      // Pour l'instant : gravite uniforme g => phi(s) = - x scalaire g
      const DoubleTab& gravite = milieu().gravite().valeurs();
      if (gravite.nb_dim() != 2 || gravite.dimension(1) != dimension)
        {
          Cerr << "Error for calculer_champ_forces_superficielles\n";
          Cerr << " gravite.dimension(1) != Objet_U::dimension" << finl;
          Process::exit();
        }
      for (int j = 0; j < dimension; j++)
        g[j] =  gravite(0,j);

      // On multiplie par les volumes entrelaces et on applique ensuite le solveur masse
      //  (traitement special des CL de Dirichlet)
      if (nbdim1)
        {
          const IntVect& orientation = ref_cast(Zone_VDF, zone_dis().valeur()).orientation();
          for (int face=0; face<n; face++)
            gravite_face(face)=volumes_entrelaces(face)*g(orientation[face]);
        }
      else
        {
          for (int face=0; face<n; face++)
            for (int dim=0; dim<m; dim++)
              gravite_face(face,dim)=volumes_entrelaces(face)*g(dim);
        }
      solveur_masse.appliquer(gravite_face);
    }
  else
    {
      // Pas de gravite :
      gravite_face = 0.;
    }

  IntTab flag_gradP(n);
  IntTab coef_TSF(n);
  coef_TSF = 1;
  int flag_diff;
  DoubleTab& gradP = variables_internes().gradient_pression.valeurs();
  if (schema_temps().diffusion_implicite())
    {
      //on calcule gradP (pour le qdm)
      gradient.calculer(la_pression.valeur().valeurs(), gradP);
      solveur_masse.appliquer(gradP);
      // si variables_internes().is_penalized=1 alors variables_internes().is_explicite=0
      if( variables_internes().is_penalized )
        flag_gradP = 0 ; //(grad P raide)
      else
        flag_gradP = 1 ;
      // on ajoute pas la diffusion cela sera fait par Gradient_conjugue_diff_impl
      // sauf si !is_explicite (terme forcage vitesse implicite; resolution conjointe forcage diffusion)
      // car dans ce cas la vitesse a imposer a besoin d'une bonne approximation de vpoint
      if( !variables_internes().is_explicite )
        flag_diff = 1;
      else
        flag_diff = 0;
    }
  else
    {
      flag_gradP = 0;
      flag_diff = 1;
    }

  bool interf_vitesse_imposee_ok = false;
  int nb_eqs = variables_internes().ref_eq_interf_vitesse_imposee.size();
  int nb_eq_non_nul = 0;
  for (int k=0; k<nb_eqs; k++)
    {
      REF(Transport_Interfaces_FT_Disc) & refeq_transport =
        variables_internes().ref_eq_interf_vitesse_imposee[k];

      if (refeq_transport.non_nul()) nb_eq_non_nul +=1;
    }
  DoubleTab terme_mul;
  if ( nb_eq_non_nul == nb_eqs && nb_eqs != 0 )
    {
      interf_vitesse_imposee_ok = true;
      terme_mul.copy(champ_rho_faces_.valeur().valeurs(), Array_base::COPY_INIT);
      terme_mul = 0.;
    }

  REF(Transport_Interfaces_FT_Disc) & refeq_transport_2pha =
    variables_internes().ref_eq_interf_proprietes_fluide;
  if (refeq_transport_2pha.non_nul() && interf_vitesse_imposee_ok && variables_internes().is_penalized)
    {
      const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
      const Zone_VF& zone_vdf = ref_cast(Zone_VDF, zone_dis().valeur());
      const IntTab& face_voisins = zone_vf.face_voisins();
      const IntTab& elem_faces = zone_vf.elem_faces();
      const IntVect& orientation = zone_vdf.orientation();
      const int   nb_faces_elem = elem_faces.dimension(1);
      for (int k=0; k<nb_eqs; k++)
        {
          REF(Transport_Interfaces_FT_Disc) & refeq_transport =
            variables_internes().ref_eq_interf_vitesse_imposee[k];
          const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_compute_indicatrice_faces().valeurs();
          for (int i = 0; i < face_voisins.dimension(0) ; i++)
            {
              if (indicatrice_faces(i) > 0.)
                {
                  flag_gradP(i) = 0;
                  coef_TSF(i) = 0; //annulation local terme source interf
                  const int ori=orientation(i);
                  const int voisin0 = face_voisins(i,0);
                  if (voisin0 >= 0)
                    {
                      int face_visavi = elem_faces(voisin0, ori) + elem_faces(voisin0, ori+Objet_U::dimension) - i;
                      for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                        {
                          const int face = elem_faces(voisin0, i_face);
                          if (indicatrice_faces(face) == 0. && face == face_visavi)
                            {
                              flag_gradP(face) = 0;
                              coef_TSF(face) = 0; //annulation local terme source interf
                            }
                        }
                    }
                  const int voisin1 = face_voisins(i,1);
                  if (voisin1 >= 0)
                    {
                      int face_visavi = elem_faces(voisin1, ori) + elem_faces(voisin1, ori+Objet_U::dimension) - i;
                      for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                        {
                          const int face = elem_faces(voisin1, i_face);
                          if (indicatrice_faces(face) == 0. && face == face_visavi)
                            {
                              flag_gradP(face) = 0;
                              coef_TSF(face) = 0; //annulation local terme source interf
                            }
                        }
                    }
                }
            }
        }
    }

  // Ajout des differentes contributions a vpoint :
  for (int i = 0; i < n; i++)
    {
      const double rho_face = tab_rho_faces(i);
      if (nbdim1)
        {
          vpoint(i) = ( - flag_gradP(i) * gradP(i) + flag_diff * tab_diffusion(i) + coef_TSF(i) * termes_sources_interf(i) + terme_source_collisions(i) ) / rho_face
                      + tab_convection(i) + termes_sources(i) + gravite_face(i);

        }
      else
        {
          for (int j = 0; j < m; j++)
            vpoint(i, j) = ( - flag_gradP(i) * gradP(i,j) + flag_diff * tab_diffusion(i,j) + coef_TSF(i) * termes_sources_interf(i,j) ) / rho_face
                           + tab_convection(i,j) + termes_sources(i,j) + gravite_face(i,j);
        }
    }
  vpoint.echange_espace_virtuel();

  //  si terme forcage vitesse explicite => J'ai tout, je peux resoudre (Gradient_conjugue_diff_impl)
  if (schema_temps().diffusion_implicite() && variables_internes().is_explicite )
    {
      DoubleTab derivee(la_vitesse.valeurs());
      // on indique au solveur masse de diviser par rho en plus du volume car l'operateur de diffusion renvoit qqqc en rho d u/dt
      solveur_masse->set_name_of_coefficient_temporel(champ_rho_faces_.valeur().le_nom());

      DoubleTrav tt(vpoint);
      tt=vpoint;
      derivee=inconnue().valeurs();
      Equation_base::Gradient_conjugue_diff_impl( tt, derivee ) ;

      solveur_masse->set_name_of_coefficient_temporel("no_coeff");

      vpoint=derivee;
      // on retire le gradient si on ne penalise pas:
      for (int i = 0; i < n; i++)
        {
          const double rho_face = tab_rho_faces(i);
          if (nbdim1)
            vpoint(i) += gradP(i) / rho_face ;
          else
            for (int j = 0; j < m; j++)
              vpoint(i, j) += gradP(i,j) / rho_face ;
        }
    }

  const int nfaces = vpoint.dimension_tot(0);

  if( interf_vitesse_imposee_ok )
    {
      int compteur_vimp_regul =0;
      DoubleTrav vpoint0(vpoint) ;
      vpoint0 = vpoint ;
      terme_mul = 1.0;
      // S'il y a une equation de transport avec vitesse imposee, on impose:
      DoubleTrav forces_tot(vpoint);
      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k=0; k<nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) & refeq_transport =
            variables_internes().ref_eq_interf_vitesse_imposee[k];

          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

          const DoubleTab& inco_val = inconnue().valeur().valeurs();
          const DoubleTab& rho_faces = champ_rho_faces_.valeur().valeurs();
          DoubleTab& source_val = variables_internes().terme_source.valeur().valeurs();
          const double temps = schema_temps().temps_courant();
          const double dt = schema_temps().pas_de_temps();

          //On ajoute un terme source a vpoint pour imposer au fluide la vitesse de l interface
          //source_val est rempli (peut etre postraite)
          eq_transport.modifier_vpoint_pour_imposer_vit(inco_val,vpoint0,vpoint,rho_faces,source_val,temps,dt,variables_internes().is_explicite,variables_internes().eta);
          source_val.echange_espace_virtuel();
          forces_tot += source_val;

          // Afin de savoir s il existe des interfaces IBC/fluide regularisees.
          // Si oui, on ne penalise que pour indic=1.
          if (eq_transport.get_vimp_regul()) compteur_vimp_regul++;
          //calcul du terme 1 + somme Xs / eta
          if ( !variables_internes().is_explicite )
            {
              const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_indicatrice_faces().valeurs();
              for (int i = 0; i < nfaces; i++)
                {
                  if ((eq_transport.get_vimp_regul()==0 && indicatrice_faces(i) > 0.) ||
                      (eq_transport.get_vimp_regul()==1 && indicatrice_faces(i)== 1.)) terme_mul(i) += 1. / variables_internes().eta;
                }
            }
        }
      terme_mul.echange_espace_virtuel();
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco terme_mul:",terme_mul);

      if (schema_temps().diffusion_implicite())
        {
          // si !variables_internes().is_explicite (terme forcage implicite) on retire la diffusion explicite
          // de vpoint (implicitee dans Gradient_conjugue_diff_impl_IBC)
          if( !variables_internes().is_explicite )
            {
              const DoubleTab& rho_faces = champ_rho_faces_.valeur().valeurs();
              const DoubleTab& diffusion = variables_internes().terme_diffusion.valeur().valeurs();

              for (int i = 0; i < vpoint.dimension(0); i++)
                {
                  const double rho_face = rho_faces(i);
                  if (vpoint.nb_dim() == 1)
                    {
                      vpoint(i) -= diffusion(i) / rho_face ;
                    }
                  else
                    {
                      for (int j = 0; j < vpoint.dimension(1); j++)
                        {
                          vpoint(i,j) -= diffusion(i,j) / rho_face ;
                        }
                    }
                }
              vpoint.echange_espace_virtuel() ;
              DoubleTab derivee(la_vitesse.valeurs());
              // on indique au solveur masse de diviser par rho en plus du volume car l'operateur de diffusion renvoit qqqc en rho d u/dt
              {
                solveur_masse->set_name_of_coefficient_temporel(champ_rho_faces_.valeur().le_nom());

                DoubleTrav tt(vpoint);
                tt=vpoint;
                Equation_base::Gradient_conjugue_diff_impl( tt, derivee, terme_mul ) ;

                solveur_masse->set_name_of_coefficient_temporel("no_coeff");
              }
              vpoint=derivee;

              // on retire le gradient
              const int nbis = vpoint.dimension(0);
              const int nbdim1bis = (vpoint.nb_dim() == 1);
              const int mbis = (nbdim1bis ? 0 : vpoint.dimension(1));
              for (int i = 0; i < nbis; i++)
                {
                  const double rho_face = rho_faces(i);
                  if (nbdim1bis)
                    {
                      vpoint(i) += ( flag_gradP(i)*gradP(i) ) / rho_face ;
                    }
                  else
                    {
                      for (int j = 0; j < mbis; j++)
                        {
                          vpoint(i, j) += ( flag_gradP(i)*gradP(i,j)  ) / rho_face ;
                        }
                    }
                }
            }
          vpoint.echange_espace_virtuel() ;
        }
      else
        {
          // si on implicite le calcul de vpoint
          // On divise vpoint par le terme multiplicatif calcule avant
          if ( !variables_internes().is_explicite )
            {
              const int nbdim1bis = (vpoint.nb_dim() == 1);
              const int mbis = (nbdim1bis ? 0 : vpoint.dimension(1));
              // calcul de vpoint : vpoint / (1 + somme Xs/eta )
              for (int i = 0; i < nfaces; i++)
                {
                  if (nbdim1bis)
                    {
                      vpoint(i) /= terme_mul(i);
                    }
                  else
                    {
                      for (int j = 0; j < mbis; j++)
                        {
                          vpoint(i,j) /= terme_mul(i);
                        }
                    }
                }
              vpoint.echange_espace_virtuel();
            }

        }
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco vpoint:",vpoint);

      // Dans le cas penalize + vitesse imposee regularisee,
      // seuls les ddl tels que indic_faces=1 sont penalisees, les
      // autres sont forces avec un DF :
      if( !variables_internes().is_explicite && compteur_vimp_regul)
        {
          vpoint0=vpoint;
          // S'il y a une equation de transport avec vitesse imposee, on impose:
          DoubleTrav forces_totbis(vpoint);
          for (int k=0; k<nb_eqs_bis; k++)
            {
              REF(Transport_Interfaces_FT_Disc) & refeq_transport =
                variables_internes().ref_eq_interf_vitesse_imposee[k];

              Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

              const DoubleTab& inco_val = inconnue().valeur().valeurs();
              const DoubleTab& rho_faces = champ_rho_faces_.valeur().valeurs();
              DoubleTab& source_val = variables_internes().terme_source.valeur().valeurs();
              const double temps = schema_temps().temps_courant();
              const double dt = schema_temps().pas_de_temps();

              //On ajoute un terme source a vpoint pour imposer au fluide la vitesse de l interface
              //source_val est rempli (peut etre postraite)
              eq_transport.modifier_vpoint_pour_imposer_vit(inco_val,vpoint0,vpoint,rho_faces,source_val,temps,dt,/* is_explicite */ 1, /* eta */ 1.);
              source_val.echange_espace_virtuel();
              forces_totbis += source_val;
            }
        }
      vpoint.echange_espace_virtuel() ;

      //Calcul des efforts exerces par le fluide sur chaque interface
      //Attention valable si les ibc ne se chevauchent pas
      if (limpr())
        {
          const DoubleTab& tab_rho_facesbis = champ_rho_faces_.valeur().valeurs();
          for (int k=0; k<nb_eqs_bis; k++)
            {
              REF(Transport_Interfaces_FT_Disc) & refeq_transport =
                variables_internes().ref_eq_interf_vitesse_imposee[k];
              Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
              eq_transport.calcul_effort_fluide_interface(vpoint,tab_rho_facesbis,forces_tot,variables_internes().is_explicite,variables_internes().eta);
              forces_tot.echange_espace_virtuel();
            }
        }
    }

  // Assemblage de la matrice INTEGRALE            ( div(1/rho * grad(P)) ) :
  //                          (volume de controle pression)
  // Si l'option "matrice_pression_invariante" est activee, on ne recalcule
  // pas la matrice :
  if ( !interf_vitesse_imposee_ok )
    {
      if ( !variables_internes().matrice_pression_invariante )
        {
          assembleur_pression().valeur().assembler_rho_variable(matrice_pression_,
                                                                champ_rho_faces_.valeur());
          // On a modifie la matrice, il faut reinitialiser le solveur
          //  (recalcul des preconditionnement, factorisation pour Cholesky, etc...)
          solveur_pression().valeur().reinit();
        }
    }

  // ====================================================================
  // Methode de projection :
  // Deuxieme etape : projection du champ de vitesse sur le sous-espace
  //  a divergence nulle.
  // Trouver la pression P telle que
  //   d(u)/dt = vpoint - 1/rho * grad(P)
  //   div( d(u)/dt ) = 0
  // Soit :
  //   div(1/rho * grad(P)) = div(vpoint)

  // Calcul du second membre :
  //  div(vpoint) a l'interieur du domaine,
  //  prise en compte des conditions aux limites de pression/vitesse

  DoubleTab& secmem = variables_internes().second_membre_projection.valeurs();
  const double& dt = schema_temps().pas_de_temps();
  const DoubleTab& inco = inconnue().valeur().valeurs();
  // secmem = div(U/dt+vpoint) = div(U(n+1)/dt)
  DoubleTab du(inco);
  du /= dt;
  du += vpoint;
  divergence.calculer(du, secmem);
  secmem *= -1;

  // Prise en compte des sources de constituant (terme source dans une equation concentration)
  // Pour ne pas faire figurer la source deux fois, je la mets uniquement dans l'equation
  // de concentration... elle produit alors un terme source de div_u dans Navier_Stokes:
  const Noms& noms_eq = variables_internes().equations_concentration_source_fluide_;
  for (int i_eq = 0; i_eq < noms_eq.size(); i_eq++)
    {
      const Equation_base& eq = probleme().get_equation_by_name(noms_eq[i_eq]);
      for (int i_source = 0; i_source < eq.sources().size(); i_source++)
        {
          if (!sub_type(Terme_Source_Constituant_Vortex_VEF_Face, eq.sources()[i_source].valeur()))
            continue;

          const Terme_Source_Constituant_Vortex_VEF_Face& src = ref_cast(Terme_Source_Constituant_Vortex_VEF_Face, eq.sources()[i_source].valeur());
          src.ajouter_terme_div_u(secmem, schema_temps().pas_de_temps());
        }

    }
  secmem.echange_espace_virtuel();

  // Prise en compte du terme source div(u) du changement de phase
  if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
    {
      // GB2016 : Le calcul de mpoint ci-dessous me semble inutile car il est fait au debut de calculer_delta_u_interface:
      // GB2016 : Mais si je ne le fais pas, j'ai parfois : 'vx.get_md_vector() == md' dans calculer_delta_u_interface
      if (variables_internes().ref_equation_mpoint_.non_nul())
        variables_internes().ref_equation_mpoint_.valeur().calculer_mpoint(variables_internes().mpoint.valeur());
      if (variables_internes().ref_equation_mpoint_vap_.non_nul())
        variables_internes().ref_equation_mpoint_vap_.valeur().calculer_mpoint(variables_internes().mpoint_vap.valeur());

      calculer_delta_u_interface(variables_internes().delta_u_interface, -1 /* vitesse de l'interface */, variables_internes().correction_courbure_ordre_ /* ordre de la correction en courbure */);

      const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
      const IntTab& face_voisins = zone_vf.face_voisins();
      const IntTab& elem_faces = zone_vf.elem_faces();
      const int   nb_faces_elem = elem_faces.dimension(1);
      REF(Transport_Interfaces_FT_Disc) & refeq_transport =
        variables_internes().ref_eq_interf_proprietes_fluide;
      const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
      // Distance a l'interface discretisee aux elements:
      const DoubleTab& distance = eq_transport.get_update_distance_interface().valeurs();
      DoubleTab secmem2(secmem);
      divergence.calculer(variables_internes().delta_u_interface, secmem2);
      // On ne conserve que la divergence des elements traverses par l'interface
      const int nb_elem = secmem2.dimension(0);
      for (int elem = 0; elem < nb_elem; elem++)
        {
          const double dist = distance(elem);
          int i_face = -1;
          if (dist < -1e20)
            {
              // Distance invalide: on est loin de l'interface
              i_face = nb_faces_elem;
            }
          else
            {
              // Y a-t-il un voisin pour lequel la distance est de signe oppose
              for (i_face = 0; i_face < nb_faces_elem; i_face++)
                {
                  const int face = elem_faces(elem, i_face);
                  const int voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                  if (voisin >= 0)
                    {
                      const double d = distance(voisin);
                      if (d > -1e20 && d * dist < 0.)
                        break; // Changement de signe
                    }
                }
            }
          if (i_face == nb_faces_elem)
            {
              // Tous les voisins sont du meme cote de l'interface
              secmem2(elem) = 0.;
            }
        }
      secmem2 /= schema_temps().pas_de_temps();
      secmem += secmem2;
      secmem.echange_espace_virtuel();
    }

  Champ_Fonc champ_rho_faces_modifie(champ_rho_faces_);
  DoubleTab& rho_faces_modifie = champ_rho_faces_modifie.valeur().valeurs();

  if ( interf_vitesse_imposee_ok )
    {

      // On verifie s il existe des interfaces IBC/fluide regularisees.
      int compteur_vimp_regul =0;
      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k=0; k<nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) & refeq_transport =
            variables_internes().ref_eq_interf_vitesse_imposee[k];
          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();


          if (eq_transport.get_vimp_regul()==1)
            compteur_vimp_regul++;
        }

      // Creation d'un nouveau rho qui prend en compte les vitesses imposees dans
      // le terme de forcage
      // Si des interfaces IBC/fluide sont regularisees, la projection devient classique
      // (a ameliorer potentiellement).
      // Mais dans le cas ou on a une interface diphasique, on bloque toutes les vitesses impossees
      // en PDF (regularise ou non)
      int modif_rho_true = 0;
      if (variables_internes().is_penalized && (compteur_vimp_regul==0 || refeq_transport_2pha.non_nul()))
        {
          for (int i = 0; i < nfaces; i++)
            {
              if (compteur_vimp_regul!=0)
                {
                  for (int k = 0; k < nb_eqs ; k++)
                    {
                      REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
                      Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
                      const DoubleTab& indicatrice_faces = eq_transport.get_indicatrice_faces().valeurs();
                      if (eq_transport.get_vimp_regul()==1 &&
                          indicatrice_faces(i) > 0.0 && indicatrice_faces(i)!= 1.) terme_mul(i) += 1.0/variables_internes().eta;
                    }
                }
              rho_faces_modifie(i) *= terme_mul(i);
            }
          rho_faces_modifie.echange_espace_virtuel();
          modif_rho_true = 1;
        }
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco rho_faces_modifie:",rho_faces_modifie);

      // Assemblage de la matrice INTEGRALE            ( div(1/rho * grad(P)) ) :
      //                          (volume de controle pression)

      assembleur_pression().valeur().assembler_rho_variable(matrice_pression_,champ_rho_faces_modifie.valeur());

      //Penalization L2 de la pression si necessaire
      if (modif_rho_true == 1 && variables_internes().p_ref_pena != -1.e40)
        {
          // On se base sur le nombre de composantes par faces pour la discretisation
          Matrice_Morse_Sym& matrice_valeurs = (vpoint.nb_dim() == 1
                                                ? ref_cast(Matrice_Morse_Sym, (ref_cast(Matrice_Bloc, matrice_pression_.valeur())).get_bloc(0,0).valeur()) // VDF (1)
                                                : ref_cast(Matrice_Morse_Sym, matrice_pression_.valeur())) ;                                              // VEF (>1)
          DoubleTab& pressu = la_pression.valeurs();
          const int nb_elem = champ_rho_elem_.valeur().valeurs().dimension(0);
          const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
          const Zone& ma_zone = zone_dis().zone ();
          const IntTab& elem_faces = zone_vf.elem_faces();
          const int   nb_faces_elem = elem_faces.dimension(1);
          const int   nb_sommet = ma_zone.nb_som_elem();
          int numero_global_som, ligne_mat;
          int point_fluide_dirichlet=-1;
          if (variables_internes().is_pfl_flottant)
            {
              if( Objet_U::dimension == 3 )
                {
                  point_fluide_dirichlet = ma_zone.chercher_elements( variables_internes().x_pfl_imp , variables_internes().y_pfl_imp , variables_internes().z_pfl_imp );
                }
              else
                {
                  point_fluide_dirichlet = ma_zone.chercher_elements( variables_internes().x_pfl_imp , variables_internes().y_pfl_imp );
                }
              if (mp_max(point_fluide_dirichlet)==-1)
                {
                  Cerr << "Point de reference pression fluide situe en dehors du domaine !" << finl;
                  Process::exit();
                }
            }
          for (int e = 0; e < nb_elem; e++)
            {
              int nbfglob = 0;
              int nbfpena = 0;
              for ( int f = 0; f < nb_faces_elem; f++ )
                {
                  int fglob = elem_faces(e,f);
                  if (fglob >= 0)
                    {
                      nbfglob += 1;
                      int dejafait =0 ;
                      for (int k = 0; k < nb_eqs_bis && dejafait == 0 ; k++)
                        {
                          REF(Transport_Interfaces_FT_Disc) & refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
                          const DoubleTab& indicatrice_faces = refeq_transport.valeur().get_indicatrice_faces().valeurs();
                          if (indicatrice_faces(fglob) > 0.0)
                            {
                              dejafait = 1 ;
                              nbfpena += 1 ;
                            }
                        }
                    }
                }
              if (nbfpena == nbfglob && nbfpena != 0)
                {
                  matrice_valeurs(e,e) += 1.0 / variables_internes().eta;
                  secmem(e) +=  variables_internes().p_ref_pena / variables_internes().eta;
                  pressu(e) = variables_internes().p_ref_pena ;
                  if (vpoint.nb_dim() != 1) // VEF
                    {
                      for (int somloc = 0; somloc < nb_sommet; somloc ++)
                        {
                          numero_global_som = ma_zone.sommet_elem(e,somloc);
                          ligne_mat = nb_elem+numero_global_som;
                          //                         matrice_valeurs(ligne_mat,ligne_mat) += 1.0 / variables_internes().eta;
                          pressu(ligne_mat) = 0.0 ;
                        }
                    }
                }
              if (point_fluide_dirichlet == e)
                {
                  if (nbfpena != nbfglob && nbfpena != 0)
                    {
                      // On impose une reference de pression (p_ref)sur un bord
                      secmem(e) +=  matrice_valeurs(e,e) * variables_internes().p_ref_pena / (float(nbfglob-nbfpena));
                      matrice_valeurs(e,e) *= float(nbfglob-nbfpena+1)/(float(nbfglob-nbfpena));
                    }
                  else
                    {
                      Cerr<<"Point de reference pression fluide non situe dans une cellule fluide voisin d'une IBC !"<<finl;
                      Cerr<<"Nb faces IBC = "<<nbfpena<<finl;
                      Process::exit();
                    }
                }
            }
          secmem.echange_espace_virtuel();
          pressu.echange_espace_virtuel();
        }

      // On a modifie la matrice, il faut reinitialiser le solveur
      //  (recalcul des preconditionnement, factorisation pour Cholesky, etc...)
      solveur_pression().valeur().reinit();

    }

  assembleur_pression_.valeur().modifier_secmem(secmem);

  // Resolution du systeme en pression : calcul de la_pression
  solveur_pression_.resoudre_systeme(matrice_pression_.valeur(),
                                     secmem,
                                     la_pression.valeurs()
                                    );
  assembleur_pression_.modifier_solution(la_pression);
  // Calcul d(u)/dt = vpoint + 1/rho*grad(P)
  gradient.calculer(la_pression.valeur().valeurs(), gradP);
  solveur_masse.appliquer(gradP);

  // Correction de vpoint :
  if (projection_a_faire()) // Temporaire pour permettre de ne pas resoudre NS avec mettant operateurs nuls et projection_initiale 0
    {
      const int nbis = vpoint.dimension(0);
      if (vpoint.nb_dim() == 1)
        {
          for (int i = 0; i < nbis; i++)
            vpoint(i) -= gradP(i) / rho_faces_modifie(i);
        }
      else
        {
          const int mbis = vpoint.dimension(1);
          for (int i = 0; i < nbis; i++)
            {
              double x = rho_faces_modifie(i);
              for (int j = 0; j < mbis; j++)
                vpoint(i,j) -= gradP(i,j) / x;
            }
        }
      vpoint.echange_espace_virtuel();
    }

  // Calcul des efforts sur les IBCs et impression dans un fichier
  if (interf_vitesse_imposee_ok && limpr())
    {
      const DoubleTab& rho = champ_rho_faces_.valeur().valeurs();

      DoubleTrav forces_tot_2(vpoint) ;
      DoubleTrav pressure_part(vpoint) ;
      DoubleTrav diffusion_part(vpoint) ;

      //Calcul de la vitesse au temps n+1
      DoubleTab vv(vpoint) ;
      vv *= schema_temps().pas_de_temps() ;
      vv += inconnue().valeur().valeurs() ;

      // Terme de diffusion
      terme_diffusif.calculer(vv, variables_internes().terme_diffusion.valeur().valeurs());
      variables_internes().terme_diffusion.valeur().valeurs().echange_espace_virtuel() ;
      solveur_masse.appliquer(variables_internes().terme_diffusion.valeur().valeurs());
      const DoubleTab& diffusion = variables_internes().terme_diffusion.valeur().valeurs();
      // Terme de convection
      DoubleTrav trav(variables_internes().terme_convection.valeur().valeurs());
      derivee_en_temps_conv( trav, la_vitesse.valeurs());
      variables_internes().terme_convection.valeur().valeurs()=trav;
      solveur_masse.appliquer(variables_internes().terme_convection.valeur().valeurs());
      const DoubleTab& convection = variables_internes().terme_convection.valeur().valeurs();

      const int nbis = vpoint.dimension(0);
      const int nbdim1bis = (vpoint.nb_dim() == 1);
      const int mbis = (nbdim1bis ? 0 : vpoint.dimension(1));
      const Zone_dis_base& la_zone_disbis = zone_dis().valeur();
      int is_VDF=0;
      if (la_zone_disbis.que_suis_je().debute_par("Zone_VDF"))
        is_VDF=1;
      else if (la_zone_disbis.que_suis_je().debute_par("Zone_VEF"))
        is_VDF=0;
      else
        {
          Cerr<<"La zone discretisee "<<la_zone_disbis.que_suis_je()<<" n est pas reconnue"<<finl;
          exit();
        }
      for (int i = 0; i<nbis; i++)
        {
          if (is_VDF)
            {
              pressure_part(i) = gradP(i) ;
              diffusion_part(i) = -rho(i) * convection(i) - diffusion(i) ;
              forces_tot_2(i) = pressure_part(i) + diffusion_part(i) ;
            }
          else
            {
              //is_VEF
              for (int j = 0; j < mbis; j++)
                {
                  pressure_part(i,j) = gradP(i,j) ;
                  diffusion_part(i,j) = -rho(i) * convection(i,j) - diffusion(i,j) ;
                  forces_tot_2(i,j) = pressure_part(i,j) + diffusion_part(i,j) ;
                }
            }
        }
      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k=0; k<nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) & refeq_transport =
            variables_internes().ref_eq_interf_vitesse_imposee[k];
          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
          // Impression des efforts
          eq_transport.impr_effort_fluide_interface( forces_tot_2, pressure_part, diffusion_part ) ;
        }
    }

  return vpoint;
}


const Probleme_FT_Disc_gen& Navier_Stokes_FT_Disc::probleme_ft() const
{
  return probleme_ft_.valeur();
}

Probleme_FT_Disc_gen& Navier_Stokes_FT_Disc::probleme_ft()
{
  return probleme_ft_.valeur();
}

// Description:
// In Front Tracking, pression is in Pa and so pression_pa field <=> pression field
void Navier_Stokes_FT_Disc::calculer_la_pression_en_pa()
{
  la_pression_en_pa.valeurs()=la_pression.valeurs();
}

const Navier_Stokes_FT_Disc_interne& Navier_Stokes_FT_Disc::variables_internes() const
{
  return *variables_internes_;
}
Navier_Stokes_FT_Disc_interne& Navier_Stokes_FT_Disc::variables_internes()
{
  return *variables_internes_;
}

// Description: Si le champ de vitesse est discontinu (calcul avec changement de phase),
//  renvoie un pointeur vers le champ delta_v de "discontinuite", tel que
//  inconnue - delta_v = vitesse de deplacement des interfaces
//  (voir Transport_Interfaces_FT_Disc::deplacer_maillage_v_fluide())
//  Si pas de changement de phase, renvoie un pointeur nul.
const Champ_base *  Navier_Stokes_FT_Disc::get_delta_vitesse_interface() const
{
  if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
    return & (variables_internes().delta_u_interface.valeur());
  else
    return 0;
}
// Description : calcul de div(n) (la courbure discretisee sur le maillage volumique)
//  Faudrait deplacer cette methode dans transport interfaces...
const Champ_base& Navier_Stokes_FT_Disc::calculer_div_normale_interface()
{
  REF(Transport_Interfaces_FT_Disc) & refeq_transport =
    variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  // Distance a l'interface discretisee aux elements:
  const DoubleTab& dist = eq_transport.get_update_distance_interface().valeurs();

  DoubleTab& phi = variables_internes().laplacien_d.valeur().valeurs();
  DoubleTab u0(inconnue().valeurs());

  //  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_div_normale", 0);
  //  statistiques().begin_count(count);

  phi = dist;

  const int n = phi.dimension(0);
  for (int i = 0; i < n; i++)
    {
      if (phi(i) < -1e20)
        phi(i) = 0;
    }
  phi.echange_espace_virtuel();
  //  static const Stat_Counter_Id count2 = statistiques().new_counter(1, "calculer_gradient", 0);
  //  statistiques().begin_count(count2);
  gradient.calculer(phi, u0);
  //  statistiques().end_count(count2);
  //  static const Stat_Counter_Id count4 = statistiques().new_counter(1, "calculer_solveur_masse", 0);
  //  statistiques().begin_count(count4);
  solveur_masse.appliquer(u0);
  //  statistiques().end_count(count4);

  // Calcul de integrale(div(u0)) sur les mailles:
  //  static const Stat_Counter_Id count3 = statistiques().new_counter(1, "calculer_div", 0);
  //  statistiques().begin_count(count3);
  divergence.calculer(u0, phi);
  //  statistiques().end_count(count3);
  // Division par le volume des mailles:
  const DoubleVect& volumes = ref_cast(Zone_VF, zone_dis().valeur()).volumes();
  for (int i = 0; i < n; i++)
    {
      const double p = phi(i);
      if (p != 0.)
        {
          const double v = volumes[i];
          phi(i) = p / v;
        }
    }

  //  statistiques().end_count(count);

  return variables_internes().laplacien_d.valeur();
}

const Champ_Fonc& Navier_Stokes_FT_Disc::champ_rho_faces() const
{
  return champ_rho_faces_;
}

//Renvoie 1 si l option GRAVITE_RHO_G est activee 0 sinon
int Navier_Stokes_FT_Disc::is_terme_gravite_rhog() const
{
  if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G)
    return 1;
  else
    return 0;
}

const Champ_Fonc& Navier_Stokes_FT_Disc::get_num_compo() const
{
  return variables_internes().num_compo;
}

void Navier_Stokes_FT_Disc::reprendre_num_compo(Entree& is)
{
  variables_internes().num_compo.reprendre(is);
}
