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
// File:        Marching_Cubes_data.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////

static const int mcubes_def_aretes_vdf_2d[][2] =
{
  {0, 2},
  {0, 1},
  {1, 3},
  {2, 3}
};

// En 2D une face d'une element n'a qu'une arete, le tableau est trivial.
static const int mcubes_def_aretes_faces_vdf_2d[][2] =
{
  {0, 1} // Les aretes d'une face vont du sommet 0 de la face au sommet 1
};

// mcubes_nb_facettes : nombre de facettes a creer pour chaque cas
// et indice dans le tableau de description des facettes.
static const int mcubes_nb_facettes_vdf_2d[16] =
{
  // signe sommet: 3,2,1,0 (bit de poids fort pour le dernier sommet)
  0, //       0,0,0,0
  1, //       0,0,0,1
  1, //       0,0,1,0
  1, //       0,0,1,1
  1, //       0,1,0,0
  1, //       0,1,0,1
  2, //       0,1,1,0
  1,  //       0,1,1,1
  1,
  2,
  1,
  1,
  1,
  1,
  1,
  0
};

// Description des facettes a creer dans chaque cas.
// Les facettes sont orientees de sorte que les sommets de signe negatif
// (liquide) se trouvent a droite dans le sens de parcours.
// Dans les cas 8-15, il faut retourner les facettes.
//GET : inversion des orientations des facettes : la normale doit pointer vers les signes 1
static const int mcubes_facettes_vdf_2d[] =
{
  // cas 0 : pas de facette
  1,0, // cas 1: 1 facette entre les faces 0 et 1
  2,1, // cas 2: 1 facette entre les faces 1 et 2
  2,0, // 3
  0,3, // 4
  1,3, // 5
  0,3,2,1,  // 6
  2,3, // 7
  3,2, // 8
  3,0,1,2, // 9
  3,1, // 10
  3,0, // 11
  0,2, // 12
  1,2, // 13
  0,1, // 14
  // 15 : pas de facette
  -1   // Signature de fin pour test...
};

// Convention sur l'ordre des aretes (2D VEF)
//    sommets     aretes
//       2          *
//      / \        1 0
//     0---1      *-2-*

static const int mcubes_def_aretes_vef_2d[][2] =
{
  {1, 2},
  {2, 0},
  {0, 1}
};

// En 2D une face d'une element n'a qu'une arete, le tableau est trivial.
static const int mcubes_def_aretes_faces_vef_2d[][2] =
{
  {0, 1} // Les aretes d'une face vont du sommet 0 de la face au sommet 1
};

static const int mcubes_nb_facettes_vef_2d[8] =
{
  // signe sommet:   2,1,0 (bit de poids fort pour le dernier sommet)
  0,      //         0,0,0
  1,      //         0,0,1
  1,      //         0,1,0
  1,       //         0,1,1
  1,
  1,
  1,
  0
};

static const int mcubes_facettes_vef_2d[] =
{
  // cas 0 : pas de facette
  2,1, // cas 1: 1 facette entre les faces 0 et 1
  0,2, // cas 2: 1 facette entre les faces 1 et 2
  0,1,
  1,0,
  2,0,
  1,2,
  // cas 7 : pas de facette
  -1   // Signature de fin de tableau
};

// Convention sur l'ordre des aretes (3D VEF)
//    sommets          aretes
//       3-_               *-_
//      / \ -_            / \ 5_
//     /   \ .-2         3   4 .-*
//    / .   \ /         / . 2 \ 1
//   0-------1         *---0---*
//
static const int mcubes_def_aretes_vef_3d[][2] =
{
  {0, 1},
  {1, 2},
  {2, 0},
  {0, 3},
  {1, 3},
  {2, 3}
};

static const int mcubes_def_aretes_faces_vef_3d[][2] =
{
  {0, 1},
  {1, 2},
  {2, 0}
};

static const int mcubes_facettes_vef_3d[] =
{
  // cas 0 : pas de facette
  0,3,2,          //0001
  0,1,4,          //0010
  2,4,3, 1,4,2,   //0011
  1,2,5,          //0100
  0,3,5, 0,5,1,   //0101
  0,2,4, 2,5,4,   //0110
  3,5,4,          //0111
  3,4,5,          //1000
  0,4,2, 2,4,5,   //1001
  0,5,3, 0,1,5,   //1010
  1,5,2,          //1011
  2,3,4, 1,2,4,   //1100
  0,4,1,          //1101
  0,2,3,          //1110
  // cas 15: pas de facette
  -1   // Signature de fin de tableau
};

static const int mcubes_nb_facettes_vef_3d[16] =
{
  // signe sommet:   2,1,0 (bit de poids fort pour le dernier sommet)
  0,
  1,
  1,
  2,
  1,
  2,
  2,
  1,
  1,
  2,
  2,
  1,
  2,
  1,
  1,
  0
};

// Convention sur l'ordre des aretes (3D VDF)
//    sommets         aretes
//      6------7            *--3---*
//     /|     /|          10|    11|
//    2------3 |          *---1--* 7
//    | |    | |          | 6    | |
//    | 4----|-5          4 *--2-5-*
//    |/     |/           |8     |9
//    0------1            *---0--*

static const int mcubes_def_aretes_vdf_3d[][2] =
{
  { 0, 1 },
  { 2, 3 },
  { 4, 5 },
  { 6, 7 },
  { 0, 2 },
  { 1, 3 },
  { 4, 6 },
  { 5, 7 },
  { 0, 4 },
  { 1, 5 },
  { 2, 6 },
  { 3, 7 }
};

static const int mcubes_def_aretes_faces_vdf_3d[][2] =
{
  {0, 1},
  {2, 3},
  {0, 2},
  {1, 3}
};

static const int mcubes_facettes_vdf_3d[] =
{
  /*   0 */
  /*   1 */     0, 8, 4,
  /*   2 */     0, 5, 9,
  /*   3 */     5, 8, 4,    5, 9, 8,
  /*   4 */     4,10, 1,
  /*   5 */     0,10, 1,    0, 8,10,
  /*   6 */     0, 5, 9,    4,10, 1,
  /*   7 */     5,10, 1,    5, 8,10,    5, 9, 8,
  /*   8 */     5, 1,11,
  /*   9 */     5, 1,11,    0, 8, 4,
  /*  10 */     0,11, 9,    0, 1,11,
  /*  11 */     4, 9, 8,    4,11, 9,    4, 1,11,
  /*  12 */     5,10,11,    5, 4,10,
  /*  13 */     0,11, 5,    0,10,11,    0, 8,10,
  /*  14 */     0,11, 9,    0,10,11,    0, 4,10,
  /*  15 */     8,11, 9,    8,10,11,
  /*  16 */     2, 6, 8,
  /*  17 */     0, 6, 4,    0, 2, 6,
  /*  18 */     2, 6, 8,    0, 5, 9,
  /*  19 */     5, 6, 4,    5, 2, 6,    5, 9, 2,
  /*  20 */     2, 6, 8,    4,10, 1,
  /*  21 */     0,10, 1,    0, 6,10,    0, 2, 6,
  /*  22 */     2, 6, 8,    4,10, 1,    0, 5, 9,
  /*  23 */     5,10, 1,    5, 6,10,    5, 2, 6,    5, 9, 2,
  /*  24 */     2, 6, 8,    5, 1,11,
  /*  25 */     0, 6, 4,    0, 2, 6,    5, 1,11,
  /*  26 */     2, 6, 8,    0,11, 9,    0, 1,11,
  /*  27 */     4, 2, 6,    4, 9, 2,    4,11, 9,    4, 1,11,
  /*  28 */     2, 6, 8,    5,10,11,    5, 4,10,
  /*  29 */     0,11, 5,    0,10,11,    0, 6,10,    0, 2, 6,
  /*  30 */     2, 6, 8,    0,11, 9,    0,10,11,    0, 4,10,
  /*  31 */     2,11, 9,    2,10,11,    2, 6,10,
  /*  32 */     2, 9, 7,
  /*  33 */     2, 9, 7,    0, 8, 4,
  /*  34 */     0, 7, 2,    0, 5, 7,
  /*  35 */     5, 8, 4,    5, 2, 8,    5, 7, 2,
  /*  36 */     4,10, 1,    2, 9, 7,
  /*  37 */     0,10, 1,    0, 8,10,    2, 9, 7,
  /*  38 */     0, 7, 2,    0, 5, 7,    4,10, 1,
  /*  39 */     5,10, 1,    5, 8,10,    5, 2, 8,    5, 7, 2,
  /*  40 */     2, 9, 7,    5, 1,11,
  /*  41 */     2, 9, 7,    5, 1,11,    0, 8, 4,
  /*  42 */     0, 7, 2,    0,11, 7,    0, 1,11,
  /*  43 */     4, 2, 8,    4, 7, 2,    4,11, 7,    4, 1,11,
  /*  44 */     2, 9, 7,    5,10,11,    5, 4,10,
  /*  45 */     2, 9, 7,    0,11, 5,    0,10,11,    0, 8,10,
  /*  46 */     0, 7, 2,    0,11, 7,    0,10,11,    0, 4,10,
  /*  47 */     2,11, 7,    2,10,11,    2, 8,10,
  /*  48 */     7, 8, 9,    7, 6, 8,
  /*  49 */     0, 6, 4,    0, 7, 6,    0, 9, 7,
  /*  50 */     0, 6, 8,    0, 7, 6,    0, 5, 7,
  /*  51 */     5, 6, 4,    5, 7, 6,
  /*  52 */     7, 8, 9,    7, 6, 8,    4,10, 1,
  /*  53 */     0,10, 1,    0, 6,10,    0, 7, 6,    0, 9, 7,
  /*  54 */     4,10, 1,    0, 6, 8,    0, 7, 6,    0, 5, 7,
  /*  55 */     5,10, 1,    5, 6,10,    5, 7, 6,
  /*  56 */     5, 1,11,    7, 8, 9,    7, 6, 8,
  /*  57 */     5, 1,11,    0, 6, 4,    0, 7, 6,    0, 9, 7,
  /*  58 */     0, 6, 8,    0, 7, 6,    0,11, 7,    0, 1,11,
  /*  59 */     4, 7, 6,    4,11, 7,    4, 1,11,
  /*  60 */     7, 8, 9,    7, 6, 8,    5,10,11,    5, 4,10,
  /*  61 */     0,10, 5,    0, 6,10,    0, 9, 6,    9, 7, 6,    5,10,11,
  /*  62 */     0, 7, 8,    0,11, 7,    0, 4,11,    4,10,11,    8, 7, 6,
  /*  63 */     7,10,11,    7, 6,10,
  /*  64 */     6, 3,10,
  /*  65 */     6, 3,10,    0, 8, 4,
  /*  66 */     6, 3,10,    0, 5, 9,
  /*  67 */     6, 3,10,    5, 8, 4,    5, 9, 8,
  /*  68 */     4, 3, 1,    4, 6, 3,
  /*  69 */     0, 3, 1,    0, 6, 3,    0, 8, 6,
  /*  70 */     0, 5, 9,    4, 3, 1,    4, 6, 3,
  /*  71 */     5, 3, 1,    5, 6, 3,    5, 8, 6,    5, 9, 8,
  /*  72 */     6, 3,10,    5, 1,11,
  /*  73 */     6, 3,10,    5, 1,11,    0, 8, 4,
  /*  74 */     6, 3,10,    0,11, 9,    0, 1,11,
  /*  75 */     6, 3,10,    4, 9, 8,    4,11, 9,    4, 1,11,
  /*  76 */     5, 3,11,    5, 6, 3,    5, 4, 6,
  /*  77 */     0,11, 5,    0, 3,11,    0, 6, 3,    0, 8, 6,
  /*  78 */     0,11, 9,    0, 3,11,    0, 6, 3,    0, 4, 6,
  /*  79 */     3, 8, 6,    3, 9, 8,    3,11, 9,
  /*  80 */     2,10, 8,    2, 3,10,
  /*  81 */     0,10, 4,    0, 3,10,    0, 2, 3,
  /*  82 */     2,10, 8,    2, 3,10,    0, 5, 9,
  /*  83 */     5,10, 4,    5, 3,10,    5, 2, 3,    5, 9, 2,
  /*  84 */     4, 3, 1,    4, 2, 3,    4, 8, 2,
  /*  85 */     0, 3, 1,    0, 2, 3,
  /*  86 */     4, 3, 1,    4, 2, 3,    4, 8, 2,    0, 5, 9,
  /*  87 */     5, 3, 1,    5, 2, 3,    5, 9, 2,
  /*  88 */     2,10, 8,    2, 3,10,    5, 1,11,
  /*  89 */     0,10, 4,    0, 3,10,    0, 2, 3,    5, 1,11,
  /*  90 */     2,10, 8,    2, 3,10,    0,11, 9,    0, 1,11,
  /*  91 */     4, 2,10,    4, 9, 2,    4, 1, 9,    1,11, 9,   10, 2, 3,
  /*  92 */     5, 3,11,    5, 2, 3,    5, 8, 2,    5, 4, 8,
  /*  93 */     0,11, 5,    0, 3,11,    0, 2, 3,
  /*  94 */     4,11, 0,    4, 3,11,    4, 8, 3,    8, 2, 3,    0,11, 9,
  /*  95 */     2,11, 9,    2, 3,11,
  /*  96 */     6, 3,10,    2, 9, 7,
  /*  97 */     6, 3,10,    2, 9, 7,    0, 8, 4,
  /*  98 */     6, 3,10,    0, 7, 2,    0, 5, 7,
  /*  99 */     6, 3,10,    5, 8, 4,    5, 2, 8,    5, 7, 2,
  /* 100 */     4, 3, 1,    4, 6, 3,    2, 9, 7,
  /* 101 */     0, 3, 1,    0, 6, 3,    0, 8, 6,    2, 9, 7,
  /* 102 */     0, 7, 2,    0, 5, 7,    4, 3, 1,    4, 6, 3,
  /* 103 */     8, 5, 2,    8, 1, 5,    8, 6, 1,    6, 3, 1,    2, 5, 7,
  /* 104 */     6, 3,10,    2, 9, 7,    5, 1,11,
  /* 105 */     6, 3,10,    2, 9, 7,    5, 1,11,    0, 8, 4,
  /* 106 */     6, 3,10,    0, 7, 2,    0,11, 7,    0, 1,11,
  /* 107 */     6, 3,10,    4, 2, 8,    4, 7, 2,    4,11, 7,    4, 1,11,
  /* 108 */     2, 9, 7,    5, 3,11,    5, 6, 3,    5, 4, 6,
  /* 109 */     2, 9, 7,    0,11, 5,    0, 3,11,    0, 6, 3,    0, 8, 6,
  /* 110 */    11, 4, 3,   11, 0, 4,   11, 7, 0,    7, 2, 0,    3, 4, 6,
  /* 111 */    11, 6, 3,   11, 8, 6,   11, 2, 8,   11, 7, 2,
  /* 112 */     7, 8, 9,    7,10, 8,    7, 3,10,
  /* 113 */     0,10, 4,    0, 3,10,    0, 7, 3,    0, 9, 7,
  /* 114 */     0,10, 8,    0, 3,10,    0, 7, 3,    0, 5, 7,
  /* 115 */     5,10, 4,    5, 3,10,    5, 7, 3,
  /* 116 */     4, 3, 1,    4, 7, 3,    4, 9, 7,    4, 8, 9,
  /* 117 */     0, 3, 1,    0, 7, 3,    0, 9, 7,
  /* 118 */     8, 3, 4,    8, 7, 3,    8, 0, 7,    0, 5, 7,    4, 3, 1,
  /* 119 */     5, 3, 1,    5, 7, 3,
  /* 120 */     5, 1,11,    7, 8, 9,    7,10, 8,    7, 3,10,
  /* 121 */     5, 1,11,    0,10, 4,    0, 3,10,    0, 7, 3,    0, 9, 7,
  /* 122 */     7, 0,11,    7, 8, 0,    7, 3, 8,    3,10, 8,   11, 0, 1,
  /* 123 */     4, 3,10,    4, 7, 3,    4,11, 7,    4, 1,11,
  /* 124 */     3, 8, 7,    3, 4, 8,    3,11, 4,   11, 5, 4,    7, 8, 9,
  /* 125 */     0,11, 5,    0, 3,11,    0, 7, 3,    0, 9, 7,
  /* 126 */     7, 3,11,    0, 4, 8,
  /* 127 */     7, 3,11,
  /* 128 */     7,11, 3,
  /* 129 */     7,11, 3,    0, 8, 4,
  /* 130 */     7,11, 3,    0, 5, 9,
  /* 131 */     7,11, 3,    5, 8, 4,    5, 9, 8,
  /* 132 */     7,11, 3,    4,10, 1,
  /* 133 */     7,11, 3,    0,10, 1,    0, 8,10,
  /* 134 */     7,11, 3,    0, 5, 9,    4,10, 1,
  /* 135 */     7,11, 3,    5,10, 1,    5, 8,10,    5, 9, 8,
  /* 136 */     5, 3, 7,    5, 1, 3,
  /* 137 */     0, 8, 4,    5, 3, 7,    5, 1, 3,
  /* 138 */     0, 7, 9,    0, 3, 7,    0, 1, 3,
  /* 139 */     4, 9, 8,    4, 7, 9,    4, 3, 7,    4, 1, 3,
  /* 140 */     5, 3, 7,    5,10, 3,    5, 4,10,
  /* 141 */     0, 7, 5,    0, 3, 7,    0,10, 3,    0, 8,10,
  /* 142 */     0, 7, 9,    0, 3, 7,    0,10, 3,    0, 4,10,
  /* 143 */     7,10, 3,    7, 8,10,    7, 9, 8,
  /* 144 */     7,11, 3,    2, 6, 8,
  /* 145 */     7,11, 3,    0, 6, 4,    0, 2, 6,
  /* 146 */     7,11, 3,    2, 6, 8,    0, 5, 9,
  /* 147 */     7,11, 3,    5, 6, 4,    5, 2, 6,    5, 9, 2,
  /* 148 */     7,11, 3,    2, 6, 8,    4,10, 1,
  /* 149 */     7,11, 3,    0,10, 1,    0, 6,10,    0, 2, 6,
  /* 150 */     7,11, 3,    2, 6, 8,    4,10, 1,    0, 5, 9,
  /* 151 */     7,11, 3,    5,10, 1,    5, 6,10,    5, 2, 6,    5, 9, 2,
  /* 152 */     2, 6, 8,    5, 3, 7,    5, 1, 3,
  /* 153 */     0, 6, 4,    0, 2, 6,    5, 3, 7,    5, 1, 3,
  /* 154 */     2, 6, 8,    0, 7, 9,    0, 3, 7,    0, 1, 3,
  /* 155 */     9, 1, 7,    9, 4, 1,    9, 2, 4,    2, 6, 4,    7, 1, 3,
  /* 156 */     2, 6, 8,    5, 3, 7,    5,10, 3,    5, 4,10,
  /* 157 */    10, 0, 6,   10, 5, 0,   10, 3, 5,    3, 7, 5,    6, 0, 2,
  /* 158 */     2, 6, 8,    0, 7, 9,    0, 3, 7,    0,10, 3,    0, 4,10,
  /* 159 */     9, 3, 7,    9,10, 3,    9, 6,10,    9, 2, 6,
  /* 160 */     2,11, 3,    2, 9,11,
  /* 161 */     2,11, 3,    2, 9,11,    0, 8, 4,
  /* 162 */     0, 3, 2,    0,11, 3,    0, 5,11,
  /* 163 */     5, 8, 4,    5, 2, 8,    5, 3, 2,    5,11, 3,
  /* 164 */     2,11, 3,    2, 9,11,    4,10, 1,
  /* 165 */     2,11, 3,    2, 9,11,    0,10, 1,    0, 8,10,
  /* 166 */     4,10, 1,    0, 3, 2,    0,11, 3,    0, 5,11,
  /* 167 */     5, 8, 1,    5, 2, 8,    5,11, 2,   11, 3, 2,    1, 8,10,
  /* 168 */     5, 2, 9,    5, 3, 2,    5, 1, 3,
  /* 169 */     5, 2, 9,    5, 3, 2,    5, 1, 3,    0, 8, 4,
  /* 170 */     0, 3, 2,    0, 1, 3,
  /* 171 */     4, 2, 8,    4, 3, 2,    4, 1, 3,
  /* 172 */     5, 2, 9,    5, 3, 2,    5,10, 3,    5, 4,10,
  /* 173 */     5, 3, 9,    5,10, 3,    5, 0,10,    0, 8,10,    9, 3, 2,
  /* 174 */     0, 3, 2,    0,10, 3,    0, 4,10,
  /* 175 */     2,10, 3,    2, 8,10,
  /* 176 */     6,11, 3,    6, 9,11,    6, 8, 9,
  /* 177 */     0, 6, 4,    0, 3, 6,    0,11, 3,    0, 9,11,
  /* 178 */     0, 6, 8,    0, 3, 6,    0,11, 3,    0, 5,11,
  /* 179 */     5, 6, 4,    5, 3, 6,    5,11, 3,
  /* 180 */     4,10, 1,    6,11, 3,    6, 9,11,    6, 8, 9,
  /* 181 */     6, 9, 3,    6, 0, 9,    6,10, 0,   10, 1, 0,    3, 9,11,
  /* 182 */     4,10, 1,    0, 6, 8,    0, 3, 6,    0,11, 3,    0, 5,11,
  /* 183 */     5,10, 1,    5, 6,10,    5, 3, 6,    5,11, 3,
  /* 184 */     5, 8, 9,    5, 6, 8,    5, 3, 6,    5, 1, 3,
  /* 185 */     9, 6, 0,    9, 3, 6,    9, 5, 3,    5, 1, 3,    0, 6, 4,
  /* 186 */     0, 6, 8,    0, 3, 6,    0, 1, 3,
  /* 187 */     4, 3, 6,    4, 1, 3,
  /* 188 */     3, 5,10,    3, 9, 5,    3, 6, 9,    6, 8, 9,   10, 5, 4,
  /* 189 */     6,10, 3,    0, 9, 5,
  /* 190 */     0, 6, 8,    0, 3, 6,    0,10, 3,    0, 4,10,
  /* 191 */     6,10, 3,
  /* 192 */     7,10, 6,    7,11,10,
  /* 193 */     7,10, 6,    7,11,10,    0, 8, 4,
  /* 194 */     7,10, 6,    7,11,10,    0, 5, 9,
  /* 195 */     7,10, 6,    7,11,10,    5, 8, 4,    5, 9, 8,
  /* 196 */     4,11, 1,    4, 7,11,    4, 6, 7,
  /* 197 */     0,11, 1,    0, 7,11,    0, 6, 7,    0, 8, 6,
  /* 198 */     0, 5, 9,    4,11, 1,    4, 7,11,    4, 6, 7,
  /* 199 */     1, 6,11,    1, 8, 6,    1, 5, 8,    5, 9, 8,   11, 6, 7,
  /* 200 */     5, 6, 7,    5,10, 6,    5, 1,10,
  /* 201 */     0, 8, 4,    7, 5, 1,    7, 1,10,    7,10, 6,
  /*  Le cas 201 a ete corrige a la main par rapport a OpenDx 4.2 */
  /*  car il produit une surface non fermee */
  /*  Donnees d'origine issues d'opendx : */
  /*          0,7,5,      0,6,7,      0,8,6,      4,1,10, */
  /* 202 */     0, 7, 9,    0, 6, 7,    0,10, 6,    0, 1,10,
  /* 203 */     1, 9, 4,    1, 7, 9,    1,10, 7,   10, 6, 7,    4, 9, 8,
  /* 204 */     5, 6, 7,    5, 4, 6,
  /* 205 */     0, 7, 5,    0, 6, 7,    0, 8, 6,
  /* 206 */     0, 7, 9,    0, 6, 7,    0, 4, 6,
  /* 207 */     7, 8, 6,    7, 9, 8,
  /* 208 */     2,10, 8,    2,11,10,    2, 7,11,
  /* 209 */     0,10, 4,    0,11,10,    0, 7,11,    0, 2, 7,
  /* 210 */     0, 5, 9,    2,10, 8,    2,11,10,    2, 7,11,
  /* 211 */     2, 4, 9,    2,10, 4,    2, 7,10,    7,11,10,    9, 4, 5,
  /* 212 */     4,11, 1,    4, 7,11,    4, 2, 7,    4, 8, 2,
  /* 213 */     0,11, 1,    0, 7,11,    0, 2, 7,
  /* 214 */     0, 5, 9,    4,11, 1,    4, 7,11,    4, 2, 7,    4, 8, 2,
  /* 215 */     1, 7,11,    1, 2, 7,    1, 9, 2,    1, 5, 9,
  /* 216 */     5, 2, 7,    5, 8, 2,    5,10, 8,    5, 1,10,
  /* 217 */    10, 7, 1,   10, 2, 7,   10, 4, 2,    4, 0, 2,    1, 7, 5,
  /* 218 */     7,10, 2,    7, 1,10,    7, 9, 1,    9, 0, 1,    2,10, 8,
  /* 219 */     2, 7, 9,    4, 1,10,
  /* 220 */     5, 2, 7,    5, 8, 2,    5, 4, 8,
  /* 221 */     0, 7, 5,    0, 2, 7,
  /* 222 */     7, 8, 2,    7, 4, 8,    7, 0, 4,    7, 9, 0,
  /* 223 */     2, 7, 9,
  /* 224 */     2,10, 6,    2,11,10,    2, 9,11,
  /* 225 */     0, 8, 4,    2,10, 6,    2,11,10,    2, 9,11,
  /* 226 */     0, 6, 2,    0,10, 6,    0,11,10,    0, 5,11,
  /* 227 */     2,11, 6,    2, 5,11,    2, 8, 5,    8, 4, 5,    6,11,10,
  /* 228 */     4,11, 1,    4, 9,11,    4, 2, 9,    4, 6, 2,
  /* 229 */     6, 1, 8,    6,11, 1,    6, 2,11,    2, 9,11,    8, 1, 0,
  /* 230 */    11, 2, 5,   11, 6, 2,   11, 1, 6,    1, 4, 6,    5, 2, 0,
  /* 231 */     5,11, 1,    2, 8, 6,
  /* 232 */     5, 2, 9,    5, 6, 2,    5,10, 6,    5, 1,10,
  /* 233 */     0, 8, 4,    5, 2, 9,    5, 6, 2,    5,10, 6,    5, 1,10,
  /* 234 */     0, 6, 2,    0,10, 6,    0, 1,10,
  /* 235 */     2,10, 6,    2, 1,10,    2, 4, 1,    2, 8, 4,
  /* 236 */     5, 2, 9,    5, 6, 2,    5, 4, 6,
  /* 237 */     5, 2, 9,    5, 6, 2,    5, 8, 6,    5, 0, 8,
  /* 238 */     0, 6, 2,    0, 4, 6,
  /* 239 */     2, 8, 6,
  /* 240 */     8,11,10,    8, 9,11,
  /* 241 */     0,10, 4,    0,11,10,    0, 9,11,
  /* 242 */     0,10, 8,    0,11,10,    0, 5,11,
  /* 243 */     5,10, 4,    5,11,10,
  /* 244 */     4,11, 1,    4, 9,11,    4, 8, 9,
  /* 245 */     0,11, 1,    0, 9,11,
  /* 246 */     8, 1, 4,    8,11, 1,    8, 5,11,    8, 0, 5,
  /* 247 */     5,11, 1,
  /* 248 */     5, 8, 9,    5,10, 8,    5, 1,10,
  /* 249 */    10, 5, 1,   10, 9, 5,   10, 0, 9,   10, 4, 0,
  /* 250 */     0,10, 8,    0, 1,10,
  /* 251 */     4, 1,10,
  /* 252 */     5, 8, 9,    5, 4, 8,
  /* 253 */     0, 9, 5,
  /* 254 */     0, 4, 8,
  /* 255 */
  -1   // Signature de fin de tableau
};

static const int mcubes_nb_facettes_vdf_3d[] =
{
  0, /*   0 */
  1, /*   1 */
  1, /*   2 */
  2, /*   3 */
  1, /*   4 */
  2, /*   5 */
  2, /*   6 */
  3, /*   7 */
  1, /*   8 */
  2, /*   9 */
  2, /*  10 */
  3, /*  11 */
  2, /*  12 */
  3, /*  13 */
  3, /*  14 */
  2, /*  15 */
  1, /*  16 */
  2, /*  17 */
  2, /*  18 */
  3, /*  19 */
  2, /*  20 */
  3, /*  21 */
  3, /*  22 */
  4, /*  23 */
  2, /*  24 */
  3, /*  25 */
  3, /*  26 */
  4, /*  27 */
  3, /*  28 */
  4, /*  29 */
  4, /*  30 */
  3, /*  31 */
  1, /*  32 */
  2, /*  33 */
  2, /*  34 */
  3, /*  35 */
  2, /*  36 */
  3, /*  37 */
  3, /*  38 */
  4, /*  39 */
  2, /*  40 */
  3, /*  41 */
  3, /*  42 */
  4, /*  43 */
  3, /*  44 */
  4, /*  45 */
  4, /*  46 */
  3, /*  47 */
  2, /*  48 */
  3, /*  49 */
  3, /*  50 */
  2, /*  51 */
  3, /*  52 */
  4, /*  53 */
  4, /*  54 */
  3, /*  55 */
  3, /*  56 */
  4, /*  57 */
  4, /*  58 */
  3, /*  59 */
  4, /*  60 */
  5, /*  61 */
  5, /*  62 */
  2, /*  63 */
  1, /*  64 */
  2, /*  65 */
  2, /*  66 */
  3, /*  67 */
  2, /*  68 */
  3, /*  69 */
  3, /*  70 */
  4, /*  71 */
  2, /*  72 */
  3, /*  73 */
  3, /*  74 */
  4, /*  75 */
  3, /*  76 */
  4, /*  77 */
  4, /*  78 */
  3, /*  79 */
  2, /*  80 */
  3, /*  81 */
  3, /*  82 */
  4, /*  83 */
  3, /*  84 */
  2, /*  85 */
  4, /*  86 */
  3, /*  87 */
  3, /*  88 */
  4, /*  89 */
  4, /*  90 */
  5, /*  91 */
  4, /*  92 */
  3, /*  93 */
  5, /*  94 */
  2, /*  95 */
  2, /*  96 */
  3, /*  97 */
  3, /*  98 */
  4, /*  99 */
  3, /* 100 */
  4, /* 101 */
  4, /* 102 */
  5, /* 103 */
  3, /* 104 */
  4, /* 105 */
  4, /* 106 */
  5, /* 107 */
  4, /* 108 */
  5, /* 109 */
  5, /* 110 */
  4, /* 111 */
  3, /* 112 */
  4, /* 113 */
  4, /* 114 */
  3, /* 115 */
  4, /* 116 */
  3, /* 117 */
  5, /* 118 */
  2, /* 119 */
  4, /* 120 */
  5, /* 121 */
  5, /* 122 */
  4, /* 123 */
  5, /* 124 */
  4, /* 125 */
  2, /* 126 */
  1, /* 127 */
  1, /* 128 */
  2, /* 129 */
  2, /* 130 */
  3, /* 131 */
  2, /* 132 */
  3, /* 133 */
  3, /* 134 */
  4, /* 135 */
  2, /* 136 */
  3, /* 137 */
  3, /* 138 */
  4, /* 139 */
  3, /* 140 */
  4, /* 141 */
  4, /* 142 */
  3, /* 143 */
  2, /* 144 */
  3, /* 145 */
  3, /* 146 */
  4, /* 147 */
  3, /* 148 */
  4, /* 149 */
  4, /* 150 */
  5, /* 151 */
  3, /* 152 */
  4, /* 153 */
  4, /* 154 */
  5, /* 155 */
  4, /* 156 */
  5, /* 157 */
  5, /* 158 */
  4, /* 159 */
  2, /* 160 */
  3, /* 161 */
  3, /* 162 */
  4, /* 163 */
  3, /* 164 */
  4, /* 165 */
  4, /* 166 */
  5, /* 167 */
  3, /* 168 */
  4, /* 169 */
  2, /* 170 */
  3, /* 171 */
  4, /* 172 */
  5, /* 173 */
  3, /* 174 */
  2, /* 175 */
  3, /* 176 */
  4, /* 177 */
  4, /* 178 */
  3, /* 179 */
  4, /* 180 */
  5, /* 181 */
  5, /* 182 */
  4, /* 183 */
  4, /* 184 */
  5, /* 185 */
  3, /* 186 */
  2, /* 187 */
  5, /* 188 */
  2, /* 189 */
  4, /* 190 */
  1, /* 191 */
  2, /* 192 */
  3, /* 193 */
  3, /* 194 */
  4, /* 195 */
  3, /* 196 */
  4, /* 197 */
  4, /* 198 */
  5, /* 199 */
  3, /* 200 */
  4, /* 201 */
  4, /* 202 */
  5, /* 203 */
  2, /* 204 */
  3, /* 205 */
  3, /* 206 */
  2, /* 207 */
  3, /* 208 */
  4, /* 209 */
  4, /* 210 */
  5, /* 211 */
  4, /* 212 */
  3, /* 213 */
  5, /* 214 */
  4, /* 215 */
  4, /* 216 */
  5, /* 217 */
  5, /* 218 */
  2, /* 219 */
  3, /* 220 */
  2, /* 221 */
  4, /* 222 */
  1, /* 223 */
  3, /* 224 */
  4, /* 225 */
  4, /* 226 */
  5, /* 227 */
  4, /* 228 */
  5, /* 229 */
  5, /* 230 */
  2, /* 231 */
  4, /* 232 */
  5, /* 233 */
  3, /* 234 */
  4, /* 235 */
  3, /* 236 */
  4, /* 237 */
  2, /* 238 */
  1, /* 239 */
  2, /* 240 */
  3, /* 241 */
  3, /* 242 */
  2, /* 243 */
  3, /* 244 */
  2, /* 245 */
  4, /* 246 */
  1, /* 247 */
  3, /* 248 */
  4, /* 249 */
  2, /* 250 */
  1, /* 251 */
  2, /* 252 */
  1, /* 253 */
  1, /* 254 */
  0 /* 255 */
};
