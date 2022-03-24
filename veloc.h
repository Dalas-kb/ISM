#ifndef _VELOCITY_H_
#define _VELOCITY_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "simul.h"



// composant de moment struct
typedef struct moment_s
{
  double mx;
  double my;
  double mz;

}moment_t;

// résultat de la structure moment et température
typedef struct cinet_s
{
    moment_t *mi;
    int Nl;
    double Ec;
    double Tc;
    
}cinet_t;

// fonction vélocité verlet ;
cinet_t compute_velocity_verlet(particles_t *p, vec_t *tv, int N);
// moment d'initialisation de la fonction
cinet_t init_moment_cinetique(particles_t *p);
// énergie cinétique
cinet_t compute_cinetique_energie(cinet_t cn, particles_t *p);
// fonction recalibrée
cinet_t compute_first_recalibrated(cinet_t cn, particles_t *p);
// deuxième recalibrage
cinet_t compute_second_recalibrated(cinet_t cn, particles_t *p);

// thermostat berendsen
cinet_t compute_mc_berendsen(cinet_t cn , particles_t *p);
// mémoire libre
void free_cinetique(cinet_t cn);
#endif
