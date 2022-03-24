#ifndef _LENNARDJONES_H_
#define _LENNARDJONES_H_
#include <string.h>
#include <stdio.h>
#include <stdlib.h>




// définit un paramètre de valeur différente
#define sigma 3.0
#define eps 0.2
#define L   30
#define R_cut 10.0
#define N_sym 27
#define dt  1
#define CONVERSION_FORCE  0.0001*4.186
#define M_i 18
#define CONSTANT_R  0.00199
#define T0  300
#define gama   0.01

// structure pour les différentes forces
typedef struct force_s
{
  double fx;
  double fy;
  double fz;
}force_t;

// structure pour lennard jones
typedef struct lennard
{
  /* structure pour lennard jones */
  
  double Uj;
  force_t **frc;
  force_t *som_frc;
  
}lennard_jones;



/* stocker la coordonnée xyz */
typedef struct coord_s
{
  double x;
  double y;
  double z;
}coord_t;

// structure pour stock les différentes valeurs pour les particules
typedef struct particles_s
{
  int N_particles_total;
  int N_particles_local;
  int type;
  coord_t coord;
}particles_t;


// vecteur traducteur
typedef struct vec_s
{
  double x;
  double y;
  double z;

}vec_t;

lennard_jones compute_Lennard_Jones(particles_t *p, double **r);
void lennard_verify(lennard_jones jon, particles_t *p);
void print_forces(lennard_jones jon, particles_t *p);
void print_sum_forces(lennard_jones jon, particles_t *p);
lennard_jones compute_Lennard_Jones_periodic(particles_t *p, double **r, int N);
void free_lennard(lennard_jones jon, particles_t *p);





// fonction pour retourner les particules de structure
particles_t *read_data(char *fname);
// fonction pour imprimer les données des particules
void print_data(particles_t *p);
// fonction pour calculer la distance rij
double **compute_distance(particles_t *p);
// vecteur de transalteur de fonction
vec_t *translator_vector_init(int N);
// distance de mise à jour de la fonction
particles_t *update_particle_data(particles_t *p, vec_t *tv, int N);
// libère de la mémoire après utilisation
void free_vector_translate( vec_t *tv);
void free_distance(double **r, particles_t *p);
void free_particle(particles_t *p);


#endif
