#include <stdio.h>
#include <stdlib.h>
#include "veloc.h"
#include "simul.h"



int main(int argc, char **argv)
{

  // vérifie les arguments
  if(argc < 2)
  {
    printf("Usage : [bin] [file] [STEP] [Berendsen]\n");
    exit(0);
  }

  char *fname = argv[1];
  ///int N_step = atoi(argv[2]);
  //int test = atoi(argv[3]);
  
  particles_t *p = NULL;
  //lennard_t lj;
  vec_t *tv;
  cinet_t Cn;
  
  //double **rt ;

  // lit les données dans le nom du fichier
  p = read_data(fname);

 
  
  
  printf(".............Molécule Compute Dynamics.............\n");

  // calcule la vitesse du verlet
  Cn = compute_velocity_verlet(p, tv, N_sym);

  // simulation_dynamique_molecule(p, cn, N_sym, N_step);
  
  free_cinetique(Cn);
  free_particle(p);
  return 0;
}
