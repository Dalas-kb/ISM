#include "simul.h"
#include <math.h>


lennard_jones compute_Lennard_Jones(particles_t *p, double **r)
{
  
  double somfx = 0.0;
  double somfy = 0.0;
  double somfz = 0.0;

  lennard_jones lj;
  double tm = 0.0;

  lj.frc=(force_t**)malloc(sizeof(force_t*)*p->N_particles_total);
  lj.som_frc=(force_t*)malloc(sizeof(force_t)*p->N_particles_total);
   
  for (int i = 0; i < p->N_particles_total; i++)
  {
      lj.frc[i]=(force_t*)malloc(sizeof(force_t)*p->N_particles_total);
  }
    
  //partie initialisation des forces   
  for (int i = 0; i < p->N_particles_total; i++)
  {
    for (int j = 0; j < p->N_particles_total; j++)
    {
      //initialisation des forces 
      lj.frc[i][j].fx = 0.0;
      lj.frc[i][j].fy = 0.0;
      lj.frc[i][j].fz = 0.0;
    }

    //initialisation des sommes des forces 
    lj.som_frc[i].fx = 0.0;
    lj.som_frc[i].fy = 0.0;
    lj.som_frc[i].fz = 0.0;
  }

    for (int i = 0; i < p->N_particles_total; i++)
    {
      for (int j = i+1; j < p->N_particles_total; j++)
       {
        
           tm +=((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
           /(r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
           - 2*((sigma*sigma*sigma*sigma*sigma*sigma)/(r[i][j] * r[i][j] *r[i][j]));

           //calculs les forces avec la formule 
           lj.frc[i][j].fx = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.x - p[j].coord.x);
           lj.frc[i][j].fy = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.y - p[j].coord.y);
           lj.frc[i][j].fz = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.z - p[j].coord.z);

          //mise a jour de la force j avec i 
          lj.frc[j][i].fx =-lj.frc[i][j].fx;
          lj.frc[j][i].fy =-lj.frc[i][j].fy;
          lj.frc[j][i].fz =-lj.frc[i][j].fz;
         
         //mise a jour de la force i avec j 
        lj.frc[i][j].fx =lj.frc[i][j].fx;
        lj.frc[i][j].fy =lj.frc[i][j].fy;
        lj.frc[i][j].fz =lj.frc[i][j].fz;

     
          somfx += lj.frc[i][j].fx;
          somfy += lj.frc[i][j].fy;
          somfz += lj.frc[i][j].fz;
       }

        //affecter les somme des forces 
        lj.som_frc[i].fx = somfx;
        lj.som_frc[i].fy = somfy;
        lj.som_frc[i].fz = somfz;
      }
    
    //calulc de lennard jones terme avec la formule du cours 
    lj.Uj = 4*eps*tm;

  return ljlj;
}


particles_t *update_particle_data(particles_t *p, vec_t *tv, int N)
{

  particles_t *tm = (particles_t*)malloc(sizeof(particles_t)*p->N_particles_total);
  tm->N_particles_total = p->N_particles_total;

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < tm->N_particles_total; j++)
    {
      tm[j].coord.x = p[j].coord.x + tv[i].x;
      tm[j].coord.y = p[j].coord.y + tv[i].y;
      tm[j].coord.z = p[j].coord.z + tv[i].z;
    }
  }
  
  return tm;
}


lennard_jones compute_Lennard_Jones_periodic(particles_t *p, double **r, int N)
{
  /* calcul periodique de lennard jones */
  double somfx = 0.0;
  double somfy = 0.0;
  double somfz = 0.0;
  double k = 0.5;

  lennard_jones lj;
  double tm = 0.0;

  lj.frc=(force_t**)malloc(sizeof(force_t*)*p->N_particles_total);
  lj.som_frc=(force_t*)malloc(sizeof(force_t)*p->N_particles_total);

  for (int i = 0; i < p->N_particles_total; i++)
  {
      lj.frc[i]=(force_t*)malloc(sizeof(force_t)*p->N_particles_total);
  }
  
  
  for (int i = 0; i < p->N_particles_total; i++)
  {
    for (int j = 0; j < p->N_particles_total; j++)
    {
      /*initialisation des forces*/
      lj.frc[i][j].fx = 0.0;
      lj.frc[i][j].fy = 0.0;
      lj.frc[i][j].fz = 0.0;
    }

    /*initialisation des sommes des forces */
    lj.som_frc[i].fx = 0.0;
    lj.som_frc[i].fy = 0.0;
    lj.som_frc[i].fz = 0.0;
  }
  
  for (int i_sym = 0; i_sym < N; i_sym++)
  {
    for (int i = 0; i < p->N_particles_total; i++)
      {
        for (int j = i+1; j < p->N_particles_total; j++)
        {
          
          if (r[i][j] > R_cut*R_cut)
              continue;
          
            tm +=((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
           /(r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
           - 2*((sigma*sigma*sigma*sigma*sigma*sigma)/(r[i][j] * r[i][j] *r[i][j]));

            /*calculs les forces */
           lj.frc[i][j].fx = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.x - p[j].coord.x);
           lj.frc[i][j].fy = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.y - p[j].coord.y);
           lj.frc[i][j].fz = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.z - p[j].coord.z);

          /*mise a jour de la force j et i */
          lj.frc[j][i].fx =-lj.frc[i][j].fx;
          lj.frc[j][i].fy =-lj.frc[i][j].fy;
          lj.frc[j][i].fz =-lj.frc[i][j].fz;
         
          /*mise a jour de la force i et j*/ 
          lj.frc[i][j].fx =lj.frc[i][j].fx;
          lj.frc[i][j].fy =lj.frc[i][j].fy;
          lj.frc[i][j].fz =lj.frc[i][j].fz;

       
          somfx += lj.frc[i][j].fx;
          somfy += lj.frc[i][j].fy;
          somfz += lj.frc[i][j].fz;
       }

        /*affcter les sommes des forces*/ 
        lj.som_frc[i].fx = somfx;
        lj.som_frc[i].fy = somfy;
        lj.som_frc[i].fz = somfz;
       }
  }
    
     if(N == 1)
        k = 1;

  /*calulc de lennard jones terme */ 
  lj.Uj = 4*k*eps*tm;

  return lj;
}

//print forces for each particles
void print_forces(lennard_jones lj, particles_t *p)
{ 
  /* affichage des forces */
  for (int i = 0; i < p->N_particles_total; i++)
      {
        for (int j = 0; j < p->N_particles_total; j++)
        {
           printf("%lf  \t  %lf  \t  %lf", lj.frc[i][j].fx, lj.frc[i][j].fy, lj.frc[j][i].fz) ;
          
        }
        printf("\n");
       }
}


void print_sum_forces(lennard_jones lj, particles_t *p)
{ 
  /*affichage de la somme des forces */
  
  int i=0;
  printf("\a\a\a ## Somme des forces calculées ##\a\a\a \n");
  while(i < p->N_particles_total)
  {
       printf("Fx[%d] = %lf \t Fy[%d] = %lf  \t  Fz[%d] = %lf \n", i, lj.som_frc[i].fx , i, lj.som_frc[i].fy , i, lj.som_frc[i].fz); 
       i++;
  }
  
}

void lennard_verify(lennard_jones lj, particles_t *p)
{ 
  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;
  int i=0;
  int j=0;
  
  while (i < p->N_particles_total)
  {
    while(j < p->N_particles_total)
    {
      Fx=Fx+lj.frc[i][j].fx;
      Fy=Fy+lj.frc[i][j].fy;
      Fz=Fz+lj.frc[i][j].fz;
      j++;
    }
    i++;
  }

  printf("Fx = %lf  Fy = %lf  Fz = %lf\n", Fx,Fy,Fz);
}


void free_lennard(lennard_jones lj, particles_t *p)
{

  int i=0;
  while (i<p->N_particles_total)
  {
    	free(lj.frc[i]);
  	i++;
  }
  
  free(lj.frc);
  free(lj.som_frc);

}


/* lire les données pour les particules coord
  fichier param : nom du fichier
  retour : table des particules de coordonnées
*/
particles_t *read_data(char *fname)
{

    // vérifie si le fichier est nul
  if (fname == NULL)
  {
    printf("Error read data %s\n",__func__);
    exit(0);
  }
  FILE *file=fopen(fname, "r+");
    // vérifier le fichier
  if(!file)
    exit(1);
  int nb_elmt =0;
  int c;
  
  // vérifie le nombre d'éléments dans le fichier
  do
  {
    if(c=='\n')
      nb_elmt++;
  } while((c=fgetc(file))!=EOF);

  particles_t *p = (particles_t*)malloc(sizeof(particles_t)*nb_elmt);
 

  rewind(file);
  fscanf(file,"%d",&p->N_particles_total);

  // vérifie si nb_particles_total nombre de coordonnées différent
  if (p->N_particles_total != nb_elmt-1)
  {
    printf("Number patricle is different of the number coordinates\n" );
    exit(0);
  }

  for (int i = 0; i < p->N_particles_total; i++)
  {
    fscanf(file,"  %d   %lf   %lf   %lf",&p[i].type,&(p[i].coord.x), &(p[i].coord.y),&(p[i].coord.z));
    
  }
  return p;
}

/* Fonction pour les données d'impression
paramètre d'entrée : particules
*/
void print_data(particles_t *p)
{
  //check if pointer is NULL
  if(p == NULL)
  {
    printf("Error pointer is null\n" );
    exit(0);
  }
  printf("P  \t| Type  \t |X \t\t |Y \t\t |Z\n");

  for (int i = 0; i < p->N_particles_total; i++)
  {
    printf("%d \t %d \t %lf \t %lf \t %lf\n",i,p[i].type,p[i].coord.x, p[i].coord.y,p[i].coord.z);
  }

  printf("=======================================================\n");
}

/* calcule la distance des différentes particules
*params : particules
*distance de retour pour les différentes particules
*/
double **compute_distance(particles_t *p)
{
  // vérifie si le pointeur est NULL
  if(p == NULL)
  {
    printf("Error pointer is null\n" );
    exit(0);
  }

   double **r = NULL;
   r=(double**)malloc(p->N_particles_total*sizeof(double*));
   for (int i = 0; i < p->N_particles_total; i++)
     r[i] =(double*)malloc(p->N_particles_total*sizeof(double));


  for (int i = 0; i < p->N_particles_total; i++)
  {
    for (int j = i+1; j < p->N_particles_total; j++)
    {
        //compute distance
        r[i][j] = ((p[i].coord.x - p[j].coord.x)*(p[i].coord.x - p[j].coord.x))+
                 ((p[i].coord.y - p[j].coord.y)*(p[i].coord.y - p[j].coord.y))+
                 ((p[i].coord.z - p[j].coord.z)*(p[i].coord.z - p[j].coord.z));
        //update  distance
        r[j][i] = r[i][j];
    }
  }

  return r;
}

/* initialisation du vecteur de traduction
*/
vec_t *translator_vector_init(int N)
{
  vec_t *tv = (vec_t*)malloc(sizeof(vec_t)*N);

  for (int i = 0; i < N; i++)
  {
    tv[i].x = (double)((i / 9) - 1)*L;
    tv[i].y = (double)((i / 3)% 3 -1)*L;
    tv[i].z = (double)((i%3) - 1) *L;
  }
  
  return tv;
}

// mémoire libre
void free_vector_translate( vec_t *tv)
{
  free(tv);
}

void free_distance(double **r, particles_t *p)
{
  for (int i = 0; i < p->N_particles_total; i++)
              free(r[i]);

  free(r);
}

void free_particle(particles_t *p)
{
  free(p);
}
