#include "veloc.h"
#include<time.h>
#include<math.h>

// x if y >= 0.0, -x else
#define sign_function(x, y) (y < 0.0 ? -x : x)



cinet_t compute_velocity_verlet(particles_t *p, vec_t *tv, int N)
{

    double **rt;

    lennard_jones lj;

    cinet_t ct;

    ct.mi = (moment_t*)malloc(sizeof(moment_t)*p->N_particles_total);
    particles_t *ps = (particles_t*)malloc(sizeof(particles_t)*p->N_particles_total);

    //met à jour les particules
    p = update_particle_data(p, tv, N_sym);
    //distance
    rt = compute_distance(p);
    printf("rt10 = %lf\n",rt[1][0]);
    
    //Partie pour calculer le périodique Lennard Jones 
    lj = compute_Lennard_Jones_periodic(p, rt, N_sym);

    printf("Ulj = %lf\n", lj.Ulj);

    printf("lj.som_frc[1].fx = %lf\n", lj.som_frc[1].fx);
    // moment d'initialisation
    for (int i = 0; i < p->N_particles_total; i++)
    {
        ct.mi[i].mx = dt * CONVERSION_FORCE * lj.som_frc[i].fx /2.0;
        ct.mi[i].my = dt * CONVERSION_FORCE * lj.som_frc[i].fy / 2.0;
        ct.mi[i].mz = dt * CONVERSION_FORCE * lj.som_frc[i].fz / 2.0;
    }

    printf("ct.mi[1].mx = %lf\n", ct.mi[1].mx);
    
    ps->N_particles_total = p->N_particles_total;
    //mise à jour de la position
    for (int i = 0; i < ps->N_particles_total ; i++)
    {
        ps[i].coord.x = dt * ct.mi[i].mx / M_i;
        ps[i].coord.y = dt * ct.mi[i].my / M_i;
        ps[i].coord.z = dt * ct.mi[i].mz / M_i;
    }
    //met à jour le calcul de la distance
    rt = compute_distance(ps);
    printf("rt10 = %lf\n",rt[1][0]);
    //uptdate force 
    lj = compute_Lennard_Jones_periodic(ps, r, N_sym);

    printf("Ulj = %lf\n", lj.Ulj);

    //mettre à jour les instants
    for (int i = 0; i < p->N_particles_total; i++)
    {
        ct.mi[i].mx = dt * CONVERSION_FORCE * lj.som_frc[i].fx ;
        ct.mi[i].my = dt * CONVERSION_FORCE * lj.som_frc[i].fy / 2.0;
        ct.mi[i].mz = dt * CONVERSION_FORCE * lj.som_frc[i].fz / 2.0;
    }

    printf("ct.mi[1].mx = %lf\n", ct.mi[1].mx);
    //init un autre paramètre
    ct.Nl = 0;
    ct.Ec = 0.0;
    ct.Tc = 0.0;
    //gratuit
    free_distance(rt,p);
    free_lennard(lj,p);

    return ct;
}

/* moment d'initialisation cinetique
*params : particules p
*/
cinet_t init_moment_cinetique(particles_t *p)
{
    cinet_t ct;
    double c, s;

    ct.mi = (moment_t*)malloc(sizeof(moment_t) * p->N_particles_total);

    srand(time(NULL));

    for (int i = 0; i < p->N_particles_total; i++)
    {
        //composant dans x
        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        ct.mi[i].mx = sign_function(1.0, 0.5 - s) * c;

        //composant dans y
        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        ct.mi[i].my = sign_function(1.0, 0.5 - s) * c;

        //composant dans z
        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        ct.mi[i].mz = sign_function(1.0, 0.5 - s) * c;
    }
    
    return ct;
}

/*
*calculer l'énergie cinétique
*params : dans le moment
*/
cinet_t compute_cinetique_energie(cinet_t ct, particles_t *p)
{
    double tp = 0;

    ct.Nl = 3 * p->N_particles_total -3;

    for (int i = 0; i < p->N_particles_total; i++)
    {
        tp += (ct.mi[i].mx * ct.mi[i].mx + ct.mi[i].my * ct.mi[i].my + ct.mi[i].mz * ct.mi[i].mz)/M_i;
    }
    
    ct.Ec = tp / (2* CONVERSION_FORCE);
    
    // Température cinétique
    ct.Tc = ct.Ec / (ct.Nl * CONSTANT_R);

    return ct;
}

/*
*premier recalibrage
*/
cinet_t compute_first_recalibrated(cinet_t ct, particles_t *p)
{
    double Re;

    Re = ct.Nl * CONSTANT_R * T0 / ct.Ec;

    for (int i = 0; i < p->N_particles_total; i++)
    {
        //composant dans x
        ct.mi[i].mx = Re * ct.mi[i].mx;

        //composant dans y
        ct.mi[i].my = Re * ct.mi[i].my;

        //composant dans z
        ct.mi[i].mz = Re * ct.mi[i].mz;
    }
    
    return ct;
}
/*
*deuxième recalibré
*/
cinet_t compute_second_recalibrated(cinet_t ct, particles_t *p)
{
    double x = 0, y = 0, z =0;
    //correction
    for (int i = 0; i < p->N_particles_total; i++)
    {
        x += ct.mi[i].mx;
        y += ct.mi[i].my;
        z += ct.mi[i].mz;
    }
    
    for (int i = 0; i < p->N_particles_total; i++)
    {
        //composant dans x
        ct.mi[i].mx = ct.mi[i].mx - x;

        //composant dans y
        ct.mi[i].my = ct.mi[i].my - y;

        //composant dans z
        ct.mi[i].mz = ct.mi[i].mz - z;
    }

    printf("x = %lf, y = %lf, z = %lf \n", x, y, z);

    return ct;
}

/*
*compute thermostat berendsen
*params: moments particles
*/
cinet_t compute_mc_berendsen(cinet_t ct, particles_t *p)
{
    for (int i = 0; i < p->N_particles_total; i++)
    {
        //composant dans x
        ct.mi[i].mx += gama*((T0/ct.Tc) - 1)*ct.mi[i].mx;

        //composant dans y
        ct.mi[i].my = gama*((T0/ct.Tc) - 1)*ct.mi[i].my;

        //composant dans z
        ct.mi[i].mz = gama*((T0/ct.Tc) - 1)*ct.mi[i].mz;
    }

    return ct;
}
// libère de la mémoire après utilisation
void free_cinetique(cinet_t ct)
{
    free(ct.mi);
}
