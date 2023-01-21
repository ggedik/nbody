//*SOA*//
#include <omp.h>
#include <math.h> //we are using only 
#include <stdio.h>
#include <stdlib.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
typedef struct particle_s {

  f32 *x, *y, *z;      //float rather than double 
  f32 *vx, *vy, *vz;   // 6*4 = 24 bytes total => padding vectors? 
  
} particle_t;

//
void init(particle_t p, u64 n)   //nothing to do here 
{
  for (u64 i = 0; i < n; i++)
    {
      //
      u64 r1   =  (u64)rand();
      u64 r2   =  (u64)rand();
      f32 sign =  (r1 > r2) ? 1 : -1;
      //
      //1;//
      p.x[i] =     sign * (f32)rand() / (f32)RAND_MAX;
      p.y[i] =     (f32)rand() / (f32)RAND_MAX;
      p.z[i] =     sign * (f32)rand() / (f32)RAND_MAX;
   //
      //1;//
      p.vx[i] =   (f32)rand() / (f32)RAND_MAX;
      p.vy[i] =   sign * (f32)rand() / (f32)RAND_MAX;
      p.vz[i] =   (f32)rand() / (f32)RAND_MAX;
    }
}

//
void move_particles(particle_t p, //a struct
                    const f32 dt, 
                    u64 n) //just this part to translate 
{
  //
  const f32 softening = 1e-20;  //what does it do? 
  
  //
  for (u64 i = 0; i < n; i++){ //for every particle 
      //initializing forces 
      f32 fx = 0.0; //do it in a smart way 
      f32 fy = 0.0;
      f32 fz = 0.0;
    
      //23 floating-point operations
      for (u64 j = 0; j < n; j++){
	       //Newton's law
	      const f32 dx  = p.x[j] - p.x[i]; //1 (sub) //how much a particle changed position
	      const f32 dy  = p.y[j] - p.y[i]; //2 (sub)
	      const f32 dz  = p.z[j] - p.z[i]; //3 (sub)
	      const f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; //9 (mul, add)
	      const f32 d_3_over_2 = pow(d_2, 3.0 / 2.0); //11 (pow, div) //WHY??? 3/2
                                                    //supposed to be euclidian distance right? why it's not root to cube?

        //div = 1 /d_3_over_2; //force it to use fma ports with 
	      //Net force
	      fx += dx / d_3_over_2; //13 (add, div)    //how to reduce divisions? 
	      fy += dy / d_3_over_2; //15 (add, div)
	      fz += dz / d_3_over_2; //17 (add, div)
	    }
        //neyse okay we'll understand the mathematics later ==> youtube 
      //
      p.vx[i] += dt * fx; //19 (mul, add)   //fma 
      p.vy[i] += dt * fy; //21 (mul, add)
      p.vz[i] += dt * fz; //23 (mul, add)
    }

  //for the ones above fx/y/z should be an array how to handle that?? did i misunderstood something?
  //ASK: unroll multiple variable such as below; 
  //3 floating-point operations

  for (u64 i = 0; i < n; i++){ //vectorize them with a filler SIMD? 
      p.x[i] += dt * p.vx[i]; 
      p.y[i] += dt * p.vy[i];
      p.z[i] += dt * p.vz[i];
  }
}

//
int main(int argc, char **argv)
{
  //
  const u64 n     = (argc > 1) ? atoll(argv[1]) : 800; //number of particles
  const u64 steps = 10;         
  const f32 dt    = 0.01;
//n tane partikulun bir birlerine gore olan konumlarini/interactionlarini 10 step icinde .... //hesapliyoruz ama dt tam olarak ne simdi 
  //
  f64 rate = 0.0, drate = 0.0;

  //Steps to skip for warm up
  const u64 warmup = 3;
  
  //
  particle_t p;
  p.x = (f32*) malloc(sizeof(f32) * n); //array of our structs : n tane element 
  p.y = (f32*) malloc(sizeof(f32) * n);
  p.z = (f32*) malloc(sizeof(f32) * n);
  p.vx= (f32*) malloc(sizeof(f32) * n);
  p.vy= (f32*) malloc(sizeof(f32) * n); 
  p.vz= (f32*) malloc(sizeof(f32) * n);
  //
  init(p, n);

  // for (int i = 0; i < n; i++){
  //  printf("p.x[i] = %f\n", p.x[i]);
  //}

  const u64 s = sizeof(particle_t) * n; //total gb stays the smea since we just change the placements
  
  printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
  
  //
  printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); 
  fflush(stdout);
  
  //
  for (u64 i = 0; i < steps; i++)
    {
      //Measure
      const f64 start = omp_get_wtime(); //why not clock get time??
                                          //

      move_particles(p, dt, n);

      const f64 end = omp_get_wtime();

      //Number of interactions/iterations
      const f32 h1 = (f32)(n) * (f32)(n - 1);

      //GFLOPS
      const f32 h2 = (23.0 * h1 + 3.0 * (f32)n) * 1e-9;
      
      if (i >= warmup)
	{
	  rate += h2 / (end - start);
	  drate += (h2 * h2) / ((end - start) * (end - start));
	}
      
      //       step  time   average number of interactions, glops
      printf("%5llu %10.3e %10.3e %8.1f %s\n",
	     i,
	     (end - start), //THE DIFFERENCE OF TIME SO PERIVIOUS ONE CAN BE FASTER THAN THE CURRENT ONE THEREFORE HAVE A LESSER VALUE
	     h1 / (end - start),
	     h2 / (end - start),
	     (i < warmup) ? "*" : "");
      
      fflush(stdout); //fflush (file stream) but why????
    }
  
  //
  rate /= (f64)(steps - warmup);
  drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));
  
  printf("-----------------------------------------------------\n");
  printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
	 "Average performance:", "", rate, drate);
  printf("-----------------------------------------------------\n");
  
  //for (int i = 0; i < n; i++){
  //  printf("p.x[i] = %f\n", p.x[i]);
  //}
  //
  free(p.x);
  free(p.y);
  free(p.z);
  free(p.vx);
  free(p.vy);
  free(p.vz);
  
  
  //
  return 0;
}
