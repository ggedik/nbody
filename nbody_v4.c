//SOA + ALLIGN + UNROLL
//ALIGNMENT IS NOT POWER OF 2 BUT IDK WHAT TO DO
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
      u64 r1   = (u64)rand();
      u64 r2   = (u64)rand();
      f32 sign = (r1 > r2) ? 1 : -1;
      
      //
      p.x[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p.y[i] = (f32)rand() / (f32)RAND_MAX;
      p.z[i] = sign * (f32)rand() / (f32)RAND_MAX;
    
      p.vx[i] = (f32)rand() / (f32)RAND_MAX;
      p.vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p.vz[i] = (f32)rand() / (f32)RAND_MAX;
    }
}

//
void move_particles(particle_t p, //a struct
                    const f32 dt, 
                    u64 n) //just this part to translate 
{
  //
  const f32 softening = 1e-20;  //what does it do? 
  
      f32 *d_2 = (f32*)aligned_alloc(16,sizeof(f32) * 4);
      f32 *d_3_over_2 = (f32*) aligned_alloc(16, sizeof(f32) * 4); //there is a tradeoff here with unrolling and malloctime we dont know if it's worth
      f32 *dx = aligned_alloc(16,sizeof(f32)* 4);  
      f32 *dy = aligned_alloc(16,sizeof(f32)* 4);
      f32 *dz = aligned_alloc(16,sizeof(f32)* 4);  
  for (u64 i = 0; i < n; i++){ //for every particle 
      //initializing forces 
      f32 fx = 0; 
      f32 fy = 0;
      f32 fz = 0;
      //23 floating-point operations
      for (u64 j = 0; j < n; j=j+4){
	       //Newton's law
	      dx[0] = p.x[j]   - p.x[i]; //1 (sub) //how much a particle changed position
	      dx[1] = p.x[j+1] - p.x[i];
        dx[2] = p.x[j+2] - p.x[i];
        dx[3] = p.x[j+3] - p.x[i];
      
        dy[0] = p.y[j]   - p.y[i]; //2 (sub)
        dy[1] = p.y[j+1] - p.y[i];
        dy[2] = p.y[j+2] - p.y[i];
        dy[3] = p.y[j+3] - p.y[i];

        dz[0] = p.z[j]   - p.z[i];
        dz[1] = p.z[j+1] - p.z[i];
        dz[2] = p.z[j+2] - p.z[i];
        dz[3] = p.z[j+3] - p.z[i]; 

        d_2[0] = (dx[0] * dx[0]) + (dy[0] * dy[0]) + (dz[0] * dz[0]) + softening; 
        d_2[1] = (dx[1] * dx[1]) + (dy[1] * dy[1]) + (dz[1] * dz[1]) + softening;
        d_2[2] = (dx[2] * dx[2]) + (dy[2] * dy[2]) + (dz[2] * dz[2]) + softening;
        d_2[3] = (dx[3] * dx[3]) + (dy[3] * dy[3]) + (dz[3] * dz[3]) + softening;//9 (mul, add)
	      
        d_3_over_2[0] = pow(d_2[0], 3.0 / 2.0); //11 (pow, div) //euclidian distance 
        d_3_over_2[1] = pow(d_2[1], 3.0 / 2.0);
        d_3_over_2[2] = pow(d_2[2], 3.0 / 2.0);
        d_3_over_2[3] = pow(d_2[3], 3.0 / 2.0);
	     // Net force //16 byte = > scalar instructions
	      fx += dx[0] / d_3_over_2[0] + dx[1] / d_3_over_2[1] + dx[2] / d_3_over_2[2] + dx[3] / d_3_over_2[3];
	      fy += dy[0] / d_3_over_2[0] + dy[1] / d_3_over_2[1] + dy[2] / d_3_over_2[2] + dy[3] / d_3_over_2[3];
	      fz += dz[0] / d_3_over_2[0] + dz[1] / d_3_over_2[1] + dz[2] / d_3_over_2[2] +  dz[3] / d_3_over_2[3];
	    }
        //neyse okay we'll understand the mathematics later ==> youtube 
      //
      p.vx[i]   += dt * fx;
      p.vy[i]   += dt * fy;
      p.vz[i]   += dt * fz;

    }

  //3 floating-point operations (unrolling of this party is okay)
  for (u64 i = 0; i < n; i=i+4){ 
      p.x[i]   += dt * p.vx[i]; 
      p.x[i+1] += dt * p.vx[i+1];
      p.x[i+2] += dt * p.vx[i+2];
      p.x[i+3] += dt * p.vx[i+3];
      
      p.y[i]   += dt * p.vy[i];
      p.y[i+1] += dt * p.vy[i+1];
      p.y[i+2] += dt * p.vy[i+2];
      p.y[i+3] += dt * p.vy[i+3];
      
      p.z[i]   += dt * p.vz[i];
      p.z[i+1] += dt * p.vz[i+1];
      p.z[i+2] += dt * p.vz[i+2];  
      p.z[i+3] += dt * p.vz[i+3]; 
    }
}


//
int main(int argc, char **argv)
{
  //
  const u64 n     = (argc > 1) ? atoll(argv[1]) : 800; //number of particles
  const u64 steps = 10;         
  const f32 dt    = 0.01;
//n tane partikulun bir birlerine gore olan konumlarini/interactionlarini 10 step icinde 
  //
  f64 rate = 0.0, drate = 0.0;

  //Steps to skip for warm up
  const u64 warmup = 3;
  
  //
  particle_t p;
   p.x = aligned_alloc(64,sizeof(f32) * n); //array of our structs n tnae 
   p.y = aligned_alloc(64,sizeof(f32) * n);
   p.z = aligned_alloc(64,sizeof(f32) * n);
   p.vx= aligned_alloc(64,sizeof(f32) * n);
   p.vy= aligned_alloc(64,sizeof(f32) * n);
   p.vz= aligned_alloc(64,sizeof(f32) * n);
  //
  init(p, n);
//
// for (int i = 0; i < n; i++){
//  printf("p.x[i] = %f\n", p.x[i]);
//}
  const u64 s = sizeof(particle_t) * n; //soa aos this doesnt change
  
  printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
  //s >> 20 = s / 2^20 smart way, shoft to right 20 times so make divide it 20 times by 2
  //how much does a sighle bitshift take in cycles? 
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
//  
//    for (int i = 0; i < n; i++){
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
