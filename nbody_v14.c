#include <omp.h>
#include <math.h> //we are using only 
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <pmmintrin.h>

typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

typedef struct particle_s {

  f32 *x, *y, *z;      //float rather than double 
  f32 *vx, *vy, *vz;   // 6*4 = 24 bytes total => padding vectors? 
  
} particle_t;

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

/***
I should have read the manual better and not be scared of unknown.
let's write everything in detail to remember later and learn better.
_m256 data type can hold 8 floats
<= why bitwise inverse of the value?
**her sey fonksiyon gibi atanmis 
_mm_set_pd(2.0, 1.0); => set packed double 2 and 1
assumption=> how do we do if we want to set arrays? => think like sse instructoions => set one portion of it and then shift the array pointer?
why the f they didnt put a simple example and explain veryline in detail with compilation etc 
_mm256_load_ps(float const *a) => VMOVAPS load aligned packed 8 single instruction from memory location
_mm256_loadu_ps(b+i) => i yi 8 arttir : i+=8 => one register takes 8 
why the f i get stressed out that much it's just code it's something that i do her gun yani streslenmeye gerek yok
keyif al stres yapma
of biktim ben kendimden
_mm256_fmadd_ps(num1, num2, num3) => num1 +num2 =+ num3
_mm256_set1_ps() : a float32 value to be initialized into the 256-bit vector
nefs al
resource: https://github.com/srinathv/ImproveHpc/blob/master/intel/2015-compilerSamples/C%2B%2B/intrinsic_samples/intrin_dot_sample.c
http://kayaogz.github.io/teaching/programmation-parallele-distribuee-2019/cours/simd.pdf
oguz abim be kral, gene budum hemserimi
The compiler defaults to using the VFMADD213PS
instruction and uses the other forms VFMADD132PS
or VFMADD231PS
only if a low level optimization decides it is useful or necessary. For example, the compiler could change the default 
if it finds that another instruction form saves a register or eliminates a move. 
ne pas oublier aligner tous, we are using aligned instructions here 
***/




void move_particles( particle_t p, const f32 dt, u64 n){
    float dtt = dt;
    float a ;
    __m256 dt_vec  = _mm256_set1_ps(dtt);
     f32 softening = 1e-20;
    __m256 soft = _mm256_set1_ps(softening);
    __m256 bir  = _mm256_set1_ps(1);
    __m256 tmp  = _mm256_set1_ps(0.0);
    __m128 top, bot ;
    __m128 top1, bot1;
    __m128 top2, bot2;
    for (u64 i = 0; i < n; i++){ //for every particle 
      //initializing forces 

      f32 fx = 0.0; 
      f32 fy = 0.0;
      f32 fz = 0.0;
    __m256 fxx = _mm256_set1_ps(0.0);
    __m256 fyy = _mm256_set1_ps(0.0);
    __m256 fzz = _mm256_set1_ps(0.0);
    __m256 px_temp = _mm256_set1_ps(p.x[i]);
    __m256 py_temp = _mm256_set1_ps(p.y[i]);
    __m256 pz_temp = _mm256_set1_ps(p.z[i]); 

      for (u64 j = 0; j < n; j+=8){
        __m256 bir  = _mm256_set1_ps(1.0);
          /*const f32 dx  = p.x[j] - p.x[i]; //1 (sub) //how much a particle changed position
	      const f32 dy  = p.y[j] - p.y[i]; //2 (sub)
	      const f32 dz  = p.z[j] - p.z[i]; //3 (sub)
	      const f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; //9 (mul, add)
	      const f32 d_3_over_2 = pow(d_2, 3.0 / 2.0); //11 (pow, div) //WHY??? 3/2
                                                    //supposed to be euclidian distance right? why it's not root to cube?

        //div = 1 /d_3_over_2; //force it to use fma ports with 
	      //Net force
	      fx += dx / d_3_over_2; //13 (add, div)    //how to reduce divisions? 
	      fy += dy / d_3_over_2; //15 (add, div)
	      fz += dz / d_3_over_2; //17 (add, div)*/
	       //Newton's law
           __m256 px_j = _mm256_load_ps(p.x +j);
           __m256 py_j = _mm256_load_ps(p.y +j);
           __m256 pz_j = _mm256_load_ps(p.z +j); 
           
//how to use mm_dp_ps? didnt get it
           __m256 dxx = _mm256_sub_ps(px_j, px_temp); 
           __m256 dyy = _mm256_sub_ps(py_j, py_temp); 
           __m256 dzz = _mm256_sub_ps(pz_j, pz_temp);
//
           __m256 d2x =_mm256_mul_ps(dxx,dxx);
           __m256 d2y =_mm256_mul_ps(dyy,dyy);
           __m256 d2z =_mm256_mul_ps(dzz,dzz);
//
           __m256 d22 = _mm256_add_ps(d2x, d2y);
                  d22 = _mm256_add_ps(d22,d2z);
                  d22 = _mm256_add_ps(d22, soft);
           __m256 d3o2 =_mm256_set1_ps(1);
           // _mm256_mul_ps(d22,d22); WHY THIS DOESNT WORK HOFFF SABIR YA SABIR
               d3o2 = _mm256_mul_ps(d22,d3o2);
                d3o2 = _mm256_rsqrt_ps(d3o2);
               // d3o2 = _mm256_div_ps(bir, d3o2);
//
           fxx = _mm256_fmadd_ps(dxx, d3o2, fxx);
           fyy = _mm256_fmadd_ps(dyy, d3o2, fyy);
           fzz = _mm256_fmadd_ps(dzz, d3o2, fzz);

	    }
        fxx = _mm256_hadd_ps(fxx,fxx); 
        top = _mm256_extractf128_ps(fxx,1);
        bot = _mm256_extractf128_ps(fxx,0);
        top = _mm_add_ps(top, bot);
        top = _mm_hadd_ps(top,top);
        _mm_store_ss(&fx,top);
      fyy = _mm256_hadd_ps(fyy,fyy); 
        top1 = _mm256_extractf128_ps(fyy,1);
        bot1 = _mm256_extractf128_ps(fyy,0);
        top1 = _mm_add_ps(top1, bot1);
        top1 = _mm_hadd_ps(top1,top1);
        _mm_store_ss(&fy,top1);

        fzz = _mm256_hadd_ps(fzz,fzz); 
        top2 = _mm256_extractf128_ps(fzz,1);
        bot2 = _mm256_extractf128_ps(fzz,0);
        top2 = _mm_add_ps(top2, bot2);
        top2 = _mm_hadd_ps(top2,top2);
        _mm_store_ss(&fz,top2);
          /*tmp = _mm256_fmadd_ps(dt_vec, fxx, tmp);
         __m128 hhalf = _mm256_extractf128_ps(tmp, 1);
         __m128 lhalf = _mm256_extractf128_ps(tmp, 0);
         __m128 total = _mm_add_ps(hhalf, lhalf); 
                //total = _mm_hadd_ps(total, total);
                a = _mm_cvtss_f32(hhalf);*/
    //FFFFF THIS PART IS HARD, BREATHHHHHH hffff
      //p.vx[i] = _mm_fmadd_ss(dt, fx, p.vx[i]);
         printf("fx, %f", fz);
      p.vx[i] += dt * fx; //19 (mul, add)   //fma 
      p.vy[i] += dt * fy; //21 (mul, add)
      p.vz[i] += dt * fz; //23 (mul, add)
    }
    for(u64 i =0; i < n; i+=16){ //unrolled by 2
        //__m256 dt_vec  = _mm256_set1_ps(dtt);
        __m256 px_vec  = _mm256_load_ps(p.x + i);
        __m256 py_vec  = _mm256_load_ps(p.y + i);
        __m256 pz_vec  = _mm256_load_ps(p.z + i);
        __m256 pvx_vec = _mm256_load_ps(p.vx + i); 
        __m256 pvy_vec = _mm256_load_ps(p.vy + i);
        __m256 pvz_vec = _mm256_load_ps(p.vz + i);
        
        px_vec = _mm256_fmadd_ps(dt_vec, pvx_vec, px_vec);
        py_vec = _mm256_fmadd_ps(dt_vec, pvy_vec, py_vec);
        pz_vec = _mm256_fmadd_ps(dt_vec, pvz_vec, pz_vec);
        
        _mm256_store_ps(p.x+i, px_vec);
        _mm256_store_ps(p.y+i, py_vec);
        _mm256_store_ps(p.z+i, pz_vec);

        __m256 px_vec2  = _mm256_load_ps(p.x + i+8);
        __m256 py_vec2  = _mm256_load_ps(p.y + i+8);
        __m256 pz_vec2  = _mm256_load_ps(p.z + i+8);
        __m256 pvx_vec2 = _mm256_load_ps(p.vx + i+8);
        __m256 pvy_vec2 = _mm256_load_ps(p.vy + i+8);
        __m256 pvz_vec2 = _mm256_load_ps(p.vz + i+8);
         
        px_vec2 = _mm256_fmadd_ps(dt_vec, pvx_vec2, px_vec2);
        py_vec2 = _mm256_fmadd_ps(dt_vec, pvy_vec2, py_vec2);
        pz_vec2 = _mm256_fmadd_ps(dt_vec, pvz_vec2, pz_vec2);

        _mm256_store_ps(p.x+i+8, px_vec2);
        _mm256_store_ps(p.y+i+8, py_vec2);
        _mm256_store_ps(p.z+i+8, pz_vec2);

    }


}


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
  p.x = (f32*) aligned_alloc(32,sizeof(f32) * n); //array of our structs : n tane element 
  p.y = (f32*) aligned_alloc(32,sizeof(f32) * n);
  p.z = (f32*) aligned_alloc(32,sizeof(f32) * n);
  p.vx= (f32*) aligned_alloc(32,sizeof(f32) * n);
  p.vy= (f32*) aligned_alloc(32,sizeof(f32) * n); 
  p.vz= (f32*) aligned_alloc(32,sizeof(f32) * n);
  //
  init(p, n);
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
  
  for (int i = 0; i < n; i++){
    printf("p.x[i] = %f\n", p.x[i]);
  }
  
  free(p.x);
  free(p.y);
  free(p.z);
  free(p.vx);
  free(p.vy);
  free(p.vz);
    //
  return 0;
}
