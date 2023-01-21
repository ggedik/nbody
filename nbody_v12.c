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
      u64 r1   = (u64)rand();
      u64 r2   = (u64)rand();
      f32 sign = (r1 > r2) ? 1 : -1;
      //
      p.x[i] =   sign * (f32)rand() / (f32)RAND_MAX;
      p.y[i] =   (f32)rand() / (f32)RAND_MAX;
      p.z[i] =   sign * (f32)rand() / (f32)RAND_MAX;
      //
      p.vx[i] =  (f32)rand() / (f32)RAND_MAX;
      p.vy[i] =  sign * (f32)rand() / (f32)RAND_MAX;
      p.vz[i] =  (f32)rand() / (f32)RAND_MAX;
    }
}

//
void move_particles(particle_t p, //a struct
                    const f32 dt, 
                    u64 n) //just this part to translate 
{
  //
 // const f32 softening = 1e-20;
  f32 softening[8] = {1e-20,1e-20,1e-20,1e-20,1e-20,1e-20,1e-20,1e-20};  //what does it do? 
   //we'll stock 8 elements of p.x[i] in there 
   //although i dont know how to do thaat with the same problem with p.y[i]
  
   //ça mange un peu memoire ici, à penser 
  float* px = p.x;
  float* py = p.y;
  float* pz = p.z;
  float* pvx= p.vx;
  float* pvy= p.vy;
  float* pvz= p.vz;
  float* dtt= (float*) aligned_alloc(32,sizeof(float)*8);
  dtt[0] = dt;
  dtt[1] = dt;
  dtt[2] = dt;
  dtt[3] = dt;
  dtt[4] = dt; //i didnt want to include whole strging.h library just to memset this
  dtt[5] = dt;
  dtt[6] = dt;
  dtt[7] = dt; //so lost omg 
  //unsigned long long int i; 
  /*for (u64 i = 0; i < n; i++){ //for every particle 
      //initializing forces 
      f32 fx = 0.0; //do it in a smart way 
      f32 fy = 0.0;
      f32 fz = 0.0;
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
*/
  //for the ones above fx/y/z should be an array how to handle that?? did i misunderstood something?
  //ASK: unroll multiple variable such as below; 
  //3 floating-point operations


      //23 floating-point operations
     
  __asm__ volatile(
    "xor %%rax, %%rax;\n"
    "xor %%rcx, %%rcx;\n"
    "xor %%rdx, %%rdx;\n"
    "xor %%rbx, %%rbx;\n"
    
    "loop1:\n"
    "vmovaps (%[_px], %%rbx), %%ymm0;\n" //a
    "vmovaps (%[_py], %%rbx), %%ymm1;\n" //b
    "vmovaps (%[_pz], %%rbx), %%ymm2;\n" //c  //store first 8 elements in these registers  //c
   
    "vxorps %%ymm10, %%ymm10, %%ymm10;\n" //fx
    "vxorps %%ymm11, %%ymm11, %%ymm11;\n" //fy
    "vxorps %%ymm12, %%ymm12, %%ymm12;\n" //fz

    //"inc (%[_l]);\n"
    "xor %%rax, %%rax;\n"
    "xor %%rdx, %%rdx;\n"
    "cmp %[_n], %%rcx;\n"
    "jne loop2;\n"

    "loop2:\n"
    
    "vmovaps (%[_px], %%rax), %%ymm3;\n"
    "vmovaps (%[_py], %%rax), %%ymm4;\n"
    "vmovaps (%[_pz], %%rax), %%ymm5;\n"
    
    "vsubps %%ymm0, %%ymm3, %%ymm3;\n"          //ymm3 is dx now
    "vsubps %%ymm1, %%ymm4, %%ymm4;\n"          //ymm4 -ymm1 ==> ymm4 //dy
    "vsubps %%ymm2, %%ymm5, %%ymm5;\n"          //dz

    "vmulps %%ymm3, %%ymm3, %%ymm7;\n"          //dx*dx to store in ymm7
    "vmulps %%ymm4, %%ymm4, %%ymm8;\n"          //dy * dy in ymm8
    "vmulps %%ymm5, %%ymm5, %%ymm9;\n"          //dz * dz in ymm9

    "vaddps %%ymm7, %%ymm8, %%ymm13;\n"          //dy*dy + dx*dx stored in 7
    "vaddps %%ymm13, %%ymm9, %%ymm13;\n"          //ymm7 is d_2
    "vaddps (%[_softening]), %%ymm13, %%ymm13;\n" //ymm13 = d_2

    "vmulps  %%ymm13, %%ymm13, %%ymm6;\n"         //d_2*d_2 store in ymm6
    "vmulps  %%ymm6, %%ymm13, %%ymm13;\n"         //pow(d_2, 3)
    "vsqrtps %%ymm13, %%ymm13;\n"   //sqrt how to use less of it, ymm6 is d_3/2 
  
  //"vdivps ymm1, ymm2, ymm3 intel synthax felix courtier: divide ymm2 by ymm3 "
  //so i guess it will be: vdivps ymm2, ymm3, ymm1 : ymm2/ymm3 => ymm1
   
    "vdivps %%ymm3, %%ymm13, %%ymm3;\n"          //fx // allxorped in outer loop=>loop1
    "vdivps %%ymm4, %%ymm13, %%ymm4;\n"          //fy
    "vdivps %%ymm5, %%ymm13, %%ymm5;\n"          //fz
    //fx +    dx 
    "vaddps %%ymm10, %%ymm3, %%ymm10;\n"        //ymm10 = fx
    "vaddps %%ymm11, %%ymm4, %%ymm11;\n"        //ymm11 = fy
    "vaddps %%ymm12, %%ymm5, %%ymm12;\n"        //ymm12 = fz
   
    //"inc (%[_k]);\n"
    "add $32, %%rax;\n"
    "add $8, %%rdx;\n"
    "cmp %[_n], %%rdx;\n"                       //8 elements at a time
    "jne loop2;\n"
    "je con;\n"
    
    "con:"
    "vmovaps (%[_pvx], %%rbx), %%ymm3;\n"     //we are done with 3,4,5 when we get out of the loop 2 so we can use them here 
    "vmovaps (%[_pvy], %%rbx), %%ymm4;\n"   
    "vmovaps (%[_pvz], %%rbx), %%ymm5;\n"     //takes by packets
    "vmovaps (%[_dtt]), %%ymm13;\n"            //takes the fisrt eight at  once 

   
    "vmulps %%ymm10, %%ymm13, %%ymm7;\n"       //fx *dt
    "vmulps %%ymm11, %%ymm13, %%ymm8;\n"       //since we are done with 7,8,9 we can us ethem here and they are the destination operand so no need to zero them
    "vmulps %%ymm12, %%ymm13, %%ymm9;\n"

    "vaddps %%ymm7, %%ymm3, %%ymm3;\n"
    "vaddps %%ymm8, %%ymm4, %%ymm4;\n"
    "vaddps %%ymm9, %%ymm5, %%ymm5;\n"

    //"vmovaps %%ymm3, (%[_pvx], %%rbx);\n"
    //"vmovaps %%ymm4, (%[_pvy], %%rbx);\n"
    //"vmovaps %%ymm5, (%[_pvz], %%rbx);\n"
    //"inc (%[_j]);\n"
   
    "add $32, %%rbx;\n"
    "add $8, %%rcx;\n"
    "cmp %[_n], %%rcx;"
    "jne loop1;\n"
    "je end;\n"
    
    "end:"
    
    //"inc (%[_i]);\n" //finito
    
    ://outputs
    ://inputs
    //[_l] "r" (&l),
   // [_k] "r" (&k),
    //[_j] "r" (&j),
    //[_i]   "r" (&i),
    [_px]  "r" (px),
    [_py]  "r" (py),
    [_n]   "r" (n),
    [_pz]  "r" (pz),
    [_pvx] "r" (pvx),
    [_pvy] "r" (pvy),
    [_pvz] "r" (pvz),
    [_dtt] "r" (dtt),
    [_softening] "r" (softening)

    ://clobbers
    "cc", "memory", "rax", "rcx", "rdx", "rbx", "ymm1", "ymm2", "ymm3", "ymm4", "ymm5", "ymm6", "ymm7", "ymm8", "ymm9", "ymm10", "ymm11", "ymm12","ymm13"

  );
  

/*will merge it later */
  __asm__ volatile(
    "xor %%rcx, %%rcx;\n"//counter = i 
    "xor %%rdx, %%rdx;\n" //at first adding this line was causing segfault on _py but now it's fine???? 

    "loop:\n"
    "vmovaps (%[_px], %%rdx), %%ymm0;\n"        //load the pointed values by packets of 8
    "vmovaps (%[_pvx],%%rdx), %%ymm1;\n"        //load by packets of 8 p.vx
    "vmovaps (%[_dtt]), %%ymm2;\n"
    "vfmadd231ps %%ymm2, %%ymm1, %%ymm0;\n"          //dtt * pvx and store it in ymm3
                 //add ymm3 to ymm0 and store it in ymm0 =>felix dit qu'il peut y avoir trois registres?
   
    "vmovaps (%[_py], %%rdx),    %%ymm3;\n"        //load by packets of 8
    "vmovaps (%[_pvy], %%rdx),   %%ymm4;\n"       //load by packets of 8 p.vx //why xoring rdx causes segmentation fault when i include this line?
    "vfmadd231ps %%ymm2, %%ymm4, %%ymm3;\n"          //dtt * pvx and store it in ymm3
            //add ymm3 to ymm0 and store it in ymm0 =>felix dit qu'il peut y avoir trois registres?
    
    "vmovaps (%[_pz], %%rdx),  %%ymm5;\n"       //load by packets of 8
    "vmovaps (%[_pvz], %%rdx), %%ymm6;\n"      //load by packets of 8 p.vx
    "vfmadd231ps %%ymm2, %%ymm6, %%ymm5;\n"        //dtt * pvx and store it in ymm3

    "vmovaps %%ymm0, (%[_px], %%rdx);\n"   
    "vmovaps %%ymm3, (%[_py], %%rdx);\n" 
    "vmovaps %%ymm5, (%[_pz], %%rdx);\n"   //WHY THIS LINE CHANGES THE RESULTS??
        
    "add $32, %%rdx;\n"                       //add 64 should give me segmentation fault but it doesnt??
    //"inc (%[_i]);\n" //loop works
    "add $8, %%rcx;\n" //add 8 to rcx and store it in rcx => because we'll load our elements by packets of 8 floats => avx2 registers are 32bytes
    "cmp %[_n], %%rcx;\n"
    "jne loop;\n"
   
    ://outputs

    ://inputs
    //[_i] "r" (&i),
    [_px]  "r" (px),
    [_pvx] "r" (pvx),
    [_dtt] "r" (dtt),
    [_py]  "r" (py),
    [_pvy] "r" (pvy),
    [_pz]  "r" (pz),
    [_pvz] "r" (pvz),
    [_n]   "r" (n)
    
    ://clobbers
    "cc", "memory", "rcx", "ymm0", "ymm1", "ymm2", "ymm3", "ymm4", "ymm5", "ymm6", "rdx"

  );
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
  p.x = (f32*) aligned_alloc(32,sizeof(f32) * n); //array of our structs : n tane element 
  p.y = (f32*) aligned_alloc(32,sizeof(f32) * n);
  p.z = (f32*) aligned_alloc(32,sizeof(f32) * n);
  p.vx= (f32*) aligned_alloc(32,sizeof(f32) * n);
  p.vy= (f32*) aligned_alloc(32,sizeof(f32) * n); 
  p.vz= (f32*) aligned_alloc(32,sizeof(f32) * n);
  //
  init(p, n);
 for (int i = 0; i < n; i++){
    printf("p.x[i] = %f\n", p.x[i]);
  }
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
