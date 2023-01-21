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
      p.x[i] =    sign * (f32)rand() / (f32)RAND_MAX;
      p.y[i] =   (f32)rand() / (f32)RAND_MAX;
      p.z[i] =    sign * (f32)rand() / (f32)RAND_MAX;
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
  //const f32 softening = 1e-20; 
  //
  f32 softening[8] = {1e-20,1e-20,1e-20,1e-20,1e-20,1e-20,1e-20,1e-20};  //what does it do? 
   //we'll stock 8 elements of p.x[i] in there 
   //although i dont know how to do thaat with the same problem with p.y[i]
  long long unsigned int i = 0;
  long long unsigned int j = 0;
  long long unsigned int k = 0;
  long long unsigned int l = 0;
   //ça mange un peu memoire ici, à penser 
  float* px = p.x;//(float*) aligned_alloc(32, sizeof(float) * 8);//p.x;
  float* py = p.y;//(float*) aligned_alloc(32, sizeof(float) * 8);//p.y;
  float* pz = p.z;//(float*) aligned_alloc(32, sizeof(float) * 8);//p.z;
  float* pvx= p.vx;//(float*) aligned_alloc(32, sizeof(float) * 8);//p.vx;
  float* pvy= p.vy;//(float*) aligned_alloc(32, sizeof(float) * 8);//p.vy;
  float* pvz= p.vz;//(float*) aligned_alloc(32, sizeof(float) * 8);//p.vz;
  float* dtt= (float*) aligned_alloc(32, sizeof(float) * 8);
  dtt[0] = dt;
  dtt[1] = dt;
  dtt[2] = dt;
  dtt[3] = dt;
  dtt[4] = dt; //i didnt want to include whole strging.h library just to memset this
  dtt[5] = dt;
  dtt[6] = dt;
  dtt[7] = dt; 
  float* dd = (float*) aligned_alloc(32,sizeof(float)*8);
  dd[0] = 1;
  dd[1] = 1;
  dd[2] = 1;
  dd[3] = 1;
  dd[4] = 1;
  dd[5] = 1;
  dd[6] = 1;
  dd[7] = 1; 
  float* tmpx =(float*) aligned_alloc(32,sizeof(f32)* 8);
  float* tmpy =(float*) aligned_alloc(32,sizeof(f32)* 8);
  float* tmpz =(float*) aligned_alloc(32,sizeof(f32)* 8);
  float* fxx = (float*) aligned_alloc(32,sizeof(f32)* 8);
  float* fyy = (float*) aligned_alloc(32,sizeof(f32)* 8);
  float* fzz = (float*) aligned_alloc(32,sizeof(f32)* 8);
  f32 *d_2 = (f32*)aligned_alloc(32,sizeof(f32) * 8);
  f32 *d_3_over_2 = (f32*) aligned_alloc(16, sizeof(f32) * 8); //there is a tradeoff here with unrolling and malloctime we dont know if it's worth
  f32 *dx = aligned_alloc(32,sizeof(f32)* 8);  
  f32 *dy = aligned_alloc(32,sizeof(f32)* 8);
  f32 *dz = aligned_alloc(32,sizeof(f32)* 8);  
  for (u64 i = 0; i < n; i++){ //for every particle 
      //initializing forces 
      f32 fx = 0; 
      f32 fy = 0;
      f32 fz = 0;
      tmpx[0] = p.x[i];
      tmpx[1] = p.x[i];
      tmpx[2] = p.x[i];
      tmpx[3] = p.x[i];
      tmpx[4] = p.x[i];
      tmpx[5] = p.x[i];
      tmpx[6] = p.x[i];      
      tmpx[7] = p.x[i];      //23 floating-point operations
      tmpx[0] = p.x[i];
      tmpy[0] = p.y[i];
      tmpy[1] = p.y[i];
      tmpy[2] = p.y[i];
      tmpy[3] = p.y[i];
      tmpy[4] = p.y[i];
      tmpy[5] = p.y[i];
      tmpy[6] = p.y[i];
      tmpy[7] = p.y[i];
      tmpz[0] = p.z[i];
      tmpz[1] = p.z[i];
      tmpz[2] = p.z[i];
      tmpz[3] = p.z[i];
      tmpz[4] = p.z[i];
      tmpz[5] = p.z[i];
      tmpz[6] = p.z[i];
      tmpz[7] = p.z[i]; //idk the smart way to do this

fxx[0] = 0;
fyy[0]= 0;
fzz[0] = 0;
fxx[1] = 0;
fyy[1]= 0;
fzz[1] = 0;
fxx[2] = 0;
fyy[2]= 0;
fzz[2] = 0;
fxx[3] = 0;
fyy[3]= 0;
fzz[3] = 0;
fxx[4] = 0;
fyy[4]= 0;
fzz[4] = 0;
fxx[5] = 0;
fyy[5]= 0;
fzz[5] = 0;
fxx[6] = 0;
fyy[6]= 0;
fzz[6] = 0;    
fxx[7] = 0;
fyy[7]= 0;
fzz[7] = 0;  
      __asm__ volatile( 
      "xor %%rcx, %%rcx;\n"
      "xor %%rdx, %%rdx;\n"
      
      "loop1:"
      "vmovaps (%[_tmpx]), %%ymm1;\n" //p.x[i]
      "vmovaps (%[_tmpy]), %%ymm2;\n"
      "vmovaps (%[_tmpz]), %%ymm3;\n"

      "vmovaps (%[_px], %%rdx), %%ymm4;\n" //p.x[j]
      "vmovaps (%[_py], %%rdx), %%ymm5;\n"
      "vmovaps (%[_pz], %%rdx), %%ymm6;\n"

      "vsubps %%ymm4, %%ymm1, %%ymm7;\n" //ymm7 is dx
      "vsubps %%ymm5, %%ymm2, %%ymm8;\n" //ymm8 = dy
      "vsubps %%ymm6, %%ymm3, %%ymm9;\n"  //ymm9 = dz

      "vmulps %%ymm7, %%ymm7, %%ymm10;\n" //dx[0..7] x dx[0..7] element by element
      "vmulps %%ymm8, %%ymm8, %%ymm11;\n"
      "vmulps %%ymm9, %%ymm9, %%ymm12;\n"
      "vaddps %%ymm10, %%ymm11, %%ymm13;\n"
      "vaddps %%ymm13, %%ymm12, %%ymm13;\n" //so d_2 = ymm13 eight elements
      "vaddps %[_softening], %%ymm13, %%ymm13;\n"

      "vmulps %%ymm13, %%ymm13, %%ymm10;\n"
      "vmulps %%ymm13, %%ymm10, %%ymm10;\n"
      "vsqrtps %%ymm10, %%ymm10;\n"   //ymm10 = d32

      "vdivps %%ymm7, %%ymm10, %%ymm13;\n" 
      "vdivps %%ymm8, %%ymm10, %%ymm11;\n"
      "vdivps %%ymm9, %%ymm10, %%ymm12;\n"

      //"vaddps %%ymm13, (%[_fxx]), %%ymm13;\n"
      //"vaddps %%ymm11, (%[_fyy]), %%ymm11;\n"
      //"vaddps %%ymm12, (%[_fzz]), %%ymm12;\n"

      //"vmovaps %%ymm13, (%[_fxx]);\n"
      //"vmovaps %%ymm12, (%[_fyy]);\n"
      //"vmovaps %%ymm11, (%[_fzz]);\n"
      
      "add $32, %%rdx;\n"
      "add $8, %%rcx;\n"
      "cmp %[_n], %%rcx;\n"
      "jne loop1;\n"
    ://outputs
    ://inputs
    [_px]  "r" (px),
    [_py]  "r" (py),
    [_n]   "r" (n),
    [_pz]  "r" (pz),
    [_tmpx] "r" (tmpx),
    [_tmpy] "r" (tmpy),
    [_tmpz] "r" (tmpz),
    [_softening] "r" (softening)
    //[_fxx]  "r" (fxx),
    //[_fyy] "r" (fyy),
    //[_fzz] "r" (fzz)

    ://clobbers
    "cc", "memory", "rax", "rcx", "rdx", "rbx", "ymm1", "ymm2", "ymm3", "ymm4", "ymm5", "ymm6", "ymm7", "ymm8", "ymm9", "ymm10", "ymm11", "ymm12", "ymm13"
      );

    fx = fxx[0] + fxx[1] + fxx[2] + fxx[3] + fxx[4] + fxx[5] + fxx[6] +fxx[7]; 
    fy = fyy[0] + fyy[1] + fyy[2] + fyy[3] + fyy[4] + fyy[5] + fyy[6] +fyy[7];
    fz = fzz[0] + fzz[1] + fzz[2] + fzz[3] + fzz[4] + fzz[5] + fzz[6] +fzz[7];
    p.vx[i]   += dt * fx;
    p.vy[i]   += dt * fy;
    p.vz[i]   += dt * fz;
    }
  
  
  /*printf("ii %llu\n", i);
  printf("jj %llu\n", j);
  printf("kk %llu\n", k);*/
  //printf("ll %llu\n", l);
  
  
//3 floating-point operations
//scale vector 
//elemanlari sekizli isliyoruz
/*will merge it later */
  __asm__ volatile( /*this part does its job flawlessly */
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
    "inc (%[_i]);\n" //loop works
    "add $8, %%rcx;\n" //add 8 to rcx and store it in rcx => because we'll load our elements by packets of 8 floats => avx2 registers are 32bytes
    "cmp %[_n], %%rcx;\n"
    "jne loop;\n"
   
    ://outputs

    ://inputs
    [_i] "r" (&i),
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
  //printf("%llu", i); //loop works 
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
