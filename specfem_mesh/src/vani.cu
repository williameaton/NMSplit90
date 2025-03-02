#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include "device_launch_parameters.h"



//double *h_eta1 = nullptr;
double *d_xcoord = nullptr;
double *d_ycoord = nullptr;
double *d_zcoord = nullptr;
double *d_m3d    = nullptr; // flattened 3D model

double *d_eta1   = nullptr;
double *d_eta2   = nullptr;

double *d_Cxyz   = nullptr;
int    *d_LUT    = nullptr;

double *d_wgll   = nullptr;

double *d_allstrain_r   = nullptr;
double *d_allstrain_i   = nullptr;

double *d_vani_real = nullptr;
double *d_vani_imag = nullptr;

  __constant__ int Vcont[9] = {0, 1, 2, 3, 3, 4, 4, 5, 5};

#define GET_SWGLL(i, j) d_wgll[(i) * nspec + (j)]


extern "C" {
int cudatestfunc(int size){

  return size*2 + 1;
}
}



extern "C" {
int allocate_eta_arrays(int size){
  // Allocates the eta arrays
  cudaMalloc(&d_eta1, size*sizeof(double));
  cudaMalloc(&d_eta2, size*sizeof(double));
  return 0;
}
}


extern "C" {
int allocate_Cxyz_array(int size){
  // Allocates the eta arrays
  cudaMalloc(&d_Cxyz, size*sizeof(double));
  return 0;
}
}


extern "C" {
int allocate_Vani_arrays(int size){
  // Allocates the Vani arrays
  cudaMalloc(&d_vani_real, size*sizeof(double));
  cudaMalloc(&d_vani_imag, size*sizeof(double));
  return 0;
}
}



extern "C" {
int copy_LUT_array(int *hloc, int size){
  // Allocates the eta arrays
  cudaMalloc(&d_LUT, size*sizeof(int));
  cudaMemcpy(d_LUT, hloc, size*sizeof(int), cudaMemcpyHostToDevice);
  return 0;
}
}



extern "C" {
int copy_wgll_array(double *hloc, int size){
  // Allocates the eta arrays
  cudaMalloc(&d_wgll, size*sizeof(double));
  cudaMemcpy(d_wgll, hloc, size*sizeof(double), cudaMemcpyHostToDevice);
  return 0;
}
}


extern "C" {
int copy_allstrains(double *hloc_r, double *hloc_i, int size){
  // Allocates the eta arrays
  cudaMalloc(&d_allstrain_r, size*sizeof(double));
  cudaMemcpy( d_allstrain_r, hloc_r, size*sizeof(double), cudaMemcpyHostToDevice);

  // Allocates the eta arrays
  cudaMalloc(&d_allstrain_i, size*sizeof(double));
  cudaMemcpy( d_allstrain_i, hloc_i, size*sizeof(double), cudaMemcpyHostToDevice);

  return 0;
}
}





__global__ void vanikernel_allstrains_allmodes(int maxtl1, int nspec, int ngll_per_loop, int nelem_in_block, int ngll_in_block, 
                                               int ngll, int maxnn1,
                                               int *d_LUT, double *d_allstrain_r, double *d_allstrain_i, 
                                               double *d_wgll, double *d_Cxyz, double *d_vani_real, double *d_vani_imag){ 
  // Kernel is launched with dimensions: 
  // <<< dim3(nblocks_for_all_elems, nn1_total, 81), dim3(nelem_in_block, ngll_in_block, 1) >>>
  int startelem, myspec, ispec, endelem, igllstart, igllend, imode, m1, m2, p, q, utripos, cxyzindex;
  double cont_r, cont_i; 

  extern __shared__ double shared_mem[];  // Dynamic shared memory
  // Assign pointers to shared memory
  double* sreal1 = shared_mem;
  double* sreal2 = sreal1 + (125 * 32);
  double* simag1 = sreal2 + (125 * 32);
  double* simag2 = simag1 + (125 * 32);
  //double* swgll  = simag2 + (125 * 32);
  
  // Block wise reduction
  double* scont_r = simag2 + (125 * 32);
  double* scont_i = scont_r + ngll_in_block*nelem_in_block;  // Offset for imaginary part

  // 
  startelem = blockIdx.x * nelem_in_block;
  myspec    = threadIdx.x;                   // 0 - 31;
  ispec     = startelem + myspec;


  if(ispec >= nspec)return; // safeguard note >= because ispec goes to 0, nspec-1 


  endelem = startelem + 31;
  if (endelem >= nspec )endelem = nspec -1 ;

  int nloc_el = endelem - startelem + 1;


  // Each warp is responsible for ngll_per_loop gll points
  // Will go from igllstart to igllstart + ngll_per_loop
  igllstart = threadIdx.y * ngll_per_loop ;
  igllend   = igllstart   + ngll_per_loop ;

  if(igllend > ngll*ngll*ngll) igllend = ngll*ngll*ngll; //safeguard




  // Which mode am i and which m1, m2 value am I solving? 
  imode    = d_LUT[blockIdx.y*4    ];  // The number of this mode
  utripos  = d_LUT[blockIdx.y*4 + 1];  // The position in the Vani matrix (upper triangular)

  m1       = d_LUT[blockIdx.y*4 + 2];  // The row of the Vani element
  m2       = d_LUT[blockIdx.y*4 + 3];  // The column of the Vani element


  // Which of the 81 contractions am i solving? 
  p = blockIdx.z/9     ;
  q = blockIdx.z - p*9 ;

 

  // Copy over 125 * 32 points
  // d_allstrain_r order is nspec, ngll,
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    int ictr = 0;
    int myind_1 =  imode*6*maxtl1*125*nspec + Vcont[p]*maxtl1*125*nspec + m1*125*nspec;
    int myind_2 =  imode*6*maxtl1*125*nspec + Vcont[q]*maxtl1*125*nspec + m2*125*nspec;

    for (int i = 0; i < 125; ++i){
      for (int j = startelem; j <= endelem; ++j){

        sreal1[ictr] = d_allstrain_r[myind_1 + (i*nspec) + j]; 
        sreal2[ictr] = d_allstrain_r[myind_2 + (i*nspec) + j]; 

        simag1[ictr] = d_allstrain_i[myind_1 + (i*nspec) + j]; 
        simag2[ictr] = d_allstrain_i[myind_2 + (i*nspec) + j]; 

        ///swgll[ictr] = d_wgll[i*nspec + j];

        ictr = ictr + 1;
      }
    } 
  }


__syncthreads();

cont_r = 0.0;
cont_i = 0.0;


// What cxyz point am i: 
for (int igll = igllstart; igll < igllend; ++igll){

        // If you are in the last block then this wont be 32 

        // Cxyz index: 
        cxyzindex = (Vcont[p]*(nspec*125)*6)  +  (Vcont[q]*nspec*125) + (igll*nspec) + ispec;        

        // Real part 
        cont_r = cont_r  +  (sreal1[nloc_el*igll + myspec] * sreal2[nloc_el*igll + myspec]  +  
                             simag1[nloc_el*igll + myspec] * simag2[nloc_el*igll + myspec]) * 
                             d_Cxyz[cxyzindex] *  d_wgll[igll*nspec + ispec];  //GET_SWGLL(i, j);//swgll[nloc_el*igll + myspec];

        
        // Imag part 
        cont_i = cont_i  +  (sreal1[nloc_el*igll + myspec] * simag2[nloc_el*igll + myspec]  -  
                             sreal2[nloc_el*igll + myspec] * simag1[nloc_el*igll + myspec]) * 
                             d_Cxyz[cxyzindex] *  d_wgll[igll*nspec + ispec];
} 


  // It needs to be this way around because for the final block_x 
  // there are (probably) not 32 elements left 
  int tid = (threadIdx.x * ngll_in_block) + threadIdx.y;  

  scont_r[tid] = cont_r;
  scont_i[tid] = cont_i;
  __syncthreads();



  /// TWO VERSIONS OF THE ATOMIC ADD - 2nd is more parallel and a bit faster
  /// but I think the requirement to sync threads lots slows it down 
  /// its a small reduction so not much difference in the speed 
  /// Currently using 1st version because the 2nd seems to give slightly
  /// wrong answer and dont want to debug it rn. 

  if (threadIdx.x == 0 && threadIdx.y == 0) {
    for (int i = 1; i < ngll_in_block*nloc_el; ++i){
      scont_r[0] += scont_r[i] ;
      scont_i[0] += scont_i[i] ;
    } 
      atomicAdd(&d_vani_real[maxnn1 * imode + utripos], scont_r[0]);
      atomicAdd(&d_vani_imag[maxnn1 * imode + utripos], scont_i[0]);
  }


  // // Perform parallel reduction in shared memory
  //  int total_threads = ngll_in_block * nloc_el;    // Effective number of threads
  // for (int stride = total_threads / 2; stride > 0; stride >>= 1) {
  //     if (tid < stride) {
  //         scont_r[tid] += scont_r[tid + stride];
  //         scont_i[tid] += scont_i[tid + stride];
  //     }
  //     __syncthreads();
  // }
  // // Final atomic add by one thread per block
  // if (tid == 0) {
  //     atomicAdd(&d_vani_real[maxnn1 * imode + utripos], scont_r[0]);
  //     atomicAdd(&d_vani_imag[maxnn1 * imode + utripos], scont_i[0]);
  // }

}




extern "C" {
int launch_vanikernel(int ngll, int nspec, int nn1_total, int maxnn1, int maxtl1){
  
  // Local variables
  cudaFuncAttributes attrib;
  cudaError_t ierr; 
  size_t sharedMemSize = 160 * 1024;

  int nelem_in_block        = 32;
  int ngll_in_block         = 16;
  int nblocks_for_all_elems = ceil(float(nspec)/float(nelem_in_block));
  int ngll_per_loop         = ceil(125.0/float(ngll_in_block));

  // CUDA event timers
  cudaEvent_t start, stop;
  float milliseconds = 0;

  // Create events
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Max out the dynamic shared memory for A100
  cudaFuncSetAttribute(vanikernel_allstrains_allmodes, cudaFuncAttributeMaxDynamicSharedMemorySize, 163840);
  cudaFuncGetAttributes(&attrib, vanikernel_allstrains_allmodes);
  cudaFuncSetAttribute(vanikernel_allstrains_allmodes, cudaFuncAttributePreferredSharedMemoryCarveout, 100);

  if(ngll == 5){

        // Start timing
        cudaEventRecord(start);

        vanikernel_allstrains_allmodes<<<dim3(nblocks_for_all_elems, nn1_total, 81), 
                                         dim3(nelem_in_block, ngll_in_block, 1),
                                         sharedMemSize>>>
                                         (maxtl1, nspec, ngll_per_loop, nelem_in_block, ngll_in_block, ngll, maxnn1, d_LUT,
                                         d_allstrain_r, d_allstrain_i, d_wgll, d_Cxyz, d_vani_real, d_vani_imag);        

        // Stop timing
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milliseconds, start, stop);

  } else {
        printf("Error. hardcoded for ngll = 5\n");
        abort();
  } // ngll=5

  cudaDeviceSynchronize();

  ierr = cudaGetLastError();
  printf("Error code:    :   %i \n", ierr);

  if (ierr != cudaSuccess) {
      printf("CUDA kernel launch error: %s\n", cudaGetErrorString(ierr));
  } else {
      printf("Kernel execution time: %f ms\n", milliseconds);
  }

  // Clean up CUDA events
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  return 0;
}
}





extern "C" {
int copythisarraytodevice(double *hloc, int size, int varid){

  
  switch (varid) {
    case 1:
      cudaMalloc(&d_xcoord  , size*sizeof(double));
      cudaMemcpy(d_xcoord, hloc, size*sizeof(double), cudaMemcpyHostToDevice);

      break;
    case 2:
      cudaMalloc(&d_ycoord  , size*sizeof(double));
      cudaMemcpy(d_ycoord, hloc, size*sizeof(double), cudaMemcpyHostToDevice);

      break;
    case 3:
      cudaMalloc(&d_zcoord  , size*sizeof(double));
      cudaMemcpy(d_zcoord, hloc, size*sizeof(double), cudaMemcpyHostToDevice);

      break;
    case 4:
      cudaMalloc(&d_m3d  , size*sizeof(double));
      cudaMemcpy(d_m3d, hloc, size*sizeof(double), cudaMemcpyHostToDevice);
      break;
    default:
      printf("Invalid VarID entered, must be 1-3 or 21 but was %i \n", varid);
      abort();
  }

  return 0;
}
}






// WHEN WE LAUNCH THIS IT SHOULD BE: 
// <<<dim3(nspec, 1, 1), dim3(32, ngll**3 // 32, 1)>>>
// Ie a warp does 32 gll points so need ~ ngll**3 // 32 = 4 warps for the element
__global__ void project_eta_to_gll(int npoints, int ngll, int nspec,
                                   double *d_xcoord, double *d_ycoord, double *d_zcoord,  
                                   double *d_m3d,    double *d_eta1,   double *d_eta2, 
                                   double modelA, double modelC, double modelL, double modelN, double modelF,
                                   double *d_Cxyz){

    // npoints = number of points in 3D model to test
    // ngll = number of GLL in one direction

    // internal variables
    int  mypt, ispec, myigll, myindex;
    double  myx, myy, myz, dist, last_dist, vp, vs, rho, rad2;
    double  myA, myC, myL, myN, myF;
    double  n1, n2, c1, c2, s1, s2, r11, r12, r13, r21, r22, r23, r31, r32, r33;

    double Cnat[6][6], Q[6][6], Crot[6][6];

    
    myigll = (threadIdx.y)*32 + (threadIdx.x);
    ispec  = blockIdx.x * ngll*ngll*ngll;

    if(myigll >= 125) return;

    // These arrays are flattened as follows - i then j, then k, then ispec
    // hence there are 125 values in a row for each ispec
    myx = d_xcoord[ispec + myigll];
    myy = d_ycoord[ispec + myigll];
    myz = d_zcoord[ispec + myigll];

    // For each possible point let us evaluate the distance
    last_dist = 1000.0;
    for (int ipt = 0; ipt < npoints; ++ipt){
        dist = std::pow((std::pow(d_m3d[ipt*5    ] - myx, 2.0) + 
                         std::pow(d_m3d[ipt*5 + 1] - myy, 2.0) + 
                         std::pow(d_m3d[ipt*5 + 2] - myz, 2.0)), 
                         0.5);
        // If distance is smaller than last then save this index 
        if (dist < last_dist) {
            last_dist = dist;
            mypt      = ipt ;
        };
    };


    // // Ultimately we set this points eta 1 and eta 2 to the final mypt 
    // // that survived: 
    n1 = d_m3d[mypt*5  + 3];
    n2 = d_m3d[mypt*5  + 4];
    
    // // Store these for benchmarking
    //d_eta1[ispec + myigll] = n1;
    //d_eta2[ispec + myigll] = n2;

    // // Get the PREM related values for ACLNF at this point: 
    // r^2 normalised
    rad2 = std::pow(myx, 2.0) + std::pow(myy, 2.0) + std::pow(myz, 2.0);

    // // Originally g/cm^3 --> kg/m^3 -> nondimensionalised
    rho = (13.088500000 - 8.838100000*rad2)*0.1813466804490;   // 1000.d0/RHOAV
    // Originally km/s --> m/s  
    vp  = (11.262200000 - 6.364000000*rad2)*0.1459938230354;   // 1000.d0/SCALE_V
    vs  = (3.667800000 - 4.447500000*rad2)*0.1459938230354;    // 1000.d0/SCALE_V

    myA = rho * vp * vp * modelA;
    myC = rho * vp * vp * modelC;
    myL = rho * vs * vs * modelL;
    myN = rho * vs * vs * modelN;
    // eta  = 1 for PREM core and F = eta * (A - 2L)
    myF = ((rho * vp * vp) - 2.0*(rho * vs * vs)) * modelF;


    // ! Create natural Stiffness matrix:
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
          Cnat[i][j] = 0.0;
      }
    }

    Cnat[0][0] = myA;
    Cnat[1][1] = myA;
    Cnat[2][2] = myC;
    Cnat[3][3] = myL;
    Cnat[4][4] = myL;
    Cnat[5][5] = myN;


    Cnat[0][1] = myA - 2.0 * myN;
    Cnat[1][0] = myA - 2.0 * myN;

    Cnat[0][2] = myF;
    Cnat[2][0] = myF;
    Cnat[1][2] = myF;
    Cnat[2][1] = myF;


    //Create bond matrix: 
    c1 = cos(n1);
    c2 = cos(n2);
    s1 = sin(n1);
    s2 = sin(n2);

    // ! Eqn 5 of Brett 2024
    r11 = c1 * c2;
    r12 = - s1;
    r13 = c1*s2;
    
    r21 = s1*c2;
    r22 = c1;
    r23 = s1*s2;

    r31 = -s2;
    r32 = 0.0;
    r33 = c2;


    // ! Eqn 8 of Brett 2024
    // r11 = c1 * c2 ;
    // r12 = s1 * c2;
    // r13 = -s2; 

    // r21 = -s1 ;
    // r22 = c1 ;
    // r23 = 0.0 ;

    // r31 = c1*s2;
    // r32 = s1*s2;
    // r33 = c2;

    Q[0][0] = r11 * r11; 
    Q[1][0] = r21 * r21; 
    Q[2][0] = r31 * r31;
    Q[3][0] = r21 * r31;
    Q[4][0] = r31 * r11;
    Q[5][0] = r11 * r21;

    Q[0][1] = r12 * r12;
    Q[1][1] = r22 * r22;
    Q[2][1] = r32 * r32;
    Q[3][1] = r22 * r32;
    Q[4][1] = r32 * r12;
    Q[5][1] = r12 * r22;

    Q[0][2] = r13 * r13;
    Q[1][2] = r23 * r23;
    Q[2][2] = r33 * r33;
    Q[3][2] = r23 * r33;
    Q[4][2] = r33 * r13;
    Q[5][2] = r13 * r23;

    Q[0][3] = 2.0 * r12 * r13;
    Q[1][3] = 2.0 * r22 * r23;
    Q[2][3] = 2.0 * r32 * r33;
    Q[3][3] = r22 * r33 + r32 * r23;
    Q[4][3] = r12 * r33 + r13 * r32;
    Q[5][3] = r12 * r23 + r13 * r22;

    Q[0][4] = 2.0 * r11 * r13;
    Q[1][4] = 2.0 * r21 * r23;
    Q[2][4] = 2.0 * r33 * r31;
    Q[3][4] = r23 * r31 + r21 * r33;
    Q[4][4] = r33 * r11 + r13 * r31;
    Q[5][4] = r13 * r21 + r23 * r11;

    Q[0][5] = 2.0 * r11 * r12;
    Q[1][5] = 2.0 * r21 * r22;
    Q[2][5] = 2.0 * r31 * r32;
    Q[3][5] = r21 * r32 + r31 * r22;
    Q[4][5] = r31 * r12 + r11 * r32;
    Q[5][5] = r11 * r22 + r21 * r12;


    // ! Compute Bond * C * Bond^T
    // ! C'_il =  B_ij C_jk B_lk  note transpose on last matrix
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
          Crot[i][j] = 0.0;
      }
    }

    for (int i = 0; i < 6; ++i){
      for (int l = 0; l < 6; ++l){
        for (int j = 0; j < 6; ++j){
          for (int k = 0; k < 6; ++k){
            Crot[i][l] = Crot[i][l] +  (Q[i][j] * Cnat[j][k] * Q[l][k]);
          } 
        } 
      }
    } 

    // order of operations: 
    // nspec, ngll, col, row 
    for (int i = 0; i < 6; ++i){
      for (int l = 0; l < 6; ++l){
            myindex = (i * (nspec*(ngll*ngll*ngll))*6)  +  (l*nspec*ngll*ngll*ngll) + (myigll*nspec) + blockIdx.x;

    
        d_Cxyz[myindex] = Crot[i][l];
      }
    }
} 


extern "C" {
  // Function that launches the CUDA cpp project to eta kernel
 int cpp_project_eta_to_gll(int model_npoints, int nspec, int ngll, double modelv[5]){

  // local variables: 
  int ngll_per_warp, ngllwarps_per_block;
  
  cudaFuncAttributes attrib;

  ngll_per_warp       = 32;
  ngllwarps_per_block = ceil(float(ngll*ngll*ngll)/float(ngll_per_warp)); 




  // Allocate 160 kB of dynamic shared memory
  // To go over the 48 kB default it needs to be dynamic
  cudaFuncSetAttribute(project_eta_to_gll, cudaFuncAttributeMaxDynamicSharedMemorySize, 163840);
  cudaFuncGetAttributes(&attrib, project_eta_to_gll);


  // printf("maxDynamicSharedSizeBytes: %i \n", attrib.maxDynamicSharedSizeBytes);
  // printf("preferredShmemCarveout   :    %i \n", attrib.preferredShmemCarveout);
  // printf("Shared memory size       : %zu bytes \n", attrib.sharedSizeBytes);


  project_eta_to_gll<<<dim3(nspec, 1, 1), 
                       dim3(ngll_per_warp, ngllwarps_per_block, 1)>>>
                       (model_npoints, ngll, nspec, d_xcoord, d_ycoord, d_zcoord, d_m3d, 
                       d_eta1, d_eta2, modelv[0], modelv[1], modelv[2], modelv[3], modelv[4],
                       d_Cxyz);
  cudaDeviceSynchronize();
  return 0;
 }
}









extern "C" {
int copyfromdevice(double *hloc, int size, int varid){
  
  switch (varid) {
    case 1:
      cudaMemcpy(hloc, d_xcoord, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    case 2:
      cudaMemcpy(hloc, d_ycoord, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    case 3:
      cudaMemcpy(hloc, d_zcoord, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    case 4:
      cudaMemcpy(hloc, d_m3d, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    case 5:
      cudaMemcpy(hloc, d_eta1, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    case 6:
      cudaMemcpy(hloc, d_eta2, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    case 7:
      cudaMemcpy(hloc, d_Cxyz, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    case 8:
      cudaMemcpy(hloc, d_vani_real, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    case 9:
      cudaMemcpy(hloc, d_vani_imag, size*sizeof(double), cudaMemcpyDeviceToHost);
      break;
    default:
      printf("Invalid VarID entered, must be 1-3 or 21 but was %i \n", varid);
      abort();
  }

  return 0;
}
}




