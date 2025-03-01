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


double *d_allstrain_r   = nullptr;
double *d_allstrain_i   = nullptr;


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
int copy_LUT_array(int *hloc, int size){
  // Allocates the eta arrays
  cudaMalloc(&d_LUT, size*sizeof(int));
  cudaMemcpy(d_LUT, hloc, size*sizeof(int), cudaMemcpyHostToDevice);
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





__global__ void vanikernel_allstrains_allmodes(int nspec, int ngll_per_loop, int nelem_in_block, int ngll,
                                               int *d_LUT, double *d_allstrain_r, double *d_allstrain_i){ 
  // Kernel is launched with dimensions: 
  // <<< dim3(nblocks_for_all_elems, nn1_total, 81), dim3(nelem_in_block, ngll_in_block, 1) >>>
  int startelem, myspec, ispec, endelem, igllstart, igllend, imode, m1, m2, p, q, utripos;
  double cont_r, cont_i; 

  extern __shared__ double shared_mem[];  // Dynamic shared memory

  // Assign pointers to shared memory
  double* sreal1 = shared_mem;
  double* sreal2 = sreal1 + (125 * 32);
  double* simag1 = sreal2 + (125 * 32);
  double* simag2 = simag1 + (125 * 32);
  

  // 
  startelem = blockIdx.x * nelem_in_block;
  myspec    = threadIdx.x;                   // 0 - 31;
  ispec     = startelem + myspec;


  if(ispec >= nspec)return; // safeguard note >= because ispec goes to 0, nspec-1 


  endelem = startelem + 31;
  if (endelem >= nspec )endelem = nspec -1 ;


  // Each warp is responsible for ngll_per_loop gll points
  // Will go from igllstart to igllstart + ngll_per_loop
  igllstart = threadIdx.y * ngll_per_loop ;
  igllend   = igllstart + ngll_per_loop - 1;

  if(igllend >= ngll*ngll*ngll)igllend = ngll*ngll*ngll; //safeguard


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
  for (int i = 0; i < 125*32; ++i){
    //index = ispec + 
    sreal1[i] = 2.0;  // d_allstrain_r[startelem:endelem, :, m1, Vcont[p], imode];
    sreal2[i] = 3.0;  // d_allstrain_r[startelem:endelem, :, m2, Vcont[p], imode];

    simag1[i] = 4.0;  // d_allstrain_i[startelem:endelem, :, m1, Vcont[q], imode];
    simag2[i] = 5.0;  // d_allstrain_i[startelem:endelem, :, m2, Vcont[q], imode];
  } 

  cont_r = 0.0;
  cont_i = 0.0;


   for (int igll = igllstart; igll < igllend; ++igll){
            // Real part 
            cont_r = cont_r  +  (sreal1[32*igll + myspec] * sreal2[32*igll + myspec]  +  
                                 simag1[32*igll + myspec] * simag2[32*igll + myspec]); // * & 
                                 // NEED TO ADD THIS d_Cxyz_g[10] * swgll(myspec, igll)

            // Imag part 
            cont_i = cont_i  +  (sreal1[32*igll + myspec] * simag2[32*igll + myspec]  -  
                                 sreal2[32*igll + myspec] * simag1[32*igll + myspec]);  //* & 
                                 //d_Cxyz_g(myspec, igll, Vcont(p), Vcont(q)) * swgll(myspec, igll)
   } 

}




extern "C" {
int launch_vanikernel(int ngll, int nspec, int nn1_total, int nn1max){
  
  // Local variables
  cudaFuncAttributes attrib;
  cudaError_t ierr; 
  size_t sharedMemSize = 160 * 1024;

  int nelem_in_block        = 32;
  int ngll_in_block         = 16;
  int nblocks_for_all_elems = ceil(float(nspec)/float(nelem_in_block));
  int ngll_per_loop         = ceil(125.0/float(ngll_in_block));


  // Max out the dynamic shared memory for A100
  cudaFuncSetAttribute(vanikernel_allstrains_allmodes, cudaFuncAttributeMaxDynamicSharedMemorySize, 163840);
  cudaFuncGetAttributes(&attrib, vanikernel_allstrains_allmodes);

  cudaFuncSetAttribute(vanikernel_allstrains_allmodes, cudaFuncAttributePreferredSharedMemoryCarveout, 100);


  printf("Dimensions for calculation: \n");
  printf("    - Nspec                     :   %i\n\n", nspec);
  printf("Grid dimensions: \n");
  printf("    - To cover all elements     :   %i\n ", nblocks_for_all_elems);
  printf("   - To cover all nn1_total    :   %i\n ", nn1_total);
  printf("   - To cover all contractions :   %i\n\n ", 9*9);
  printf("Block dimensions: \n");
  printf("    - Elements in a block       :   %i\n", nelem_in_block);
  printf("    - Ngll covered in the block :   %i\n", ngll_in_block);
  printf("    - Ngll in loop              :   %i\n\n", ngll_per_loop);
  printf("Shared memory: \n");
  printf("maxDynamicSharedSizeBytes       :   %i  \n", attrib.maxDynamicSharedSizeBytes);
  printf("preferredShmemCarveout          :   %i  \n", attrib.preferredShmemCarveout);
  printf("Static(?) shared memory size    :   %zu \n", attrib.sharedSizeBytes);


  if(ngll == 5){
        vanikernel_allstrains_allmodes<<<dim3(nblocks_for_all_elems, nn1_total, 81), 
                                         dim3(nelem_in_block, ngll_in_block, 1),
                                         sharedMemSize>>>
                                         (nspec, ngll_per_loop, nelem_in_block, ngll, d_LUT,
                                         d_allstrain_r, d_allstrain_i);        
  } else {
        printf("Error. hardcoded for ngll = 5\n");
        abort();
  } // ngll=5


  cudaDeviceSynchronize();


  ierr = cudaGetLastError();
  printf("Error code:    :   %i \n", ierr);

  if (ierr != cudaSuccess) {
      printf("CUDA kernel launch error: %s\n", cudaGetErrorString(ierr));
  }

  return 0;
}
}





extern "C" {
int copythisarraytodevice(double *hloc, int size, int varid){

  
  switch (varid) {
    case 1:
      cudaMalloc(&d_xcoord  , size*sizeof(double));
      cudaMemcpy(d_xcoord, hloc, size*sizeof(double), cudaMemcpyHostToDevice);
      printf("Copying xcoord \n");
      break;
    case 2:
      cudaMalloc(&d_ycoord  , size*sizeof(double));
      cudaMemcpy(d_ycoord, hloc, size*sizeof(double), cudaMemcpyHostToDevice);
      printf("Copying ycoord \n");
      break;
    case 3:
      cudaMalloc(&d_zcoord  , size*sizeof(double));
      cudaMemcpy(d_zcoord, hloc, size*sizeof(double), cudaMemcpyHostToDevice);
      printf("Copying zcoord \n");
      break;
    case 4:
      cudaMalloc(&d_m3d  , size*sizeof(double));
      cudaMemcpy(d_m3d, hloc, size*sizeof(double), cudaMemcpyHostToDevice);
      printf("Copying 3D model \n");
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
    // ! r11 = c1 * c2 
    // ! r12 = s1 * c2
    // ! r13 = -s2 

    // ! r21 = -s1 
    // ! r22 = c1 
    // ! r23 = 0.0 

    // ! r31 = c1*s2
    // ! r32 = s1*s2
    // ! r33 = c2

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


  printf("maxDynamicSharedSizeBytes: %i \n", attrib.maxDynamicSharedSizeBytes);
  printf("preferredShmemCarveout   :    %i \n", attrib.preferredShmemCarveout);
  printf("Shared memory size       : %zu bytes \n", attrib.sharedSizeBytes);


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
      printf("Copying xcoord to host \n");
      break;
    case 2:
      cudaMemcpy(hloc, d_ycoord, size*sizeof(double), cudaMemcpyDeviceToHost);
      printf("Copying ycoord to host \n");
      break;
    case 3:
      cudaMemcpy(hloc, d_zcoord, size*sizeof(double), cudaMemcpyDeviceToHost);
      printf("Copying zcoord to host \n");
      break;
    case 4:
      cudaMemcpy(hloc, d_m3d, size*sizeof(double), cudaMemcpyDeviceToHost);
      printf("Copying 3D model to host \n");
      break;
    case 5:
      cudaMemcpy(hloc, d_eta1, size*sizeof(double), cudaMemcpyDeviceToHost);
      printf("Copying eta1 to host \n");
      break;
    case 6:
      cudaMemcpy(hloc, d_eta2, size*sizeof(double), cudaMemcpyDeviceToHost);
      printf("Copying eta2 to host \n");
      break;
    case 7:
      cudaMemcpy(hloc, d_Cxyz, size*sizeof(double), cudaMemcpyDeviceToHost);
      printf("Copying d_Cxyz to host \n");
      break;
    default:
      printf("Invalid VarID entered, must be 1-3 or 21 but was %i \n", varid);
      abort();
  }

  return 0;
}
}




