#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include "device_launch_parameters.h"
#include "driver_types.h"
#include "precisioncpp.h"
#include <nvml.h>

// Debugging: 
__device__ int d_errorFlag = 0;



cudaGraph_t Vanigraph;
cudaGraphExec_t VanigraphExec;
cudaStream_t Vanistream;
cudaGraphNode_t       VanikernelNode;
cudaKernelNodeParams  kernelParams = {};
bool graphConstructed = false;  // Flag to track whether the graph is created or not


CPPCUSTOM_REAL *d_xcoord = nullptr;
CPPCUSTOM_REAL *d_ycoord = nullptr;
CPPCUSTOM_REAL *d_zcoord = nullptr;
CPPCUSTOM_REAL *d_m3d    = nullptr; // flattened 3D model


CPPCUSTOM_REAL *d_Cxyz   = nullptr;
int            *d_LUT    = nullptr;

CPPCUSTOM_REAL *d_wgll   = nullptr;

CPPCUSTOM_REAL *d_allstrain_r   = nullptr;
CPPCUSTOM_REAL *d_allstrain_i   = nullptr;

CPPCUSTOM_REAL *d_vani_real = nullptr;
CPPCUSTOM_REAL *d_vani_imag = nullptr;

// There are 81 contractions to be done but only 
// 36 of them are unique. I.e. originally a 9 x 9
// where vcont = [0, 1, 2, 3, 3, 4, 4, 5, 5]
// Here I list [p1, q1, occurrences, p1, q1, occurrences ...]
// for the 36 unique combinations
__constant__ int Vcont[108] = {0, 0, 1, 0, 1, 1, 0, 2, 1, 0, 3, 2, 0, 4, 2, 
                               0, 5, 2, 1, 0, 1, 1, 1, 1, 1, 2, 1, 1, 3, 2, 
                               1, 4, 2, 1, 5, 2, 2, 0, 1, 2, 1, 1, 2, 2, 1, 
                               2, 3, 2, 2, 4, 2, 2, 5, 2, 3, 0, 2, 3, 1, 2, 
                               3, 2, 2, 3, 3, 4, 3, 4, 4, 3, 5, 4, 4, 0, 2, 
                               4, 1, 2, 4, 2, 2, 4, 3, 4, 4, 4, 4, 4, 5, 4, 
                               5, 0, 2, 5, 1, 2, 5, 2, 2, 5, 3, 4, 5, 4, 4, 
                               5, 5, 4};


__constant__ int Lvals[27] = {8, 7,   6, 5,  5,  5, 5,  5,  4,  4,   4,  3,  3, 3,  3, 3, 3, 2 , 2, 2,  2, 2,  1,  1, 1,  1,  1};

// extern "C" {
// int allocate_eta_arrays(int size){
//   // Allocates the eta arrays
//   int ierr = cudaMalloc(&d_eta1, size*sizeof(CPPCUSTOM_REAL));
//   if(ierr != 0)return -1 ;

//   ierr = cudaMalloc(&d_eta2, size*sizeof(CPPCUSTOM_REAL));
//   return ierr;
// }
// }



extern "C" {
void check_gpu_utilization(int gpu_id) {
    nvmlDevice_t device;
    nvmlUtilization_t utilization;
    
    nvmlInit();
    nvmlDeviceGetHandleByIndex(gpu_id, &device);
    nvmlDeviceGetUtilizationRates(device, &utilization);
    
    printf("GPU %d: Compute=%d%%, Memory=%d%%\n", 
           gpu_id, utilization.gpu, utilization.memory);
    
    nvmlShutdown();
}
}




extern "C" 
int fh_cuda_preinit()
{
    int cnt = 0;
    if (cudaGetDeviceCount(&cnt) != cudaSuccess || cnt <= 0) return -1;
    // Use device 0 under CUDA_VISIBLE_DEVICES; mpirun/Slurm can remap per rank.
    cudaSetDevice(0);
    // Force context creation so UCX finds a valid context during MPI_Init
    cudaFree(0);
    return 0;
}



extern "C" {
int get_cpp_precision(){
  return sizeof(CPPCUSTOM_REAL);
}
}


extern "C" {
int allocate_Cxyz_array(int size){
  // Allocates the eta arrays
  int ierr = cudaMalloc(&d_Cxyz, size*sizeof(CPPCUSTOM_REAL));
  return ierr;
}
}


extern "C" {
int allocate_Vani_arrays(int size){
  // Allocates the Vani arrays
  int ierr = cudaMalloc(&d_vani_real, size*sizeof(CPPCUSTOM_REAL));
  if(ierr != 0)return -1 ;

  ierr = cudaMalloc(&d_vani_imag, size*sizeof(CPPCUSTOM_REAL));
  return ierr;
}
}


extern "C" {
int allocate_M3D_array(int size){
  // Allocates the Vani arrays
  int ierr = cudaMalloc(&d_m3d  , size*sizeof(CPPCUSTOM_REAL));
  return ierr;
}
}





extern "C" {
int copy_M3D_array(CPPCUSTOM_REAL *hloc, int size){
  // Allocates the eta arrays
  int ierr = cudaMemcpy(d_m3d, hloc, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);

  return ierr;
  }
}





extern "C" {
int copy_LUT_array(int *hloc, int size){
  // Allocates the eta arrays
  int ierr = cudaMalloc(&d_LUT, size*sizeof(int));
  if(ierr != 0)return -1 ;

  ierr = cudaMemcpy(d_LUT, hloc, size*sizeof(int), cudaMemcpyHostToDevice);
  return ierr;
}
}



extern "C" {
int copy_wgll_array(CPPCUSTOM_REAL *hloc, int size){
  // Allocates the eta arrays
  int ierr = cudaMalloc(&d_wgll, size*sizeof(CPPCUSTOM_REAL));
  if (ierr !=0 ){
    printf("Error allocating d_wgll...");
    return -1 ;
  }
  ierr = cudaMemcpy(d_wgll, hloc, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);
    if (ierr !=0 ){
    printf("Error copying d_wgll...");
    return -1 ;
  }
  return 0;
}
}




extern "C" {
int copy_allstrains(CPPCUSTOM_REAL *hloc_r, CPPCUSTOM_REAL *hloc_i, int64_t size){
  // Original function to copy all strains from CPU to device 
  // For PT we use the rank0 version with IPC 

  int ierr;

  ierr = cudaMalloc(&d_allstrain_r, size*sizeof(CPPCUSTOM_REAL));
  if (ierr !=0 ){
    printf("Error allocating strain real...");
    return -1 ;
  }

  ierr = cudaMemcpy( d_allstrain_r, hloc_r, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);
  if (ierr !=0 ){
    printf("Error copying strain real...");
    return -1 ;
  }

  // Allocates the eta arrays
  ierr = cudaMalloc(&d_allstrain_i, size*sizeof(CPPCUSTOM_REAL));
    if (ierr !=0 ){
    printf("Error allocating strain imag...");
    return -1 ;
  }
  ierr = cudaMemcpy( d_allstrain_i, hloc_i, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);
  if (ierr !=0 ){
    printf("Error copying strain imag...");
    return -1 ;
  }
  return 0;
}
}




extern "C" {
int copy_allstrains_fromrank0(CPPCUSTOM_REAL *hloc_r, CPPCUSTOM_REAL *hloc_i, int64_t size, 
                              int myf90rank, char* ipc_handle_buffer){
  int ierr, gpu_id;

  // Safety check the device 
    ierr = cudaGetDevice(&gpu_id);
    if (ierr != cudaSuccess) {
      printf("Error: Unable to get current device!\n");
      return -1;
    }

  if (myf90rank != gpu_id){
      printf("ERROR: myf90rank (%i) != gpu_id (%i) ", myf90rank, gpu_id);
      return -1;
  }

  // Allocate real strain memory 
  ierr = cudaMalloc(&d_allstrain_r, size*sizeof(CPPCUSTOM_REAL));
  if (ierr !=0 ){
        printf("Error allocating strain real on GPU %d\n", gpu_id);
    return -1 ;
  }

  // Copy over the real strain eigenfunctions
  ierr = cudaMemcpy(d_allstrain_r, hloc_r, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);
  if (ierr !=0 ){
        printf("Error allocating strain real on GPU %d\n", gpu_id);
    return -1 ;
  }

  // Allocates the imaginary strains 
  ierr = cudaMalloc(&d_allstrain_i, size*sizeof(CPPCUSTOM_REAL));
    if (ierr !=0 ){
      printf("Error allocating strain real on GPU %d\n", gpu_id);
    return -1 ;
  }
  
  // Copy imaginary strains to the device
  ierr = cudaMemcpy( d_allstrain_i, hloc_i, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);
  if (ierr !=0 ){
    printf("Error copying strain imag...");
    return -1 ;
  }


  // Create IPC handles: 
  cudaIpcMemHandle_t handle_r, handle_i;

  ierr = cudaIpcGetMemHandle(&handle_r, d_allstrain_r);
  if (ierr != 0){
      printf("Error getting IPC handle for strain_r on GPU %d\n", gpu_id);
      return -1;
  }

  ierr = cudaIpcGetMemHandle(&handle_i, d_allstrain_i);
  if (ierr != 0) {
      printf("Error getting IPC handle for strain_r on GPU %d\n", gpu_id);
      return -1;
  }


  // Pack both handles into output buffer (128 bytes total: 64 + 64)
  memcpy(ipc_handle_buffer,      &handle_r, sizeof(cudaIpcMemHandle_t));
  memcpy(ipc_handle_buffer + 64, &handle_i, sizeof(cudaIpcMemHandle_t));


  //printf("My GPU: %d just set strains from rank %d \n", gpu_id, myf90rank);


  return 0;
}
}


extern "C" {
int copy_allstrains_higherranks(int myf90rank, char* ipc_handle_in) {
    int ierr, gpu_id;
    // Assigns the correct pointer 
    
    // Safety check the device 
    ierr = cudaGetDevice(&gpu_id);
    if (ierr != cudaSuccess) {
      printf("Error: Unable to get current device!\n");
      return -1;
    }

    if (myf90rank != gpu_id){
      printf("ERROR in copy_allstrains_higherranks: myf90rank (%i) != gpu_id (%i) ", myf90rank, gpu_id);
      return -1;
    }
    
    // Unpack handles from the buffer
    cudaIpcMemHandle_t handle_r, handle_i;
    memcpy(&handle_r, ipc_handle_in, sizeof(cudaIpcMemHandle_t));
    memcpy(&handle_i, ipc_handle_in + 64, sizeof(cudaIpcMemHandle_t));
    

    // Open shared memory for real strains
    ierr = cudaIpcOpenMemHandle((void**)&d_allstrain_r, handle_r, 
                                cudaIpcMemLazyEnablePeerAccess);
    if (ierr != 0) {
        printf("Error opening IPC handle for strain_r on GPU %d\n", gpu_id);
        return -1;
    }
    
    // Open shared memory for imaginary strains
    ierr = cudaIpcOpenMemHandle((void**)&d_allstrain_i, handle_i, 
                                cudaIpcMemLazyEnablePeerAccess);
    if (ierr != 0) {
        printf("Error opening IPC handle for strain_i on GPU %d\n", gpu_id);
        return -1;
    }
    
    //printf("My GPU: %d accessed from rank %d \n", gpu_id, myf90rank);

    return 0;
}
}


extern "C"{
int assign_proc_to_device(int nprocs, int myrank){
    int err, devcount, mydevice, mydev , nsets_per_gpu;

    // Number of devices
    err = cudaGetDeviceCount(&devcount);
    if (err != cudaSuccess) {
      printf("Error: Unable to get device count!\n");
      return -1;
    }

    //devcount = 1;

    if(myrank == 1){
      printf("Number of GPU devices:  %i\n", devcount);
    }
   if(myrank == 1){
      printf("Number of procs and my rank:  %i  %i\n", nprocs, myrank);
    }
    nsets_per_gpu = (nprocs + devcount - 1) / devcount;  

   if(myrank == 1){
      printf("Number of nsets_per_gpu:  %i\n", nsets_per_gpu);
    }

    mydev = myrank / nsets_per_gpu;  // This mimics FLOOR(myrank / nsets_per_gpu)

    if(myrank == 1){
      printf("myrank and dev:  %i %i\n", myrank, mydev);
    }
    

    // ensure we don't go out of bounds
    if(mydev >= devcount){
      printf("Error mydev is >= devcount: mydev: %i,   devcount: %i\n", mydev, devcount);
      return -1; 
    } 
    if(mydev < 0) {
      printf("Error mydev is < 0: %i\n", mydev);
      return -1;
    }

    err = cudaSetDevice(mydev);
    if (err != cudaSuccess) {
      printf("Error: Unable to set device!\n");
      return -1 ;
    }

    err = cudaGetDevice(&mydevice);
    if (err != cudaSuccess) {
      printf("Error: Unable to get current device!\n");
      return -1;
    }

    printf("myrank = %i -- on device %i\n", myrank, mydevice);
    return 0;
}
}




extern "C"{
int force_proc_to_device(int nprocs, int myrank){
    int err, devcount, mydevice, mydev;

    // A less cutesy setup when we know which procs we want where
    err = cudaGetDeviceCount(&devcount);
    if (err != cudaSuccess) {
      printf("Error: Unable to get device count!\n");
      return -1;
    }

    if(nprocs != devcount){
      printf("Error nprocs != device count: procs %i devices %i", nprocs, devcount );
      return -1; 
    } 

    mydev = myrank;

    // ensure we don't go out of bounds
    if(mydev >= devcount){
      printf("Error mydev is >= devcount: mydev: %i,   devcount: %i\n", mydev, devcount);
      return -1; 
    } 
    if(mydev < 0) {
      printf("Error mydev is < 0: %i\n", mydev);
      return -1;
    }

    err = cudaSetDevice(mydev);
    if (err != cudaSuccess) {
      printf("Error: Unable to set device!\n");
      return -1 ;
    }

    err = cudaGetDevice(&mydevice);
    if (err != cudaSuccess) {
      printf("Error: Unable to get current device!\n");
      return -1;
    }
  
    //printf("myrank = %i -- on device %i\n", myrank, mydevice);
    return 0;
}
}




__global__ void vanikernel_allstrains_allmodes(int maxtl1, int nspec, int ngll_per_loop, int nelem_in_block, int ngll_in_block, 
                                               int ngll, int maxnn1, int elemperthread,
                                               int * __restrict__ d_LUT, 
                                               CPPCUSTOM_REAL * __restrict__ d_allstrain_r, 
                                               CPPCUSTOM_REAL * __restrict__ d_allstrain_i, 
                                               CPPCUSTOM_REAL * __restrict__ d_Cxyz,
                                               CPPCUSTOM_REAL * __restrict__ d_vani_real, 
                                               CPPCUSTOM_REAL * __restrict__ d_vani_imag){ 
  // Kernel is launched with dimensions: 
  // <<< dim3(nblocks_for_all_elems, nn1_total, 81), dim3(nelem_in_block, ngll_in_block, 1) >>>
  int startelem, myspec, ispec, endelem, igllstart, igllend, imode, 
      p, utripos, endispec, row, col, lval;

  CPPCUSTOM_REAL cont_r_m, cont_i_m, cont_r_p, cont_i_p; 

  extern __shared__ CPPCUSTOM_REAL shared_mem[];  // Dynamic shared memory
  
  
  // Block wise reduction
  CPPCUSTOM_REAL* scont_r_m = shared_mem;
  CPPCUSTOM_REAL* scont_i_m = scont_r_m + ngll_in_block*32;  

  CPPCUSTOM_REAL* scont_r_p = scont_i_m + ngll_in_block*32;  
  CPPCUSTOM_REAL* scont_i_p = scont_r_p + ngll_in_block*32;  


  startelem = blockIdx.x  * nelem_in_block;  // First element in the block
  myspec    = threadIdx.x * elemperthread;   // threadid * number of elements in block (0-63 )              
  ispec     = startelem + myspec;            // global element number
  endispec  = ispec + elemperthread -1 ;     // End element for this thread
 
  // for individual thread
  if(ispec >= nspec)return;                // if thread is out of the nspec
  if(endispec >= nspec)endispec= nspec-1;  // if thread would go over 


  endelem = startelem + nelem_in_block - 1;
  if (endelem >= nspec )endelem = nspec -1 ;
  int nloc_el = endelem - startelem + 1;

  // Each warp is responsible for ngll_per_loop gll points
  // Will go from igllstart to igllstart + ngll_per_loop
  igllstart = threadIdx.y * ngll_per_loop ;
  igllend   = igllstart   + ngll_per_loop ;

  if(igllend > ngll*ngll*ngll) igllend = ngll*ngll*ngll; //safeguard
  
  // Which mode am i and which m1, m2 value am I solving? 
  int4 lut_values = __ldg(reinterpret_cast<const int4*>(&d_LUT[blockIdx.y * 4]));
  imode   = lut_values.x;  // The number of this mode
  utripos = lut_values.y;  // Index of the 
  row     = lut_values.z;  // The row of the Vani element
  col     = lut_values.w;  // The column of the Vani element

  lval    = Lvals[imode];


  // reflected (positive m1) position
  int ltripos = utripos + 2*(lval-row)*(lval+1) - (lval-row)*(lval-row+1)/2;



  // if(blockIdx.x == 0 && blockIdx.z == 0 && threadIdx.x == 0 && threadIdx.y == 0 && imode==0){
  // printf("%i %i %i %i %i %i \n", blockIdx.y, lval, row, col, utripos, ltripos);
  // }


  // -1 ** |m1| ie the positive m1 
  int m1 = row-lval;

  //float sign = ((lval-row) % 2 == 0) ? 1.0f : -1.0f;


  // Which of the 81 contractions am i solving? 
  // Any time that we use vp or vq we use it as 
  // vp * 125 * nspec 
  // Except in the if statement below but in that case this will still work
  // so lets premultiply: 
  // There are 81 contractions but only 36 are unique so we launch 
  // Block dim z of 36 and then work out the occurrences
  p = blockIdx.z;
  int VP125nspec  = Vcont[3*p]    * 125*nspec;
  int VQ125nspec  = Vcont[3*p +1] * 125*nspec;
  int occurrences = Vcont[3*p +2];
  float occurrences_float = static_cast<float>(occurrences);  // Convert to float

  cont_r_m = 0.0;  // negative m2
  cont_i_m = 0.0;  // negative m2

  cont_r_p = 0.0;  // positive m2
  cont_i_p = 0.0;  // positive m2

  int myind_1   =  imode*6*maxtl1*125*nspec + VP125nspec*maxtl1 + row*125*nspec  + ispec;
  int myind_2   =  imode*6*maxtl1*125*nspec + VQ125nspec*maxtl1 + col*125*nspec  + ispec;
  int startcxyz = (VP125nspec*6)+  VQ125nspec + ispec ;


  for (int igll = igllstart; igll < igllend; ++igll){
    // m1 <= 0 
    float4  real_1     = __ldg(reinterpret_cast<const float4*>(&d_allstrain_r[myind_1 + (igll * nspec)]));
    float4  imag_1     = __ldg(reinterpret_cast<const float4*>(&d_allstrain_i[myind_1 + (igll * nspec)]));

    // m2 >= 0 
    float4  real_2     = __ldg(reinterpret_cast<const float4*>(&d_allstrain_r[myind_2 + (igll * nspec)]));
    float4  imag_2     = __ldg(reinterpret_cast<const float4*>(&d_allstrain_i[myind_2 + (igll * nspec)]));
    float4  cxyz       = __ldg(reinterpret_cast<const float4*>(&d_Cxyz[startcxyz  + (igll*nspec)]));


    // Case for upper right corner -- m1 <= 0 and m2 >= 0
    cont_r_m  +=  (real_1.x * real_2.x  + imag_1.x * imag_2.x) * cxyz.x +
                  (real_1.y * real_2.y  + imag_1.y * imag_2.y) * cxyz.y +
                  (real_1.z * real_2.z  + imag_1.z * imag_2.z) * cxyz.z +  
                  (real_1.w * real_2.w  + imag_1.w * imag_2.w) * cxyz.w;  

    cont_i_m  +=  (real_1.x * imag_2.x  - real_2.x * imag_1.x) * cxyz.x +
                  (real_1.y * imag_2.y  - real_2.y * imag_1.y) * cxyz.y +
                  (real_1.z * imag_2.z  - real_2.z * imag_1.z) * cxyz.z +
                  (real_1.w * imag_2.w  - real_2.w * imag_1.w) * cxyz.w;


  // Now we can use the fact that 
  // E_{-m} = (-1)**m E_m conjugate 
  // We can therefore compute the equivalent for m1 > 0 as well

      cont_r_p +=  (real_1.x * real_2.x  - imag_1.x * imag_2.x) * cxyz.x +
                  (real_1.y * real_2.y  - imag_1.y * imag_2.y) * cxyz.y +
                  (real_1.z * real_2.z  - imag_1.z * imag_2.z) * cxyz.z +  
                  (real_1.w * real_2.w  - imag_1.w * imag_2.w) * cxyz.w;  

      cont_i_p +=  (real_1.x * imag_2.x  + real_2.x * imag_1.x) * cxyz.x +
                  (real_1.y * imag_2.y  + real_2.y * imag_1.y) * cxyz.y +
                  (real_1.z * imag_2.z  + real_2.z * imag_1.z) * cxyz.z +
                  (real_1.w * imag_2.w  + real_2.w * imag_1.w) * cxyz.w;

  } // loop gll  


  if(m1%2!=0){
    cont_r_p = cont_r_p * -1.0;
    cont_i_p = cont_i_p * -1.0;
  }

  // ADD TO THE GLOBAL MATRIX USING A 2-STEP REDUCTIOn 
  // It needs to be this way around because for the final block_x 
  // there are (probably) not 32 elements left 
  int tid = (threadIdx.x * ngll_in_block) + threadIdx.y;  

  scont_r_m[tid] = cont_r_m * occurrences_float;
  scont_i_m[tid] = cont_i_m * occurrences_float;

  scont_r_p[tid] = cont_r_p * occurrences_float;
  scont_i_p[tid] = cont_i_p * occurrences_float;
  
  __syncthreads();

  if (nloc_el < nelem_in_block){ 
    if (tid == 0) {
      for (int i = 1; i < ngll_in_block*nloc_el/elemperthread; ++i){
        scont_r_m[0] += scont_r_m[i] ;
        scont_i_m[0] += scont_i_m[i] ;

        scont_r_p[0] += scont_r_p[i] ;
        scont_i_p[0] += scont_i_p[i] ;

      } 

      // Add value in neg m1, pos m2
      atomicAdd(&d_vani_real[maxnn1 * imode + utripos], scont_r_m[0]);
      atomicAdd(&d_vani_imag[maxnn1 * imode + utripos], scont_i_m[0]);

      // Add reflected value in positive m1, m2
      atomicAdd(&d_vani_real[maxnn1 * imode + ltripos], scont_r_p[0]);
      atomicAdd(&d_vani_imag[maxnn1 * imode + ltripos], scont_i_p[0]);

    }
  } else  {
    // Perform parallel reduction in shared memory
    for (int stride = ngll_in_block * 16; stride > 0; stride >>= 1) {
      if (tid < stride) {
          scont_r_m[tid] += scont_r_m[tid + stride];
          scont_i_m[tid] += scont_i_m[tid + stride];

          scont_r_p[tid] += scont_r_p[tid + stride];
          scont_i_p[tid] += scont_i_p[tid + stride];


      }
      __syncthreads();
    }

    if (tid == 0) {
        int gidx1 = (maxnn1 * imode + utripos);

        atomicAdd(&d_vani_real[gidx1], scont_r_m[0]);
        atomicAdd(&d_vani_imag[gidx1], scont_i_m[0]);
    
    
      if(col -lval >= lval - row && row < lval && col > lval){
        
        int gidx2 = (maxnn1 * imode + ltripos);

        atomicAdd(&d_vani_real[gidx2], scont_r_p[0]);
        atomicAdd(&d_vani_imag[gidx2], scont_i_p[0]);
      }

    }
  } // nloc_el < nelem_in_block


}



extern "C" {
int launch_vanikernel(int ngll, int nspec, int nn1_total, int maxnn1, 
                      int maxtl1, int nmodes, int myrank){
  
  // Local variables
  cudaError_t ierr; 

  cudaEvent_t startEvent, stopEvent;
  float elapsedTime;

  cudaEventCreate(&startEvent);
  cudaEventCreate(&stopEvent);



  // First time we launch the kernel we create a graph which 
  // Stores the launch parameters 
  if (!graphConstructed) {

    int elem_per_thread       = 4 ; // Each thread is responsible for 2 elements
    int nelem_in_block        = 32*elem_per_thread;
    int ngll_in_block         = 4;
    int nblocks_for_all_elems = ceil(float(nspec)/float(nelem_in_block)) ;
    int ngll_per_loop         = ceil(125.0/float(ngll_in_block));
    size_t sharedMemSize      = 4 * sizeof(CPPCUSTOM_REAL) * ngll_in_block * 32;


    cudaStreamCreate(&Vanistream);
    cudaGraphCreate(&Vanigraph, 0);

    void* kernelArgs[] = {
          &maxtl1, &nspec, &ngll_per_loop, &nelem_in_block, &ngll_in_block, &ngll,
          &maxnn1, &elem_per_thread, &d_LUT, &d_allstrain_r, &d_allstrain_i,
          &d_Cxyz, &d_vani_real, &d_vani_imag
      };

    //printf(" Launching: %i x %i x %i blocks \n", nblocks_for_all_elems, nn1_total, 36);

    // Prepare kernel launch parameters
    kernelParams.func           = (void*)vanikernel_allstrains_allmodes;
    kernelParams.gridDim        = dim3(nblocks_for_all_elems, nn1_total, 36);
    kernelParams.blockDim       = dim3(32, ngll_in_block, 1);
    kernelParams.sharedMemBytes = sharedMemSize;
    kernelParams.kernelParams   = kernelArgs;

    //printf("Launch params: (%i %i %i ) x (%i %i %i )", nblocks_for_all_elems, nn1_total, 81, 32, ngll_in_block, 1 );

    // Add the kernel launch to the graph
    cudaGraphAddKernelNode(&VanikernelNode, Vanigraph, nullptr, 0, &kernelParams);
    
    // Instantiate the graph
    cudaGraphInstantiate(&VanigraphExec, Vanigraph, nullptr, nullptr, 0);
    
    graphConstructed = true;  // Set the flag to true once the graph is constructed
  } 

  cudaEventRecord(startEvent, Vanistream);

  // Reset d_vani_real and d_vani_imag to zero before the kernel launch
  // I tried adding this to the graph but hit a brick wall...
  cudaMemsetAsync(d_vani_real, 0, nmodes*maxnn1 * sizeof(CPPCUSTOM_REAL), Vanistream);
  cudaMemsetAsync(d_vani_imag, 0, nmodes*maxnn1 * sizeof(CPPCUSTOM_REAL), Vanistream);

  // Synchronize to ensure memory is cleared before launching the graph
  cudaStreamSynchronize(Vanistream);


  cudaGraphLaunch(VanigraphExec, Vanistream);

  ierr = cudaGetLastError();
    if (ierr != cudaSuccess) {
        printf("CUDA kernel launch error: %s\n", cudaGetErrorString(ierr));
    }


  cudaEventRecord(stopEvent, Vanistream);
  cudaEventSynchronize(stopEvent);

  cudaEventElapsedTime(&elapsedTime, startEvent, stopEvent);
  // if(myrank==0)printf(" Kernel time         : %f ms\n", elapsedTime);
  cudaEventDestroy(startEvent);
  cudaEventDestroy(stopEvent);

  return 0;
}
}



extern "C" {
int copythisarraytodevice(CPPCUSTOM_REAL *hloc, int size, int varid){

  
  switch (varid) {
    case 1:
      cudaMalloc(&d_xcoord  , size*sizeof(CPPCUSTOM_REAL));
      cudaMemcpy(d_xcoord, hloc, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);

      break;
    case 2:
      cudaMalloc(&d_ycoord  , size*sizeof(CPPCUSTOM_REAL));
      cudaMemcpy(d_ycoord, hloc, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);

      break;
    case 3:
      cudaMalloc(&d_zcoord  , size*sizeof(CPPCUSTOM_REAL));
      cudaMemcpy(d_zcoord, hloc, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyHostToDevice);
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
                                   CPPCUSTOM_REAL *d_xcoord, 
                                   CPPCUSTOM_REAL *d_ycoord, 
                                   CPPCUSTOM_REAL *d_zcoord,  
                                   CPPCUSTOM_REAL *d_m3d,   
                                   CPPCUSTOM_REAL modelA, 
                                   CPPCUSTOM_REAL modelC, 
                                   CPPCUSTOM_REAL modelL, 
                                   CPPCUSTOM_REAL modelN, 
                                   CPPCUSTOM_REAL modelF,
                                   CPPCUSTOM_REAL *d_Cxyz){

    // npoints = number of points in 3D model to test
    // ngll = number of GLL in one direction

    // internal variables
    int  mypt, ispec, myigll, myindex;
    CPPCUSTOM_REAL  myx, myy, myz, dist, last_dist, vp, vs, rho, rad2; 
    CPPCUSTOM_REAL  myA, myC, myL, myN, myF;
    CPPCUSTOM_REAL  n1, n2, c1, c2, s1, s2, r11, r12, r13, r21, r22, r23, r31, r32, r33;

    CPPCUSTOM_REAL Cnat[6][6], Q[6][6], Crot[6][6];


    myigll = (threadIdx.y)*32 + (threadIdx.x);
    ispec  = blockIdx.x ;

    if(myigll >= 125)return; 


    // // These arrays are flattened as follows - i then j, then k, then ispec
    // // hence there are 125 values in a row for each ispec
    myx = d_xcoord[myigll*nspec + ispec];
    myy = d_ycoord[myigll*nspec + ispec];
    myz = d_zcoord[myigll*nspec + ispec];

    if (!isfinite(myx)) {
          printf("ERROR: Invalid myx calculation at point %d: dist_squared=%f\n", 
                  myigll, (double)myx);
    }
        if (!isfinite(myy)) {
          printf("ERROR: Invalid myy calculation at point %d: dist_squared=%f\n", 
                  myigll, (double)myy);
    }
    if (!isfinite(myz)) {
          printf("ERROR: Invalid myz calculation at point %d: dist_squared=%f\n", 
                  myigll, (double)myz);
    }
    
    last_dist = 100000.0;

    // Test 6: Check for potential integer overflow in model data access
    for (int ipt = 0; ipt < npoints; ++ipt){
        
        int model_index = ipt*5;
        if (model_index + 4 >= npoints * 5) {
            printf("ERROR: Model index %d out of bounds\n", model_index + 4);
            continue;
        }

        // Test 7: Check model data validity before distance calculation
        CPPCUSTOM_REAL model_x = d_m3d[model_index];
        CPPCUSTOM_REAL model_y = d_m3d[model_index + 1];
        CPPCUSTOM_REAL model_z = d_m3d[model_index + 2];

        // Calculate distance with overflow protection
        CPPCUSTOM_REAL dx = model_x - myx;
        CPPCUSTOM_REAL dy = model_y - myy;
        CPPCUSTOM_REAL dz = model_z - myz;
        
        if (!isfinite(model_x)) {
            printf("ERROR: Invalid model_x calculation at point %d: dist_squared=%f\n", 
                   ipt, (double)model_x);
            continue;
        }
        if (!isfinite(myx)) {
            printf("ERROR: Invalid dy calculation at point %d: dist_squared=%f\n", 
                   ipt, (double)dy);
            continue;
        }
 



        // Test 8: Check for potential overflow in distance calculation
        CPPCUSTOM_REAL dist_squared = dx*dx + dy*dy + dz*dz;
        if (!isfinite(dist_squared) || dist_squared < 0) {
            printf("ERROR: Invalid distance calculation at point %d: dist_squared=%f\n", 
                   ipt, (double)dist_squared);
            continue;
        }

        dist = sqrt(dist_squared);
        
        if (!isfinite(dist)) {
            printf("ERROR: Invalid distance at point %d: dist=%f\n", ipt, (double)dist);
            continue;
        }

        if (dist < last_dist) {
            last_dist = dist;
            mypt = ipt;
        }
    }
    

    if (mypt < 0 || mypt >= npoints ){
      if(myigll==1){
        printf("ERROR: in mypoint");
      }
      atomicExch(&d_errorFlag, 1);
    }
   
    // // Ultimately we set this points eta 1 and eta 2 to the final mypt 
    // // that survived: 
    n1 = d_m3d[mypt*5  + 3];
    n2 = d_m3d[mypt*5  + 4];

    if (!isfinite(n1)) {
        printf("ERROR: Invalid n1 calculation at point %d: n1=%f\n", 
                mypt, (double)n1);
    }
    if (!isfinite(n2)) {
        printf("ERROR: Invalid n2 calculation at point %d: n2=%f\n", 
                mypt, (double)n2);
    }

    if (n1 < -3.14159265359 || n1 > 3.14159265359 ){
      printf("ERROR: Invalid n1 at %d  eta1 = %f:\n",  mypt, (double)n1);
    }
    if (n2 < 0 || n2 > 3.14159265359/2.0 ){
      printf("ERROR: Invalid n2 at %d  eta1 = %f:\n",  mypt, (double)n2);
    }
    



    // // Get the PREM related values for ACLNF at this point: 
    // // r^2 normalised
    rad2 = myx*myx + myy*myy + myz*myz;

    if (!isfinite(rad2)) {
        printf("ERROR: Invalid rad2=%f\n", 
               (double)rad2);
    }

    if (rad2 < 0 || rad2>1221.5/6371.0){
      printf("ERROR: Invalid rad2 rad2=%f\n", 
        (double)rad2);
    }

    // // printf("%f %f %f %f \n", rad2, rho, vp, vs);
    // printf("%f %f %f %f\n", 
    //    (double)myx, 
    //    (double)myy, 
    //    (double)myz, 
    //    (double)rad2);

    // // Originally g/cm^3 --> kg/m^3 -> nondimensionalised
    rho = (13.088500000 - 8.838100000*rad2)*0.1813466804490;   // 1000.d0/RHOAV
    // Originally km/s --> m/s  
    vp  = (11.262200000 - 6.364000000*rad2)*0.1459938230354;   // 1000.d0/SCALE_V
    vs  = (3.667800000  - 4.447500000*rad2)*0.1459938230354;    // 1000.d0/SCALE_V

    myA = rho * vp * vp * modelA;
    myC = rho * vp * vp * modelC;
    myL = rho * vs * vs * modelL;
    myN = rho * vs * vs * modelN;

    // // eta  = 1 for PREM core and F = eta * (A - 2L)
    myF = ((rho * vp * vp) - 2.0*(rho * vs * vs)) * modelF;

    // Create natural Stiffness matrix:
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


    // Eqn 5 - DO NOT USE EXCEPT OLD BENCHMARKS
    // This is not even the correct equations for the real eqn 5
    // The commented ones are the real eqn 5    
      // r11 = c1 * c2;  // r11 = c1 * c2 ;
      // r12 = - s1;     // r12 = s1 * c2;
      // r13 = c1*s2;    // r13 = -s2; 
      
      // r21 = s1*c2;    // r21 = -s1 ;
      // r22 = c1;       // r22 = c1 ;
      // r23 = s1*s2;    // r23 = 0.0 ;

      // r31 = -s2;      // r31 = c1*s2;
      // r32 = 0.0;      // r32 = s1*s2;
      // r33 = c2;       // r33 = c2;


    // Eqn 8 of Brett 2024
    r11 = c1 * c2 ;
    r12 = -s1;
    r13 = s2*c1; 

    r21 = s1*c2 ;
    r22 = c1 ;
    r23 = s1*s2 ;

    r31 = -s2;
    r32 = 0.0;
    r33 = c2;

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
  cudaError_t ierr;
  CPPCUSTOM_REAL  A, C, L, N, F;

  // Timing: 

  ngll_per_warp       = 32;
  ngllwarps_per_block = ceil(float(ngll*ngll*ngll)/float(ngll_per_warp)); 

  A = modelv[0];
  C = modelv[1];
  L = modelv[2];
  N = modelv[3];
  F = modelv[4];



  project_eta_to_gll<<<dim3(nspec, 1, 1), 
                       dim3(ngll_per_warp, ngllwarps_per_block, 1)>>>
                       (model_npoints, ngll, nspec, d_xcoord, d_ycoord, d_zcoord, d_m3d, 
                       A, C, L, N, F,
                       d_Cxyz);

  cudaDeviceSynchronize();

  int h_errorFlag;
  cudaMemcpyFromSymbol(&h_errorFlag, d_errorFlag, sizeof(int), 0, cudaMemcpyDeviceToHost);
  if (h_errorFlag==1) {
      printf("Error: invalid index detected in kernel.\n");
      cudaDeviceReset();

      exit(1);
  } 

  cudaDeviceSynchronize();
  ierr = cudaGetLastError();
    if (ierr != cudaSuccess) {
        printf("Error running project_eta_to_gll: %s\n", cudaGetErrorString(ierr));
    }

  return 0;
 }
}









extern "C" {
int copyfromdevice(CPPCUSTOM_REAL *hloc, int size, int varid){
  
  switch (varid) {
    case 1:
      cudaMemcpy(hloc, d_xcoord, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyDeviceToHost);
      break;
    case 2:
      cudaMemcpy(hloc, d_ycoord, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyDeviceToHost);
      break;
    case 3:
      cudaMemcpy(hloc, d_zcoord, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyDeviceToHost);
      break;
    case 4:
      cudaMemcpy(hloc, d_m3d, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyDeviceToHost);
      break;
    case 7:
      cudaMemcpy(hloc, d_Cxyz, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyDeviceToHost);
      break;
    case 8:
      cudaMemcpy(hloc, d_vani_real, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyDeviceToHost);
      break;
    case 9:
      cudaMemcpy(hloc, d_vani_imag, size*sizeof(CPPCUSTOM_REAL), cudaMemcpyDeviceToHost);
      break;
    default:
      printf("Invalid VarID entered, must be 1-3 or 21 but was %i \n", varid);
      abort();
  }

  return 0;
}
}




