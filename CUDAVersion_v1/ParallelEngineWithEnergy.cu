//
// CUDA Sagan Simulator - a CUDA version of the sagan Simulator Software
// version 1.0 - Alpha 01
// by TTYMan
//
// para compilar:  nvcc ParallelEngineWithEnergy.cu -o ParEngEnergy -I/home/ttyman/NVIDIA_GPU_Computing_SDK/C/common/inc -arch sm_21

#include <stdio.h>
#include <string.h>
#include <cutil.h>
//cutil.h ==> /home/ttyman/NVIDIA_GPU_Computing_SDK/C/common/inc

#define BSIZE 256
#define BUFFER_MAX 2048      // buffer size for reading each line of initial conditions file.
#define offset      4        // this offset space is used to define particle position information (float4 CUDA datatype):
                             //   0 (x) - Pos X ==> array[index*offset+0]
                             //   1 (y) - Pos Y ==> array[index*offset+1]
                             //   2 (z) - Pos Z
                             //   3 (w) - Mass
                             // The velocity information is stored ina float3 CUDA datatype, as follows:
                             //   0 (x) - Vel X
                             //   1 (y) - Vel Y
                             //   2 (z) - Vel Z

                             // this approuch is useful to take advantage of CUDA float3 and float4 datatypes.

// declaracoes memoria device

// device memory for particles array
__device__ float4*  d_ParticleArrayPosition;
__device__ float3*  d_ParticleArrayVelocity;
int quantidade;


// Kernels
__global__ void UpdateParticleData(float4* d_ParticleArrayPosition, float3* d_ParticleArrayVelocity)
{
        d_ParticleArrayVelocity[threadIdx.x].x = 10.0f;
        d_ParticleArrayVelocity[threadIdx.x].z = (float)threadIdx.x;
}

__device__ float3 ParticleToParticleAccelerationWithPotentialEnergy ( float4 ParticleIPosition,
                                                                      float4 ParticleJPosition,
                                                                      float3 Acceleration,
                                                                      double  SofteningFactor,
                                                                      double* EPot )
{

    double XDist;
    double YDist;
    double ZDist;
    double XT;
    double YT;
    double ZT;
    double DistanceSQR;
    double Rinv;
    double Rinv3;
    double PotentialEnergy = 0.0;

    //float eps = 4.0f/N; // softening factor
	// SofteningFactor = eps

	XDist = ParticleIPosition.x - ParticleJPosition.x;
    XT = XDist * XDist;

    YDist = ParticleIPosition.y - ParticleJPosition.y;
    YT = YDist * YDist;

    ZDist = ParticleIPosition.z - ParticleJPosition.z;
    ZT = ZDist * ZDist;

    DistanceSQR = XT + YT + ZT;

    Rinv = 1.0 / sqrtf( DistanceSQR + (SofteningFactor * SofteningFactor) );

    Rinv3 = (Rinv * Rinv * Rinv);

    Acceleration.x -= ParticleJPosition.w * XDist * Rinv3;
	Acceleration.y -= ParticleJPosition.w * YDist * Rinv3;
	Acceleration.z -= ParticleJPosition.w * ZDist * Rinv3;

    // Potential Energy calculation
    PotentialEnergy = (ParticleIPosition.w * ParticleJPosition.w) / sqrt( DistanceSQR + (SofteningFactor * SofteningFactor) );
    //PotentialEnergy = (ParticleIPosition.w * ParticleJPosition.w) / sqrt( DistanceSQR + (SofteningFactor * SofteningFactor) );
    *EPot = PotentialEnergy;

    return Acceleration;

}


__device__ float3 ParticleToParticleAcceleration(float4 ParticleIPosition,
                                                 float4 ParticleJPosition,
                                                 float3 Acceleration,
                                                 double  SofteningFactor)
{

    double XDist;
    double YDist;
    double ZDist;
    double XT;
    double YT;
    double ZT;
    double DistanceSQR;
    double Rinv;
    double Rinv3;

    //float eps = 4.0f/N; // softening factor
	// SofteningFactor = eps

	XDist = ParticleIPosition.x - ParticleJPosition.x;
    XT = XDist * XDist;

    YDist = ParticleIPosition.y - ParticleJPosition.y;
    YT = YDist * YDist;

    ZDist = ParticleIPosition.z - ParticleJPosition.z;
    ZT = ZDist * ZDist;

    DistanceSQR = XT + YT + ZT;

    Rinv = 1.0 / sqrtf( DistanceSQR + (SofteningFactor * SofteningFactor) );

    Rinv3 = (Rinv * Rinv * Rinv);

    Acceleration.x -= ParticleJPosition.w * XDist * Rinv3;
	Acceleration.y -= ParticleJPosition.w * YDist * Rinv3;
	Acceleration.z -= ParticleJPosition.w * ZDist * Rinv3;

    return Acceleration;

}


__global__ void UpdateParticlesAcceleration(float4* d_ParticleArrayPosition,
                                            float3* d_ParticleArrayVelocity,
                                            int     N                      ,
                                            float   SofteningFactor        ,
                                            float   TimeStep               )
{

//    // set initial positions
//    d_ParticleArrayPosition[threadIdx.x].x += ( d_ParticleArrayVelocity[threadIdx.x].x * TimeStep ) / 2.0f;
//    d_ParticleArrayPosition[threadIdx.x].y += ( d_ParticleArrayVelocity[threadIdx.x].y * TimeStep ) / 2.0f;
//    d_ParticleArrayPosition[threadIdx.x].z += ( d_ParticleArrayVelocity[threadIdx.x].z * TimeStep ) / 2.0f;
//
//    __syncthreads();

    float3 Acceleration = {0.0f, 0.0f, 0.0f };

    for(int j=0; j<N; j++)
    {
        if(j!=threadIdx.x) // avoid comparition between the same particle!
        {
                Acceleration = ParticleToParticleAcceleration( d_ParticleArrayPosition[threadIdx.x],
                                                               d_ParticleArrayPosition[j],
                                                               Acceleration,
                                                               SofteningFactor);
        }
    }
    __syncthreads();

        // updating the velocities
        d_ParticleArrayVelocity[threadIdx.x].x   += Acceleration.x * TimeStep;
        d_ParticleArrayVelocity[threadIdx.x].y   += Acceleration.y * TimeStep;
        d_ParticleArrayVelocity[threadIdx.x].z   += Acceleration.z * TimeStep;

        // updating the coordinates
        d_ParticleArrayPosition[threadIdx.x].x += (d_ParticleArrayVelocity[threadIdx.x].x * TimeStep)/2.0f;
        d_ParticleArrayPosition[threadIdx.x].y += (d_ParticleArrayVelocity[threadIdx.x].y * TimeStep)/2.0f;
        d_ParticleArrayPosition[threadIdx.x].z += (d_ParticleArrayVelocity[threadIdx.x].z * TimeStep)/2.0f;

//    __syncthreads();
}

__global__ void UpdateParticlesAcceleration2(float4* d_ParticleArrayPosition,
                                            float3* d_ParticleArrayVelocity,
                                            int     N                      ,
                                            float   SofteningFactor        ,
                                            float   TimeStep               )
{

    // set initial positions
    d_ParticleArrayPosition[threadIdx.x].x += ( d_ParticleArrayVelocity[threadIdx.x].x * TimeStep ) / 2.0f;
    d_ParticleArrayPosition[threadIdx.x].y += ( d_ParticleArrayVelocity[threadIdx.x].y * TimeStep ) / 2.0f;
    d_ParticleArrayPosition[threadIdx.x].z += ( d_ParticleArrayVelocity[threadIdx.x].z * TimeStep ) / 2.0f;

    __syncthreads();

    unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;
    if(index>N) return;

    float3 Acceleration = {0.0f, 0.0f, 0.0f };

    for(int j=0; j<N; j++)
    {
        if(j!=index) // avoid comparition between the same particle!
        {
                Acceleration = ParticleToParticleAcceleration( d_ParticleArrayPosition[index],
                                                               d_ParticleArrayPosition[j],
                                                               Acceleration,
                                                               SofteningFactor);
        }
    }
    __syncthreads();

        // updating the velocities
        d_ParticleArrayVelocity[index].x   += Acceleration.x * TimeStep;
        d_ParticleArrayVelocity[index].y   += Acceleration.y * TimeStep;
        d_ParticleArrayVelocity[index].z   += Acceleration.z * TimeStep;

        // updating the coordinates
        d_ParticleArrayPosition[index].x += (d_ParticleArrayVelocity[index].x * TimeStep)/2.0f;
        d_ParticleArrayPosition[index].y += (d_ParticleArrayVelocity[index].y * TimeStep)/2.0f;
        d_ParticleArrayPosition[index].z += (d_ParticleArrayVelocity[index].z * TimeStep)/2.0f;

//    __syncthreads();
}


extern __shared__ double2 SharedMemExternArray[]; //.x for KE and .y for PE!
__global__ void UpdateParticlesAccelerationAndEnergyCalculation(float4* d_ParticleArrayPosition,
                                                                float3* d_ParticleArrayVelocity,
                                                                int     N                      ,
                                                                double   SofteningFactor        ,
                                                                double   TimeStep               ,
                                                                double   *d_KineticEnergy        ,
                                                                double   *d_PotentialEnergy      )
{

    // array alocated on shared memory to save kinetic and potential energy values for each particle.
    // SharedMemArray[index].x holds the Kinetic Energy.
    // SharedMemArray[index].y holds the Potential Energy.
    double2* SharedMemArray = (double2*)SharedMemExternArray;

    // set initial positions
    d_ParticleArrayPosition[threadIdx.x].x += ( d_ParticleArrayVelocity[threadIdx.x].x * TimeStep ) / 2.0;
    d_ParticleArrayPosition[threadIdx.x].y += ( d_ParticleArrayVelocity[threadIdx.x].y * TimeStep ) / 2.0;
    d_ParticleArrayPosition[threadIdx.x].z += ( d_ParticleArrayVelocity[threadIdx.x].z * TimeStep ) / 2.0;

    __syncthreads();

    unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;
    if(index>N) return;

    // Kinetic Energy calculation - teste teste teste
    //SharedMemArray[index].x = 1.0f;	// KE
    //SharedMemArray[index].y = 2.2f;	// PE


    float3 Acceleration    = {0.0f, 0.0f, 0.0f };
    double PotentialEnergy = 0.0;
    double PotentialEnergyACUM = 0.0;

    for(int j=0; j<N; j++)
    {
        if(j!=index) // avoid comparition between the same particle!
        {
                Acceleration = ParticleToParticleAccelerationWithPotentialEnergy( d_ParticleArrayPosition[index],
                                                                                  d_ParticleArrayPosition[j],
                                                                                  Acceleration,
                                                                                  SofteningFactor,
                                                                                  &PotentialEnergy );
                PotentialEnergyACUM += PotentialEnergy;
        }
    }
    __syncthreads();

        // updating the velocities
        d_ParticleArrayVelocity[index].x   += Acceleration.x * TimeStep;
        d_ParticleArrayVelocity[index].y   += Acceleration.y * TimeStep;
        d_ParticleArrayVelocity[index].z   += Acceleration.z * TimeStep;

        // updating the coordinates
        d_ParticleArrayPosition[index].x += (d_ParticleArrayVelocity[index].x * TimeStep)/2.0;
        d_ParticleArrayPosition[index].y += (d_ParticleArrayVelocity[index].y * TimeStep)/2.0;
        d_ParticleArrayPosition[index].z += (d_ParticleArrayVelocity[index].z * TimeStep)/2.0;

        // Kinetic Energy calculation
        SharedMemArray[index].x = (( d_ParticleArrayVelocity[threadIdx.x].x * d_ParticleArrayVelocity[threadIdx.x].x )  +
                                   ( d_ParticleArrayVelocity[threadIdx.x].y * d_ParticleArrayVelocity[threadIdx.x].y )  +
                                   ( d_ParticleArrayVelocity[threadIdx.x].z * d_ParticleArrayVelocity[threadIdx.x].z )) *
                                     d_ParticleArrayPosition[threadIdx.x].w;

        // Potential Energy saving
        SharedMemArray[index].y = PotentialEnergyACUM;

    __syncthreads();

    if(threadIdx.x == 0)
    {
        double KEsum = 0.0;   *d_KineticEnergy    = 0.0;
        double PEsum = 0.0;   *d_PotentialEnergy  = 0.0;
        for(int i=0; i < N; i++)
        {
            // adds the kinetic energy for all particles.
            KEsum += SharedMemArray[i].x;

            // sums the potential energy for each particle.
            PEsum -= SharedMemArray[i].y;
        }
        //Kinetic Energy totalization
        KEsum = KEsum * 0.5; // from the physical definition of the problem!
        //atomicAdd(d_KineticEnergy, KEsum);
        *d_KineticEnergy += KEsum;

        //Potential Energy totalization
        PEsum = PEsum / 2.0; // to avoid doubled values!
        //atomicAdd(d_PotentialEnergy, PEsum);
        *d_PotentialEnergy += PEsum;

    }
}

__global__ void EnergyCalculation(float4* d_ParticleArrayPosition,
                                  float3* d_ParticleArrayVelocity,
                                  int     N                      ,
                                  double   SofteningFactor       ,
                                  double   TimeStep              ,
                                  double   *d_KineticEnergy      ,
                                  double   *d_PotentialEnergy     )
{

    // This is the same routine "UpdateParticlesAccelerationAndEnergyCalculation", but instead to
    // update particles acceleration and position data, only calculates the KE and PE!

    // array alocated on shared memory to save kinetic and potential energy values for each particle.
    // SharedMemArray[index].x holds the Kinetic Energy.
    // SharedMemArray[index].y holds the Potential Energy.
    float2* SharedMemArray = (float2*)SharedMemExternArray;

    // set initial positions
    //d_ParticleArrayPosition[threadIdx.x].x += ( d_ParticleArrayVelocity[threadIdx.x].x * TimeStep ) / 2.0;
    //d_ParticleArrayPosition[threadIdx.x].y += ( d_ParticleArrayVelocity[threadIdx.x].y * TimeStep ) / 2.0;
    //d_ParticleArrayPosition[threadIdx.x].z += ( d_ParticleArrayVelocity[threadIdx.x].z * TimeStep ) / 2.0;

    //__syncthreads();

    unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;
    if(index>N) return;

    // Kinetic Energy calculation - teste teste teste
    //SharedMemArray[index].x = 1.0f;	// KE
    //SharedMemArray[index].y = 2.2f;	// PE


    float3 Acceleration    = {0.0f, 0.0f, 0.0f };
    double PotentialEnergy = 0.0;
    double PotentialEnergyACUM = 0.0;

    for(int j=0; j<N; j++)
    {
        if(j!=index) // avoid comparition between the same particle!
        {
                Acceleration = ParticleToParticleAccelerationWithPotentialEnergy( d_ParticleArrayPosition[index],
                                                                                  d_ParticleArrayPosition[j],
                                                                                  Acceleration,
                                                                                  SofteningFactor,
                                                                                  &PotentialEnergy );
                PotentialEnergyACUM += PotentialEnergy;
        }
    }
    __syncthreads();

        // updating the velocities
        //d_ParticleArrayVelocity[index].x   += Acceleration.x * TimeStep;
        //d_ParticleArrayVelocity[index].y   += Acceleration.y * TimeStep;
        //d_ParticleArrayVelocity[index].z   += Acceleration.z * TimeStep;

        // updating the coordinates
        //d_ParticleArrayPosition[index].x += (d_ParticleArrayVelocity[index].x * TimeStep)/2.0;
        //d_ParticleArrayPosition[index].y += (d_ParticleArrayVelocity[index].y * TimeStep)/2.0;
        //d_ParticleArrayPosition[index].z += (d_ParticleArrayVelocity[index].z * TimeStep)/2.0;

        // Kinetic Energy calculation
        SharedMemArray[index].x = (( d_ParticleArrayVelocity[threadIdx.x].x * d_ParticleArrayVelocity[threadIdx.x].x )  +
                                   ( d_ParticleArrayVelocity[threadIdx.x].y * d_ParticleArrayVelocity[threadIdx.x].y )  +
                                   ( d_ParticleArrayVelocity[threadIdx.x].z * d_ParticleArrayVelocity[threadIdx.x].z )) *
                                     d_ParticleArrayPosition[threadIdx.x].w;

        // Potential Energy saving
        SharedMemArray[index].y = PotentialEnergyACUM;

    __syncthreads();

    if(threadIdx.x == 0)
    {
        double KEsum = 0.0;   *d_KineticEnergy    = 0.0;
        double PEsum = 0.0;   *d_PotentialEnergy  = 0.0;
        for(int i=0; i < N; i++)
        {
            // adds the kinetic energy for all particles.
            KEsum += SharedMemArray[i].x;

            // sums the potential energy for each particle.
            PEsum -= SharedMemArray[i].y;
        }
        //Kinetic Energy totalization
        KEsum = KEsum * 0.5; // from the physical definition of the problem!
        //atomicAdd(d_KineticEnergy, KEsum);
        *d_KineticEnergy += KEsum;

        //Potential Energy totalization
        PEsum = PEsum / 2.0; // to avoid doubled values!
        //atomicAdd(d_PotentialEnergy, PEsum);
        *d_PotentialEnergy += PEsum;

    }
}



//Auxiliary functions

//void EnergyLog(char* FileName, float time, float energy, float energyerr, char* fmode)
//{
//
//    //
//    // used to generated the Energy Log.
//    //
//
//    FILE *file_p;
//
//    file_p = fopen(FileName, fmode);
//
//	fprintf(file_p, "%e %e %e\n", time, energy, energyerr);
//
//    fclose(file_p);
//
//}
void EnergyLog(char* FileName, double time, double Etot, double Ek, double Ep, double energyerr, char* fmode)
{

    //
    // used to generated the Energy Log.
    //

    FILE *file_p;

    file_p = fopen(FileName, fmode);

	fprintf(file_p, "%e %e %e %e %e\n", time, Etot, Ek, Ep, energyerr);

    fclose(file_p);

}


void LoadInitialConditionsData(char* FileName, int* N, float4** ParticleArrayPosition, float3** ParticleArrayVelocity)
{
    FILE* fptr;

    if ((fptr = fopen(FileName, "r")))
    {
        char    LineBuffer[BUFFER_MAX];
        int     LineCount = 0;

        // The first line of the initial conditions file should contain the amount of particles information.
        // So, read this value and use it for allocate memory for the particles arrays.
        fgets (LineBuffer, BUFFER_MAX, fptr);
        sscanf(LineBuffer, "%d ", &LineCount);
        *N = LineCount; // 'N' value finally discovered!

        // allocate memory for the Particle arrays.
        *ParticleArrayPosition = (float4*)malloc( sizeof(float4) * LineCount );
        *ParticleArrayVelocity = (float3*)malloc( sizeof(float3) * LineCount );

        // read the file line by line, saving position and velocity into data arrays.
        // each file line represents only one particle.
        for (int i=0; i<LineCount; i++)
        {
            fgets (LineBuffer, BUFFER_MAX, fptr);
            sscanf(LineBuffer, "%f %f %f %f %f %f %f " , &(*ParticleArrayPosition)[i].x,      // PosX ==> &(((*ParticleArrayPosition)+i)->x),
                                                         &(*ParticleArrayPosition)[i].y,      // PosY
                                                         &(*ParticleArrayPosition)[i].z,      // PosZ
                                                         &(*ParticleArrayVelocity)[i].x,      // VelX
                                                         &(*ParticleArrayVelocity)[i].y,      // VelY
                                                         &(*ParticleArrayVelocity)[i].z,      // VelZ
                                                         &(*ParticleArrayPosition)[i].w);     // Mass
        }

        //close the file and release ay used resource
        fclose(fptr);

        return;

    }
    else
    {
        //Error! Couldn't open the file!
        printf("\nError opening file %s!\n", FileName);
        return;
    }
}

void GenerateResultFile(char* FileName, float4* ParticleArrayPosition, float3* ParticleArrayVelocity, int N)
{

    //
    // used to generated the data file containing the updated particles positions.
    //

    FILE *file_p;

    file_p = fopen(FileName, "w");

    //first line must cointain the number of particles
    //fprintf(file_p, "%d\n", N);

    for(int i=0; i<N; i++)
    {
        fprintf(file_p, "%f %f %f %f %f %f %f\n", ParticleArrayPosition[i].x,
                                                  ParticleArrayPosition[i].y,
                                                  ParticleArrayPosition[i].z,
                                                  ParticleArrayVelocity[i].x,
                                                  ParticleArrayVelocity[i].y,
                                                  ParticleArrayVelocity[i].z,
                                                  ParticleArrayPosition[i].w);

    }

    fclose(file_p);

}



int main(int argc, char **argv)
{
    printf("\nCUDASaganSimulator version 1.0 Beta\nStarting...\n:: Getting information about available GPU...\n");


	// Getting information about installed GPU
	// this information will be used to calculate the size of execution kernel blocks.
    int devID;
    int MAX_NUMBER_THREADS_PER_BLOCK;
    int NUMBER_OF_BLOCKS;
    int THREADS_PER_BLOCK;
    cudaDeviceProp props;
    CUDA_SAFE_CALL( cudaGetDevice(&devID) );
    CUDA_SAFE_CALL( cudaGetDeviceProperties(&props, devID) );
    if (props.major == 9999 && props.minor == 9999)
    {
        // there is not CUDA devices available!!!
        printf("\n ==>>> There is no CUDA device available! Exiting now...\n");
        return 0;
    }
    MAX_NUMBER_THREADS_PER_BLOCK = props.maxThreadsPerBlock;

    printf(":::: GPU found: %s\n", props.name);
    printf(":::: Maximum Number os threads per block allowed: %d\n", MAX_NUMBER_THREADS_PER_BLOCK);




    // discrete variables
    char    Arquivo[] = "./data/figura8.ini";
    int     N = 0;

    // host memory for particle arrays.
    float4* ParticleArrayPosition;
    float3* ParticleArrayVelocity;


    // Load 'N' value and particle information from initial conditions file.
    printf(":: Loading Initial Conditions File (%s)...\n", Arquivo);
    LoadInitialConditionsData(Arquivo, &N, &ParticleArrayPosition, &ParticleArrayVelocity);

    // Once 'N' is now knowed, it's possible to define the size of the execution block
    printf(":: Initial Conditions File loaded. %d particles were found.\n", N);
    if(N>MAX_NUMBER_THREADS_PER_BLOCK)
    {
        NUMBER_OF_BLOCKS =  ((N % MAX_NUMBER_THREADS_PER_BLOCK)>0 ? 1 : 0) + (N / MAX_NUMBER_THREADS_PER_BLOCK);
        THREADS_PER_BLOCK = MAX_NUMBER_THREADS_PER_BLOCK;
    }
    else
    {
        NUMBER_OF_BLOCKS = 1;
        THREADS_PER_BLOCK = N;
    }

    printf(":: Number of execution blocks = %d\n", NUMBER_OF_BLOCKS);
    printf(":: Number of threads per block = %d\n", THREADS_PER_BLOCK);

    //return 0;

    // device memory allocation for particle data arrays.
    CUDA_SAFE_CALL ( cudaMalloc((void**)&d_ParticleArrayPosition, sizeof(float4)*N) );
    CUDA_SAFE_CALL ( cudaMalloc((void**)&d_ParticleArrayVelocity, sizeof(float3)*N) );


    // copy particle information already loaded from initial conditions file to data arrays on device memory.
    CUDA_SAFE_CALL( cudaMemcpy(d_ParticleArrayPosition, ParticleArrayPosition, sizeof(float4)*N, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy(d_ParticleArrayVelocity, ParticleArrayVelocity, sizeof(float3)*N, cudaMemcpyHostToDevice) );

    // Run the Kernel => no blocks, only threads for now :-)
    //UpdateParticleData<<<1, N>>>(d_ParticleArrayPosition, d_ParticleArrayVelocity);

    int i=0;    // we need these sequential numbers to name the output files.
    double t = 0.0;        // controls the amount of time steps already processed for the simulation.
    double dt = 0.0001;    // indicates the contribution of each scan for the simulation evolution.
    double dtout = 0.01;   // indicates the value of the step where a result file has to be generated.
    double tout = t+dtout;
    int    timeSteps = 10; //WARNING!! Don't forget to update this value with the desired time steps amount!
    double epson = 0.0;
    if(N>100)epson = 4.0/N; // softening factor

    char ResultFileName[100];

    // Energy calculation
    double* KineticEnergy;
    double* PotentialEnergy;
    double* d_KineticEnergy;
    double* d_PotentialEnergy;
    KineticEnergy   = (double*)malloc(sizeof(double));
    PotentialEnergy = (double*)malloc(sizeof(double));
    cudaMalloc( (void**)&d_KineticEnergy,   sizeof(double));
    cudaMalloc( (void**)&d_PotentialEnergy, sizeof(double));

    // this first execution is used to calculate the initial energy of the system. No updates are done.
    size_t SharedMemorySize = sizeof(double2) * N;
    EnergyCalculation<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, SharedMemorySize>>>( d_ParticleArrayPosition,
                                                                                  d_ParticleArrayVelocity,
                                                                                  N,
                                                                                  epson,
                                                                                  dt,
                                                                                  d_KineticEnergy,
                                                                                  d_PotentialEnergy );

    // copy calculated energy values from device memory back to host memory.
    CUDA_SAFE_CALL( cudaMemcpy(KineticEnergy,         d_KineticEnergy,         sizeof(double),    cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy(PotentialEnergy,       d_PotentialEnergy,       sizeof(double),    cudaMemcpyDeviceToHost) );

    double InitialKE = *KineticEnergy;
    double InitialPE = *PotentialEnergy;
    double InitialTotalEnergy = InitialKE + InitialPE;
    EnergyLog("./data/energy.log", 0, InitialTotalEnergy , InitialKE, InitialPE, 0.0, "w");


    while (t<=timeSteps)
    {

        // run the Kernel
        //UpdateParticlesAcceleration<<<1, N>>>( d_ParticleArrayPosition, d_ParticleArrayVelocity, N, 4.0f/395.0f, dt );
        //UpdateParticlesAcceleration2<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>( d_ParticleArrayPosition, d_ParticleArrayVelocity, N, epson, dt );
        size_t SharedMemorySize = sizeof(double2) * N;
        UpdateParticlesAccelerationAndEnergyCalculation<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, SharedMemorySize>>>( d_ParticleArrayPosition,
                                                                                                                    d_ParticleArrayVelocity,
                                                                                                                    N,
                                                                                                                    epson,
                                                                                                                    dt,
                                                                                                                    d_KineticEnergy,
                                                                                                                    d_PotentialEnergy );


        t += dt;

        if(t>tout)
        {
                // copy the updated data arrays from device global memory back to host memory.
                CUDA_SAFE_CALL( cudaMemcpy(ParticleArrayPosition, d_ParticleArrayPosition, sizeof(float4)*N, cudaMemcpyDeviceToHost) );
                CUDA_SAFE_CALL( cudaMemcpy(ParticleArrayVelocity, d_ParticleArrayVelocity, sizeof(float3)*N, cudaMemcpyDeviceToHost) );
                CUDA_SAFE_CALL( cudaMemcpy(KineticEnergy,         d_KineticEnergy,         sizeof(double),    cudaMemcpyDeviceToHost) );
                CUDA_SAFE_CALL( cudaMemcpy(PotentialEnergy,       d_PotentialEnergy,       sizeof(double),    cudaMemcpyDeviceToHost) );

                tout = t + dtout;
                sprintf(ResultFileName, "./data/%03d.dat", i++);
                printf("\n::: gerando arquivo de dados %s...\n", ResultFileName);
                GenerateResultFile(ResultFileName, ParticleArrayPosition, ParticleArrayVelocity, N);

                // Energy Logging
                double EK = *KineticEnergy;
                double EP = *PotentialEnergy;
                double EnergiaTotal  = EK + EP;
                double EnergyError = ( EnergiaTotal - InitialTotalEnergy) / fabs(InitialTotalEnergy);
                EnergyLog("./data/energy.log", t, EnergiaTotal , EK, EP, EnergyError, "a");
                //EnergyLog("./data/energy.log", EnergiaTotal , 0.0f, 0.0f, "a");
        }

    }


    // copy the updated data arrays from device global memory back to host memory.
    CUDA_SAFE_CALL( cudaMemcpy(ParticleArrayPosition, d_ParticleArrayPosition, sizeof(float4)*N, cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy(ParticleArrayVelocity, d_ParticleArrayVelocity, sizeof(float3)*N, cudaMemcpyDeviceToHost) );

    cudaFree(d_ParticleArrayPosition);
    cudaFree(d_ParticleArrayVelocity);
    free(ParticleArrayPosition);
    free(ParticleArrayVelocity);


//    // compara o valor dos dois arrays
//    printf("\n Listando valores\n=====================\n");
//    for(int i =0; i<N; i++)
//    {
//	printf("Posicao [%d]::: X: %f Y: %f Z: %f M: %f Vx: %f Vy: %f Vz: %f \n", i,
//	                                                                          ParticleArrayPosition[i].x,
//	                                                                          ParticleArrayPosition[i].y,
//	                                                                          ParticleArrayPosition[i].z,
//	                                                                          ParticleArrayPosition[i].w,
//	                                                                          ParticleArrayVelocity[i].x,
//	                                                                          ParticleArrayVelocity[i].y,
//	                                                                          ParticleArrayVelocity[i].z);
//    }









    printf("\nValor de N: %d\n", N);

    // // teste teste
    // for(int i=0; i<N; i++)
            // printf("main i=%d: %f %f %f %f %f %f %f \n" , i,          ParticleArray[i*offset+0],      /* PosX     */
                                                         // ParticleArray[i*offset+1],      /* PosY     */
                                                         // ParticleArray[i*offset+2],      /* PosZ     */
                                                         // ParticleArray[i*offset+3],      /* Mass     */
                                                         // ParticleArray[i*offset+4],      /* VelX     */
                                                         // ParticleArray[i*offset+5],      /* VelY     */
                                                         // ParticleArray[i*offset+6]);     /* VelZ     */

    printf("\n\nYou ran me by typing: %s\n", argv[0]);
    printf("\nPlaca: %s\n", props.name);
    return (0);
}
