// Load and instantiate global namespaces with relevant data and working space.
#include "solver.h"



__device__ void load_data(Params *params, float* Hi) {
  // In this function, load all problem instance data.
  for(int i = 0; i < 15; i++){
	params->Hi[i] = Hi[i];
  }
}

__device__ void use_solution(Vars vars, Params params, int tid) {
  // In this function, use the optimization result.
  for(int i = 0; i < 15; i++){
	   printf("Ui[%d], %f, (%d) \n", i, vars.Ui[i], tid);
  }

  for(int i = 0; i < 15; i++){
	   printf("Hi[%d], %f, (%d) \n", i, params.Hi[i], tid);
  }
}

__global__ void call_solver(float** H){
    int tid = threadIdx.x;
    printf("Threadid %d", tid);
    int num_iters;

    Vars vars;
    Params params;
    Workspace work;
    Settings settings;

    set_defaults(&settings);  // Set basic algorithm parameters.
    setup_pointers(&work, &vars);


    load_data(&params, H[tid]);


   // Solve our problem at high speed!
    num_iters = solve(&work, &params, &settings, &vars);


    use_solution(vars, params, tid);

}

int main(int argc, char **argv){
    dim3 grid(1, 1);
    dim3 block(2, 1);

    float** H;
    int n = 4;
    int numcentroids = 15;

    cudaMallocManaged(&H, n * sizeof(float*));
    for(int i = 0; i < n; i++){
      cudaMallocManaged(&H[i], (numcentroids * sizeof(float)));
    }

    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 15; j++){
		    H[i][j] = (i+1)*(j+1) * .01f;
      }
	  }


	  call_solver<<<grid,block>>>(H);

     cudaDeviceSynchronize();

     //use_solution(vars);

}
