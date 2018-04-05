// Load and instantiate global namespaces with relevant data and working space.
#include "solver.h"



__device__ void load_data(Params *params, float* Hi) {
  // In this function, load all problem instance data.
  for(int i = 0; i < 15; i++){
	params->Hi[i] = Hi[i];
  }
}

__device__ void use_solution(Vars vars, Params params) {
  // In this function, use the optimization result.
  for(int i = 0; i < 15; i++){
	printf("Ui[%d], %f \n", i, vars.Ui[i]);
  }

  for(int i = 0; i < 15; i++){
	printf("Hi[%d], %f \n", i, params.Hi[i]);
  }
}

__global__ void call_solver(float* Hi){
	set_defaults();  // Set basic algorithm parameters.
	setup_indexing();
	int num_iters;

    __shared__ extern Vars vars;
    __shared__ extern Params params;
    __shared__ extern Workspace work;
    __shared__ extern Settings settings;

    load_data(&params, Hi);

   // Solve our problem at high speed!
    num_iters = solve();

}

int main(int argc, char **argv){
    dim3 grid(1, 1);
    dim3 block(1, 1);

    float* Hi;

    cudaMallocManaged(&Hi, 15 * sizeof(float));

    for(int i = 0; i < 15; i++){
		Hi[i] = i * .01f;
	}


	call_solver<<<grid,block>>>(Hi);

    //use_solution(vars);

}
