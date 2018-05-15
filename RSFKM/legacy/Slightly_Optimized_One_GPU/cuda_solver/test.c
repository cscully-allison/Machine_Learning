// Load and instantiate global namespaces with relevant data and working space.
#include "solver.h"
Vars vars;
Params params;
Workspace work;
Settings settings;

void load_data(Params *params, float* Hi) {
  // In this function, load all problem instance data.
  for(int i = 0; i < 15; i++){
	params->Hi[i] = Hi[i];
  }
}

void use_solution(Vars vars) {
  // In this function, use the optimization result.
  for(int i = 0; i < 15; i++){
	printf("Ui[%d], %f \n", i, vars.Ui[i]);
  }
  
  for(int i = 0; i < 15; i++){
	printf("Hi[%d], %f \n", i, params.Hi[i]);
  }
}

void call_solver(float* Hi){
	set_defaults();  // Set basic algorithm parameters.
	setup_indexing();
	int num_iters;
	Vars vars;
	Params params;

    load_data(&params, Hi);

   // Solve our problem at high speed!
    num_iters = solve();
	
    //use_solution(vars);
}

int main(int argc, char **argv) {
	float* Hi = malloc(15*sizeof(float));
	for(int i = 0; i < 15; i++){
		Hi[i] = i * .01f;
	}
	
	call_solver(Hi);
	
	
}