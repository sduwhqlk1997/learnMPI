/*
Jing-Yuan 2025.02.16 UM-CAM E11-2024, based on ksp/tutorial/ex29.c.
Inhomogeneous Laplacian in 2D. Modeled by the partial differential equation

   -div grad u = f,  0 < x,y < 1,

with forcing function

   f = 2 * \sin \pi x * \cos \pi y

with Dirichlet boundary conditions

   u = \sin \pi x * \cos \pi y, for x = 0, x = 1, y = 0, y = 1,
   
the exact solution is 
   u = \sin \pi x * \cos \pi y.
*/

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx);
PetscErrorCode ComputeMatrix(KSP ksp, Mat J, Mat jac, void *ctx);
PetscErrorCode EvaluateAnalyticSolution(DM da, Vec sol);

int main(int argc, char **argv)
{
  KSP         ksp;
  DM          da;
  Vec         x;
  PetscLogDouble	time1,time2;
  PetscFunctionBeginUser;

  //--------------------------Preparation--------------------------------------

  // Inintialize PETSc envirionment, read the command line instructions.
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));			
  // Create a parallel KSP (linear solver) on all the processors.
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));			
  
  PetscInt	N = 4;
  // Load the number of mesh points in one direction in command line.
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-mesh_points",&N,NULL));
  N = N + 1;
  PetscScalar   Hx, Hy;  
  Hx = 1.0 / (PetscReal)(N - 1);
  Hy = 1.0 / (PetscReal)(N - 1);

  // Create a parallel structured mesh with 4 points on x and y direction. Results in a 4 by 4 grid. The dof and stencil width is 1.
  PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, N, N, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da));   
  // Read instructions in command line that is related to DM.
  PetscCall(DMSetFromOptions(da));							
  PetscCall(DMSetUp(da));
  //Set the range of compute domain, [0,1] * [0,1]
  PetscCall(DMDASetUniformCoordinates(da, 0, 1, 0, 1, 0, 0));
  PetscCall(DMDASetFieldName(da, 0, "Potential"));

  //-------------------------Computation----------------------------------------
  // This part is for solving linear system Ax = b in parallel using KSP (linear solver)

  // Evaluate the right hand side b using user defined function ComputeRHS.
  PetscCall(KSPSetComputeRHS(ksp, ComputeRHS, NULL));
  // Evaluate the linear operator A using user defind function ComputeMatrix.
  PetscCall(KSPSetComputeOperators(ksp, ComputeMatrix, NULL));
  // Set the context of da into KSP such that we can get it out in user defined function.
  PetscCall(KSPSetDM(ksp, da));
  PC	pc;
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCNONE));
  // Read instructions in command line that is related to KSP.
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSetUp(ksp));
  PetscTime(&time1);

  // Solve the linear system.
  PetscCall(KSPSolve(ksp, NULL, NULL));
  PetscTime(&time2);
  PetscPrintf(PETSC_COMM_WORLD,"Computing time: %f PETSc seconds\n",time2-time1);
  // Get the solution x out.
  PetscCall(KSPGetSolution(ksp,&x));

  //-------------------------Error Estimation-----------------------------------
  // Compute the difference between the numerical solution and analytic solution
  
  Vec	sol;
  PetscReal	 error_norm_L2, error_norm_Linf;
  PetscCall(DMCreateGlobalVector(da, &sol));
  // Compute the analytic solution
  PetscCall(EvaluateAnalyticSolution(da, sol));
  // Compute the difference by sol = sol - x;
  PetscCall(VecAXPY(sol, -1.0, x));
  // Compute the L2 norm of the error
  PetscCall(VecNorm(sol, NORM_2, &error_norm_L2));
  error_norm_L2 = error_norm_L2 * sqrt(Hx * Hy); 
  // Compute the L_inf norm of the error
  PetscCall(VecNorm(sol, NORM_INFINITY, &error_norm_Linf));
  error_norm_Linf = error_norm_Linf;
  // Output result
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "The L2 error is %.3e, the L inf error is %.3e\n",error_norm_L2, error_norm_Linf));
  // Free memory
  PetscCall(DMDestroy(&da));
  PetscCall(KSPDestroy(&ksp));
  PetscCall(PetscFinalize());
  return 0;
}

PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx)
{
  /******************************************************************************
  *	This function is to evaluate the right hand side of the linear system used  *
  *	by KSP. In this case, for the Dirichlet boundary point, the value is 		*
  *		u = \sin \pi x * \cos \pi y, for x = 0, x = 1, y = 0, y = 1				*
  *	For the inner points, the value is the discretization of function f			*
  *				f = -2 \pi^2 \sin \pi x * \cos \pi y							*
  ******************************************************************************/

  PetscInt      i, j, mx, my, xm, ym, xs, ys;
  PetscScalar   Hx, Hy;
  PetscScalar **array;
  DM            da;

  PetscFunctionBeginUser;
  // Get the context of da out to get the domain information on current processor.
  PetscCall(KSPGetDM(ksp, &da));
  // Get the global dimension in the first and the second direction.  
  PetscCall(DMDAGetInfo(da, NULL, &mx, &my, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
  // Get the mesh size Hx and Hy.
  Hx = 1.0 / (PetscReal)(mx - 1);
  Hy = 1.0 / (PetscReal)(my - 1);
  // Get the global indices (xs, ys) of the lower left corner and the size (xm, ym) of local region 
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));
  // Get the multiple dimension array that stores the data in Vector b associated with da.
  PetscCall(DMDAVecGetArray(da, b, &array));
  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      if (i == 0 || j == 0 || i == mx - 1 || j == my - 1){
	   // Dirichlet boundary condition.
       array[j][i] = PetscSinScalar((PetscReal)i * Hx * PETSC_PI) * PetscCosScalar((PetscReal)j * Hy * PETSC_PI);
      }
      else{
	   // Source term for inner points.
       array[j][i] =PetscSinScalar((PetscReal)i * Hx * PETSC_PI) * PetscCosScalar((PetscReal)j * Hy * PETSC_PI)*Hx*Hy+ 2 * pow(PETSC_PI,2) * PetscSinScalar((PetscReal)i * Hx * PETSC_PI) * PetscCosScalar((PetscReal)j * Hy * PETSC_PI) * Hx * Hy;
      }
    }
  }
  // Restore array after each get array for the purpose of resource management and parallel operation.
  PetscCall(DMDAVecRestoreArray(da, b, &array));
  // Assemble vector to make it ready to use.
  PetscCall(VecAssemblyBegin(b));
  PetscCall(VecAssemblyEnd(b));

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ComputeMatrix(KSP ksp, Mat J, Mat jac, void *ctx)
{
  /******************************************************************************
  *	This function is to evaluate the linear operator of the linear system used  *
  *	by KSP. In this case, for the Dirichlet boundary point, we set the diagonal *
  *	of the matrix to 1. For the inner point, we use finite difference method to *
  *	discretize it.																*
  ******************************************************************************/
  PetscInt     i, j, mx, my, xm, ym, xs, ys;
  PetscScalar  v[5];
  PetscReal    Hx, Hy, HydHx, HxdHy;
  MatStencil   row, col[5];
  DM           da;

  PetscFunctionBeginUser;
  // Get the context of da out to get the domain information on current processor.
  PetscCall(KSPGetDM(ksp, &da));
  // Get the global dimension in the first and the second direction.  
  PetscCall(DMDAGetInfo(da, NULL, &mx, &my, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
  // Get the mesh size Hx and Hy.
  Hx    = 1.0 / (PetscReal)(mx - 1);
  Hy    = 1.0 / (PetscReal)(my - 1);
  HxdHy = Hx / Hy;
  HydHx = Hy / Hx;
  // Get the global indices (xs, ys) of the lower left corner and the size (xm, ym) of local region 
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));
  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      row.i = i;      row.j = j;
	  // Dirichlet boundary condition.
      if (i == 0 || j == 0 || i == mx - 1 || j == my - 1) {
        v[0] = 1.0;
        PetscCall(MatSetValuesStencil(jac, 1, &row, 1, &row, v, INSERT_VALUES));
      } else {
		// Inner points
        col[0].i = i;            col[0].j = j - 1;      v[0] = -HxdHy;       
        col[1].i = i - 1;        col[1].j = j;          v[1] = -HydHx;       
        col[2].i = i;            col[2].j = j;          v[2] = 2.0 * (HxdHy + HydHx) + Hx*Hy;        
        col[3].i = i + 1;        col[3].j = j;          v[3] = -HydHx;        
        col[4].i = i;            col[4].j = j + 1;      v[4] = -HxdHy;       
		// Set the values in the matrix to specific location using information in col[]. 
        PetscCall(MatSetValuesStencil(jac, 1, &row, 5, col, v, INSERT_VALUES));
      }
    }
  }
  // Assemble matrix to make it ready to use.
  PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(MatViewFromOptions(jac, NULL, "-view_mat"));
  
  PetscFunctionReturn(PETSC_SUCCESS);
} 


PetscErrorCode EvaluateAnalyticSolution(DM da, Vec sol)
{
  /******************************************************************************
  *	This function is to evaluate the analytic solution of the problem           *
  *	            	u = \sin \pi x * \cos \pi y,                     			*
  ******************************************************************************/

  PetscInt      i, j, mx, my, xm, ym, xs, ys;
  PetscScalar   Hx, Hy;
  PetscScalar **array;

  PetscFunctionBeginUser;
  // Get the global dimension in the first and the second direction.  
  PetscCall(DMDAGetInfo(da, NULL, &mx, &my, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL));
  // Get the mesh size Hx and Hy.
  Hx = 1.0 / (PetscReal)(mx - 1);
  Hy = 1.0 / (PetscReal)(my - 1);
  // Get the global indices (xs, ys) of the lower left corner and the size (xm, ym) of local region 
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));
  // Get the multiple dimension array that stores the data in Vector b associated with da.
  PetscCall(DMDAVecGetArray(da, sol, &array));
  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
	   // Analytic solution.
       array[j][i] = PetscSinScalar((PetscReal)i * Hx * PETSC_PI) * PetscCosScalar((PetscReal)j * Hy * PETSC_PI);
    }
  }
  // Restore array after each get array for the purpose of resource management and parallel operation.
  PetscCall(DMDAVecRestoreArray(da, sol, &array));
  // Assemble vector to make it ready to use.
  PetscCall(VecAssemblyBegin(sol));
  PetscCall(VecAssemblyEnd(sol));

  PetscFunctionReturn(PETSC_SUCCESS);
}

