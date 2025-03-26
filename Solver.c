/* C libraries */
#include <time.h>

/* PETSc libraries */
#include <petsc.h>

/* PETGEM functions */
#include "solver.h"
#include "constants.h"
#include "common.h"

// =============================================================================
// Function: solveSystem
// Description: Solves the linear system using ASM domain decomposition
// =============================================================================
PetscErrorCode solveSystem(DM dm, Mat A, Vec b, Vec x) {
    PetscFunctionBeginUser;

    /* Create KSP and PC objects */
    KSP ksp;
    PC  pc;
    PetscLogDouble start_time, end_time;
    PetscInt its;
    PetscReal norm;
    PetscBool assembled;
    PetscMPIInt size;

    /* Check input parameters */
    if (!dm || !A || !b || !x) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "Input parameters cannot be NULL");
    }

    /* Check if matrix is assembled */
    PetscCall(MatAssembled(A, &assembled));
    if (!assembled) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Matrix must be assembled");
    }

    /* Get MPI size */
    PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));

    /* Start performance monitoring */
    PetscCall(PetscTime(&start_time));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Solution of linear system:\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Solver status     = started\n"));

    /* Create solver */
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));

    /* Set operators first */
    PetscCall(KSPSetOperators(ksp, A, A));

    /* Get the preconditioner context */
    PetscCall(KSPGetPC(ksp, &pc));

    /* Set up ASM preconditioner */
    PetscCall(PCSetType(pc, PCASM));

    /* Configure ASM for electromagnetic problems */
    PetscInt dim;
    PetscCall(DMGetDimension(dm, &dim));

    if (dim == 3) {
        /* For 3D electromagnetic problems, use overlapping subdomains */
        PetscCall(PCASMSetOverlap(pc, 1));  // One layer overlap
        PetscCall(PCASMSetType(pc, PC_ASM_BASIC));
    }

    /* Set up the Krylov solver */
    PetscCall(KSPSetType(ksp, KSPCG));  // Use CG for symmetric positive definite systems

    /* Configure solver tolerances and monitoring */
    PetscCall(KSPSetTolerances(ksp, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 1000));

    /* Print problem info */
    PetscInt M, N;
    PetscCall(MatGetSize(A, &M, &N));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Matrix size       = %d x %d\n", M, N));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Number of procs   = %d\n", size));

    /* Setup KSP */
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSetUp(ksp));

    /* Solve the system */
    PetscCall(KSPSolve(ksp, b, x));

    /* Get performance metrics */
    PetscCall(KSPGetIterationNumber(ksp, &its));
    PetscCall(KSPGetResidualNorm(ksp, &norm));

    /* End performance monitoring */
    PetscCall(PetscTime(&end_time));

    /* Print performance statistics */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n Performance Statistics:\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Total time        = %g seconds\n", end_time - start_time));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Iterations        = %d\n", its));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Residual norm     = %g\n", norm));

    /* Clean up */
    PetscCall(KSPDestroy(&ksp));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Solver status     = completed\n"));

    PetscFunctionReturn(PETSC_SUCCESS);
}
