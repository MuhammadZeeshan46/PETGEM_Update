/*
 * Filename: hvfem.c
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-12
 *
 * Description:
 * This file contains functions for high-order vector finite element method (HVFEM) computations. 
  *
 * List of Functions:
 * 
*/

/* C libraries */ 


/* PETSc libraries */
#include <petsc.h>
#include <petscsys.h> 

/* PETGEM functions */ 
#include "constants.h"

// =============================================================================
// Function: computeJacobian
// =============================================================================
PetscErrorCode tetrahedronXYZToXiEtaZeta(PetscScalar *cellCoords, PetscReal point[NUM_DIMENSIONS], PetscReal XiEtaZeta[NUM_DIMENSIONS]){
    /*Compute the reference tetrahedron coordinates from xyz global tetrahedron coordinates.

    :param ndarray cellCoords: spatial coordinates of the nodes with dimensions = (4,3)
    :param ndarray points: xyz points coordinates to be transformed
    :return: xietazeta points coordinates
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    PetscReal J, xi, eta, zeta, cellCoordsReal[NUM_VERTICES_PER_ELEMENT*NUM_DIMENSIONS];

    for (PetscInt i=0; i < NUM_VERTICES_PER_ELEMENT*NUM_DIMENSIONS; i++){
        cellCoordsReal[i] = PetscRealPart(cellCoords[i]);
    }

    J = cellCoordsReal[5] * ( cellCoordsReal[0] * (cellCoordsReal[10] - cellCoordsReal[7])
        + cellCoordsReal[6] * (cellCoordsReal[1] - cellCoordsReal[10])
        + cellCoordsReal[9] * (cellCoordsReal[7] - cellCoordsReal[1]) )
        + cellCoordsReal[2] * ( cellCoordsReal[3] * (cellCoordsReal[7] - cellCoordsReal[10])
        + cellCoordsReal[6] * (cellCoordsReal[10] - cellCoordsReal[4])
        + cellCoordsReal[9] * (cellCoordsReal[4] - cellCoordsReal[7]) )
        + cellCoordsReal[8] * ( cellCoordsReal[3] * (cellCoordsReal[10] - cellCoordsReal[1])
        + cellCoordsReal[0] * (cellCoordsReal[4] - cellCoordsReal[10])
        + cellCoordsReal[9] * (cellCoordsReal[1] - cellCoordsReal[4]) )
        + cellCoordsReal[11] * ( cellCoordsReal[3] * (cellCoordsReal[1] - cellCoordsReal[7])
        + cellCoordsReal[0] * (cellCoordsReal[7] - cellCoordsReal[4])
        + cellCoordsReal[6] * (cellCoordsReal[4] - cellCoordsReal[1]) );

    // Compute affine transformation for xi
    xi = ( cellCoordsReal[11] * (cellCoordsReal[7] - cellCoordsReal[4]) + cellCoordsReal[5] * (cellCoordsReal[10] - cellCoordsReal[7])
         + cellCoordsReal[8] * (cellCoordsReal[4] - cellCoordsReal[10]) ) / J * point[0] +
         ( cellCoordsReal[5] * (cellCoordsReal[6] - cellCoordsReal[9]) + cellCoordsReal[11] * (cellCoordsReal[3] - cellCoordsReal[6])
         + cellCoordsReal[8] * (cellCoordsReal[9] - cellCoordsReal[3]) ) / J * point[1] +
         ( cellCoordsReal[3] * (cellCoordsReal[7] - cellCoordsReal[10]) + cellCoordsReal[9] * (cellCoordsReal[4] - cellCoordsReal[7])
         + cellCoordsReal[6] * (cellCoordsReal[10] - cellCoordsReal[4]) ) / J * point[2] +
         ( cellCoordsReal[8] * (cellCoordsReal[3] * cellCoordsReal[10] - cellCoordsReal[9] * cellCoordsReal[4])
         + cellCoordsReal[5] * (cellCoordsReal[9] * cellCoordsReal[7] - cellCoordsReal[6] * cellCoordsReal[10])
         + cellCoordsReal[11] * (cellCoordsReal[6] * cellCoordsReal[4] - cellCoordsReal[3] * cellCoordsReal[7]) ) / J;
        
    // Compute affine transformation for eta
    eta = ( cellCoordsReal[2] * (cellCoordsReal[10] - cellCoordsReal[4]) + cellCoordsReal[11] * (cellCoordsReal[4] - cellCoordsReal[1])
          + cellCoordsReal[5] * (cellCoordsReal[1] - cellCoordsReal[10]) ) / J * point[0] +
          ( cellCoordsReal[2] * (cellCoordsReal[3] - cellCoordsReal[9]) + cellCoordsReal[5] * (cellCoordsReal[9] - cellCoordsReal[0])
          + cellCoordsReal[11] * (cellCoordsReal[0] - cellCoordsReal[3]) ) / J * point[1] +
          ( cellCoordsReal[0] * (cellCoordsReal[4] - cellCoordsReal[10]) + cellCoordsReal[3] * (cellCoordsReal[10] - cellCoordsReal[1])
          + cellCoordsReal[9] * (cellCoordsReal[1] - cellCoordsReal[4]) ) / J * point[2] +
          ( cellCoordsReal[2] * (cellCoordsReal[9] * cellCoordsReal[4] - cellCoordsReal[3] * cellCoordsReal[10])
          + cellCoordsReal[5] * (cellCoordsReal[0] * cellCoordsReal[10] - cellCoordsReal[9] * cellCoordsReal[1])
          + cellCoordsReal[11] * (cellCoordsReal[3] * cellCoordsReal[1] - cellCoordsReal[0] * cellCoordsReal[4]) ) / J;
        
    // Compute affine transformation for zeta
    zeta = ( cellCoordsReal[5] * (cellCoordsReal[7] - cellCoordsReal[1]) + cellCoordsReal[8] * (cellCoordsReal[1] - cellCoordsReal[4])
           + cellCoordsReal[2] * (cellCoordsReal[4] - cellCoordsReal[7]) ) / J * point[0] +
           ( cellCoordsReal[8] * (cellCoordsReal[3] - cellCoordsReal[0]) + cellCoordsReal[5] * (cellCoordsReal[0] - cellCoordsReal[6])
           + cellCoordsReal[2] * (cellCoordsReal[6] - cellCoordsReal[3]) ) / J * point[1] +
           ( cellCoordsReal[3] * (cellCoordsReal[1] - cellCoordsReal[7]) + cellCoordsReal[0] * (cellCoordsReal[7] - cellCoordsReal[4])
           + cellCoordsReal[6] * (cellCoordsReal[4] - cellCoordsReal[1]) ) / J * point[2] +
           ( cellCoordsReal[5] * ( cellCoordsReal[9] * (cellCoordsReal[1] - cellCoordsReal[7])
           + cellCoordsReal[10] * (cellCoordsReal[6] - cellCoordsReal[0]) )
           + cellCoordsReal[8] * ( cellCoordsReal[9] * (cellCoordsReal[4] - cellCoordsReal[1])
           + cellCoordsReal[10] * (cellCoordsReal[0] - cellCoordsReal[3]) )
           + cellCoordsReal[2] * ( cellCoordsReal[9] * (cellCoordsReal[7] - cellCoordsReal[4])
           + cellCoordsReal[10] * (cellCoordsReal[3] - cellCoordsReal[6]) )
           + cellCoordsReal[11] * ( cellCoordsReal[0] * (cellCoordsReal[4] - cellCoordsReal[7])
           + cellCoordsReal[3] * (cellCoordsReal[7] - cellCoordsReal[1])
           + cellCoordsReal[6] * (cellCoordsReal[1] - cellCoordsReal[4]) ) + J ) / J;
            
    XiEtaZeta[0] = xi;
    XiEtaZeta[1] = eta;
    XiEtaZeta[2] = zeta;
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode vectorRotation(PetscReal azimuth, PetscReal dip, PetscReal rotationVector[NUM_DIMENSIONS]){
    /*Compute the weigths vector for source rotation in the xyz plane.

    :param float azimuth: degrees for x-y plane rotation
    :param float dip: degrees for x-z plane rotation
    :return: weigths for source rotation
    :rtype: ndarray.
    */

    PetscFunctionBeginUser;

    PetscReal base_vector[NUM_DIMENSIONS] = {1., 0., 0.};

    // ---------------------------------------------------------------
    // Compute vector for source rotation
    // ---------------------------------------------------------------
    // Convert degrees to radians for rotation
    PetscReal alpha = azimuth * PETSC_PI / 180.;    // x-y plane
    PetscReal beta  = dip * PETSC_PI / 180.;        // x-z plane
    PetscReal tetha = 0.0 * PETSC_PI / 180.;        // y-z plane

    // Define rotation matrices for each plane
    // x-y plane
    PetscReal M1[NUM_DIMENSIONS][NUM_DIMENSIONS] = {{PetscCosReal(alpha), -PetscSinReal(alpha),   0.},
                                                    {PetscSinReal(alpha),  PetscCosReal(alpha),   0.},
                                                    {                 0.,                   0.,   1.}};

    // x-z plane
    PetscReal M2[NUM_DIMENSIONS][NUM_DIMENSIONS] = {{PetscCosReal(beta),  0.,  -PetscSinReal(beta)},
                                                    {                0.,  1.,                   0.},
                                                    {PetscSinReal(beta),  0.,   PetscCosReal(beta)}};

     // y-z plane
    PetscReal M3[NUM_DIMENSIONS][NUM_DIMENSIONS] = {{1.,   0.,                                     0.},
                                                    {0.,   PetscCosReal(tetha),  -PetscSinReal(tetha)},
                                                    {0.,   PetscSinReal(tetha),   PetscCosReal(tetha)}};
    
    PetscReal temp1[NUM_DIMENSIONS][NUM_DIMENSIONS], temp2[NUM_DIMENSIONS][NUM_DIMENSIONS];

    // Perform matrix multiplications
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
            PetscReal sum = 0.0;
            for (PetscInt k = 0; k < NUM_DIMENSIONS; k++) {
                sum += M1[i][k] * M2[k][j];
            }
            temp1[i][j] = sum;
        }
    }

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
            PetscReal sum = 0.0;
            for (PetscInt k = 0; k < NUM_DIMENSIONS; k++) {
                sum += temp1[i][k] * M3[k][j];
            }
            temp2[i][j] = sum;
        }
    }

    // Apply rotation
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        rotationVector[i] = 0.0;
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
            rotationVector[i] += temp2[i][j] * base_vector[j];
        }
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode crossProduct(PetscReal vector1[NUM_DIMENSIONS], PetscReal vector2[NUM_DIMENSIONS], PetscReal result[NUM_DIMENSIONS]){
    PetscFunctionBeginUser;

    result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1]; // Compute x component
    result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2]; // Compute y component
    result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0]; // Compute z component

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode matrixVectorProduct(PetscReal vector[NUM_DIMENSIONS], PetscReal matrix[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal result[NUM_DIMENSIONS]){
    PetscFunctionBeginUser;

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        result[i] = 0.0;
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++){
            result[i] += matrix[i][j] * vector[j];
        }
    }    

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode vectorMatrixProduct(PetscReal vector[NUM_DIMENSIONS], PetscReal matrix[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal result[NUM_DIMENSIONS]){
    PetscFunctionBeginUser;

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        result[i] = 0.0;
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++){
            result[i] += matrix[j][i] * vector[j];
        }
    }    

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode dotProduct(PetscReal vector1[NUM_DIMENSIONS], PetscReal vector2[NUM_DIMENSIONS], PetscReal *result){
    PetscFunctionBeginUser;

    *result = 0.0;
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        *result += vector1[i] * vector2[i];
    }    

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode computeJacobian(PetscScalar *cellCoords, PetscReal jacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal invJacobian[NUM_DIMENSIONS][NUM_DIMENSIONS]) {
    PetscFunctionBeginUser;

    /* Reset matrices */
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
            jacobian[i][j] = 0.0;
            invJacobian[i][j] = 0.0;
        }
    }

    /* Compute jacobian */ 
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        jacobian[0][i] = PetscRealPart(cellCoords[i])   - PetscRealPart(cellCoords[3+i]); 
        jacobian[1][i] = PetscRealPart(cellCoords[6+i]) - PetscRealPart(cellCoords[3+i]); 
        jacobian[2][i] = PetscRealPart(cellCoords[9+i]) - PetscRealPart(cellCoords[3+i]);
    }

    /* Compute determinant */
    PetscReal determinant; 

    determinant = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1]) -
                  jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0]) +
                  jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);

    /* Compute cofactor matrix */
    PetscReal coFactorMatrix[NUM_DIMENSIONS][NUM_DIMENSIONS];

    coFactorMatrix[0][0] =   jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1];
    coFactorMatrix[0][1] = -(jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0]);
    coFactorMatrix[0][2] =   jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0];
    coFactorMatrix[1][0] = -(jacobian[0][1] * jacobian[2][2] - jacobian[0][2] * jacobian[2][1]);
    coFactorMatrix[1][1] =   jacobian[0][0] * jacobian[2][2] - jacobian[0][2] * jacobian[2][0];
    coFactorMatrix[1][2] = -(jacobian[0][0] * jacobian[2][1] - jacobian[0][1] * jacobian[2][0]);
    coFactorMatrix[2][0] =   jacobian[0][1] * jacobian[1][2] - jacobian[0][2] * jacobian[1][1];
    coFactorMatrix[2][1] = -(jacobian[0][0] * jacobian[1][2] - jacobian[0][2] * jacobian[1][0]);
    coFactorMatrix[2][2] =   jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];

    /* Compute adjugate matrix */ 
    PetscReal adjugateMatrix[NUM_DIMENSIONS][NUM_DIMENSIONS];

    for (PetscInt i = 0; i < 3; i++) {
        for (PetscInt j = 0; j < 3; j++) {
            adjugateMatrix[i][j] = coFactorMatrix[j][i];
        }
    }

    /* Compute inverse of jacobian */
    PetscReal invDeterminant;

    invDeterminant = 1.0 / determinant;
            
    for (PetscInt i = 0; i < 3; i++) {
        for (PetscInt j = 0; j < 3; j++) {
            invJacobian[i][j] = invDeterminant * adjugateMatrix[i][j];
        }
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}
    

PetscErrorCode computeNumGaussPoints3D(PetscInt nord, PetscInt *numGaussPoints){
    PetscFunctionBeginUser;

    PetscInt gaussOrder;
                
    /* Compute gauss order */ 
    gaussOrder = 2*nord;
            
    PetscCheck(gaussOrder >= 1, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Error: Orders lower than 1 are not supported.\n");
    PetscCheck(gaussOrder <= 12, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Error: Orders higher than 6 are not supported by PETGEM.\n");
    
    switch (gaussOrder) {
        case 1:  *numGaussPoints = 1; break;
        case 2:  *numGaussPoints = 4; break;
        case 3:  *numGaussPoints = 5; break;
        case 4:  *numGaussPoints = 11; break;
        case 5:  *numGaussPoints = 14; break;
        case 6:  *numGaussPoints = 24; break;
        case 7:  *numGaussPoints = 31; break;
        case 8:  *numGaussPoints = 43; break;
        case 9:  *numGaussPoints = 53; break;
        case 10: *numGaussPoints = 126; break;
        case 11: *numGaussPoints = 126; break;
        case 12: *numGaussPoints = 210; break;
        default: break;    
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode renormalization3DGaussPoints(PetscInt numPoints, const PetscReal (*gaussPoints)[4], PetscReal** points, PetscReal* weights){
    PetscFunctionBeginUser;

    for(PetscInt i = 0; i < numPoints; i++){
        weights[i] = gaussPoints[i][3]/8;
        points[i][0] = (1 + gaussPoints[i][1])/2;
        points[i][1] = -(1 + gaussPoints[i][0] + gaussPoints[i][1] + gaussPoints[i][2])/2;
        points[i][2] = (1 + gaussPoints[i][0])/2;
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode computeGaussPoints3D(PetscInt numPoints, PetscReal **points, PetscReal *weights){
    PetscFunctionBeginUser;

    const PetscReal nord1_3DGaussPoints[1][4] = {{-0.500000000000000, -0.500000000000000, -0.500000000000000, 1.333333333333333}};

    const PetscReal nord2_3DGaussPoints[4][4] = {{-0.723606797749979, -0.723606797749979, -0.723606797749979, 0.333333333333333},
                                                 { 0.170820393249937, -0.723606797749979, -0.723606797749979, 0.333333333333333},
                                                 {-0.723606797749979,  0.170820393249937, -0.723606797749979, 0.333333333333333},
                                                 {-0.723606797749979, -0.723606797749979,  0.170820393249937, 0.333333333333333}};

    const PetscReal nord3_3DGaussPoints[5][4] = {{-0.500000000000000, -0.500000000000000, -0.500000000000000, -1.066666666666667},
                                                 {-0.666666666666667, -0.666666666666667, -0.666666666666667,  0.600000000000000},
                                                 {-0.666666666666667, -0.666666666666667,  0.000000000000000,  0.600000000000000},
                                                 {-0.666666666666667,  0.000000000000000, -0.666666666666667,  0.600000000000000},
                                                 { 0.000000000000000, -0.666666666666667, -0.666666666666667,  0.600000000000000}};

    const PetscReal nord4_3DGaussPoints[11][4] = {{-0.500000000000000, -0.500000000000000, -0.500000000000000, -0.105244444444444},
                                                  {-0.857142857142857, -0.857142857142857, -0.857142857142857,  0.060977777777778},
                                                  {-0.857142857142857, -0.857142857142857,  0.571428571428571,  0.060977777777778},
                                                  {-0.857142857142857,  0.571428571428571, -0.857142857142857,  0.060977777777778},
                                                  { 0.571428571428571, -0.857142857142857, -0.857142857142857,  0.060977777777778},
                                                  {-0.201192847666402, -0.201192847666402, -0.798807152333598,  0.199111111111111},
                                                  {-0.201192847666402, -0.798807152333598, -0.201192847666402,  0.199111111111111},
                                                  {-0.798807152333598, -0.201192847666402, -0.201192847666402,  0.199111111111111},
                                                  {-0.201192847666402, -0.798807152333598, -0.798807152333598,  0.199111111111111},
                                                  {-0.798807152333598, -0.201192847666402, -0.798807152333598,  0.199111111111111},
                                                  {-0.798807152333598, -0.798807152333598, -0.201192847666402,  0.199111111111111}};

    const PetscReal nord5_3DGaussPoints[14][4] = {{-0.814529499378218, -0.814529499378218, -0.814529499378218,  0.097990724155149},
                                                  { 0.443588498134653, -0.814529499378218, -0.814529499378218,  0.097990724155149},
                                                  {-0.814529499378218,  0.443588498134653, -0.814529499378218,  0.097990724155149},
                                                  {-0.814529499378218, -0.814529499378218,  0.443588498134653,  0.097990724155149},
                                                  {-0.378228161473399, -0.378228161473399, -0.378228161473399,  0.150250567624021},
                                                  {-0.865315515579804, -0.378228161473399, -0.378228161473399,  0.150250567624021},
                                                  {-0.378228161473399, -0.865315515579804, -0.378228161473399,  0.150250567624021},
                                                  {-0.378228161473399, -0.378228161473399, -0.865315515579804,  0.150250567624021},
                                                  {-0.091007408251299, -0.091007408251299, -0.908992591748701,  0.056728027702775},
                                                  {-0.091007408251299, -0.908992591748701, -0.091007408251299,  0.056728027702775},
                                                  {-0.908992591748701, -0.091007408251299, -0.091007408251299,  0.056728027702775},
                                                  {-0.091007408251299, -0.908992591748701, -0.908992591748701,  0.056728027702775},
                                                  {-0.908992591748701, -0.091007408251299, -0.908992591748701,  0.056728027702775},
                                                  {-0.908992591748701, -0.908992591748701, -0.091007408251299,  0.056728027702775}};

    const PetscReal nord6_3DGaussPoints[24][4] = {{-0.570794257481696, -0.570794257481696, -0.570794257481696, 0.053230333677557},
                                                  {-0.287617227554912, -0.570794257481696, -0.570794257481696, 0.053230333677557},
                                                  {-0.570794257481696, -0.287617227554912, -0.570794257481696, 0.053230333677557},
                                                  {-0.570794257481696, -0.570794257481696, -0.287617227554912, 0.053230333677557},
                                                  {-0.918652082930777, -0.918652082930777, -0.918652082930777, 0.013436281407094},
                                                  { 0.755956248792332, -0.918652082930777, -0.918652082930777, 0.013436281407094},
                                                  {-0.918652082930777,  0.755956248792332, -0.918652082930777 ,0.013436281407094},
                                                  {-0.918652082930777, -0.918652082930777,  0.755956248792332 ,0.013436281407094},
                                                  {-0.355324219715449, -0.355324219715449, -0.355324219715449, 0.073809575391540},
                                                  {-0.934027340853653, -0.355324219715449, -0.355324219715449, 0.073809575391540},
                                                  {-0.355324219715449, -0.934027340853653, -0.355324219715449, 0.073809575391540},
                                                  {-0.355324219715449, -0.355324219715449, -0.934027340853653, 0.073809575391540},
                                                  {-0.872677996249965, -0.872677996249965, -0.460655337083368, 0.064285714285714},
                                                  {-0.872677996249965, -0.460655337083368, -0.872677996249965, 0.064285714285714},
                                                  {-0.872677996249965, -0.872677996249965,  0.206011329583298, 0.064285714285714},
                                                  {-0.872677996249965,  0.206011329583298, -0.872677996249965, 0.064285714285714},
                                                  {-0.872677996249965, -0.460655337083368,  0.206011329583298, 0.064285714285714},
                                                  {-0.872677996249965,  0.206011329583298, -0.460655337083368, 0.064285714285714},
                                                  {-0.460655337083368, -0.872677996249965, -0.872677996249965, 0.064285714285714},
                                                  {-0.460655337083368, -0.872677996249965,  0.206011329583298, 0.064285714285714},
                                                  {-0.460655337083368,  0.206011329583298, -0.872677996249965, 0.064285714285714},
                                                  { 0.206011329583298, -0.872677996249965, -0.460655337083368, 0.064285714285714},
                                                  { 0.206011329583298, -0.872677996249965, -0.872677996249965, 0.064285714285714},
                                                  { 0.206011329583298, -0.460655337083368, -0.872677996249965, 0.064285714285714}};

    const PetscReal nord7_3DGaussPoints[31][4] = {{ 0.000000000000000,   0.000000000000000,  -1.000000000000000,   0.007760141093474},
                                                  { 0.000000000000000,  -1.000000000000000,   0.000000000000000,   0.007760141093474},
                                                  {-1.000000000000000,   0.000000000000000,   0.000000000000000,   0.007760141093474},
                                                  {-1.000000000000000,  -1.000000000000000,   0.000000000000000,   0.007760141093474},
                                                  {-1.000000000000000,   0.000000000000000,  -1.000000000000000,   0.007760141093474},
                                                  { 0.000000000000000,  -1.000000000000000,  -1.000000000000000,   0.007760141093474},
                                                  {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,   0.146113787728871},
                                                  {-0.843573615339364,  -0.843573615339364,  -0.843573615339364,   0.084799532195309},
                                                  {-0.843573615339364,  -0.843573615339364,   0.530720846018092,   0.084799532195309},
                                                  {-0.843573615339364,   0.530720846018092,  -0.843573615339364,   0.084799532195309},
                                                  { 0.530720846018092,  -0.843573615339364,  -0.843573615339364,   0.084799532195309},
                                                  {-0.756313566672190,  -0.756313566672190,  -0.756313566672190,  -0.500141920914655},
                                                  {-0.756313566672190,  -0.756313566672190,   0.268940700016569,  -0.500141920914655},
                                                  {-0.756313566672190,   0.268940700016569,  -0.756313566672190,  -0.500141920914655},
                                                  { 0.268940700016569,  -0.756313566672190,  -0.756313566672190,  -0.500141920914655},
                                                  {-0.334921671107159,  -0.334921671107159,  -0.334921671107159,   0.039131402104588},
                                                  {-0.334921671107159,  -0.334921671107159,  -0.995234986678524,   0.039131402104588},
                                                  {-0.334921671107159,  -0.995234986678524,  -0.334921671107159,   0.039131402104588},
                                                  {-0.995234986678524,  -0.334921671107159,  -0.334921671107159,   0.039131402104588},
                                                  {-0.800000000000000,  -0.800000000000000,  -0.600000000000000,   0.220458553791887},
                                                  {-0.800000000000000,  -0.600000000000000,  -0.800000000000000,   0.220458553791887},
                                                  {-0.800000000000000,  -0.800000000000000,   0.200000000000000,   0.220458553791887},
                                                  {-0.800000000000000,   0.200000000000000,  -0.800000000000000,   0.220458553791887},
                                                  {-0.800000000000000,  -0.600000000000000,   0.200000000000000,   0.220458553791887},
                                                  {-0.800000000000000,   0.200000000000000,  -0.600000000000000,   0.220458553791887},
                                                  {-0.600000000000000,  -0.800000000000000,  -0.800000000000000,   0.220458553791887},
                                                  {-0.600000000000000,  -0.800000000000000,   0.200000000000000,   0.220458553791887},
                                                  {-0.600000000000000,   0.200000000000000,  -0.800000000000000,   0.220458553791887},
                                                  { 0.200000000000000,  -0.800000000000000,  -0.600000000000000,   0.220458553791887},
                                                  { 0.200000000000000,  -0.800000000000000,  -0.800000000000000,   0.220458553791887},
                                                  { 0.200000000000000,  -0.600000000000000,  -0.800000000000000,   0.220458553791887}};

    const PetscReal nord8_3DGaussPoints[43][4] = {{-0.500000000000000,  -0.500000000000000,  -0.500000000000000, -0.164001509269119},
                                                  {-0.586340136778654,  -0.586340136778654,  -0.586340136778654,  0.114002446582935},
                                                  {-0.586340136778654,  -0.586340136778654,  -0.240979589664039,  0.114002446582935},
                                                  {-0.586340136778654,  -0.240979589664039,  -0.586340136778654,  0.114002446582935},
                                                  {-0.240979589664039,  -0.586340136778654,  -0.586340136778654,  0.114002446582935},
                                                  {-0.835792823378907,  -0.835792823378907,  -0.835792823378907,  0.015736266505071},
                                                  {-0.835792823378907,  -0.835792823378907,   0.507378470136720,  0.015736266505071},
                                                  {-0.835792823378907,   0.507378470136720,  -0.835792823378907,  0.015736266505071},
                                                  { 0.507378470136720,  -0.835792823378907,  -0.835792823378907,  0.015736266505071},
                                                  {-0.988436098989604,  -0.988436098989604,  -0.988436098989604,  0.001358672872743},
                                                  {-0.988436098989604,  -0.988436098989604,   0.965308296968812,  0.001358672872743},
                                                  {-0.988436098989604,   0.965308296968812,  -0.988436098989604,  0.001358672872743},
                                                  { 0.965308296968812,  -0.988436098989604,  -0.988436098989604,  0.001358672872743},
                                                  {-0.898934519962212,  -0.898934519962212,  -0.101065480037788,  0.036637470595738},
                                                  {-0.898934519962212,  -0.101065480037788,  -0.898934519962212,  0.036637470595738},
                                                  {-0.101065480037788,  -0.898934519962212,  -0.898934519962212,  0.036637470595738},
                                                  {-0.898934519962212,  -0.101065480037788,  -0.101065480037788,  0.036637470595738},
                                                  {-0.101065480037788,  -0.898934519962212,  -0.101065480037788,  0.036637470595738},
                                                  {-0.101065480037788,  -0.101065480037788,  -0.898934519962212,  0.036637470595738},
                                                  {-0.541866927766378,  -0.541866927766378,  -0.928720834422932,  0.045635886469455},
                                                  {-0.541866927766378,  -0.928720834422932,  -0.541866927766378,  0.045635886469455},
                                                  {-0.541866927766378,  -0.541866927766378,   0.012454689955687,  0.045635886469455},
                                                  {-0.541866927766378,   0.012454689955687,  -0.541866927766378,  0.045635886469455},
                                                  {-0.541866927766378,  -0.928720834422932,   0.012454689955687,  0.045635886469455},
                                                  {-0.541866927766378,   0.012454689955687,  -0.928720834422932,  0.045635886469455},
                                                  {-0.928720834422932,  -0.541866927766378,  -0.541866927766378,  0.045635886469455},
                                                  {-0.928720834422932,  -0.541866927766378,   0.012454689955687,  0.045635886469455},
                                                  {-0.928720834422932,   0.012454689955687,  -0.541866927766378,  0.045635886469455},
                                                  { 0.012454689955687,  -0.541866927766378,  -0.928720834422932,  0.045635886469455},
                                                  { 0.012454689955687,  -0.541866927766378,  -0.541866927766378,  0.045635886469455},
                                                  { 0.012454689955687,  -0.928720834422932,  -0.541866927766378,  0.045635886469455},
                                                  {-0.926784500893605,  -0.926784500893605,  -0.619027916130733,  0.017124153129297},
                                                  {-0.926784500893605,  -0.619027916130733,  -0.926784500893605,  0.017124153129297},
                                                  {-0.926784500893605,  -0.926784500893605,   0.472596917917943,  0.017124153129297},
                                                  {-0.926784500893605,   0.472596917917943,  -0.926784500893605,  0.017124153129297},
                                                  {-0.926784500893605,  -0.619027916130733,   0.472596917917943,  0.017124153129297},
                                                  {-0.926784500893605,   0.472596917917943,  -0.619027916130733,  0.017124153129297},
                                                  {-0.619027916130733,  -0.926784500893605,  -0.926784500893605,  0.017124153129297},
                                                  {-0.619027916130733,  -0.926784500893605,   0.472596917917943,  0.017124153129297},
                                                  {-0.619027916130733,   0.472596917917943,  -0.926784500893605,  0.017124153129297},
                                                  { 0.472596917917943,  -0.926784500893605,  -0.619027916130733,  0.017124153129297},
                                                  { 0.472596917917943,  -0.926784500893605,  -0.926784500893605,  0.017124153129297},
                                                  { 0.472596917917943,  -0.619027916130733,  -0.926784500893605,  0.017124153129297}};

    const PetscReal nord9_3DGaussPoints[53][4] = {{-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  -1.102392306608869},
                                                  {-0.903297922900526,  -0.903297922900526,  -0.903297922900526,   0.014922692552682},
                                                  {-0.903297922900526,  -0.903297922900526,   0.709893768701580,   0.014922692552682},
                                                  {-0.903297922900526,   0.709893768701580,  -0.903297922900526,   0.014922692552682},
                                                  { 0.709893768701580,  -0.903297922900526,  -0.903297922900526,   0.014922692552682},
                                                  {-0.350841439764235,  -0.350841439764235,  -0.350841439764235,   0.034475391755947},
                                                  {-0.350841439764235,  -0.350841439764235,  -0.947475680707294,   0.034475391755947},
                                                  {-0.350841439764235,  -0.947475680707294,  -0.350841439764235,   0.034475391755947},
                                                  {-0.947475680707294,  -0.350841439764235,  -0.350841439764235,   0.034475391755947},
                                                  {-0.770766919552010,  -0.770766919552010,  -0.770766919552010,  -0.721478131849612},
                                                  {-0.770766919552010,  -0.770766919552010,   0.312300758656029,  -0.721478131849612},
                                                  {-0.770766919552010,   0.312300758656029,  -0.770766919552010,  -0.721478131849612},
                                                  { 0.312300758656029,  -0.770766919552010,  -0.770766919552010,  -0.721478131849612},
                                                  {-0.549020096176972,  -0.549020096176972,  -0.549020096176972,   0.357380609620092},
                                                  {-0.549020096176972,  -0.549020096176972,  -0.352939711469084,   0.357380609620092},
                                                  {-0.549020096176972,  -0.352939711469084,  -0.549020096176972,   0.357380609620092},
                                                  {-0.352939711469084,  -0.549020096176972,  -0.549020096176972,   0.357380609620092},
                                                  {-0.736744381506260,  -0.736744381506260,  -0.832670596765630,   0.277603247076406},
                                                  {-0.736744381506260,  -0.832670596765630,  -0.736744381506260,   0.277603247076406},
                                                  {-0.736744381506260,  -0.736744381506260,   0.306159359778151,   0.277603247076406},
                                                  {-0.736744381506260,   0.306159359778151,  -0.736744381506260,   0.277603247076406},
                                                  {-0.736744381506260,  -0.832670596765630,   0.306159359778151,   0.277603247076406},
                                                  {-0.736744381506260,   0.306159359778151,  -0.832670596765630,   0.277603247076406},
                                                  {-0.832670596765630,  -0.736744381506260,  -0.736744381506260,   0.277603247076406},
                                                  {-0.832670596765630,  -0.736744381506260,   0.306159359778151,   0.277603247076406},
                                                  {-0.832670596765630,   0.306159359778151,  -0.736744381506260,   0.277603247076406},
                                                  { 0.306159359778151,  -0.736744381506260,  -0.832670596765630,   0.277603247076406},
                                                  { 0.306159359778151,  -0.736744381506260,  -0.736744381506260,   0.277603247076406},
                                                  { 0.306159359778151,  -0.832670596765630,  -0.736744381506260,   0.277603247076406},
                                                  {-0.132097077177186,  -0.132097077177186,  -0.784460280901143,   0.026820671221285},
                                                  {-0.132097077177186,  -0.784460280901143,  -0.132097077177186,   0.026820671221285},
                                                  {-0.132097077177186,  -0.132097077177186,  -0.951345564744484,   0.026820671221285},
                                                  {-0.132097077177186,  -0.951345564744484,  -0.132097077177186,   0.026820671221285},
                                                  {-0.132097077177186,  -0.784460280901143,  -0.951345564744484,   0.026820671221285},
                                                  {-0.132097077177186,  -0.951345564744484,  -0.784460280901143,   0.026820671221285},
                                                  {-0.784460280901143,  -0.132097077177186,  -0.132097077177186,   0.026820671221285},
                                                  {-0.784460280901143,  -0.132097077177186,  -0.951345564744484,   0.026820671221285},
                                                  {-0.784460280901143,  -0.951345564744484,  -0.132097077177186,   0.026820671221285},
                                                  {-0.951345564744484,  -0.132097077177186,  -0.784460280901143,   0.026820671221285},
                                                  {-0.951345564744484,  -0.132097077177186,  -0.132097077177186,   0.026820671221285},
                                                  {-0.951345564744484,  -0.784460280901143,  -0.132097077177186,   0.026820671221285},
                                                  {-1.002752554636276,  -1.002752554636276,  -0.446893054726385,   0.003453031004456},
                                                  {-1.002752554636276,  -0.446893054726385,  -1.002752554636276,   0.003453031004456},
                                                  {-1.002752554636276,  -1.002752554636276,   0.452398163998938,   0.003453031004456},
                                                  {-1.002752554636276,   0.452398163998938,  -1.002752554636276,   0.003453031004456},
                                                  {-1.002752554636276,  -0.446893054726385,   0.452398163998938,   0.003453031004456},
                                                  {-1.002752554636276,   0.452398163998938,  -0.446893054726385,   0.003453031004456},
                                                  {-0.446893054726385,  -1.002752554636276,  -1.002752554636276,   0.003453031004456},
                                                  {-0.446893054726385,  -1.002752554636276,   0.452398163998938,   0.003453031004456},
                                                  {-0.446893054726385,   0.452398163998938,  -1.002752554636276,   0.003453031004456},
                                                  { 0.452398163998938,  -1.002752554636276,  -0.446893054726385,   0.003453031004456},
                                                  { 0.452398163998938,  -1.002752554636276,  -1.002752554636276,   0.003453031004456},
                                                  { 0.452398163998938,  -0.446893054726385,  -1.002752554636276,   0.003453031004456}};

    const PetscReal nord10_3DGaussPoints[126][4] = {{-0.857142857142857,  -0.857142857142857,  0.571428571428571 ,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,   0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,  -0.285714285714286,   0.000000000000000,   0.362902592520648},
                                                    {-0.857142857142857,   0.000000000000000,  -0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,   0.285714285714286,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,   0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,   0.285714285714286,   0.362902592520648},
                                                    {-0.571428571428571,  -0.571428571428571,   0.000000000000000,   0.362902592520648},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.285714285714286,   0.362902592520648},
                                                    {-0.571428571428571,   0.000000000000000,  -0.571428571428571,   0.362902592520648},
                                                    {-0.571428571428571,   0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    {-0.285714285714286,  -0.857142857142857,   0.000000000000000,   0.362902592520648},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.285714285714286,   0.362902592520648},
                                                    {-0.285714285714286,  -0.285714285714286,  -0.571428571428571,   0.362902592520648},
                                                    {-0.285714285714286,   0.000000000000000,  -0.857142857142857,   0.362902592520648},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.285714285714286,   0.362902592520648},
                                                    { 0.000000000000000,  -0.571428571428571,  -0.571428571428571,   0.362902592520648},
                                                    { 0.000000000000000,  -0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    { 0.285714285714286,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    { 0.285714285714286,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    { 0.571428571428571,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,   0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,   0.000000000000000,   0.362902592520648},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,   0.000000000000000,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,   0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,   0.000000000000000,   0.362902592520648},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.285714285714286,   0.362902592520648},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.571428571428571,   0.362902592520648},
                                                    {-0.571428571428571,   0.000000000000000,  -0.857142857142857,   0.362902592520648},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.285714285714286,   0.362902592520648},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.571428571428571,   0.362902592520648},
                                                    {-0.285714285714286,  -0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    { 0.000000000000000,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    { 0.285714285714286,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,   0.000000000000000,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,   0.000000000000000,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.285714285714286,   0.362902592520648},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.571428571428571,   0.362902592520648},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.833333333333333,  -0.833333333333333,   0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.500000000000000,   0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,   0.166666666666667,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,   0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.833333333333333,   0.166666666666667,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.166666666666667,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.500000000000000,   0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.500000000000000,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    { 0.166666666666667,  -0.833333333333333,  -0.500000000000000,  -0.932187812187812},
                                                    { 0.166666666666667,  -0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    { 0.500000000000000,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.833333333333333,   0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,   0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    { 0.166666666666667,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.800000000000000,  -0.800000000000000,   0.400000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.400000000000000,   0.000000000000000,   0.815498319838598},
                                                    {-0.800000000000000,   0.000000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.800000000000000,   0.400000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.800000000000000,   0.000000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.400000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.400000000000000,   0.000000000000000,  -0.800000000000000,   0.815498319838598},
                                                    { 0.000000000000000,  -0.800000000000000,  -0.400000000000000,   0.815498319838598},
                                                    { 0.000000000000000,  -0.400000000000000,  -0.800000000000000,   0.815498319838598},
                                                    { 0.400000000000000,  -0.800000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.800000000000000,   0.000000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.400000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.800000000000000,   0.000000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.800000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.400000000000000,  -0.800000000000000,   0.815498319838598},
                                                    { 0.000000000000000,  -0.800000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.800000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.400000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.800000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.800000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.750000000000000,  -0.750000000000000,   0.250000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,  -0.250000000000000,  -0.250000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,   0.250000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.250000000000000,  -0.750000000000000,  -0.250000000000000,  -0.280203089091978},
                                                    {-0.250000000000000,  -0.250000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    { 0.250000000000000,  -0.750000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,  -0.750000000000000,  -0.250000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,  -0.250000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.250000000000000,  -0.750000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,  -0.750000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.666666666666667,  -0.666666666666667,   0.000000000000000,   0.032544642857143},
                                                    {-0.666666666666667,   0.000000000000000,  -0.666666666666667,   0.032544642857143},
                                                    { 0.000000000000000,  -0.666666666666667,  -0.666666666666667,   0.032544642857143},
                                                    {-0.666666666666667,  -0.666666666666667,  -0.666666666666667,   0.032544642857143},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  -0.000752498530276}};

    const PetscReal nord12_3DGaussPoints[210][4] = {{-0.875000000000000,  -0.875000000000000,   0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,   0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.125000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.125000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,   0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.375000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.125000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,   0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,   0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.625000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.375000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.375000000000000,   0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.875000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.625000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.375000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    { 0.375000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.625000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,   0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.125000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.375000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,   0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.625000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.375000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.857142857142857,  -0.857142857142857,   0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,   0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,  -0.285714285714286,   0.000000000000000, -1.185481802234117},
                                                    {-0.857142857142857,   0.000000000000000,  -0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,   0.285714285714286,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,   0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,   0.285714285714286, -1.185481802234117},
                                                    {-0.571428571428571,  -0.571428571428571,   0.000000000000000, -1.185481802234117},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.285714285714286, -1.185481802234117},
                                                    {-0.571428571428571,   0.000000000000000,  -0.571428571428571, -1.185481802234117},
                                                    {-0.571428571428571,   0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    {-0.285714285714286,  -0.857142857142857,   0.000000000000000, -1.185481802234117},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.285714285714286, -1.185481802234117},
                                                    {-0.285714285714286,  -0.285714285714286,  -0.571428571428571, -1.185481802234117},
                                                    {-0.285714285714286,   0.000000000000000,  -0.857142857142857, -1.185481802234117},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.285714285714286, -1.185481802234117},
                                                    { 0.000000000000000,  -0.571428571428571,  -0.571428571428571, -1.185481802234117},
                                                    { 0.000000000000000,  -0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    { 0.285714285714286,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    { 0.285714285714286,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    { 0.571428571428571,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,   0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,   0.000000000000000, -1.185481802234117},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,   0.000000000000000,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,   0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,   0.000000000000000, -1.185481802234117},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.285714285714286, -1.185481802234117},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.571428571428571, -1.185481802234117},
                                                    {-0.571428571428571,   0.000000000000000,  -0.857142857142857, -1.185481802234117},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.285714285714286, -1.185481802234117},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.571428571428571, -1.185481802234117},
                                                    {-0.285714285714286,  -0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    { 0.000000000000000,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    { 0.285714285714286,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,   0.000000000000000, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,   0.000000000000000,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.285714285714286, -1.185481802234117},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.571428571428571, -1.185481802234117},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.833333333333333,  -0.833333333333333,   0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,  -0.500000000000000,   0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,   0.166666666666667,  -0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,   0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    {-0.500000000000000,  -0.833333333333333,   0.166666666666667,  1.198527187098616},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.166666666666667,  1.198527187098616},
                                                    {-0.500000000000000,  -0.166666666666667,  -0.500000000000000,  1.198527187098616},
                                                    {-0.500000000000000,   0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.166666666666667,  1.198527187098616},
                                                    {-0.166666666666667,  -0.500000000000000,  -0.500000000000000,  1.198527187098616},
                                                    {-0.166666666666667,  -0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    { 0.166666666666667,  -0.833333333333333,  -0.500000000000000,  1.198527187098616},
                                                    { 0.166666666666667,  -0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    { 0.500000000000000,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.833333333333333,  -0.833333333333333,   0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,   0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.166666666666667,  1.198527187098616},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  1.198527187098616},
                                                    {-0.500000000000000,  -0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.500000000000000,  1.198527187098616},
                                                    {-0.166666666666667,  -0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    { 0.166666666666667,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.500000000000000,  1.198527187098616},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.800000000000000,  -0.800000000000000,   0.400000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.400000000000000,   0.000000000000000, -0.522755333229870},
                                                    {-0.800000000000000,   0.000000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.800000000000000,   0.400000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.800000000000000,   0.000000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.400000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.400000000000000,   0.000000000000000,  -0.800000000000000, -0.522755333229870},
                                                    { 0.000000000000000,  -0.800000000000000,  -0.400000000000000, -0.522755333229870},
                                                    { 0.000000000000000,  -0.400000000000000,  -0.800000000000000, -0.522755333229870},
                                                    { 0.400000000000000,  -0.800000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.800000000000000,   0.000000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.400000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.800000000000000,   0.000000000000000 , -0.800000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.800000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.400000000000000,  -0.800000000000000, -0.522755333229870},
                                                    { 0.000000000000000,  -0.800000000000000 , -0.800000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.800000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.400000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.800000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.800000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.750000000000000,  -0.750000000000000,   0.250000000000000,  0.093401029697326},
                                                    {-0.750000000000000,  -0.250000000000000,  -0.250000000000000,  0.093401029697326},
                                                    {-0.750000000000000,   0.250000000000000,  -0.750000000000000,  0.093401029697326},
                                                    {-0.250000000000000,  -0.750000000000000,  -0.250000000000000,  0.093401029697326},
                                                    {-0.250000000000000,  -0.250000000000000,  -0.750000000000000,  0.093401029697326},
                                                    { 0.250000000000000,  -0.750000000000000 , -0.750000000000000,  0.093401029697326},
                                                    {-0.750000000000000,  -0.750000000000000,  -0.250000000000000,  0.093401029697326},
                                                    {-0.750000000000000,  -0.250000000000000,  -0.750000000000000,  0.093401029697326},
                                                    {-0.250000000000000,  -0.750000000000000,  -0.750000000000000,  0.093401029697326},
                                                    {-0.750000000000000,  -0.750000000000000,  -0.750000000000000,  0.093401029697326},
                                                    {-0.666666666666667,  -0.666666666666667,   0.000000000000000, -0.005325487012987},
                                                    {-0.666666666666667,   0.000000000000000,  -0.666666666666667, -0.005325487012987},
                                                    { 0.000000000000000,  -0.666666666666667,  -0.666666666666667, -0.005325487012987},
                                                    {-0.666666666666667,  -0.666666666666667,  -0.666666666666667, -0.005325487012987},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  0.000050166568685}};


    switch (numPoints)
    {
        case 1:
            PetscCall(renormalization3DGaussPoints(numPoints, nord1_3DGaussPoints, points, weights));
            break;
        case 4:
            PetscCall(renormalization3DGaussPoints(numPoints, nord2_3DGaussPoints, points, weights));
            break;
        case 5:
            PetscCall(renormalization3DGaussPoints(numPoints, nord3_3DGaussPoints, points, weights));
            break;
        case 11:
            PetscCall(renormalization3DGaussPoints(numPoints, nord4_3DGaussPoints, points, weights));
            break;
        case 14:
            PetscCall(renormalization3DGaussPoints(numPoints, nord5_3DGaussPoints, points, weights));
            break;
        case 24:
            PetscCall(renormalization3DGaussPoints(numPoints, nord6_3DGaussPoints, points, weights));
            break;
        case 31:
            PetscCall(renormalization3DGaussPoints(numPoints, nord7_3DGaussPoints, points, weights));
            break;
        case 43:
            PetscCall(renormalization3DGaussPoints(numPoints, nord8_3DGaussPoints, points, weights));
            break;
        case 53:
            PetscCall(renormalization3DGaussPoints(numPoints, nord9_3DGaussPoints, points, weights));
            break;
        case 126:
            PetscCall(renormalization3DGaussPoints(numPoints, nord10_3DGaussPoints, points, weights));
            break;
        case 210:
            PetscCall(renormalization3DGaussPoints(numPoints, nord12_3DGaussPoints, points, weights));
            break;
        default:
            break;
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode AffineTetrahedron(PetscReal X[NUM_DIMENSIONS], PetscReal Lam[4], PetscReal DLam[NUM_DIMENSIONS][4]){
    /*Compute affine coordinates and their gradients.

    :param ndarray X: point coordinates
    :return: affine coordinates and gradients of affine coordinates
    :rtype: ndarray

    .. note:: References:\n
       Fuentes, F., Keith, B., Demkowicz, L., & Nagaraj, S. (2015). Orientation
       embedded high order shape functions for the exact sequence elements of
       all shapes. Computers & Mathematics with applications, 70(4), 353-458.
    */
    PetscFunctionBeginUser;
    
    // Define affine coordinates
    Lam[0] = 1.-X[0]-X[1]-X[2];
    Lam[1] = X[0];
    Lam[2] = X[1];
    Lam[3] = X[2];
    
    // and their gradients
    DLam[0][0] = -1;
    DLam[0][1] =  1;
    DLam[1][0] = -1;
    DLam[1][2] =  1;
    DLam[2][0] = -1;
    DLam[2][3] =  1;

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ProjectTetE(PetscReal Lam[4], PetscReal DLam[NUM_DIMENSIONS][4], PetscReal LampE[NUM_EDGES_PER_ELEMENT][2], PetscReal DLampE[NUM_EDGES_PER_ELEMENT][NUM_DIMENSIONS][2], PetscBool* IdecE){
    /*Projection of tetrahedral edges in concordance with numbering of topological entities (vertices, edges, faces).

    :param ndarray Lam: affine coordinates
    :param ndarray DLam: gradients of affine coordinates
    :return: projection of affine coordinates on edges, projection of gradients of affine coordinates on edges
    :rtype: ndarray

    .. note:: References:\n
       Fuentes, F., Keith, B., Demkowicz, L., & Nagaraj, S. (2015). Orientation
       embedded high order shape functions for the exact sequence elements of
       all shapes. Computers & Mathematics with applications, 70(4), 353-458.
    */
    PetscFunctionBeginUser;
    
    // ---------------------------------------------------------------
    // Compute projection
    // ---------------------------------------------------------------
    // e=1 --> edge01 with local orientation v0->v1
    // PetscInt e = 0;
    // LampE[e][0] = Lam[0];
    // LampE[e][1] = Lam[1];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][0];
    //     DLampE[e][i][1] = DLam[i][1];
    // }
    
    // // e=2 --> edge12 with local orientation v1->v2
    // e = 1;
    // LampE[e][0] = Lam[1];
    // LampE[e][1] = Lam[2];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][1];
    //     DLampE[e][i][1] = DLam[i][2];
    // }
    
    // // e=3 --> edge20 with local orientation v0->v2
    // e = 2;
    // LampE[e][0] = Lam[0];
    // LampE[e][1] = Lam[2];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][0];
    //     DLampE[e][i][1] = DLam[i][2];
    // }
    
    // // e=4 --> edge03 with local orientation v0->v3
    // e = 3;
    // LampE[e][0] = Lam[0];
    // LampE[e][1] = Lam[3];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][0];
    //     DLampE[e][i][1] = DLam[i][3];
    // }
    
    // // e=5 --> edge13 with local orientation v1->v3
    // e = 4;
    // LampE[e][0] = Lam[1];
    // LampE[e][1] = Lam[3];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][1];
    //     DLampE[e][i][1] = DLam[i][3];
    // }
    
    // // e=6 --> edge23 with local orientation v2->v3
    // e = 5;
    // LampE[e][0] = Lam[2];
    // LampE[e][1] = Lam[3];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][2];
    //     DLampE[e][i][1] = DLam[i][3];
    // }

    PetscInt e = 0;
    LampE[e][0] = Lam[1];
    LampE[e][1] = Lam[0];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][1];
        DLampE[e][i][1] = DLam[i][0];
    }
    
    // e=2 --> edge12 with local orientation v1->v2
    e = 1;
    LampE[e][0] = Lam[1];
    LampE[e][1] = Lam[2];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][1];
        DLampE[e][i][1] = DLam[i][2];
    }
    
    // e=3 --> edge20 with local orientation v0->v2
    e = 2;
    LampE[e][0] = Lam[2];
    LampE[e][1] = Lam[0];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][2];
        DLampE[e][i][1] = DLam[i][0];
    }
    
    // e=4 --> edge03 with local orientation v0->v3
    e = 3;
    LampE[e][0] = Lam[0];
    LampE[e][1] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][0];
        DLampE[e][i][1] = DLam[i][3];
    }
    
    // e=5 --> edge13 with local orientation v1->v3
    e = 4;
    LampE[e][0] = Lam[3];
    LampE[e][1] = Lam[1];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][3];
        DLampE[e][i][1] = DLam[i][1];
    }
    
    // e=6 --> edge23 with local orientation v2->v3
    e = 5;
    LampE[e][0] = Lam[2];
    LampE[e][1] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][2];
        DLampE[e][i][1] = DLam[i][3];
    }

    /* Projected coordinates are Lam, so IdecE=false for all edges */
    *IdecE = PETSC_FALSE;

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode OrientE(PetscReal S[2], PetscReal DS[NUM_DIMENSIONS][2], PetscInt Nori, PetscReal GS[2], PetscReal GDS[NUM_DIMENSIONS][2]){
    /*Compute the local to global transformations of edges.

    :param ndarray S: projection of affine coordinates on edges
    :param ndarray DS: projection of gradients of affine coordinates on edges
    :param ndarray Nori: edge orientation
    :param int N: number of dimensions
    :return: global transformation of edges and global transformation of gradients of edges
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    // ---------------------------------------------------------------
    // Initialization
    // ---------------------------------------------------------------
    // Allocate
    PetscInt Or[2][2];
    
    // Nori=0 => (s0,s1)->(s0,s1)
    Or[0][0] = 0;
    Or[0][1] = 1;
    // Nori=1 => (s0,s1)->(s1,s0)
    Or[1][0] = 1;
    Or[1][1] = 0;

    // Local-to-global transformation
    GS[0] = S[Or[Nori][0]];
    GS[1] = S[Or[Nori][1]];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        GDS[i][0] = DS[i][Or[Nori][0]];
        GDS[i][1] = DS[i][Or[Nori][1]];
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode PolyLegendre(PetscReal X, PetscReal T, PetscInt Nord, PetscReal P[]){
    /*Compute values of shifted scaled Legendre polynomials.

    :param ndarray X: coordinate from [0,1]
    :param float T: scaling parameter
    :param int Nord: polynomial order
    :return: polynomial values
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    // i stands for the order of the polynomial, stored in P(i)
    // lowest order case (order 0)
    P[0] = 1.;
    
    PetscReal y;
    // First order case (order 1) if necessary
    if(Nord >= 1){
        y = 2.*X - T;
        P[1] = y;
    }
  
    if(Nord >= 2){
        PetscReal tt = pow(T,2);
        for(PetscInt i = 1; i < Nord - 1; i++){
            P[i+1] = (2.*i+1)*y*P[i] - (i)*tt*P[i-1];
            P[i+1] = P[i+1]/(PetscReal)(i+1);
        }
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode PolyJacobi(PetscReal X, PetscReal T, PetscInt Nord, PetscInt Minalpha, PetscReal **P){
    /*Compute values of shifted scaled Jacobi polynomials P**alpha-i.

    Result is a half of a matrix with each row associated to a fixed alpha.
    Alpha grows by 2 in each row.

    :param ndarray X: coordinate from [0,1]
    :param float T: scaling parameter
    :param int Nord: max polynomial order
    :param int Minalpha: first row value of alpha (integer)
    :return: polynomial values
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    PetscReal *alpha;

    // Allocate
    PetscCall(PetscCalloc1(Nord + 1, &alpha));

    // Clearly (minI,maxI)=(0,Nord), but the syntax is written as it is
    // because it reflects how the indexing is called from outside
    PetscInt minI = 0;
    PetscInt maxI = minI + Nord;
    
    for(PetscInt i = 0; i < maxI + 1; i++){
        alpha[i] = Minalpha + 2*(i - minI);
    }
    
    // Initiate first column (order 0)
    for(PetscInt i = minI; i < maxI + 1; i++){
        P[i][0] = 1.;
    }
    
    PetscReal y;
    // Initiate second column (order 1) if necessary
    if(Nord >= 1){
        y = 2*X - T;
        for(PetscInt i = minI; i < maxI; i++){
            P[i][1] = y + alpha[i]*X;
        }
    }
    
    // Fill the last columns if necessary
    if(Nord >= 2){
        PetscReal tt = pow(T, 2);
        PetscInt ni = -1;
        for(PetscInt i = 0; i < maxI - 1; i++){
            PetscReal al = alpha[i];
            PetscReal aa = pow(al, 2);
            ni += 1;
            // Use recursion in order, i, to compute P^alpha_i for i>=2
            for(PetscInt j = 2; j < Nord - ni + 1; j++){
                PetscReal ai = 2*j*(j+al)*(2*j+al-2);
                PetscReal bi = 2*j+al-1;
                PetscReal ci = (2*j+al)*(2*j+al-2);
                PetscReal di = 2*(j+al-1)*(j-1)*(2*j+al);
                
                P[i][j] = bi*(ci*y+aa*T)*P[i][j-1]-di*tt*P[i][j-2];
                P[i][j] = P[i][j]/ai;
            }
        }
    }

    // Free memory
    PetscCall(PetscFree(alpha));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode HomLegendre(PetscReal S[2], PetscInt Nord, PetscReal HomP[]){
    /*Compute values of homogenized Legendre polynomials.

    :param ndarray S: affine(like) coordinates
    :param int Nord: polynomial order
    :return: polynomial values
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    PetscCall(PolyLegendre(S[1], S[0] + S[1], Nord, HomP));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode AncEE(PetscReal S[2], PetscReal DS[NUM_DIMENSIONS][2], PetscInt Nord, PetscBool Idec, PetscReal **EE, PetscReal **CurlEE){
    /*Compute edge Hcurl ancillary functions and their curls.

    :param ndarray S: affine coordinates associated to edge
    :param ndarray DS: derivatives of S in R^N
    :param int Nord: polynomial order
    :param bool Idec: Binary flag
    :param int N: spatial dimension
    :return: edge Hcurl ancillary functions, curls of edge Hcurl ancillary functions
    :rtype: ndarray

    .. note:: References:\n
       Idec: = FALSE  s0+s1 != 1
             = TRUE   s0+s1  = 1
    */
    PetscFunctionBeginUser;

    // Local parameters
    PetscInt minI = 1;
    PetscInt maxI = Nord;
    PetscInt Ncurl = 2*NUM_DIMENSIONS-3;
    
    PetscReal *homP;
    
    if (Nord <= 1){
        PetscCall(PetscCalloc1(Nord+1, &homP));
    } else{
        PetscCall(PetscCalloc1(Nord, &homP));
    } 
    
    // Extract homogenized Legendre polyomials first
    PetscCall(HomLegendre(S, maxI, homP));
    
    // Simplified case
    if(Idec){
        for(PetscInt i = minI; i < maxI + 1; i++){
            for(PetscInt j = 0; j < NUM_DIMENSIONS; j++){
                EE[j][i-1] = homP[i-1]*DS[j][1];
            }
        }
        // No need to compute Whitney function or curl
        for (PetscInt i = 0; i < Ncurl; ++i) {
            for (PetscInt j = minI - 1; j < maxI - 1; ++j) {
                CurlEE[i][j] = 0;
            }
        }
    } else {
        // Lowest order Whitney function and its curl
        PetscReal whiE[NUM_DIMENSIONS];
        for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
            whiE[i] = S[0]*DS[i][1] - S[1]*DS[i][0];
        }
        PetscReal curlwhiE[NUM_DIMENSIONS];
        PetscReal temp1[NUM_DIMENSIONS], temp2[NUM_DIMENSIONS];
        for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
            temp1[i] = DS[i][0];
            temp2[i] = DS[i][1];
        }
        PetscCall(crossProduct(temp1, temp2, curlwhiE));
        
        // Now construct the higher order elements
        for(PetscInt i = minI; i < maxI + 1; i++){
            for(PetscInt j = 0; j < NUM_DIMENSIONS; j++){
                EE[j][i-1] = homP[i-1]*whiE[j];
            }
            for(PetscInt j = 0; j < Ncurl; j++){
                CurlEE[j][i-1] = (i+1)*homP[i-1]*curlwhiE[j];
            }
        }
    }

    PetscCall(PetscFree(homP));   
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ProjectTetF(PetscReal Lam[4], PetscReal DLam[NUM_DIMENSIONS][4], PetscReal LampF[NUM_FACES_PER_ELEMENT][NUM_DIMENSIONS], PetscReal DLampF[NUM_FACES_PER_ELEMENT][NUM_DIMENSIONS][NUM_DIMENSIONS], PetscBool* IdecF){
    /*Projection of tetrahedral faces in concordance with numbering of topological entities (vertices, edges, faces).

    :param ndarray Lam: affine coordinates
    :param ndarray DLam: gradients of affine coordinates
    :return: projection of affine coordinates on faces, projection of gradients of affine coordinates on faces
    :rtype: ndarray
    */
    PetscFunctionBeginUser;
    
    // ---------------------------------------------------------------
    // Compute projection
    // ---------------------------------------------------------------
    // f=1 --> face012 with local orientation v0->v1->v2
    PetscInt f = 0;
    LampF[f][0] = Lam[0];
    LampF[f][1] = Lam[1];
    LampF[f][2] = Lam[2];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampF[f][i][0] = DLam[i][0];
        DLampF[f][i][1] = DLam[i][1];
        DLampF[f][i][2] = DLam[i][2];
    }
    
    // f=2 --> face013 with local orientation v0->v1->v3
    f = 1;
    LampF[f][0] = Lam[0];
    LampF[f][1] = Lam[1];
    LampF[f][2] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampF[f][i][0] = DLam[i][0];
        DLampF[f][i][1] = DLam[i][1];
        DLampF[f][i][2] = DLam[i][3];
    }
    
    // f=3 --> face123 with local orientation v1->v2->v3
    f = 2;
    LampF[f][0] = Lam[1];
    LampF[f][1] = Lam[2];
    LampF[f][2] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampF[f][i][0] = DLam[i][1];
        DLampF[f][i][1] = DLam[i][2];
        DLampF[f][i][2] = DLam[i][3];
    }
    
    // f=4 --> face023 with local orientation v0->v2->v3
    f = 3;
    LampF[f][0] = Lam[0];
    LampF[f][1] = Lam[2];
    LampF[f][2] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampF[f][i][0] = DLam[i][0];
        DLampF[f][i][1] = DLam[i][2];
        DLampF[f][i][2] = DLam[i][3];
    }
    
    *IdecF = PETSC_FALSE;

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode OrientTri(PetscReal S[NUM_DIMENSIONS], PetscReal DS[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscInt Nori, PetscReal GS[NUM_DIMENSIONS], PetscReal GDS[NUM_DIMENSIONS][NUM_DIMENSIONS]){
    /*Compute the local to global transformations of edges.

    :param ndarray S: projection of affine coordinates on faces
    :param ndarray DS: projection of gradients of affine coordinates on faces
    :param ndarray Nori: face orientation
    :param int N: number of dimensions
    :return: global transformation of faces and global transformation of gradients of faces
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    PetscInt Or[NUM_DIMENSIONS*2][NUM_DIMENSIONS];
    
    // Nori=0 => (s0,s1,s2)->(s0,s1,s2)
    Or[0][0] = 0;
    Or[0][1] = 1;
    Or[0][2] = 2;
    // Nori=1 => (s0,s1,s2)->(s1,s2,s0)
    Or[1][0] = 1;
    Or[1][1] = 2;
    Or[1][2] = 0;
    // Nori=2 => (s0,s1,s2)->(s2,s0,s1)
    Or[2][0] = 2;
    Or[2][1] = 0;
    Or[2][2] = 1;
    // Nori=3 => (s0,s1,s2)->(s0,s2,s1)
    Or[3][0] = 0;
    Or[3][1] = 2;
    Or[3][2] = 1;
    // Nori=4 => (s0,s1,s2)->(s1,s0,s2)
    Or[4][0] = 1;
    Or[4][1] = 0;
    Or[4][2] = 2;
    // Nori=5 => (s0,s1,s2)->(s2,s1,s0)
    Or[5][0] = 2;
    Or[5][1] = 1;
    Or[5][2] = 0;
    
    // Local-to-global transformation
    GS[0] = S[Or[Nori][0]];
    GS[1] = S[Or[Nori][1]];
    GS[2] = S[Or[Nori][2]];
    
    for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        GDS[i][0] = DS[i][Or[Nori][0]];
        GDS[i][1] = DS[i][Or[Nori][1]];
        GDS[i][2] = DS[i][Or[Nori][2]];
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode PolyIJacobi(PetscReal X, PetscReal T, PetscInt Nord, PetscInt Minalpha, PetscBool Idec, PetscReal **L, PetscReal **P, PetscReal **R){
    /*Compute values of integrated shifted scaled Jacobi polynomials and their derivatives starting with p=1.

    Result is 'half' of a  matrix with each row  associated to a fixed alpha.
    Alpha grows by 2 in each row.

    :param ndarray X: coordinate from [0,1]
    :param ndarray T: scaling parameter
    :param int Nord: max polynomial order
    :param int Minalpha: first row value of alpha
    :param bool Idec: decision flag to compute (= FALSE polynomials with x and t derivatives, = TRUE  polynomials with x derivatives only)
    :return: polynomial values, derivatives in x (Jacobi polynomials), derivatives in t
    */
    PetscFunctionBeginUser;
    
    // clearly (minI,maxI)=(1,Nord), but the syntax is written as it is
    // because it reflects how the indexing is called from outside
    PetscInt minI = 0;
    PetscInt maxI = minI + Nord;
    PetscReal *alpha;
    PetscReal **ptemp;

    // Allocate 
    PetscCall(PetscCalloc1(Nord, &alpha));

    PetscCall(PetscCalloc1(Nord + 1, &ptemp));
    for (PetscInt i = 0; i < Nord + 1; i++){
        PetscCall(PetscCalloc1(Nord + 1, &ptemp[i]));
    }
    
    PetscCall(PolyJacobi(X, T, Nord, Minalpha, ptemp));
    
    // Define P. Note that even though P is defined at all entries,
    // because of the way Jacobi computes ptemp, only the necessary entries,
    // and those on the first subdiagonal (which are never used later)
    // are actually accurate.
    for(PetscInt i = minI; i < maxI; i++){
        for(PetscInt j = 0; j < Nord; j++){
            P[i][j] = ptemp[i][j];
        }
    }

    // Create vector alpha first
    for(PetscInt i = 0; i < maxI; i++){
        alpha[i] = Minalpha + 2*(i - minI);
    }
    
    // Initiate first column (order 1 in L)
    for(PetscInt i = minI; i < maxI; i++){
        L[i][0] = X;
    }
    
    // General case; compute R
    for(PetscInt i = minI; i < maxI; i++){
        for(PetscInt j = 0; j < Nord; j++){
            R[i][j] = 0;
        }
    }
    
    // Fill the last columns if necessary
    if(Nord >= 2){
        PetscReal tt = pow(T, 2);
        PetscInt ni = -1;
        for(PetscInt i = 0; i < maxI - 1; i++){
            PetscReal al = alpha[i];
            ni += 1;
            for(PetscInt j = 2; j < Nord - ni + 1; j++){
                PetscReal tia = j+j+al;
                PetscReal tiam1 = tia-1;
                PetscReal tiam2 = tia-2;
                PetscReal ai = (j+al)/(tiam1*tia);
                PetscReal bi = (al)/(tiam2*tia);
                PetscReal ci = (j-1)/(tiam2*tiam1);
                L[i][j-1] = ai*ptemp[i][j]+bi*T*ptemp[i][j-1]-ci*tt*ptemp[i][j-2];
                R[i][j-1] = -(j-1)*(ptemp[i][j-1]+T*ptemp[i][j-2]);
                R[i][j-1] = R[i][j-1]/tiam2;
            }
        }
    }
    
    // Free memory
    PetscCall(PetscFree(alpha));
    
    for (PetscInt i = 0; i < Nord + 1; i++){
        PetscCall(PetscFree(ptemp[i]));
    }
    PetscCall(PetscFree(ptemp));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode HomIJacobi(PetscReal S[2], PetscReal DS[NUM_DIMENSIONS][2], PetscInt Nord, PetscInt Minalpha, PetscBool Idec, PetscReal **HomL, PetscReal ***DHomL){
    /*Compute values of integrated homogenized Jacobi polynomials and their gradients.
    Result is half of a  matrix with each row  associated to a fixed alpha.
    Alpha grows by 2 in each row.

    :param ndarray S: (s0,s1) affine(like) coordinates
    :param ndarray DS: gradients of S in R(N)
    :param int Nord: max polynomial order
    :param int Minalpha: first row value of alpha (integer)
    :param bool Idec: decision flag to compute
    :return: polynomial values and derivatives in x (Jacobi polynomials)
    :rtype: ndarray
    */
    PetscFunctionBeginUser;
    
    // clearly (minI,maxI)=(1,Nord), but the syntax is written as it is
    // because it reflects how the indexing is called from outside
    PetscInt minI = 1;
    PetscInt maxI = minI+Nord-1;
    
    PetscReal **homP;   // homP[Nord][Nord]
    PetscReal **homR;   // homR[Nord][Nord] 

    // Allocate
    PetscCall(PetscCalloc1(Nord, &homP));
    for (PetscInt i = 0; i < Nord; i++){
        PetscCall(PetscCalloc1(Nord, &homP[i]));
    }
     
    PetscCall(PetscCalloc1(Nord, &homR));
    for (PetscInt i = 0; i < Nord; i++){
        PetscCall(PetscCalloc1(Nord, &homR[i]));
    }
        
    PetscInt ni = -1;
    
    if(Idec){
        PetscCall(PolyIJacobi(S[1], 1, Nord, Minalpha, Idec, HomL, homP, homR));
        for(PetscInt i = minI; i < maxI + 1; i++){
            ni += 1;
            for(PetscInt j = 1; j < Nord - ni + 1; j++){
                for(PetscInt k = 0; k < NUM_DIMENSIONS; k++){
                    DHomL[k][i-1][j-1] = homP[i-1][j-1] * DS[k][1];
                }
            }
        }
    } else {
        // If sum of S different from 1 -> Idec=.FALSE.
        
        PetscCall(PolyIJacobi(S[1], S[0] + S[1], Nord, Minalpha, Idec, HomL, homP, homR));
        
        PetscReal DS01[NUM_DIMENSIONS];
        for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
            DS01[i] = DS[i][0] + DS[i][1];
        }
        
        for(PetscInt i = minI; i < maxI + 1; i++){
            ni += 1;
            for(PetscInt j = 1; j < Nord - ni + 1; j++){
                for(PetscInt k = 0; k < NUM_DIMENSIONS; k++){
                    DHomL[k][i-1][j-1] = homP[i-1][j-1] * DS[k][1] + homR[i-1][j-1]*DS01[k];
                }
            }
        }    
    }
    
    // Free memory
    for (PetscInt i = 0; i < Nord; i++){
        PetscCall(PetscFree(homP[i]));
    }
    PetscCall(PetscFree(homP));

    for (PetscInt i = 0; i < Nord; i++){
        PetscCall(PetscFree(homR[i]));
    }
    PetscCall(PetscFree(homR));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode AncETri(PetscReal S[NUM_DIMENSIONS], PetscReal DS[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscInt Nord, PetscBool Idec, PetscReal ***ETri, PetscReal ***CurlETri){
    /*Compute triangle face Hcurl ancillary functions and their curls.

    :param ndarray S: (s0,s1,s2) affine coordinates associated to triangle face
    :param ndarray DS: derivatives of S0,S1,S2
    :param int Nord: polynomial order
    :param bool Idec: Binary flag:
    :param int N: spatial dimension
    :return: triangle Hcurl ancillary functions and curls of triangle Hcurl ancillary functions
    :rtype: ndarray
    */
    PetscFunctionBeginUser;
    
    PetscReal DsL[NUM_DIMENSIONS][2];
    PetscReal sL[2];

    // Local parameters
    PetscInt minI = 0;
    PetscInt minJ = 1;
    PetscInt maxJ = Nord-1;
    PetscInt maxIJ = Nord-1;
    PetscInt minalpha = 2*minI+1;
    PetscInt Ncurl = 2*NUM_DIMENSIONS-3;
    PetscBool IdecE = PETSC_FALSE;
    
    // get EE - this is never a simplified case (IdecE=0)
    PetscReal tempS[2] = {S[0], S[1]};
    PetscReal tempDS[NUM_DIMENSIONS][2];
    for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        tempDS[i][0] = DS[i][0];
        tempDS[i][1] = DS[i][1];
    }
    
    PetscReal **EE;     // EE[NUM_DIMENSIONS][Nord-minJ]
    PetscReal **curlEE; // curlEE[2*NUM_DIMENSIONS-3][Nord-minJ]

    // Allocate
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &EE));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(Nord-minJ, &EE[i]));
    }

    PetscCall(PetscCalloc1(2*NUM_DIMENSIONS-3, &curlEE));
    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        PetscCall(PetscCalloc1(Nord-minJ, &curlEE[i]));
    }

    PetscCall(AncEE(tempS, tempDS, Nord-minJ, IdecE, EE, curlEE));
    
    // get homogenized Integrated Jacobi polynomials, homLal, and gradients
    sL[0] = S[0]+S[1];
    sL[1] = S[2];
    for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        DsL[i][0] = DS[i][0] + DS[i][1];
        DsL[i][1] = DS[i][2];
    }
    
    PetscReal **homLal;     // homLal[maxJ][maxJ]
    PetscReal ***DhomLal;   // DhomLal[NUM_DIMENSIONS][maxJ][maxJ]

    // Allocate
    PetscCall(PetscCalloc1(maxJ, &homLal));
    for (PetscInt i = 0; i < maxJ; i++){
        PetscCall(PetscCalloc1(maxJ, &homLal[i]));
    }

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &DhomLal));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(maxJ, &DhomLal[i]));
        for (PetscInt j = 0; j < maxJ; j++){
            PetscCall(PetscCalloc1(maxJ, &DhomLal[i][j]));
        }
    }
    
    PetscCall(HomIJacobi(sL, DsL, maxJ, minalpha, Idec, homLal, DhomLal));

    // Simply complete the required information
    for(PetscInt i = 0; i < maxIJ + 1; i++){
        for(PetscInt j = minI; j < i-minJ+1; j++){
            PetscInt k = i - j;
            for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                ETri[n][j][k-1] = EE[n][j] * homLal[j][k-1];
            }
            
            PetscReal DhomLalxEE[NUM_DIMENSIONS];
            PetscReal temp1[NUM_DIMENSIONS], temp2[NUM_DIMENSIONS];
            for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                temp1[n] = DhomLal[n][j][k-1];
                temp2[n] = EE[n][j];
            }
            
            PetscCall(crossProduct(temp1, temp2, DhomLalxEE));
            
            for(PetscInt n = 0; n < Ncurl; n++){
                CurlETri[n][j][k-1] = homLal[j][k-1]*curlEE[n][j] + DhomLalxEE[n];
            }
        }
    }
    
    // Free memory
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(EE[i]));   
    }
    PetscCall(PetscFree(EE));   

    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        PetscCall(PetscFree(curlEE[i]));   
    }
    PetscCall(PetscFree(curlEE));   

    for (PetscInt i = 0; i < maxJ; i++){
        PetscCall(PetscFree(homLal[i]));   
    }
    PetscCall(PetscFree(homLal));

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < maxJ; j++){
            PetscCall(PetscFree(DhomLal[i][j]));   
            }
        PetscCall(PetscFree(DhomLal[i]));       
    }
    PetscCall(PetscFree(DhomLal));       

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode computeCellOrientation(DM dm, PetscInt cell, PetscInt cellOrientation[10]){
    
    PetscFunctionBeginUser;


    PetscInt cellFaces[NUM_FACES_PER_ELEMENT];
    PetscInt cellEdges[NUM_EDGES_PER_ELEMENT];
    PetscInt cellVertices[NUM_EDGES_PER_ELEMENT];
    PetscInt faceEdges[NUM_FACES_PER_ELEMENT][NUM_EDGES_PER_FACE];
    PetscInt faceVertices[NUM_FACES_PER_ELEMENT][NUM_VERTICES_PER_FACE];
    PetscInt edgeVertices[NUM_EDGES_PER_ELEMENT][NUM_VERTICES_PER_EDGE];

    PetscInt transitiveClosureCellSize;
    PetscInt  *transitiveClosureCellPoints = NULL;


    PetscInt transitiveClosureFaceSize;
    PetscInt *transitiveClosureFacePoints = NULL;
    const PetscInt *conePoints;
    PetscInt currentPoint;
    PetscInt currentFace;   

    
    /* Get transitive clousure for cell */
    /* Orden convention:
        - Faces indices start on position 2
        - Edges indices start on position 2 + NUM_FACES_PER_ELEMENT * 2
        - Vertices indices start on position 2 + NUM_FACES_PER_ELEMENT * 2 + NUM_EDGES_PER_ELEMENT * 2
    */
    PetscCall(DMPlexGetTransitiveClosure(dm, cell, PETSC_TRUE, &transitiveClosureCellSize, &transitiveClosureCellPoints));

    /* Get faces indices for cell, edges for each face, and vertices for each face */
    currentPoint = 2; 
    for(PetscInt i = 0; i < NUM_FACES_PER_ELEMENT; i++) {
            /* Face indices */
            cellFaces[i] = transitiveClosureCellPoints[currentPoint + i*2];
            PetscCall(DMPlexGetCone(dm, cellFaces[i], &conePoints));
            
            /* Edges for each face*/
            for(PetscInt j = 0; j < NUM_EDGES_PER_FACE; j++){
                faceEdges[i][j] = conePoints[j];
            }

            /* Vertices for each face */ 
            /* Orden convention:
                - Edges indices start on position 2
                - Vertices indices start on position 2 + NUM_EDGES_PER_FACE * 2
            */
            PetscCall(DMPlexGetTransitiveClosure(dm, cellFaces[i], PETSC_TRUE, &transitiveClosureFaceSize, &transitiveClosureFacePoints));
            
            currentFace = 8;
            for(PetscInt j = 0; j < NUM_VERTICES_PER_FACE; j++){
                faceVertices[i][j] = transitiveClosureFacePoints[currentFace + j*2];
            }
            
            PetscCall(DMPlexRestoreTransitiveClosure(dm, cellFaces[i], PETSC_TRUE, &transitiveClosureFaceSize, &transitiveClosureFacePoints));
    }

    /* Get edges indices for cell */
    currentPoint = 2 + NUM_FACES_PER_ELEMENT*2; 
    for(PetscInt i = 0; i < NUM_EDGES_PER_ELEMENT; i++) {
        cellEdges[i] = transitiveClosureCellPoints[currentPoint + i*2];
        
        /* Get edge vertices */
        PetscCall(DMPlexGetCone(dm, cellEdges[i], &conePoints));
        
        for(PetscInt j = 0; j < NUM_VERTICES_PER_EDGE; j++){
            edgeVertices[i][j] = conePoints[j];
            }
    }

    /* Get vertices indices for cell */
    currentPoint = 2 + NUM_FACES_PER_ELEMENT*2 + NUM_EDGES_PER_ELEMENT*2;
    for(PetscInt i = 0; i < NUM_VERTICES_PER_ELEMENT; i++) {
        cellVertices[i] = transitiveClosureCellPoints[currentPoint + i*2];
    }

    /* Restore transitive closure */
    PetscCall(DMPlexRestoreTransitiveClosure(dm, cell, PETSC_TRUE, &transitiveClosureCellSize, &transitiveClosureCellPoints));

    /* Compute orientation for faces */
    for(PetscInt i = 0; i < NUM_FACES_PER_ELEMENT; i++) {
        cellOrientation[i] = 0;
    }


    PetscInt localNodes[2];

    for(PetscInt i = 0; i < NUM_EDGES_PER_ELEMENT; i++) {
        switch (i) {
            case 0:
                localNodes[0] = 0; localNodes[1] = 1;
                break;
            case 1:
                localNodes[0] = 1; localNodes[1] = 2;
                break;
            case 2:
                localNodes[0] = 2; localNodes[1] = 0;
                break;
            case 3:
                localNodes[0] = 0; localNodes[1] = 3;
                break;
            case 4:
                localNodes[0] = 3; localNodes[1] = 1;
                break;
            case 5:
                localNodes[0] = 2; localNodes[1] = 3;
                break;
        }

        // Get the global node indices for this edge
        PetscInt nodesEleForThisEdge[2] = {cellVertices[localNodes[0]], cellVertices[localNodes[1]]};
        PetscInt globalNodesInEdge[2] = {edgeVertices[i][0], edgeVertices[i][1]};

        // Determine the orientation for this edge
        PetscInt orientationForThisEdge = 0;
        if ((nodesEleForThisEdge[0] == globalNodesInEdge[1]) && (nodesEleForThisEdge[1] == globalNodesInEdge[0])) {
            orientationForThisEdge = 1;  // Edge is inverted
        }
        cellOrientation[4 + i] = orientationForThisEdge;
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode shape3DETet(PetscReal X[NUM_DIMENSIONS], PetscInt nord, PetscInt cellOrientation[10], PetscReal **ShapE, PetscReal **CurlE){
    /*Compute values of 3D tetrahedron element H(curl) shape functions and their derivatives.

    :param ndarray X: master tetrahedron coordinates from (0,1)^3
    :param int Nord: polynomial order
    :param ndarray NoriE: edge orientation
    :param ndarray NoriF: face orientation
    :return: number of dof, values of the shape functions at the point, curl of the shape functions
    :rtype: ndarray.

    .. note:: References:\n
       Fuentes, F., Keith, B., Demkowicz, L., & Nagaraj, S. (2015). Orientation
       embedded high order shape functions for the exact sequence elements of
       all shapes. Computers & Mathematics with applications, 70(4), 353-458.
    */
    PetscFunctionBeginUser;
   
    // Local parameters
    PetscBool IdecB[2] = {PETSC_FALSE, PETSC_FALSE};
    PetscInt minI = 0;
    PetscInt minJ = 1;
    PetscInt minK = 1;
    PetscInt minIJ = minI + minJ;
    PetscInt minIJK = minIJ + minK;

    // Initialize counter for shape functions
    PetscInt m = 0;
    
    // Define affine coordinates and gradients
    PetscReal Lam[NUM_DIMENSIONS + 1] = {0.0};
    PetscReal DLam[NUM_DIMENSIONS][NUM_DIMENSIONS + 1] = {{0.0}};

    PetscCall(AffineTetrahedron(X, Lam, DLam));
    
    /* Shape functions over edges */
    PetscReal LampE[NUM_EDGES_PER_ELEMENT][2];
    PetscReal DLampE[NUM_EDGES_PER_ELEMENT][NUM_DIMENSIONS][2];
    PetscBool IdecE;
    
    // Compute edges projection
    PetscCall(ProjectTetE(Lam, DLam, LampE, DLampE, &IdecE));

    // Polynomial order and orientation for faces and edges 
    PetscInt NoriF[4];  // Orientation for faces
    PetscInt NoriE[6];  // Orientation for 


    // Extract orientation for faces
    for (PetscInt i = 0; i < NUM_FACES_PER_ELEMENT; ++i){
        NoriF[i] = cellOrientation[i];
    }

    /* Extract orientation for edges */
    for (PetscInt i = 0; i < NUM_EDGES_PER_ELEMENT; ++i){
        NoriE[i] = cellOrientation[i+4]; 
    }


    /* Allocate memory for shape functions on edges */
    PetscReal **EE;       // EE[NUM_DIMENSIONS][nordEdge];
    PetscReal **CurlEE;   // CurlEE[2*NUM_DIMENSIONS-3][nordEdge];

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &EE));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(nord, &EE[i]));
    }

    PetscCall(PetscCalloc1(2*NUM_DIMENSIONS-3, &CurlEE));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(nord, &CurlEE[i]));
    }


    /* Reset matrices */
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++){
            ShapE[i][j] = 0.0;
            CurlE[i][j] = 0.0;
        }
    }

    
    /* Shape functions for edges */
    PetscInt nordEdge = 0;
    PetscInt numDofEdge = 0;    
    for(PetscInt i = 0; i < NUM_EDGES_PER_ELEMENT; i++){
        // Local parameters
        nordEdge = nord;
        numDofEdge = nordEdge;
        if(numDofEdge > 0){
            // Local parameters
            PetscInt maxI = nordEdge - 1;
            // Orient 
            PetscReal GLampE[2] = {0.0};
            PetscReal GDLampE[3][2] = {{0.0}};
            PetscReal S[2] = {0.0} ;
            PetscReal D[NUM_DIMENSIONS][2] = {{0.0}};

            S[0] = LampE[i][0];
            S[1] = LampE[i][1];

            /* Extract the slice into D */
            for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
                for (PetscInt k = 0; k < 2; k++) {
                    D[j][k] = DLampE[i][j][k];
                }
            }

            PetscCall(OrientE(S, D, NoriE[i], GLampE, GDLampE));

            /* Construct the shape functions */
            PetscCall(AncEE(GLampE, GDLampE, nordEdge, IdecE, EE, CurlEE));

            for(PetscInt j = minI; j < maxI + 1; j++){
                for(PetscInt k = 0; k < NUM_DIMENSIONS; k++){
                    ShapE[k][m] = EE[k][j];
                    CurlE[k][m] = CurlEE[k][j];
                }
                m += 1;
            }
        }
    }

    // Free memory for shape functions on edges
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(EE[i]));
    }
    PetscCall(PetscFree(EE));

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(CurlEE[i]));
    }
    PetscCall(PetscFree(CurlEE));

    /* Shape functions over faces */ 
    PetscReal LampF[NUM_FACES_PER_ELEMENT][NUM_DIMENSIONS];
    PetscReal DLampF[NUM_FACES_PER_ELEMENT][NUM_DIMENSIONS][NUM_DIMENSIONS];
    PetscBool IdecF;    

    // Compute faces projection
    PetscCall(ProjectTetF(Lam, DLam, LampF, DLampF, &IdecF));

    // Allocate memory for funcions on faces
    PetscReal ***ETri;       // ETri[NUM_DIMENSIONS][nord - 1][nord - 1]
    PetscReal ***CurlETri;   // CurlETri[2*NUM_DIMENSIONS-3][nord-1][nord-1]

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &ETri));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(nord-1, &ETri[i]));
        for (PetscInt j = 0; j < nord-1; j++){
            PetscCall(PetscCalloc1(nord-1, &ETri[i][j]));
        }
    }

    PetscCall(PetscCalloc1(2*NUM_DIMENSIONS-3, &CurlETri));
    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        PetscCall(PetscCalloc1(nord-1, &CurlETri[i]));
        for (PetscInt j = 0; j < nord-1; j++){
            PetscCall(PetscCalloc1(nord-1, &CurlETri[i][j]));
        }
    }    

    // Shape functions for faces
    PetscInt nordFace = 0;
    PetscInt numDofFace = 0;

    for(PetscInt i = 0; i < NUM_FACES_PER_ELEMENT; i++){
        // Local parameters
        nordFace = nord;
        numDofFace = nordFace*(nordFace-1)/2;
        if(numDofFace > 0){
            // Local parameters (again)
            PetscInt maxIJ = nordFace - 1;

            // Orient
            PetscReal GLampF[3];
            PetscReal GDLampF[3][3];
            PetscReal tmpLampF[3];
            PetscReal tempDLampF[3][3];

            // Prepare input matrices
            for (PetscInt j = 0; j<3; j++){
                tmpLampF[j] = LampF[i][j];
                for (PetscInt k = 0; k<3; k++){
                    tempDLampF[j][k] = DLampF[i][j][k];
                }
            }

            PetscCall(OrientTri(tmpLampF, tempDLampF, NoriF[i], GLampF, GDLampF));

            // Loop over families
            PetscInt famctr = m;
            for(PetscInt j = 0; j < 2; j++){
                m = famctr + j - 1;
                PetscInt abc[3];
                for(PetscInt k = 0; k < 3; k++){
                    PetscInt pos = (k - j) % 3;
                    if (pos < 0){
                        pos += 3;
                    } 
                    abc[pos] = k;
                }
        
                PetscReal tempGLampF[3];
                PetscReal tempGDLampF[3][3];
                for(PetscInt k = 0; k < 3; k++){
                    tempGLampF[k] = GLampF[abc[k]];
                    for(PetscInt t = 0; t < NUM_DIMENSIONS; t++){
                        tempGDLampF[t][k] = GDLampF[t][abc[k]];
                    }
                }

                // Construct the shape functions
                PetscCall(AncETri(tempGLampF, tempGDLampF, nordFace, IdecF, ETri, CurlETri));

                 for(PetscInt k = minIJ; k < maxIJ + 1; k++){
                    for(PetscInt r = minI; r < k-minJ+1; r++){
                        PetscInt p = k - r;
                        m += 2;
                        for(PetscInt t = 0; t < NUM_DIMENSIONS; t++){
                            ShapE[t][m-1] = ETri[t][r][p-1];
                            CurlE[t][m-1] = CurlETri[t][r][p-1];
                        }
                    }
                }
            }
        }
    }

    // Free memory
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < nord-1; j++){
            PetscCall(PetscFree(ETri[i][j]));
        }
        PetscCall(PetscFree(ETri[i]));
    }
    PetscCall(PetscFree(ETri));   

    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        for (PetscInt j = 0; j < nord-1; j++){
            PetscCall(PetscFree(CurlETri[i][j]));
        }
        PetscCall(PetscFree(CurlETri[i]));
    }
    PetscCall(PetscFree(CurlETri));

    /* Shape functions over volume */
    PetscInt nordB = nord;
    PetscInt ndofB = nordB*(nordB-1)*(nordB-2)/6;
    PetscInt minbeta = 2*minIJ;
    PetscInt maxIJK = nordB-1;
    PetscInt maxK = maxIJK-minIJ;

    // Allocate memory for shape functions in volume
    PetscReal ***ETriV;         // ETri[NUM_DIMENSIONS][nord-minK-1][nord-minK-1]
    PetscReal ***CurlETriV;     // CurlETri[2*NUM_DIMENSIONS-3][nord-minK-1][nord-minK-1]
    PetscReal **homLbet;        // homLbet[maxK][maxK] 
    PetscReal ***DhomLbet;      // DhomLbet[NUM_DIMENSIONS][maxK][maxK]    

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &ETriV));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(nord-minK-1, &ETriV[i]));
        for (PetscInt j = 0; j < nord-minK-1; j++){
            PetscCall(PetscCalloc1(nord-minK-1, &ETriV[i][j]));
        }
    }

    PetscCall(PetscCalloc1(2*NUM_DIMENSIONS-3, &CurlETriV));
    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        PetscCall(PetscCalloc1(nord-minK-1, &CurlETriV[i]));
        for (PetscInt j = 0; j < nord-minK-1; j++){
            PetscCall(PetscCalloc1(nord-minK-1, &CurlETriV[i][j]));
        }
    }

    PetscCall(PetscCalloc1(maxK, &homLbet));
    for (PetscInt i = 0; i < maxK; i++){
        PetscCall(PetscCalloc1(maxK, &homLbet[i]));
    }

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &DhomLbet));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(maxK, &DhomLbet[i]));
        for (PetscInt j = 0; j < maxK; j++){
            PetscCall(PetscCalloc1(maxK, &DhomLbet[i][j]));
        }
    }    

    // If necessary, create bubbles
    if(ndofB > 0){
        // Local parameters (again)
        IdecB[0] = IdecF;
        IdecB[1] = PETSC_TRUE;
        
        // Loop over families
        PetscInt famctr = m;
        for(PetscInt i = 0; i < 3; i++){
            m = famctr + i - 2;
            PetscInt abcd[4];
            for(PetscInt j = 0; j < 4; j++){
                PetscInt pos = (j - i) % 4;
                if (pos < 0){
                    pos += 4;
                } 
                abcd[pos] = j;
            }
            PetscInt abc[3] = {abcd[0], abcd[1], abcd[2]};
            PetscInt d = abcd[3];

            // Now construct the shape functions (no need to orient)
            PetscReal tempLam[3];
            PetscReal tempDLam[3][3];
            for(PetscInt j = 0; j < 3; j++){
                tempLam[j] = Lam[abc[j]];
                for(PetscInt k = 0; k < NUM_DIMENSIONS; k++){
                    tempDLam[k][j] = DLam[k][abc[j]];
                }
            }
            
            PetscCall(AncETri(tempLam, tempDLam, nordB-minK, IdecB[0], ETriV, CurlETriV));

            PetscReal tmp1[2] = {1-Lam[d], Lam[d]};
            PetscReal tmp2[3][2];

            // Initialize input matrix
            for(PetscInt j = 0; j < NUM_DIMENSIONS; j++){
                tmp2[j][0] = -DLam[j][d]; 
                tmp2[j][1] = DLam[j][d]; 
            }

            PetscCall(HomIJacobi(tmp1, tmp2, maxK, minbeta, IdecB[1], homLbet, DhomLbet));

            for(PetscInt j = minIJK; j < maxIJK+1; j++){
                for(PetscInt k = minIJ; k < j-minK+1; k++){
                    for(PetscInt r = minI; r < k-minJ+1; r++){
                        PetscInt p = k - r;
                        PetscInt q = j - k;
                        m += 3;

                        for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                            ShapE[n][m-1] = ETriV[n][r][p-1]*homLbet[k-1][q-1];
                        }

                        PetscReal DhomLbetxETri[NUM_DIMENSIONS];
                        PetscReal v1[NUM_DIMENSIONS], v2[NUM_DIMENSIONS];
    
                        for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                            v1[n] = DhomLbet[n][k-1][q-1];
                            v2[n] = ETriV[n][r][p-1];
                        }
    
                        PetscCall(crossProduct(v1, v2, DhomLbetxETri));
    
                        for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                            CurlE[n][m-1] = homLbet[k-1][q-1]*CurlETriV[n][r][p-1] + DhomLbetxETri[n];
                        }
                    }
                }
            }
        }
    }

    /* Compute reordering vector (PETGEM to PETSc convention) */
    /*PetscInt tmp1[] = {1, 2, 3, 4, 5, 6};
    PetscInt tmp2[] = {13, 14, 15, 16, 19, 20, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    PetscInt tmp3[] = {43, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 37, 38, 39, 40, 41, 42, 31, 32, 33, 34, 35, 36, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    PetscInt tmp4[] = {73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
    PetscInt tmp5[] = {111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    PetscInt tmp6[] = {157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36};*/


    // PetscInt tmp1[] = {0, 1, 2, 3, 4, 5}; 
    // PetscInt tmp2[] = {12, 13, 14, 15, 18, 19, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; // Subtracted values from original array
    // PetscInt tmp3[] = {42, 43, 44, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 36, 37, 38, 39, 40, 41, 30, 31, 32, 33, 34, 35, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}; // Subtracted values from original array
    // PetscInt tmp4[] = {72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23}; // Subtracted values from original array
    // PetscInt tmp5[] = {110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29}; // Subtracted values from original array
    // PetscInt tmp6[] = {156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35}; // Subtracted values from original array

    // // Print the updated arrays
 


    // PetscInt numDofInCell; 
    // numDofInCell = nord*(nord+2)*(nord+3)/2;
    // PetscInt orderPermutation[numDofInCell]; 

    // switch (nord){
    // case 1:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp1[i];
    //     }
    //     break;
    // case 2:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp2[i];
    //     }
    //     break;
    // case 3:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp3[i];
    //     }
    //     break;
    // case 4:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp4[i];
    //     }
    //     break;
    // case 5:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp5[i];
    //     }
    //     break;
    // case 6:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp6[i];
    //     }
    //     break;        
    // }


    // /* Copy basis and curl from PETGEM order convention */
    // PetscReal tmpShapE[NUM_DIMENSIONS][numDofInCell], tmpCurlE[NUM_DIMENSIONS][numDofInCell]; 

    // for (PetscInt i = 0; i<NUM_DIMENSIONS; i++){
    //     for (PetscInt j = 0; j<numDofInCell; j++){
    //         tmpShapE[i][j] = ShapE[i][j];
    //         tmpCurlE[i][j] = CurlE[i][j];        
    //     }    
    // }

    // /* Apply PETSc ordering */
    // for (PetscInt i = 0; i<NUM_DIMENSIONS; i++){
    //     for (PetscInt j = 0; j<numDofInCell; j++){
    //         ShapE[i][j] = tmpShapE[i][orderPermutation[j]];
    //         CurlE[i][j] = tmpCurlE[i][orderPermutation[j]];
    //     }
    // }

    /* Free memory */
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < nordB-minK-1; j++){
            PetscCall(PetscFree(ETriV[i][j]));
        }
        PetscCall(PetscFree(ETriV[i]));
    }
    PetscCall(PetscFree(ETriV));   

    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        for (PetscInt j = 0; j < nordB-minK-1; j++){
            PetscCall(PetscFree(CurlETriV[i][j]));
        }
        PetscCall(PetscFree(CurlETriV[i]));
    }
    PetscCall(PetscFree(CurlETriV));
    
    for (PetscInt i = 0; i < maxK; i++){
        PetscCall(PetscFree(homLbet[i]));   
    }
    PetscCall(PetscFree(homLbet));

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < maxK; j++){
            PetscCall(PetscFree(DhomLbet[i][j]));   
        }
        PetscCall(PetscFree(DhomLbet[i]));       
    }
    PetscCall(PetscFree(DhomLbet));           

    PetscFunctionReturn(PETSC_SUCCESS);
}

        
PetscErrorCode computeElementalMatrix(PetscInt nord, PetscInt cellOrientation[10], PetscReal jacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal invJacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscInt numGaussPoints, PetscReal **gaussPoints, PetscReal *weigths, PetscReal *cellResistivity, PetscReal **Me, PetscReal **Ke){
    PetscFunctionBeginUser;

    /* Initial declarations */
    PetscReal e_r[NUM_DIMENSIONS][NUM_DIMENSIONS]  = {{0.0}}; 
    PetscReal mu_r[NUM_DIMENSIONS][NUM_DIMENSIONS] = {{0.0}}; 
    PetscReal iPoint[NUM_DIMENSIONS] = {0.0};
    PetscReal det;
    PetscReal temp1[NUM_DIMENSIONS];
    PetscReal temp2[NUM_DIMENSIONS];
    PetscReal temp3[NUM_DIMENSIONS];
    PetscReal temp4[NUM_DIMENSIONS];
    PetscReal dotResult; 

    /* Tensor for integration (Vertical transverse electric permitivity) */
    e_r[0][0] = cellResistivity[0];
    e_r[1][1] = cellResistivity[1];
    e_r[2][2] = cellResistivity[2];

    /* Tensor for integration (Constant magnetic permittivity) */
    mu_r[0][0] = 1.0;
    mu_r[1][1] = 1.0;
    mu_r[2][2] = 1.0;

    /* Compute the determinant of the Jacobian */ 
    det = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
        - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
        + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);

    /* Compute number of dofs per cell */
    PetscInt numDofInCell; 
    numDofInCell = nord*(nord+2)*(nord+3)/2;    

    /* Allocate matrices for shape functions */    
    PetscReal **ShapE;
    PetscReal **CurlE;
    PetscReal **NiReal;

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &ShapE));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &CurlE));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &NiReal));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(numDofInCell, &ShapE[i]));
        PetscCall(PetscCalloc1(numDofInCell, &CurlE[i]));
        PetscCall(PetscCalloc1(numDofInCell, &NiReal[i]));
        
    }

    /* Reset elemental matrices */
    for (PetscInt i = 0; i < numDofInCell; ++i){
        for (PetscInt j = 0; j < numDofInCell; ++j){
                Me[i][j] = 0.0;
                Ke[i][j] = 0.0;
        }
    }

    /* Compute elemental matrices (mass and stifness matrix)*/    
    for (PetscInt i = 0; i < numGaussPoints; ++i){
        // Get gauss for i point 
        iPoint[0] = gaussPoints[i][0];
        iPoint[1] = gaussPoints[i][1];
        iPoint[2] = gaussPoints[i][2]; 

        /* Compute basis function for i point */
        PetscCall(shape3DETet(iPoint, nord, cellOrientation, ShapE, CurlE));

        // NiReal = Ni in real element        
        for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
            for (PetscInt k = 0; k < numDofInCell; ++k){
                NiReal[j][k] = 0.0;
            }
        }

        for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
            for (PetscInt k = 0; k < numDofInCell; ++k){
                for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                    NiReal[j][k] += (invJacobian[j][m] * ShapE[m][k]);
                }
            }
        }

        /* Perform mass matrix integration */
        for (PetscInt j = 0; j < numDofInCell; ++j){
            for (PetscInt k = 0; k < numDofInCell; ++k){
                /* Prepare data */
                for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                    /* Extract slide for row j */
                    temp1[m] = NiReal[m][j];
                    // Extract slide for row k
                    temp2[m] = NiReal[m][k];
                }

                // Perform matrix vector multiplication
                PetscCall(matrixVectorProduct(temp1, e_r, temp3));
                 
                // Perform dot product
                PetscCall(dotProduct(temp3, temp2, &dotResult));

                // Integration
                Me[j][k] += weigths[i] * dotResult * det;
            }
        }

        // Transform curl on reference element to real element
        for (PetscInt j = 0; j < numDofInCell; ++j){
            // Prepare data
            for (PetscInt k = 0; k < NUM_DIMENSIONS; ++k){
                // Extract slide for row j
                temp1[k] = CurlE[k][j];            
            }

            // Perform vector matrix multiplication
            PetscCall(vectorMatrixProduct(temp1, jacobian, temp2));

            // Update data
            for (PetscInt k = 0; k < NUM_DIMENSIONS; ++k){
                // Update slide for row j
                CurlE[k][j] = temp2[k] / det;
            }
        }

        // Perform stiffness matrix integration
        for (PetscInt j = 0; j < numDofInCell; ++j){
            for (PetscInt k = 0; k < numDofInCell; ++k){
                // Prepare data
                for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                    // Extract slide for row j
                    temp1[m] = CurlE[m][j];
                    temp2[m] = CurlE[m][k];                                        
                }

                // Perform matrix vector multiplication
                PetscCall(matrixVectorProduct(temp1, mu_r, temp3));
                 
                // Perform point-wise multiplication
                for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                    temp4[m] = temp2[m] * det;                 
                }

                // Perform dot product 
                PetscCall(dotProduct(temp3, temp4, &dotResult));

                // Integration
                Ke[j][k] += weigths[i] * dotResult;
            }
        }
    }

    /* Free memory */ 
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(ShapE[i]));
        PetscCall(PetscFree(CurlE[i]));
        PetscCall(PetscFree(NiReal[i]));
    }
    PetscCall(PetscFree(ShapE));
    PetscCall(PetscFree(CurlE));   
    PetscCall(PetscFree(NiReal));   

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode computeBasisFunctions(PetscInt nord, PetscInt orientation[10], PetscReal jacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal invJacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal *point, PetscReal **basisFunctions, PetscReal **curlBasisFunctions){
    PetscFunctionBeginUser;

    /* Initial declarations */
    PetscReal iPoint[NUM_DIMENSIONS]  = {0.0};
    PetscReal temp1[NUM_DIMENSIONS];
    PetscReal temp2[NUM_DIMENSIONS];
    PetscReal det;

    /* Compute number of dofs per cell */
    PetscInt numDofInCell; 
    numDofInCell = nord*(nord+2)*(nord+3)/2;    

    /* Compute the determinant of the Jacobian */ 
    det = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
        - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
        + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);


    /* Allocate matrices for shape functions */    
    PetscReal **ShapE;
    PetscReal **CurlE;
    
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &ShapE));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &CurlE));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(numDofInCell, &ShapE[i]));
        PetscCall(PetscCalloc1(numDofInCell, &CurlE[i]));
    }

    // Compute basis functions for iPoint
    iPoint[0] = point[0];
    iPoint[1] = point[1];
    iPoint[2] = point[2];

    // Compute basis function for i point
    PetscCall(shape3DETet(iPoint, nord, orientation, ShapE, CurlE));    

    /* Reset basis functions and curl functions */
    for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
        for (PetscInt k = 0; k < numDofInCell; ++k){
            basisFunctions[j][k] = 0.0;
            curlBasisFunctions[j][k] = 0.0;
        }
    }

    /* NiReal = Ni in real element */        
    for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
        for (PetscInt k = 0; k < numDofInCell; ++k){
            for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                basisFunctions[j][k] += invJacobian[j][m] * ShapE[m][k];
            }
        }
    }

    /* Transform curl on reference element to real element */
    for (PetscInt i = 0; i < numDofInCell; ++i){
        // Prepare data
        for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
            // Extract slide for row j
            temp1[j] = CurlE[j][i];            
        }

        // Perform vector matrix multiplication
        PetscCall(vectorMatrixProduct(temp1, jacobian, temp2));

        // Update data
        for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
            // Update slide for row j
            curlBasisFunctions[j][i] = temp2[j] / det;
        }
    }

    /* Free memory */ 
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(ShapE[i]));
        PetscCall(PetscFree(CurlE[i]));
    }
    PetscCall(PetscFree(ShapE));
    PetscCall(PetscFree(CurlE));   

    PetscFunctionReturn(PETSC_SUCCESS);
}


root@e8cea040cfce:/workspace/PETGEM_Muhammad/src# cat hvfem.c 
/*
 * Filename: hvfem.c
 * Author: Octavio Castillo Reyes (UPC/BSC)
 * Date: 2024-06-12
 *
 * Description:
 * This file contains functions for high-order vector finite element method (HVFEM) computations. 
  *
 * List of Functions:
 * 
*/

/* C libraries */ 


/* PETSc libraries */
#include <petsc.h>
#include <petscsys.h> 

/* PETGEM functions */ 
#include "constants.h"

// =============================================================================
// Function: computeJacobian
// =============================================================================
PetscErrorCode tetrahedronXYZToXiEtaZeta(PetscScalar *cellCoords, PetscReal point[NUM_DIMENSIONS], PetscReal XiEtaZeta[NUM_DIMENSIONS]){
    /*Compute the reference tetrahedron coordinates from xyz global tetrahedron coordinates.

    :param ndarray cellCoords: spatial coordinates of the nodes with dimensions = (4,3)
    :param ndarray points: xyz points coordinates to be transformed
    :return: xietazeta points coordinates
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    PetscReal J, xi, eta, zeta, cellCoordsReal[NUM_VERTICES_PER_ELEMENT*NUM_DIMENSIONS];

    for (PetscInt i=0; i < NUM_VERTICES_PER_ELEMENT*NUM_DIMENSIONS; i++){
        cellCoordsReal[i] = PetscRealPart(cellCoords[i]);
    }

    J = cellCoordsReal[5] * ( cellCoordsReal[0] * (cellCoordsReal[10] - cellCoordsReal[7])
        + cellCoordsReal[6] * (cellCoordsReal[1] - cellCoordsReal[10])
        + cellCoordsReal[9] * (cellCoordsReal[7] - cellCoordsReal[1]) )
        + cellCoordsReal[2] * ( cellCoordsReal[3] * (cellCoordsReal[7] - cellCoordsReal[10])
        + cellCoordsReal[6] * (cellCoordsReal[10] - cellCoordsReal[4])
        + cellCoordsReal[9] * (cellCoordsReal[4] - cellCoordsReal[7]) )
        + cellCoordsReal[8] * ( cellCoordsReal[3] * (cellCoordsReal[10] - cellCoordsReal[1])
        + cellCoordsReal[0] * (cellCoordsReal[4] - cellCoordsReal[10])
        + cellCoordsReal[9] * (cellCoordsReal[1] - cellCoordsReal[4]) )
        + cellCoordsReal[11] * ( cellCoordsReal[3] * (cellCoordsReal[1] - cellCoordsReal[7])
        + cellCoordsReal[0] * (cellCoordsReal[7] - cellCoordsReal[4])
        + cellCoordsReal[6] * (cellCoordsReal[4] - cellCoordsReal[1]) );

    // Compute affine transformation for xi
    xi = ( cellCoordsReal[11] * (cellCoordsReal[7] - cellCoordsReal[4]) + cellCoordsReal[5] * (cellCoordsReal[10] - cellCoordsReal[7])
         + cellCoordsReal[8] * (cellCoordsReal[4] - cellCoordsReal[10]) ) / J * point[0] +
         ( cellCoordsReal[5] * (cellCoordsReal[6] - cellCoordsReal[9]) + cellCoordsReal[11] * (cellCoordsReal[3] - cellCoordsReal[6])
         + cellCoordsReal[8] * (cellCoordsReal[9] - cellCoordsReal[3]) ) / J * point[1] +
         ( cellCoordsReal[3] * (cellCoordsReal[7] - cellCoordsReal[10]) + cellCoordsReal[9] * (cellCoordsReal[4] - cellCoordsReal[7])
         + cellCoordsReal[6] * (cellCoordsReal[10] - cellCoordsReal[4]) ) / J * point[2] +
         ( cellCoordsReal[8] * (cellCoordsReal[3] * cellCoordsReal[10] - cellCoordsReal[9] * cellCoordsReal[4])
         + cellCoordsReal[5] * (cellCoordsReal[9] * cellCoordsReal[7] - cellCoordsReal[6] * cellCoordsReal[10])
         + cellCoordsReal[11] * (cellCoordsReal[6] * cellCoordsReal[4] - cellCoordsReal[3] * cellCoordsReal[7]) ) / J;
        
    // Compute affine transformation for eta
    eta = ( cellCoordsReal[2] * (cellCoordsReal[10] - cellCoordsReal[4]) + cellCoordsReal[11] * (cellCoordsReal[4] - cellCoordsReal[1])
          + cellCoordsReal[5] * (cellCoordsReal[1] - cellCoordsReal[10]) ) / J * point[0] +
          ( cellCoordsReal[2] * (cellCoordsReal[3] - cellCoordsReal[9]) + cellCoordsReal[5] * (cellCoordsReal[9] - cellCoordsReal[0])
          + cellCoordsReal[11] * (cellCoordsReal[0] - cellCoordsReal[3]) ) / J * point[1] +
          ( cellCoordsReal[0] * (cellCoordsReal[4] - cellCoordsReal[10]) + cellCoordsReal[3] * (cellCoordsReal[10] - cellCoordsReal[1])
          + cellCoordsReal[9] * (cellCoordsReal[1] - cellCoordsReal[4]) ) / J * point[2] +
          ( cellCoordsReal[2] * (cellCoordsReal[9] * cellCoordsReal[4] - cellCoordsReal[3] * cellCoordsReal[10])
          + cellCoordsReal[5] * (cellCoordsReal[0] * cellCoordsReal[10] - cellCoordsReal[9] * cellCoordsReal[1])
          + cellCoordsReal[11] * (cellCoordsReal[3] * cellCoordsReal[1] - cellCoordsReal[0] * cellCoordsReal[4]) ) / J;
        
    // Compute affine transformation for zeta
    zeta = ( cellCoordsReal[5] * (cellCoordsReal[7] - cellCoordsReal[1]) + cellCoordsReal[8] * (cellCoordsReal[1] - cellCoordsReal[4])
           + cellCoordsReal[2] * (cellCoordsReal[4] - cellCoordsReal[7]) ) / J * point[0] +
           ( cellCoordsReal[8] * (cellCoordsReal[3] - cellCoordsReal[0]) + cellCoordsReal[5] * (cellCoordsReal[0] - cellCoordsReal[6])
           + cellCoordsReal[2] * (cellCoordsReal[6] - cellCoordsReal[3]) ) / J * point[1] +
           ( cellCoordsReal[3] * (cellCoordsReal[1] - cellCoordsReal[7]) + cellCoordsReal[0] * (cellCoordsReal[7] - cellCoordsReal[4])
           + cellCoordsReal[6] * (cellCoordsReal[4] - cellCoordsReal[1]) ) / J * point[2] +
           ( cellCoordsReal[5] * ( cellCoordsReal[9] * (cellCoordsReal[1] - cellCoordsReal[7])
           + cellCoordsReal[10] * (cellCoordsReal[6] - cellCoordsReal[0]) )
           + cellCoordsReal[8] * ( cellCoordsReal[9] * (cellCoordsReal[4] - cellCoordsReal[1])
           + cellCoordsReal[10] * (cellCoordsReal[0] - cellCoordsReal[3]) )
           + cellCoordsReal[2] * ( cellCoordsReal[9] * (cellCoordsReal[7] - cellCoordsReal[4])
           + cellCoordsReal[10] * (cellCoordsReal[3] - cellCoordsReal[6]) )
           + cellCoordsReal[11] * ( cellCoordsReal[0] * (cellCoordsReal[4] - cellCoordsReal[7])
           + cellCoordsReal[3] * (cellCoordsReal[7] - cellCoordsReal[1])
           + cellCoordsReal[6] * (cellCoordsReal[1] - cellCoordsReal[4]) ) + J ) / J;
            
    XiEtaZeta[0] = xi;
    XiEtaZeta[1] = eta;
    XiEtaZeta[2] = zeta;
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode vectorRotation(PetscReal azimuth, PetscReal dip, PetscReal rotationVector[NUM_DIMENSIONS]){
    /*Compute the weigths vector for source rotation in the xyz plane.

    :param float azimuth: degrees for x-y plane rotation
    :param float dip: degrees for x-z plane rotation
    :return: weigths for source rotation
    :rtype: ndarray.
    */

    PetscFunctionBeginUser;

    PetscReal base_vector[NUM_DIMENSIONS] = {1., 0., 0.};

    // ---------------------------------------------------------------
    // Compute vector for source rotation
    // ---------------------------------------------------------------
    // Convert degrees to radians for rotation
    PetscReal alpha = azimuth * PETSC_PI / 180.;    // x-y plane
    PetscReal beta  = dip * PETSC_PI / 180.;        // x-z plane
    PetscReal tetha = 0.0 * PETSC_PI / 180.;        // y-z plane

    // Define rotation matrices for each plane
    // x-y plane
    PetscReal M1[NUM_DIMENSIONS][NUM_DIMENSIONS] = {{PetscCosReal(alpha), -PetscSinReal(alpha),   0.},
                                                    {PetscSinReal(alpha),  PetscCosReal(alpha),   0.},
                                                    {                 0.,                   0.,   1.}};

    // x-z plane
    PetscReal M2[NUM_DIMENSIONS][NUM_DIMENSIONS] = {{PetscCosReal(beta),  0.,  -PetscSinReal(beta)},
                                                    {                0.,  1.,                   0.},
                                                    {PetscSinReal(beta),  0.,   PetscCosReal(beta)}};

     // y-z plane
    PetscReal M3[NUM_DIMENSIONS][NUM_DIMENSIONS] = {{1.,   0.,                                     0.},
                                                    {0.,   PetscCosReal(tetha),  -PetscSinReal(tetha)},
                                                    {0.,   PetscSinReal(tetha),   PetscCosReal(tetha)}};
    
    PetscReal temp1[NUM_DIMENSIONS][NUM_DIMENSIONS], temp2[NUM_DIMENSIONS][NUM_DIMENSIONS];

    // Perform matrix multiplications
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
            PetscReal sum = 0.0;
            for (PetscInt k = 0; k < NUM_DIMENSIONS; k++) {
                sum += M1[i][k] * M2[k][j];
            }
            temp1[i][j] = sum;
        }
    }

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
            PetscReal sum = 0.0;
            for (PetscInt k = 0; k < NUM_DIMENSIONS; k++) {
                sum += temp1[i][k] * M3[k][j];
            }
            temp2[i][j] = sum;
        }
    }

    // Apply rotation
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        rotationVector[i] = 0.0;
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
            rotationVector[i] += temp2[i][j] * base_vector[j];
        }
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode crossProduct(PetscReal vector1[NUM_DIMENSIONS], PetscReal vector2[NUM_DIMENSIONS], PetscReal result[NUM_DIMENSIONS]){
    PetscFunctionBeginUser;

    result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1]; // Compute x component
    result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2]; // Compute y component
    result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0]; // Compute z component

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode matrixVectorProduct(PetscReal vector[NUM_DIMENSIONS], PetscReal matrix[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal result[NUM_DIMENSIONS]){
    PetscFunctionBeginUser;

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        result[i] = 0.0;
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++){
            result[i] += matrix[i][j] * vector[j];
        }
    }    

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode vectorMatrixProduct(PetscReal vector[NUM_DIMENSIONS], PetscReal matrix[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal result[NUM_DIMENSIONS]){
    PetscFunctionBeginUser;

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        result[i] = 0.0;
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++){
            result[i] += matrix[j][i] * vector[j];
        }
    }    

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode dotProduct(PetscReal vector1[NUM_DIMENSIONS], PetscReal vector2[NUM_DIMENSIONS], PetscReal *result){
    PetscFunctionBeginUser;

    *result = 0.0;
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        *result += vector1[i] * vector2[i];
    }    

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode computeJacobian(PetscScalar *cellCoords, PetscReal jacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal invJacobian[NUM_DIMENSIONS][NUM_DIMENSIONS]) {
    PetscFunctionBeginUser;

    /* Reset matrices */
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
            jacobian[i][j] = 0.0;
            invJacobian[i][j] = 0.0;
        }
    }

    /* Compute jacobian */ 
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        jacobian[0][i] = PetscRealPart(cellCoords[i])   - PetscRealPart(cellCoords[3+i]); 
        jacobian[1][i] = PetscRealPart(cellCoords[6+i]) - PetscRealPart(cellCoords[3+i]); 
        jacobian[2][i] = PetscRealPart(cellCoords[9+i]) - PetscRealPart(cellCoords[3+i]);
    }

    /* Compute determinant */
    PetscReal determinant; 

    determinant = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1]) -
                  jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0]) +
                  jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);

    /* Compute cofactor matrix */
    PetscReal coFactorMatrix[NUM_DIMENSIONS][NUM_DIMENSIONS];

    coFactorMatrix[0][0] =   jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1];
    coFactorMatrix[0][1] = -(jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0]);
    coFactorMatrix[0][2] =   jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0];
    coFactorMatrix[1][0] = -(jacobian[0][1] * jacobian[2][2] - jacobian[0][2] * jacobian[2][1]);
    coFactorMatrix[1][1] =   jacobian[0][0] * jacobian[2][2] - jacobian[0][2] * jacobian[2][0];
    coFactorMatrix[1][2] = -(jacobian[0][0] * jacobian[2][1] - jacobian[0][1] * jacobian[2][0]);
    coFactorMatrix[2][0] =   jacobian[0][1] * jacobian[1][2] - jacobian[0][2] * jacobian[1][1];
    coFactorMatrix[2][1] = -(jacobian[0][0] * jacobian[1][2] - jacobian[0][2] * jacobian[1][0]);
    coFactorMatrix[2][2] =   jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];

    /* Compute adjugate matrix */ 
    PetscReal adjugateMatrix[NUM_DIMENSIONS][NUM_DIMENSIONS];

    for (PetscInt i = 0; i < 3; i++) {
        for (PetscInt j = 0; j < 3; j++) {
            adjugateMatrix[i][j] = coFactorMatrix[j][i];
        }
    }

    /* Compute inverse of jacobian */
    PetscReal invDeterminant;

    invDeterminant = 1.0 / determinant;
            
    for (PetscInt i = 0; i < 3; i++) {
        for (PetscInt j = 0; j < 3; j++) {
            invJacobian[i][j] = invDeterminant * adjugateMatrix[i][j];
        }
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}
    

PetscErrorCode computeNumGaussPoints3D(PetscInt nord, PetscInt *numGaussPoints){
    PetscFunctionBeginUser;

    PetscInt gaussOrder;
                
    /* Compute gauss order */ 
    gaussOrder = 2*nord;
            
    PetscCheck(gaussOrder >= 1, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Error: Orders lower than 1 are not supported.\n");
    PetscCheck(gaussOrder <= 12, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Error: Orders higher than 6 are not supported by PETGEM.\n");
    
    switch (gaussOrder) {
        case 1:  *numGaussPoints = 1; break;
        case 2:  *numGaussPoints = 4; break;
        case 3:  *numGaussPoints = 5; break;
        case 4:  *numGaussPoints = 11; break;
        case 5:  *numGaussPoints = 14; break;
        case 6:  *numGaussPoints = 24; break;
        case 7:  *numGaussPoints = 31; break;
        case 8:  *numGaussPoints = 43; break;
        case 9:  *numGaussPoints = 53; break;
        case 10: *numGaussPoints = 126; break;
        case 11: *numGaussPoints = 126; break;
        case 12: *numGaussPoints = 210; break;
        default: break;    
    }
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode renormalization3DGaussPoints(PetscInt numPoints, const PetscReal (*gaussPoints)[4], PetscReal** points, PetscReal* weights){
    PetscFunctionBeginUser;

    for(PetscInt i = 0; i < numPoints; i++){
        weights[i] = gaussPoints[i][3]/8;
        points[i][0] = (1 + gaussPoints[i][1])/2;
        points[i][1] = -(1 + gaussPoints[i][0] + gaussPoints[i][1] + gaussPoints[i][2])/2;
        points[i][2] = (1 + gaussPoints[i][0])/2;
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode computeGaussPoints3D(PetscInt numPoints, PetscReal **points, PetscReal *weights){
    PetscFunctionBeginUser;

    const PetscReal nord1_3DGaussPoints[1][4] = {{-0.500000000000000, -0.500000000000000, -0.500000000000000, 1.333333333333333}};

    const PetscReal nord2_3DGaussPoints[4][4] = {{-0.723606797749979, -0.723606797749979, -0.723606797749979, 0.333333333333333},
                                                 { 0.170820393249937, -0.723606797749979, -0.723606797749979, 0.333333333333333},
                                                 {-0.723606797749979,  0.170820393249937, -0.723606797749979, 0.333333333333333},
                                                 {-0.723606797749979, -0.723606797749979,  0.170820393249937, 0.333333333333333}};

    const PetscReal nord3_3DGaussPoints[5][4] = {{-0.500000000000000, -0.500000000000000, -0.500000000000000, -1.066666666666667},
                                                 {-0.666666666666667, -0.666666666666667, -0.666666666666667,  0.600000000000000},
                                                 {-0.666666666666667, -0.666666666666667,  0.000000000000000,  0.600000000000000},
                                                 {-0.666666666666667,  0.000000000000000, -0.666666666666667,  0.600000000000000},
                                                 { 0.000000000000000, -0.666666666666667, -0.666666666666667,  0.600000000000000}};

    const PetscReal nord4_3DGaussPoints[11][4] = {{-0.500000000000000, -0.500000000000000, -0.500000000000000, -0.105244444444444},
                                                  {-0.857142857142857, -0.857142857142857, -0.857142857142857,  0.060977777777778},
                                                  {-0.857142857142857, -0.857142857142857,  0.571428571428571,  0.060977777777778},
                                                  {-0.857142857142857,  0.571428571428571, -0.857142857142857,  0.060977777777778},
                                                  { 0.571428571428571, -0.857142857142857, -0.857142857142857,  0.060977777777778},
                                                  {-0.201192847666402, -0.201192847666402, -0.798807152333598,  0.199111111111111},
                                                  {-0.201192847666402, -0.798807152333598, -0.201192847666402,  0.199111111111111},
                                                  {-0.798807152333598, -0.201192847666402, -0.201192847666402,  0.199111111111111},
                                                  {-0.201192847666402, -0.798807152333598, -0.798807152333598,  0.199111111111111},
                                                  {-0.798807152333598, -0.201192847666402, -0.798807152333598,  0.199111111111111},
                                                  {-0.798807152333598, -0.798807152333598, -0.201192847666402,  0.199111111111111}};

    const PetscReal nord5_3DGaussPoints[14][4] = {{-0.814529499378218, -0.814529499378218, -0.814529499378218,  0.097990724155149},
                                                  { 0.443588498134653, -0.814529499378218, -0.814529499378218,  0.097990724155149},
                                                  {-0.814529499378218,  0.443588498134653, -0.814529499378218,  0.097990724155149},
                                                  {-0.814529499378218, -0.814529499378218,  0.443588498134653,  0.097990724155149},
                                                  {-0.378228161473399, -0.378228161473399, -0.378228161473399,  0.150250567624021},
                                                  {-0.865315515579804, -0.378228161473399, -0.378228161473399,  0.150250567624021},
                                                  {-0.378228161473399, -0.865315515579804, -0.378228161473399,  0.150250567624021},
                                                  {-0.378228161473399, -0.378228161473399, -0.865315515579804,  0.150250567624021},
                                                  {-0.091007408251299, -0.091007408251299, -0.908992591748701,  0.056728027702775},
                                                  {-0.091007408251299, -0.908992591748701, -0.091007408251299,  0.056728027702775},
                                                  {-0.908992591748701, -0.091007408251299, -0.091007408251299,  0.056728027702775},
                                                  {-0.091007408251299, -0.908992591748701, -0.908992591748701,  0.056728027702775},
                                                  {-0.908992591748701, -0.091007408251299, -0.908992591748701,  0.056728027702775},
                                                  {-0.908992591748701, -0.908992591748701, -0.091007408251299,  0.056728027702775}};

    const PetscReal nord6_3DGaussPoints[24][4] = {{-0.570794257481696, -0.570794257481696, -0.570794257481696, 0.053230333677557},
                                                  {-0.287617227554912, -0.570794257481696, -0.570794257481696, 0.053230333677557},
                                                  {-0.570794257481696, -0.287617227554912, -0.570794257481696, 0.053230333677557},
                                                  {-0.570794257481696, -0.570794257481696, -0.287617227554912, 0.053230333677557},
                                                  {-0.918652082930777, -0.918652082930777, -0.918652082930777, 0.013436281407094},
                                                  { 0.755956248792332, -0.918652082930777, -0.918652082930777, 0.013436281407094},
                                                  {-0.918652082930777,  0.755956248792332, -0.918652082930777 ,0.013436281407094},
                                                  {-0.918652082930777, -0.918652082930777,  0.755956248792332 ,0.013436281407094},
                                                  {-0.355324219715449, -0.355324219715449, -0.355324219715449, 0.073809575391540},
                                                  {-0.934027340853653, -0.355324219715449, -0.355324219715449, 0.073809575391540},
                                                  {-0.355324219715449, -0.934027340853653, -0.355324219715449, 0.073809575391540},
                                                  {-0.355324219715449, -0.355324219715449, -0.934027340853653, 0.073809575391540},
                                                  {-0.872677996249965, -0.872677996249965, -0.460655337083368, 0.064285714285714},
                                                  {-0.872677996249965, -0.460655337083368, -0.872677996249965, 0.064285714285714},
                                                  {-0.872677996249965, -0.872677996249965,  0.206011329583298, 0.064285714285714},
                                                  {-0.872677996249965,  0.206011329583298, -0.872677996249965, 0.064285714285714},
                                                  {-0.872677996249965, -0.460655337083368,  0.206011329583298, 0.064285714285714},
                                                  {-0.872677996249965,  0.206011329583298, -0.460655337083368, 0.064285714285714},
                                                  {-0.460655337083368, -0.872677996249965, -0.872677996249965, 0.064285714285714},
                                                  {-0.460655337083368, -0.872677996249965,  0.206011329583298, 0.064285714285714},
                                                  {-0.460655337083368,  0.206011329583298, -0.872677996249965, 0.064285714285714},
                                                  { 0.206011329583298, -0.872677996249965, -0.460655337083368, 0.064285714285714},
                                                  { 0.206011329583298, -0.872677996249965, -0.872677996249965, 0.064285714285714},
                                                  { 0.206011329583298, -0.460655337083368, -0.872677996249965, 0.064285714285714}};

    const PetscReal nord7_3DGaussPoints[31][4] = {{ 0.000000000000000,   0.000000000000000,  -1.000000000000000,   0.007760141093474},
                                                  { 0.000000000000000,  -1.000000000000000,   0.000000000000000,   0.007760141093474},
                                                  {-1.000000000000000,   0.000000000000000,   0.000000000000000,   0.007760141093474},
                                                  {-1.000000000000000,  -1.000000000000000,   0.000000000000000,   0.007760141093474},
                                                  {-1.000000000000000,   0.000000000000000,  -1.000000000000000,   0.007760141093474},
                                                  { 0.000000000000000,  -1.000000000000000,  -1.000000000000000,   0.007760141093474},
                                                  {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,   0.146113787728871},
                                                  {-0.843573615339364,  -0.843573615339364,  -0.843573615339364,   0.084799532195309},
                                                  {-0.843573615339364,  -0.843573615339364,   0.530720846018092,   0.084799532195309},
                                                  {-0.843573615339364,   0.530720846018092,  -0.843573615339364,   0.084799532195309},
                                                  { 0.530720846018092,  -0.843573615339364,  -0.843573615339364,   0.084799532195309},
                                                  {-0.756313566672190,  -0.756313566672190,  -0.756313566672190,  -0.500141920914655},
                                                  {-0.756313566672190,  -0.756313566672190,   0.268940700016569,  -0.500141920914655},
                                                  {-0.756313566672190,   0.268940700016569,  -0.756313566672190,  -0.500141920914655},
                                                  { 0.268940700016569,  -0.756313566672190,  -0.756313566672190,  -0.500141920914655},
                                                  {-0.334921671107159,  -0.334921671107159,  -0.334921671107159,   0.039131402104588},
                                                  {-0.334921671107159,  -0.334921671107159,  -0.995234986678524,   0.039131402104588},
                                                  {-0.334921671107159,  -0.995234986678524,  -0.334921671107159,   0.039131402104588},
                                                  {-0.995234986678524,  -0.334921671107159,  -0.334921671107159,   0.039131402104588},
                                                  {-0.800000000000000,  -0.800000000000000,  -0.600000000000000,   0.220458553791887},
                                                  {-0.800000000000000,  -0.600000000000000,  -0.800000000000000,   0.220458553791887},
                                                  {-0.800000000000000,  -0.800000000000000,   0.200000000000000,   0.220458553791887},
                                                  {-0.800000000000000,   0.200000000000000,  -0.800000000000000,   0.220458553791887},
                                                  {-0.800000000000000,  -0.600000000000000,   0.200000000000000,   0.220458553791887},
                                                  {-0.800000000000000,   0.200000000000000,  -0.600000000000000,   0.220458553791887},
                                                  {-0.600000000000000,  -0.800000000000000,  -0.800000000000000,   0.220458553791887},
                                                  {-0.600000000000000,  -0.800000000000000,   0.200000000000000,   0.220458553791887},
                                                  {-0.600000000000000,   0.200000000000000,  -0.800000000000000,   0.220458553791887},
                                                  { 0.200000000000000,  -0.800000000000000,  -0.600000000000000,   0.220458553791887},
                                                  { 0.200000000000000,  -0.800000000000000,  -0.800000000000000,   0.220458553791887},
                                                  { 0.200000000000000,  -0.600000000000000,  -0.800000000000000,   0.220458553791887}};

    const PetscReal nord8_3DGaussPoints[43][4] = {{-0.500000000000000,  -0.500000000000000,  -0.500000000000000, -0.164001509269119},
                                                  {-0.586340136778654,  -0.586340136778654,  -0.586340136778654,  0.114002446582935},
                                                  {-0.586340136778654,  -0.586340136778654,  -0.240979589664039,  0.114002446582935},
                                                  {-0.586340136778654,  -0.240979589664039,  -0.586340136778654,  0.114002446582935},
                                                  {-0.240979589664039,  -0.586340136778654,  -0.586340136778654,  0.114002446582935},
                                                  {-0.835792823378907,  -0.835792823378907,  -0.835792823378907,  0.015736266505071},
                                                  {-0.835792823378907,  -0.835792823378907,   0.507378470136720,  0.015736266505071},
                                                  {-0.835792823378907,   0.507378470136720,  -0.835792823378907,  0.015736266505071},
                                                  { 0.507378470136720,  -0.835792823378907,  -0.835792823378907,  0.015736266505071},
                                                  {-0.988436098989604,  -0.988436098989604,  -0.988436098989604,  0.001358672872743},
                                                  {-0.988436098989604,  -0.988436098989604,   0.965308296968812,  0.001358672872743},
                                                  {-0.988436098989604,   0.965308296968812,  -0.988436098989604,  0.001358672872743},
                                                  { 0.965308296968812,  -0.988436098989604,  -0.988436098989604,  0.001358672872743},
                                                  {-0.898934519962212,  -0.898934519962212,  -0.101065480037788,  0.036637470595738},
                                                  {-0.898934519962212,  -0.101065480037788,  -0.898934519962212,  0.036637470595738},
                                                  {-0.101065480037788,  -0.898934519962212,  -0.898934519962212,  0.036637470595738},
                                                  {-0.898934519962212,  -0.101065480037788,  -0.101065480037788,  0.036637470595738},
                                                  {-0.101065480037788,  -0.898934519962212,  -0.101065480037788,  0.036637470595738},
                                                  {-0.101065480037788,  -0.101065480037788,  -0.898934519962212,  0.036637470595738},
                                                  {-0.541866927766378,  -0.541866927766378,  -0.928720834422932,  0.045635886469455},
                                                  {-0.541866927766378,  -0.928720834422932,  -0.541866927766378,  0.045635886469455},
                                                  {-0.541866927766378,  -0.541866927766378,   0.012454689955687,  0.045635886469455},
                                                  {-0.541866927766378,   0.012454689955687,  -0.541866927766378,  0.045635886469455},
                                                  {-0.541866927766378,  -0.928720834422932,   0.012454689955687,  0.045635886469455},
                                                  {-0.541866927766378,   0.012454689955687,  -0.928720834422932,  0.045635886469455},
                                                  {-0.928720834422932,  -0.541866927766378,  -0.541866927766378,  0.045635886469455},
                                                  {-0.928720834422932,  -0.541866927766378,   0.012454689955687,  0.045635886469455},
                                                  {-0.928720834422932,   0.012454689955687,  -0.541866927766378,  0.045635886469455},
                                                  { 0.012454689955687,  -0.541866927766378,  -0.928720834422932,  0.045635886469455},
                                                  { 0.012454689955687,  -0.541866927766378,  -0.541866927766378,  0.045635886469455},
                                                  { 0.012454689955687,  -0.928720834422932,  -0.541866927766378,  0.045635886469455},
                                                  {-0.926784500893605,  -0.926784500893605,  -0.619027916130733,  0.017124153129297},
                                                  {-0.926784500893605,  -0.619027916130733,  -0.926784500893605,  0.017124153129297},
                                                  {-0.926784500893605,  -0.926784500893605,   0.472596917917943,  0.017124153129297},
                                                  {-0.926784500893605,   0.472596917917943,  -0.926784500893605,  0.017124153129297},
                                                  {-0.926784500893605,  -0.619027916130733,   0.472596917917943,  0.017124153129297},
                                                  {-0.926784500893605,   0.472596917917943,  -0.619027916130733,  0.017124153129297},
                                                  {-0.619027916130733,  -0.926784500893605,  -0.926784500893605,  0.017124153129297},
                                                  {-0.619027916130733,  -0.926784500893605,   0.472596917917943,  0.017124153129297},
                                                  {-0.619027916130733,   0.472596917917943,  -0.926784500893605,  0.017124153129297},
                                                  { 0.472596917917943,  -0.926784500893605,  -0.619027916130733,  0.017124153129297},
                                                  { 0.472596917917943,  -0.926784500893605,  -0.926784500893605,  0.017124153129297},
                                                  { 0.472596917917943,  -0.619027916130733,  -0.926784500893605,  0.017124153129297}};

    const PetscReal nord9_3DGaussPoints[53][4] = {{-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  -1.102392306608869},
                                                  {-0.903297922900526,  -0.903297922900526,  -0.903297922900526,   0.014922692552682},
                                                  {-0.903297922900526,  -0.903297922900526,   0.709893768701580,   0.014922692552682},
                                                  {-0.903297922900526,   0.709893768701580,  -0.903297922900526,   0.014922692552682},
                                                  { 0.709893768701580,  -0.903297922900526,  -0.903297922900526,   0.014922692552682},
                                                  {-0.350841439764235,  -0.350841439764235,  -0.350841439764235,   0.034475391755947},
                                                  {-0.350841439764235,  -0.350841439764235,  -0.947475680707294,   0.034475391755947},
                                                  {-0.350841439764235,  -0.947475680707294,  -0.350841439764235,   0.034475391755947},
                                                  {-0.947475680707294,  -0.350841439764235,  -0.350841439764235,   0.034475391755947},
                                                  {-0.770766919552010,  -0.770766919552010,  -0.770766919552010,  -0.721478131849612},
                                                  {-0.770766919552010,  -0.770766919552010,   0.312300758656029,  -0.721478131849612},
                                                  {-0.770766919552010,   0.312300758656029,  -0.770766919552010,  -0.721478131849612},
                                                  { 0.312300758656029,  -0.770766919552010,  -0.770766919552010,  -0.721478131849612},
                                                  {-0.549020096176972,  -0.549020096176972,  -0.549020096176972,   0.357380609620092},
                                                  {-0.549020096176972,  -0.549020096176972,  -0.352939711469084,   0.357380609620092},
                                                  {-0.549020096176972,  -0.352939711469084,  -0.549020096176972,   0.357380609620092},
                                                  {-0.352939711469084,  -0.549020096176972,  -0.549020096176972,   0.357380609620092},
                                                  {-0.736744381506260,  -0.736744381506260,  -0.832670596765630,   0.277603247076406},
                                                  {-0.736744381506260,  -0.832670596765630,  -0.736744381506260,   0.277603247076406},
                                                  {-0.736744381506260,  -0.736744381506260,   0.306159359778151,   0.277603247076406},
                                                  {-0.736744381506260,   0.306159359778151,  -0.736744381506260,   0.277603247076406},
                                                  {-0.736744381506260,  -0.832670596765630,   0.306159359778151,   0.277603247076406},
                                                  {-0.736744381506260,   0.306159359778151,  -0.832670596765630,   0.277603247076406},
                                                  {-0.832670596765630,  -0.736744381506260,  -0.736744381506260,   0.277603247076406},
                                                  {-0.832670596765630,  -0.736744381506260,   0.306159359778151,   0.277603247076406},
                                                  {-0.832670596765630,   0.306159359778151,  -0.736744381506260,   0.277603247076406},
                                                  { 0.306159359778151,  -0.736744381506260,  -0.832670596765630,   0.277603247076406},
                                                  { 0.306159359778151,  -0.736744381506260,  -0.736744381506260,   0.277603247076406},
                                                  { 0.306159359778151,  -0.832670596765630,  -0.736744381506260,   0.277603247076406},
                                                  {-0.132097077177186,  -0.132097077177186,  -0.784460280901143,   0.026820671221285},
                                                  {-0.132097077177186,  -0.784460280901143,  -0.132097077177186,   0.026820671221285},
                                                  {-0.132097077177186,  -0.132097077177186,  -0.951345564744484,   0.026820671221285},
                                                  {-0.132097077177186,  -0.951345564744484,  -0.132097077177186,   0.026820671221285},
                                                  {-0.132097077177186,  -0.784460280901143,  -0.951345564744484,   0.026820671221285},
                                                  {-0.132097077177186,  -0.951345564744484,  -0.784460280901143,   0.026820671221285},
                                                  {-0.784460280901143,  -0.132097077177186,  -0.132097077177186,   0.026820671221285},
                                                  {-0.784460280901143,  -0.132097077177186,  -0.951345564744484,   0.026820671221285},
                                                  {-0.784460280901143,  -0.951345564744484,  -0.132097077177186,   0.026820671221285},
                                                  {-0.951345564744484,  -0.132097077177186,  -0.784460280901143,   0.026820671221285},
                                                  {-0.951345564744484,  -0.132097077177186,  -0.132097077177186,   0.026820671221285},
                                                  {-0.951345564744484,  -0.784460280901143,  -0.132097077177186,   0.026820671221285},
                                                  {-1.002752554636276,  -1.002752554636276,  -0.446893054726385,   0.003453031004456},
                                                  {-1.002752554636276,  -0.446893054726385,  -1.002752554636276,   0.003453031004456},
                                                  {-1.002752554636276,  -1.002752554636276,   0.452398163998938,   0.003453031004456},
                                                  {-1.002752554636276,   0.452398163998938,  -1.002752554636276,   0.003453031004456},
                                                  {-1.002752554636276,  -0.446893054726385,   0.452398163998938,   0.003453031004456},
                                                  {-1.002752554636276,   0.452398163998938,  -0.446893054726385,   0.003453031004456},
                                                  {-0.446893054726385,  -1.002752554636276,  -1.002752554636276,   0.003453031004456},
                                                  {-0.446893054726385,  -1.002752554636276,   0.452398163998938,   0.003453031004456},
                                                  {-0.446893054726385,   0.452398163998938,  -1.002752554636276,   0.003453031004456},
                                                  { 0.452398163998938,  -1.002752554636276,  -0.446893054726385,   0.003453031004456},
                                                  { 0.452398163998938,  -1.002752554636276,  -1.002752554636276,   0.003453031004456},
                                                  { 0.452398163998938,  -0.446893054726385,  -1.002752554636276,   0.003453031004456}};

    const PetscReal nord10_3DGaussPoints[126][4] = {{-0.857142857142857,  -0.857142857142857,  0.571428571428571 ,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,   0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,  -0.285714285714286,   0.000000000000000,   0.362902592520648},
                                                    {-0.857142857142857,   0.000000000000000,  -0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,   0.285714285714286,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,   0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,   0.285714285714286,   0.362902592520648},
                                                    {-0.571428571428571,  -0.571428571428571,   0.000000000000000,   0.362902592520648},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.285714285714286,   0.362902592520648},
                                                    {-0.571428571428571,   0.000000000000000,  -0.571428571428571,   0.362902592520648},
                                                    {-0.571428571428571,   0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    {-0.285714285714286,  -0.857142857142857,   0.000000000000000,   0.362902592520648},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.285714285714286,   0.362902592520648},
                                                    {-0.285714285714286,  -0.285714285714286,  -0.571428571428571,   0.362902592520648},
                                                    {-0.285714285714286,   0.000000000000000,  -0.857142857142857,   0.362902592520648},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.285714285714286,   0.362902592520648},
                                                    { 0.000000000000000,  -0.571428571428571,  -0.571428571428571,   0.362902592520648},
                                                    { 0.000000000000000,  -0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    { 0.285714285714286,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    { 0.285714285714286,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    { 0.571428571428571,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,   0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,   0.000000000000000,   0.362902592520648},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,   0.000000000000000,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,   0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,   0.000000000000000,   0.362902592520648},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.285714285714286,   0.362902592520648},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.571428571428571,   0.362902592520648},
                                                    {-0.571428571428571,   0.000000000000000,  -0.857142857142857,   0.362902592520648},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.285714285714286,   0.362902592520648},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.571428571428571,   0.362902592520648},
                                                    {-0.285714285714286,  -0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    { 0.000000000000000,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    { 0.285714285714286,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,   0.000000000000000,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,   0.000000000000000,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.285714285714286,   0.362902592520648},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.571428571428571,   0.362902592520648},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.285714285714286,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.571428571428571,   0.362902592520648},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.857142857142857,   0.362902592520648},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.857142857142857,   0.362902592520648},
                                                    {-0.833333333333333,  -0.833333333333333,   0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.500000000000000,   0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,   0.166666666666667,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,   0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.833333333333333,   0.166666666666667,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.166666666666667,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.500000000000000,   0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.500000000000000,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    { 0.166666666666667,  -0.833333333333333,  -0.500000000000000,  -0.932187812187812},
                                                    { 0.166666666666667,  -0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    { 0.500000000000000,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.833333333333333,   0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,   0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    { 0.166666666666667,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.166666666666667,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.500000000000000,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.833333333333333,  -0.932187812187812},
                                                    {-0.800000000000000,  -0.800000000000000,   0.400000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.400000000000000,   0.000000000000000,   0.815498319838598},
                                                    {-0.800000000000000,   0.000000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.800000000000000,   0.400000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.800000000000000,   0.000000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.400000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.400000000000000,   0.000000000000000,  -0.800000000000000,   0.815498319838598},
                                                    { 0.000000000000000,  -0.800000000000000,  -0.400000000000000,   0.815498319838598},
                                                    { 0.000000000000000,  -0.400000000000000,  -0.800000000000000,   0.815498319838598},
                                                    { 0.400000000000000,  -0.800000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.800000000000000,   0.000000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.400000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.800000000000000,   0.000000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.800000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.400000000000000,  -0.800000000000000,   0.815498319838598},
                                                    { 0.000000000000000,  -0.800000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.800000000000000,  -0.400000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.400000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.400000000000000,  -0.800000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.800000000000000,  -0.800000000000000,  -0.800000000000000,   0.815498319838598},
                                                    {-0.750000000000000,  -0.750000000000000,   0.250000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,  -0.250000000000000,  -0.250000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,   0.250000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.250000000000000,  -0.750000000000000,  -0.250000000000000,  -0.280203089091978},
                                                    {-0.250000000000000,  -0.250000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    { 0.250000000000000,  -0.750000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,  -0.750000000000000,  -0.250000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,  -0.250000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.250000000000000,  -0.750000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.750000000000000,  -0.750000000000000,  -0.750000000000000,  -0.280203089091978},
                                                    {-0.666666666666667,  -0.666666666666667,   0.000000000000000,   0.032544642857143},
                                                    {-0.666666666666667,   0.000000000000000,  -0.666666666666667,   0.032544642857143},
                                                    { 0.000000000000000,  -0.666666666666667,  -0.666666666666667,   0.032544642857143},
                                                    {-0.666666666666667,  -0.666666666666667,  -0.666666666666667,   0.032544642857143},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  -0.000752498530276}};

    const PetscReal nord12_3DGaussPoints[210][4] = {{-0.875000000000000,  -0.875000000000000,   0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,   0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.125000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.125000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,   0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.375000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.125000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,   0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,   0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.625000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.375000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.375000000000000,   0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.875000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.625000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.375000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    { 0.375000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.625000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,   0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.125000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.375000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,   0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.625000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.375000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,   0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.125000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,   0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    { 0.125000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,  -0.125000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.125000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.125000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,  -0.375000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.375000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.375000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,  -0.625000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.625000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.625000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.875000000000000,  -0.875000000000000,  -0.875000000000000,  0.420407272132140},
                                                    {-0.857142857142857,  -0.857142857142857,   0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,   0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,  -0.285714285714286,   0.000000000000000, -1.185481802234117},
                                                    {-0.857142857142857,   0.000000000000000,  -0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,   0.285714285714286,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,   0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,   0.285714285714286, -1.185481802234117},
                                                    {-0.571428571428571,  -0.571428571428571,   0.000000000000000, -1.185481802234117},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.285714285714286, -1.185481802234117},
                                                    {-0.571428571428571,   0.000000000000000,  -0.571428571428571, -1.185481802234117},
                                                    {-0.571428571428571,   0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    {-0.285714285714286,  -0.857142857142857,   0.000000000000000, -1.185481802234117},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.285714285714286, -1.185481802234117},
                                                    {-0.285714285714286,  -0.285714285714286,  -0.571428571428571, -1.185481802234117},
                                                    {-0.285714285714286,   0.000000000000000,  -0.857142857142857, -1.185481802234117},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.285714285714286, -1.185481802234117},
                                                    { 0.000000000000000,  -0.571428571428571,  -0.571428571428571, -1.185481802234117},
                                                    { 0.000000000000000,  -0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    { 0.285714285714286,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    { 0.285714285714286,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    { 0.571428571428571,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,   0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,   0.000000000000000, -1.185481802234117},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,   0.000000000000000,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,   0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,   0.000000000000000, -1.185481802234117},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.285714285714286, -1.185481802234117},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.571428571428571, -1.185481802234117},
                                                    {-0.571428571428571,   0.000000000000000,  -0.857142857142857, -1.185481802234117},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.285714285714286, -1.185481802234117},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.571428571428571, -1.185481802234117},
                                                    {-0.285714285714286,  -0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    { 0.000000000000000,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    { 0.285714285714286,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,   0.000000000000000, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,   0.000000000000000,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.285714285714286, -1.185481802234117},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.571428571428571, -1.185481802234117},
                                                    {-0.571428571428571,  -0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    {-0.285714285714286,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    { 0.000000000000000,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.285714285714286, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,  -0.285714285714286,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    {-0.571428571428571,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    {-0.285714285714286,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.571428571428571, -1.185481802234117},
                                                    {-0.857142857142857,  -0.571428571428571,  -0.857142857142857, -1.185481802234117},
                                                    {-0.571428571428571,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.857142857142857,  -0.857142857142857,  -0.857142857142857, -1.185481802234117},
                                                    {-0.833333333333333,  -0.833333333333333,   0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,  -0.500000000000000,   0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,   0.166666666666667,  -0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,   0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    {-0.500000000000000,  -0.833333333333333,   0.166666666666667,  1.198527187098616},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.166666666666667,  1.198527187098616},
                                                    {-0.500000000000000,  -0.166666666666667,  -0.500000000000000,  1.198527187098616},
                                                    {-0.500000000000000,   0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.166666666666667,  1.198527187098616},
                                                    {-0.166666666666667,  -0.500000000000000,  -0.500000000000000,  1.198527187098616},
                                                    {-0.166666666666667,  -0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    { 0.166666666666667,  -0.833333333333333,  -0.500000000000000,  1.198527187098616},
                                                    { 0.166666666666667,  -0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    { 0.500000000000000,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.833333333333333,  -0.833333333333333,   0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,   0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.166666666666667,  1.198527187098616},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  1.198527187098616},
                                                    {-0.500000000000000,  -0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.500000000000000,  1.198527187098616},
                                                    {-0.166666666666667,  -0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    { 0.166666666666667,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.166666666666667,  1.198527187098616},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,  -0.166666666666667,  -0.833333333333333,  1.198527187098616},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.500000000000000,  1.198527187098616},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    {-0.166666666666667,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.500000000000000,  1.198527187098616},
                                                    {-0.833333333333333,  -0.500000000000000,  -0.833333333333333,  1.198527187098616},
                                                    {-0.500000000000000,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.833333333333333,  -0.833333333333333,  -0.833333333333333,  1.198527187098616},
                                                    {-0.800000000000000,  -0.800000000000000,   0.400000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.400000000000000,   0.000000000000000, -0.522755333229870},
                                                    {-0.800000000000000,   0.000000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.800000000000000,   0.400000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.800000000000000,   0.000000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.400000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.400000000000000,   0.000000000000000,  -0.800000000000000, -0.522755333229870},
                                                    { 0.000000000000000,  -0.800000000000000,  -0.400000000000000, -0.522755333229870},
                                                    { 0.000000000000000,  -0.400000000000000,  -0.800000000000000, -0.522755333229870},
                                                    { 0.400000000000000,  -0.800000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.800000000000000,   0.000000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.400000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.800000000000000,   0.000000000000000 , -0.800000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.800000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.400000000000000,  -0.800000000000000, -0.522755333229870},
                                                    { 0.000000000000000,  -0.800000000000000 , -0.800000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.800000000000000,  -0.400000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.400000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.400000000000000,  -0.800000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.800000000000000,  -0.800000000000000,  -0.800000000000000, -0.522755333229870},
                                                    {-0.750000000000000,  -0.750000000000000,   0.250000000000000,  0.093401029697326},
                                                    {-0.750000000000000,  -0.250000000000000,  -0.250000000000000,  0.093401029697326},
                                                    {-0.750000000000000,   0.250000000000000,  -0.750000000000000,  0.093401029697326},
                                                    {-0.250000000000000,  -0.750000000000000,  -0.250000000000000,  0.093401029697326},
                                                    {-0.250000000000000,  -0.250000000000000,  -0.750000000000000,  0.093401029697326},
                                                    { 0.250000000000000,  -0.750000000000000 , -0.750000000000000,  0.093401029697326},
                                                    {-0.750000000000000,  -0.750000000000000,  -0.250000000000000,  0.093401029697326},
                                                    {-0.750000000000000,  -0.250000000000000,  -0.750000000000000,  0.093401029697326},
                                                    {-0.250000000000000,  -0.750000000000000,  -0.750000000000000,  0.093401029697326},
                                                    {-0.750000000000000,  -0.750000000000000,  -0.750000000000000,  0.093401029697326},
                                                    {-0.666666666666667,  -0.666666666666667,   0.000000000000000, -0.005325487012987},
                                                    {-0.666666666666667,   0.000000000000000,  -0.666666666666667, -0.005325487012987},
                                                    { 0.000000000000000,  -0.666666666666667,  -0.666666666666667, -0.005325487012987},
                                                    {-0.666666666666667,  -0.666666666666667,  -0.666666666666667, -0.005325487012987},
                                                    {-0.500000000000000,  -0.500000000000000,  -0.500000000000000,  0.000050166568685}};


    switch (numPoints)
    {
        case 1:
            PetscCall(renormalization3DGaussPoints(numPoints, nord1_3DGaussPoints, points, weights));
            break;
        case 4:
            PetscCall(renormalization3DGaussPoints(numPoints, nord2_3DGaussPoints, points, weights));
            break;
        case 5:
            PetscCall(renormalization3DGaussPoints(numPoints, nord3_3DGaussPoints, points, weights));
            break;
        case 11:
            PetscCall(renormalization3DGaussPoints(numPoints, nord4_3DGaussPoints, points, weights));
            break;
        case 14:
            PetscCall(renormalization3DGaussPoints(numPoints, nord5_3DGaussPoints, points, weights));
            break;
        case 24:
            PetscCall(renormalization3DGaussPoints(numPoints, nord6_3DGaussPoints, points, weights));
            break;
        case 31:
            PetscCall(renormalization3DGaussPoints(numPoints, nord7_3DGaussPoints, points, weights));
            break;
        case 43:
            PetscCall(renormalization3DGaussPoints(numPoints, nord8_3DGaussPoints, points, weights));
            break;
        case 53:
            PetscCall(renormalization3DGaussPoints(numPoints, nord9_3DGaussPoints, points, weights));
            break;
        case 126:
            PetscCall(renormalization3DGaussPoints(numPoints, nord10_3DGaussPoints, points, weights));
            break;
        case 210:
            PetscCall(renormalization3DGaussPoints(numPoints, nord12_3DGaussPoints, points, weights));
            break;
        default:
            break;
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode AffineTetrahedron(PetscReal X[NUM_DIMENSIONS], PetscReal Lam[4], PetscReal DLam[NUM_DIMENSIONS][4]){
    /*Compute affine coordinates and their gradients.

    :param ndarray X: point coordinates
    :return: affine coordinates and gradients of affine coordinates
    :rtype: ndarray

    .. note:: References:\n
       Fuentes, F., Keith, B., Demkowicz, L., & Nagaraj, S. (2015). Orientation
       embedded high order shape functions for the exact sequence elements of
       all shapes. Computers & Mathematics with applications, 70(4), 353-458.
    */
    PetscFunctionBeginUser;
    
    // Define affine coordinates
    Lam[0] = 1.-X[0]-X[1]-X[2];
    Lam[1] = X[0];
    Lam[2] = X[1];
    Lam[3] = X[2];
    
    // and their gradients
    DLam[0][0] = -1;
    DLam[0][1] =  1;
    DLam[1][0] = -1;
    DLam[1][2] =  1;
    DLam[2][0] = -1;
    DLam[2][3] =  1;

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ProjectTetE(PetscReal Lam[4], PetscReal DLam[NUM_DIMENSIONS][4], PetscReal LampE[NUM_EDGES_PER_ELEMENT][2], PetscReal DLampE[NUM_EDGES_PER_ELEMENT][NUM_DIMENSIONS][2], PetscBool* IdecE){
    /*Projection of tetrahedral edges in concordance with numbering of topological entities (vertices, edges, faces).

    :param ndarray Lam: affine coordinates
    :param ndarray DLam: gradients of affine coordinates
    :return: projection of affine coordinates on edges, projection of gradients of affine coordinates on edges
    :rtype: ndarray

    .. note:: References:\n
       Fuentes, F., Keith, B., Demkowicz, L., & Nagaraj, S. (2015). Orientation
       embedded high order shape functions for the exact sequence elements of
       all shapes. Computers & Mathematics with applications, 70(4), 353-458.
    */
    PetscFunctionBeginUser;
    
    // ---------------------------------------------------------------
    // Compute projection
    // ---------------------------------------------------------------
    // e=1 --> edge01 with local orientation v0->v1
    // PetscInt e = 0;
    // LampE[e][0] = Lam[0];
    // LampE[e][1] = Lam[1];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][0];
    //     DLampE[e][i][1] = DLam[i][1];
    // }
    
    // // e=2 --> edge12 with local orientation v1->v2
    // e = 1;
    // LampE[e][0] = Lam[1];
    // LampE[e][1] = Lam[2];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][1];
    //     DLampE[e][i][1] = DLam[i][2];
    // }
    
    // // e=3 --> edge20 with local orientation v0->v2
    // e = 2;
    // LampE[e][0] = Lam[0];
    // LampE[e][1] = Lam[2];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][0];
    //     DLampE[e][i][1] = DLam[i][2];
    // }
    
    // // e=4 --> edge03 with local orientation v0->v3
    // e = 3;
    // LampE[e][0] = Lam[0];
    // LampE[e][1] = Lam[3];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][0];
    //     DLampE[e][i][1] = DLam[i][3];
    // }
    
    // // e=5 --> edge13 with local orientation v1->v3
    // e = 4;
    // LampE[e][0] = Lam[1];
    // LampE[e][1] = Lam[3];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][1];
    //     DLampE[e][i][1] = DLam[i][3];
    // }
    
    // // e=6 --> edge23 with local orientation v2->v3
    // e = 5;
    // LampE[e][0] = Lam[2];
    // LampE[e][1] = Lam[3];
    // for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
    //     DLampE[e][i][0] = DLam[i][2];
    //     DLampE[e][i][1] = DLam[i][3];
    // }

    PetscInt e = 0;
    LampE[e][0] = Lam[1];
    LampE[e][1] = Lam[0];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][1];
        DLampE[e][i][1] = DLam[i][0];
    }
    
    // e=2 --> edge12 with local orientation v1->v2
    e = 1;
    LampE[e][0] = Lam[1];
    LampE[e][1] = Lam[2];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][1];
        DLampE[e][i][1] = DLam[i][2];
    }
    
    // e=3 --> edge20 with local orientation v0->v2
    e = 2;
    LampE[e][0] = Lam[2];
    LampE[e][1] = Lam[0];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][2];
        DLampE[e][i][1] = DLam[i][0];
    }
    
    // e=4 --> edge03 with local orientation v0->v3
    e = 3;
    LampE[e][0] = Lam[0];
    LampE[e][1] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][0];
        DLampE[e][i][1] = DLam[i][3];
    }
    
    // e=5 --> edge13 with local orientation v1->v3
    e = 4;
    LampE[e][0] = Lam[3];
    LampE[e][1] = Lam[1];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][3];
        DLampE[e][i][1] = DLam[i][1];
    }
    
    // e=6 --> edge23 with local orientation v2->v3
    e = 5;
    LampE[e][0] = Lam[2];
    LampE[e][1] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampE[e][i][0] = DLam[i][2];
        DLampE[e][i][1] = DLam[i][3];
    }

    /* Projected coordinates are Lam, so IdecE=false for all edges */
    *IdecE = PETSC_FALSE;

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode OrientE(PetscReal S[2], PetscReal DS[NUM_DIMENSIONS][2], PetscInt Nori, PetscReal GS[2], PetscReal GDS[NUM_DIMENSIONS][2]){
    /*Compute the local to global transformations of edges.

    :param ndarray S: projection of affine coordinates on edges
    :param ndarray DS: projection of gradients of affine coordinates on edges
    :param ndarray Nori: edge orientation
    :param int N: number of dimensions
    :return: global transformation of edges and global transformation of gradients of edges
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    // ---------------------------------------------------------------
    // Initialization
    // ---------------------------------------------------------------
    // Allocate
    PetscInt Or[2][2];
    
    // Nori=0 => (s0,s1)->(s0,s1)
    Or[0][0] = 0;
    Or[0][1] = 1;
    // Nori=1 => (s0,s1)->(s1,s0)
    Or[1][0] = 1;
    Or[1][1] = 0;

    // Local-to-global transformation
    GS[0] = S[Or[Nori][0]];
    GS[1] = S[Or[Nori][1]];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        GDS[i][0] = DS[i][Or[Nori][0]];
        GDS[i][1] = DS[i][Or[Nori][1]];
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode PolyLegendre(PetscReal X, PetscReal T, PetscInt Nord, PetscReal P[]){
    /*Compute values of shifted scaled Legendre polynomials.

    :param ndarray X: coordinate from [0,1]
    :param float T: scaling parameter
    :param int Nord: polynomial order
    :return: polynomial values
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    // i stands for the order of the polynomial, stored in P(i)
    // lowest order case (order 0)
    P[0] = 1.;
    
    PetscReal y;
    // First order case (order 1) if necessary
    if(Nord >= 1){
        y = 2.*X - T;
        P[1] = y;
    }
  
    if(Nord >= 2){
        PetscReal tt = pow(T,2);
        for(PetscInt i = 1; i < Nord - 1; i++){
            P[i+1] = (2.*i+1)*y*P[i] - (i)*tt*P[i-1];
            P[i+1] = P[i+1]/(PetscReal)(i+1);
        }
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode PolyJacobi(PetscReal X, PetscReal T, PetscInt Nord, PetscInt Minalpha, PetscReal **P){
    /*Compute values of shifted scaled Jacobi polynomials P**alpha-i.

    Result is a half of a matrix with each row associated to a fixed alpha.
    Alpha grows by 2 in each row.

    :param ndarray X: coordinate from [0,1]
    :param float T: scaling parameter
    :param int Nord: max polynomial order
    :param int Minalpha: first row value of alpha (integer)
    :return: polynomial values
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    PetscReal *alpha;

    // Allocate
    PetscCall(PetscCalloc1(Nord + 1, &alpha));

    // Clearly (minI,maxI)=(0,Nord), but the syntax is written as it is
    // because it reflects how the indexing is called from outside
    PetscInt minI = 0;
    PetscInt maxI = minI + Nord;
    
    for(PetscInt i = 0; i < maxI + 1; i++){
        alpha[i] = Minalpha + 2*(i - minI);
    }
    
    // Initiate first column (order 0)
    for(PetscInt i = minI; i < maxI + 1; i++){
        P[i][0] = 1.;
    }
    
    PetscReal y;
    // Initiate second column (order 1) if necessary
    if(Nord >= 1){
        y = 2*X - T;
        for(PetscInt i = minI; i < maxI; i++){
            P[i][1] = y + alpha[i]*X;
        }
    }
    
    // Fill the last columns if necessary
    if(Nord >= 2){
        PetscReal tt = pow(T, 2);
        PetscInt ni = -1;
        for(PetscInt i = 0; i < maxI - 1; i++){
            PetscReal al = alpha[i];
            PetscReal aa = pow(al, 2);
            ni += 1;
            // Use recursion in order, i, to compute P^alpha_i for i>=2
            for(PetscInt j = 2; j < Nord - ni + 1; j++){
                PetscReal ai = 2*j*(j+al)*(2*j+al-2);
                PetscReal bi = 2*j+al-1;
                PetscReal ci = (2*j+al)*(2*j+al-2);
                PetscReal di = 2*(j+al-1)*(j-1)*(2*j+al);
                
                P[i][j] = bi*(ci*y+aa*T)*P[i][j-1]-di*tt*P[i][j-2];
                P[i][j] = P[i][j]/ai;
            }
        }
    }

    // Free memory
    PetscCall(PetscFree(alpha));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode HomLegendre(PetscReal S[2], PetscInt Nord, PetscReal HomP[]){
    /*Compute values of homogenized Legendre polynomials.

    :param ndarray S: affine(like) coordinates
    :param int Nord: polynomial order
    :return: polynomial values
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    PetscCall(PolyLegendre(S[1], S[0] + S[1], Nord, HomP));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode AncEE(PetscReal S[2], PetscReal DS[NUM_DIMENSIONS][2], PetscInt Nord, PetscBool Idec, PetscReal **EE, PetscReal **CurlEE){
    /*Compute edge Hcurl ancillary functions and their curls.

    :param ndarray S: affine coordinates associated to edge
    :param ndarray DS: derivatives of S in R^N
    :param int Nord: polynomial order
    :param bool Idec: Binary flag
    :param int N: spatial dimension
    :return: edge Hcurl ancillary functions, curls of edge Hcurl ancillary functions
    :rtype: ndarray

    .. note:: References:\n
       Idec: = FALSE  s0+s1 != 1
             = TRUE   s0+s1  = 1
    */
    PetscFunctionBeginUser;

    // Local parameters
    PetscInt minI = 1;
    PetscInt maxI = Nord;
    PetscInt Ncurl = 2*NUM_DIMENSIONS-3;
    
    PetscReal *homP;
    
    if (Nord <= 1){
        PetscCall(PetscCalloc1(Nord+1, &homP));
    } else{
        PetscCall(PetscCalloc1(Nord, &homP));
    } 
    
    // Extract homogenized Legendre polyomials first
    PetscCall(HomLegendre(S, maxI, homP));
    
    // Simplified case
    if(Idec){
        for(PetscInt i = minI; i < maxI + 1; i++){
            for(PetscInt j = 0; j < NUM_DIMENSIONS; j++){
                EE[j][i-1] = homP[i-1]*DS[j][1];
            }
        }
        // No need to compute Whitney function or curl
        for (PetscInt i = 0; i < Ncurl; ++i) {
            for (PetscInt j = minI - 1; j < maxI - 1; ++j) {
                CurlEE[i][j] = 0;
            }
        }
    } else {
        // Lowest order Whitney function and its curl
        PetscReal whiE[NUM_DIMENSIONS];
        for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
            whiE[i] = S[0]*DS[i][1] - S[1]*DS[i][0];
        }
        PetscReal curlwhiE[NUM_DIMENSIONS];
        PetscReal temp1[NUM_DIMENSIONS], temp2[NUM_DIMENSIONS];
        for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
            temp1[i] = DS[i][0];
            temp2[i] = DS[i][1];
        }
        PetscCall(crossProduct(temp1, temp2, curlwhiE));
        
        // Now construct the higher order elements
        for(PetscInt i = minI; i < maxI + 1; i++){
            for(PetscInt j = 0; j < NUM_DIMENSIONS; j++){
                EE[j][i-1] = homP[i-1]*whiE[j];
            }
            for(PetscInt j = 0; j < Ncurl; j++){
                CurlEE[j][i-1] = (i+1)*homP[i-1]*curlwhiE[j];
            }
        }
    }

    PetscCall(PetscFree(homP));   
    
    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode ProjectTetF(PetscReal Lam[4], PetscReal DLam[NUM_DIMENSIONS][4], PetscReal LampF[NUM_FACES_PER_ELEMENT][NUM_DIMENSIONS], PetscReal DLampF[NUM_FACES_PER_ELEMENT][NUM_DIMENSIONS][NUM_DIMENSIONS], PetscBool* IdecF){
    /*Projection of tetrahedral faces in concordance with numbering of topological entities (vertices, edges, faces).

    :param ndarray Lam: affine coordinates
    :param ndarray DLam: gradients of affine coordinates
    :return: projection of affine coordinates on faces, projection of gradients of affine coordinates on faces
    :rtype: ndarray
    */
    PetscFunctionBeginUser;
    
    // ---------------------------------------------------------------
    // Compute projection
    // ---------------------------------------------------------------
    // f=1 --> face012 with local orientation v0->v1->v2
    PetscInt f = 0;
    LampF[f][0] = Lam[0];
    LampF[f][1] = Lam[1];
    LampF[f][2] = Lam[2];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampF[f][i][0] = DLam[i][0];
        DLampF[f][i][1] = DLam[i][1];
        DLampF[f][i][2] = DLam[i][2];
    }
    
    // f=2 --> face013 with local orientation v0->v1->v3
    f = 1;
    LampF[f][0] = Lam[0];
    LampF[f][1] = Lam[1];
    LampF[f][2] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampF[f][i][0] = DLam[i][0];
        DLampF[f][i][1] = DLam[i][1];
        DLampF[f][i][2] = DLam[i][3];
    }
    
    // f=3 --> face123 with local orientation v1->v2->v3
    f = 2;
    LampF[f][0] = Lam[1];
    LampF[f][1] = Lam[2];
    LampF[f][2] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampF[f][i][0] = DLam[i][1];
        DLampF[f][i][1] = DLam[i][2];
        DLampF[f][i][2] = DLam[i][3];
    }
    
    // f=4 --> face023 with local orientation v0->v2->v3
    f = 3;
    LampF[f][0] = Lam[0];
    LampF[f][1] = Lam[2];
    LampF[f][2] = Lam[3];
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++) {
        DLampF[f][i][0] = DLam[i][0];
        DLampF[f][i][1] = DLam[i][2];
        DLampF[f][i][2] = DLam[i][3];
    }
    
    *IdecF = PETSC_FALSE;

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode OrientTri(PetscReal S[NUM_DIMENSIONS], PetscReal DS[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscInt Nori, PetscReal GS[NUM_DIMENSIONS], PetscReal GDS[NUM_DIMENSIONS][NUM_DIMENSIONS]){
    /*Compute the local to global transformations of edges.

    :param ndarray S: projection of affine coordinates on faces
    :param ndarray DS: projection of gradients of affine coordinates on faces
    :param ndarray Nori: face orientation
    :param int N: number of dimensions
    :return: global transformation of faces and global transformation of gradients of faces
    :rtype: ndarray
    */
    PetscFunctionBeginUser;

    PetscInt Or[NUM_DIMENSIONS*2][NUM_DIMENSIONS];
    
    // Nori=0 => (s0,s1,s2)->(s0,s1,s2)
    Or[0][0] = 0;
    Or[0][1] = 1;
    Or[0][2] = 2;
    // Nori=1 => (s0,s1,s2)->(s1,s2,s0)
    Or[1][0] = 1;
    Or[1][1] = 2;
    Or[1][2] = 0;
    // Nori=2 => (s0,s1,s2)->(s2,s0,s1)
    Or[2][0] = 2;
    Or[2][1] = 0;
    Or[2][2] = 1;
    // Nori=3 => (s0,s1,s2)->(s0,s2,s1)
    Or[3][0] = 0;
    Or[3][1] = 2;
    Or[3][2] = 1;
    // Nori=4 => (s0,s1,s2)->(s1,s0,s2)
    Or[4][0] = 1;
    Or[4][1] = 0;
    Or[4][2] = 2;
    // Nori=5 => (s0,s1,s2)->(s2,s1,s0)
    Or[5][0] = 2;
    Or[5][1] = 1;
    Or[5][2] = 0;
    
    // Local-to-global transformation
    GS[0] = S[Or[Nori][0]];
    GS[1] = S[Or[Nori][1]];
    GS[2] = S[Or[Nori][2]];
    
    for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        GDS[i][0] = DS[i][Or[Nori][0]];
        GDS[i][1] = DS[i][Or[Nori][1]];
        GDS[i][2] = DS[i][Or[Nori][2]];
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode PolyIJacobi(PetscReal X, PetscReal T, PetscInt Nord, PetscInt Minalpha, PetscBool Idec, PetscReal **L, PetscReal **P, PetscReal **R){
    /*Compute values of integrated shifted scaled Jacobi polynomials and their derivatives starting with p=1.

    Result is 'half' of a  matrix with each row  associated to a fixed alpha.
    Alpha grows by 2 in each row.

    :param ndarray X: coordinate from [0,1]
    :param ndarray T: scaling parameter
    :param int Nord: max polynomial order
    :param int Minalpha: first row value of alpha
    :param bool Idec: decision flag to compute (= FALSE polynomials with x and t derivatives, = TRUE  polynomials with x derivatives only)
    :return: polynomial values, derivatives in x (Jacobi polynomials), derivatives in t
    */
    PetscFunctionBeginUser;
    
    // clearly (minI,maxI)=(1,Nord), but the syntax is written as it is
    // because it reflects how the indexing is called from outside
    PetscInt minI = 0;
    PetscInt maxI = minI + Nord;
    PetscReal *alpha;
    PetscReal **ptemp;

    // Allocate 
    PetscCall(PetscCalloc1(Nord, &alpha));

    PetscCall(PetscCalloc1(Nord + 1, &ptemp));
    for (PetscInt i = 0; i < Nord + 1; i++){
        PetscCall(PetscCalloc1(Nord + 1, &ptemp[i]));
    }
    
    PetscCall(PolyJacobi(X, T, Nord, Minalpha, ptemp));
    
    // Define P. Note that even though P is defined at all entries,
    // because of the way Jacobi computes ptemp, only the necessary entries,
    // and those on the first subdiagonal (which are never used later)
    // are actually accurate.
    for(PetscInt i = minI; i < maxI; i++){
        for(PetscInt j = 0; j < Nord; j++){
            P[i][j] = ptemp[i][j];
        }
    }

    // Create vector alpha first
    for(PetscInt i = 0; i < maxI; i++){
        alpha[i] = Minalpha + 2*(i - minI);
    }
    
    // Initiate first column (order 1 in L)
    for(PetscInt i = minI; i < maxI; i++){
        L[i][0] = X;
    }
    
    // General case; compute R
    for(PetscInt i = minI; i < maxI; i++){
        for(PetscInt j = 0; j < Nord; j++){
            R[i][j] = 0;
        }
    }
    
    // Fill the last columns if necessary
    if(Nord >= 2){
        PetscReal tt = pow(T, 2);
        PetscInt ni = -1;
        for(PetscInt i = 0; i < maxI - 1; i++){
            PetscReal al = alpha[i];
            ni += 1;
            for(PetscInt j = 2; j < Nord - ni + 1; j++){
                PetscReal tia = j+j+al;
                PetscReal tiam1 = tia-1;
                PetscReal tiam2 = tia-2;
                PetscReal ai = (j+al)/(tiam1*tia);
                PetscReal bi = (al)/(tiam2*tia);
                PetscReal ci = (j-1)/(tiam2*tiam1);
                L[i][j-1] = ai*ptemp[i][j]+bi*T*ptemp[i][j-1]-ci*tt*ptemp[i][j-2];
                R[i][j-1] = -(j-1)*(ptemp[i][j-1]+T*ptemp[i][j-2]);
                R[i][j-1] = R[i][j-1]/tiam2;
            }
        }
    }
    
    // Free memory
    PetscCall(PetscFree(alpha));
    
    for (PetscInt i = 0; i < Nord + 1; i++){
        PetscCall(PetscFree(ptemp[i]));
    }
    PetscCall(PetscFree(ptemp));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode HomIJacobi(PetscReal S[2], PetscReal DS[NUM_DIMENSIONS][2], PetscInt Nord, PetscInt Minalpha, PetscBool Idec, PetscReal **HomL, PetscReal ***DHomL){
    /*Compute values of integrated homogenized Jacobi polynomials and their gradients.
    Result is half of a  matrix with each row  associated to a fixed alpha.
    Alpha grows by 2 in each row.

    :param ndarray S: (s0,s1) affine(like) coordinates
    :param ndarray DS: gradients of S in R(N)
    :param int Nord: max polynomial order
    :param int Minalpha: first row value of alpha (integer)
    :param bool Idec: decision flag to compute
    :return: polynomial values and derivatives in x (Jacobi polynomials)
    :rtype: ndarray
    */
    PetscFunctionBeginUser;
    
    // clearly (minI,maxI)=(1,Nord), but the syntax is written as it is
    // because it reflects how the indexing is called from outside
    PetscInt minI = 1;
    PetscInt maxI = minI+Nord-1;
    
    PetscReal **homP;   // homP[Nord][Nord]
    PetscReal **homR;   // homR[Nord][Nord] 

    // Allocate
    PetscCall(PetscCalloc1(Nord, &homP));
    for (PetscInt i = 0; i < Nord; i++){
        PetscCall(PetscCalloc1(Nord, &homP[i]));
    }
     
    PetscCall(PetscCalloc1(Nord, &homR));
    for (PetscInt i = 0; i < Nord; i++){
        PetscCall(PetscCalloc1(Nord, &homR[i]));
    }
        
    PetscInt ni = -1;
    
    if(Idec){
        PetscCall(PolyIJacobi(S[1], 1, Nord, Minalpha, Idec, HomL, homP, homR));
        for(PetscInt i = minI; i < maxI + 1; i++){
            ni += 1;
            for(PetscInt j = 1; j < Nord - ni + 1; j++){
                for(PetscInt k = 0; k < NUM_DIMENSIONS; k++){
                    DHomL[k][i-1][j-1] = homP[i-1][j-1] * DS[k][1];
                }
            }
        }
    } else {
        // If sum of S different from 1 -> Idec=.FALSE.
        
        PetscCall(PolyIJacobi(S[1], S[0] + S[1], Nord, Minalpha, Idec, HomL, homP, homR));
        
        PetscReal DS01[NUM_DIMENSIONS];
        for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
            DS01[i] = DS[i][0] + DS[i][1];
        }
        
        for(PetscInt i = minI; i < maxI + 1; i++){
            ni += 1;
            for(PetscInt j = 1; j < Nord - ni + 1; j++){
                for(PetscInt k = 0; k < NUM_DIMENSIONS; k++){
                    DHomL[k][i-1][j-1] = homP[i-1][j-1] * DS[k][1] + homR[i-1][j-1]*DS01[k];
                }
            }
        }    
    }
    
    // Free memory
    for (PetscInt i = 0; i < Nord; i++){
        PetscCall(PetscFree(homP[i]));
    }
    PetscCall(PetscFree(homP));

    for (PetscInt i = 0; i < Nord; i++){
        PetscCall(PetscFree(homR[i]));
    }
    PetscCall(PetscFree(homR));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode AncETri(PetscReal S[NUM_DIMENSIONS], PetscReal DS[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscInt Nord, PetscBool Idec, PetscReal ***ETri, PetscReal ***CurlETri){
    /*Compute triangle face Hcurl ancillary functions and their curls.

    :param ndarray S: (s0,s1,s2) affine coordinates associated to triangle face
    :param ndarray DS: derivatives of S0,S1,S2
    :param int Nord: polynomial order
    :param bool Idec: Binary flag:
    :param int N: spatial dimension
    :return: triangle Hcurl ancillary functions and curls of triangle Hcurl ancillary functions
    :rtype: ndarray
    */
    PetscFunctionBeginUser;
    
    PetscReal DsL[NUM_DIMENSIONS][2];
    PetscReal sL[2];

    // Local parameters
    PetscInt minI = 0;
    PetscInt minJ = 1;
    PetscInt maxJ = Nord-1;
    PetscInt maxIJ = Nord-1;
    PetscInt minalpha = 2*minI+1;
    PetscInt Ncurl = 2*NUM_DIMENSIONS-3;
    PetscBool IdecE = PETSC_FALSE;
    
    // get EE - this is never a simplified case (IdecE=0)
    PetscReal tempS[2] = {S[0], S[1]};
    PetscReal tempDS[NUM_DIMENSIONS][2];
    for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        tempDS[i][0] = DS[i][0];
        tempDS[i][1] = DS[i][1];
    }
    
    PetscReal **EE;     // EE[NUM_DIMENSIONS][Nord-minJ]
    PetscReal **curlEE; // curlEE[2*NUM_DIMENSIONS-3][Nord-minJ]

    // Allocate
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &EE));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(Nord-minJ, &EE[i]));
    }

    PetscCall(PetscCalloc1(2*NUM_DIMENSIONS-3, &curlEE));
    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        PetscCall(PetscCalloc1(Nord-minJ, &curlEE[i]));
    }

    PetscCall(AncEE(tempS, tempDS, Nord-minJ, IdecE, EE, curlEE));
    
    // get homogenized Integrated Jacobi polynomials, homLal, and gradients
    sL[0] = S[0]+S[1];
    sL[1] = S[2];
    for(PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        DsL[i][0] = DS[i][0] + DS[i][1];
        DsL[i][1] = DS[i][2];
    }
    
    PetscReal **homLal;     // homLal[maxJ][maxJ]
    PetscReal ***DhomLal;   // DhomLal[NUM_DIMENSIONS][maxJ][maxJ]

    // Allocate
    PetscCall(PetscCalloc1(maxJ, &homLal));
    for (PetscInt i = 0; i < maxJ; i++){
        PetscCall(PetscCalloc1(maxJ, &homLal[i]));
    }

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &DhomLal));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(maxJ, &DhomLal[i]));
        for (PetscInt j = 0; j < maxJ; j++){
            PetscCall(PetscCalloc1(maxJ, &DhomLal[i][j]));
        }
    }
    
    PetscCall(HomIJacobi(sL, DsL, maxJ, minalpha, Idec, homLal, DhomLal));

    // Simply complete the required information
    for(PetscInt i = 0; i < maxIJ + 1; i++){
        for(PetscInt j = minI; j < i-minJ+1; j++){
            PetscInt k = i - j;
            for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                ETri[n][j][k-1] = EE[n][j] * homLal[j][k-1];
            }
            
            PetscReal DhomLalxEE[NUM_DIMENSIONS];
            PetscReal temp1[NUM_DIMENSIONS], temp2[NUM_DIMENSIONS];
            for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                temp1[n] = DhomLal[n][j][k-1];
                temp2[n] = EE[n][j];
            }
            
            PetscCall(crossProduct(temp1, temp2, DhomLalxEE));
            
            for(PetscInt n = 0; n < Ncurl; n++){
                CurlETri[n][j][k-1] = homLal[j][k-1]*curlEE[n][j] + DhomLalxEE[n];
            }
        }
    }
    
    // Free memory
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(EE[i]));   
    }
    PetscCall(PetscFree(EE));   

    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        PetscCall(PetscFree(curlEE[i]));   
    }
    PetscCall(PetscFree(curlEE));   

    for (PetscInt i = 0; i < maxJ; i++){
        PetscCall(PetscFree(homLal[i]));   
    }
    PetscCall(PetscFree(homLal));

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < maxJ; j++){
            PetscCall(PetscFree(DhomLal[i][j]));   
            }
        PetscCall(PetscFree(DhomLal[i]));       
    }
    PetscCall(PetscFree(DhomLal));       

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode computeCellOrientation(DM dm, PetscInt cell, PetscInt cellOrientation[10]){
    
    PetscFunctionBeginUser;


    PetscInt cellFaces[NUM_FACES_PER_ELEMENT];
    PetscInt cellEdges[NUM_EDGES_PER_ELEMENT];
    PetscInt cellVertices[NUM_EDGES_PER_ELEMENT];
    PetscInt faceEdges[NUM_FACES_PER_ELEMENT][NUM_EDGES_PER_FACE];
    PetscInt faceVertices[NUM_FACES_PER_ELEMENT][NUM_VERTICES_PER_FACE];
    PetscInt edgeVertices[NUM_EDGES_PER_ELEMENT][NUM_VERTICES_PER_EDGE];

    PetscInt transitiveClosureCellSize;
    PetscInt  *transitiveClosureCellPoints = NULL;


    PetscInt transitiveClosureFaceSize;
    PetscInt *transitiveClosureFacePoints = NULL;
    const PetscInt *conePoints;
    PetscInt currentPoint;
    PetscInt currentFace;   

    
    /* Get transitive clousure for cell */
    /* Orden convention:
        - Faces indices start on position 2
        - Edges indices start on position 2 + NUM_FACES_PER_ELEMENT * 2
        - Vertices indices start on position 2 + NUM_FACES_PER_ELEMENT * 2 + NUM_EDGES_PER_ELEMENT * 2
    */
    PetscCall(DMPlexGetTransitiveClosure(dm, cell, PETSC_TRUE, &transitiveClosureCellSize, &transitiveClosureCellPoints));

    /* Get faces indices for cell, edges for each face, and vertices for each face */
    currentPoint = 2; 
    for(PetscInt i = 0; i < NUM_FACES_PER_ELEMENT; i++) {
            /* Face indices */
            cellFaces[i] = transitiveClosureCellPoints[currentPoint + i*2];
            PetscCall(DMPlexGetCone(dm, cellFaces[i], &conePoints));
            
            /* Edges for each face*/
            for(PetscInt j = 0; j < NUM_EDGES_PER_FACE; j++){
                faceEdges[i][j] = conePoints[j];
            }

            /* Vertices for each face */ 
            /* Orden convention:
                - Edges indices start on position 2
                - Vertices indices start on position 2 + NUM_EDGES_PER_FACE * 2
            */
            PetscCall(DMPlexGetTransitiveClosure(dm, cellFaces[i], PETSC_TRUE, &transitiveClosureFaceSize, &transitiveClosureFacePoints));
            
            currentFace = 8;
            for(PetscInt j = 0; j < NUM_VERTICES_PER_FACE; j++){
                faceVertices[i][j] = transitiveClosureFacePoints[currentFace + j*2];
            }
            
            PetscCall(DMPlexRestoreTransitiveClosure(dm, cellFaces[i], PETSC_TRUE, &transitiveClosureFaceSize, &transitiveClosureFacePoints));
    }

    /* Get edges indices for cell */
    currentPoint = 2 + NUM_FACES_PER_ELEMENT*2; 
    for(PetscInt i = 0; i < NUM_EDGES_PER_ELEMENT; i++) {
        cellEdges[i] = transitiveClosureCellPoints[currentPoint + i*2];
        
        /* Get edge vertices */
        PetscCall(DMPlexGetCone(dm, cellEdges[i], &conePoints));
        
        for(PetscInt j = 0; j < NUM_VERTICES_PER_EDGE; j++){
            edgeVertices[i][j] = conePoints[j];
            }
    }

    /* Get vertices indices for cell */
    currentPoint = 2 + NUM_FACES_PER_ELEMENT*2 + NUM_EDGES_PER_ELEMENT*2;
    for(PetscInt i = 0; i < NUM_VERTICES_PER_ELEMENT; i++) {
        cellVertices[i] = transitiveClosureCellPoints[currentPoint + i*2];
    }

    /* Restore transitive closure */
    PetscCall(DMPlexRestoreTransitiveClosure(dm, cell, PETSC_TRUE, &transitiveClosureCellSize, &transitiveClosureCellPoints));

    /* Compute orientation for faces */
    for(PetscInt i = 0; i < NUM_FACES_PER_ELEMENT; i++) {
        cellOrientation[i] = 0;
    }


    PetscInt localNodes[2];

    for(PetscInt i = 0; i < NUM_EDGES_PER_ELEMENT; i++) {
        switch (i) {
            case 0:
                localNodes[0] = 0; localNodes[1] = 1;
                break;
            case 1:
                localNodes[0] = 1; localNodes[1] = 2;
                break;
            case 2:
                localNodes[0] = 2; localNodes[1] = 0;
                break;
            case 3:
                localNodes[0] = 0; localNodes[1] = 3;
                break;
            case 4:
                localNodes[0] = 3; localNodes[1] = 1;
                break;
            case 5:
                localNodes[0] = 2; localNodes[1] = 3;
                break;
        }

        // Get the global node indices for this edge
        PetscInt nodesEleForThisEdge[2] = {cellVertices[localNodes[0]], cellVertices[localNodes[1]]};
        PetscInt globalNodesInEdge[2] = {edgeVertices[i][0], edgeVertices[i][1]};

        // Determine the orientation for this edge
        PetscInt orientationForThisEdge = 0;
        if ((nodesEleForThisEdge[0] == globalNodesInEdge[1]) && (nodesEleForThisEdge[1] == globalNodesInEdge[0])) {
            orientationForThisEdge = 1;  // Edge is inverted
        }
        cellOrientation[4 + i] = orientationForThisEdge;
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode shape3DETet(PetscReal X[NUM_DIMENSIONS], PetscInt nord, PetscInt cellOrientation[10], PetscReal **ShapE, PetscReal **CurlE){
    /*Compute values of 3D tetrahedron element H(curl) shape functions and their derivatives.

    :param ndarray X: master tetrahedron coordinates from (0,1)^3
    :param int Nord: polynomial order
    :param ndarray NoriE: edge orientation
    :param ndarray NoriF: face orientation
    :return: number of dof, values of the shape functions at the point, curl of the shape functions
    :rtype: ndarray.

    .. note:: References:\n
       Fuentes, F., Keith, B., Demkowicz, L., & Nagaraj, S. (2015). Orientation
       embedded high order shape functions for the exact sequence elements of
       all shapes. Computers & Mathematics with applications, 70(4), 353-458.
    */
    PetscFunctionBeginUser;
   
    // Local parameters
    PetscBool IdecB[2] = {PETSC_FALSE, PETSC_FALSE};
    PetscInt minI = 0;
    PetscInt minJ = 1;
    PetscInt minK = 1;
    PetscInt minIJ = minI + minJ;
    PetscInt minIJK = minIJ + minK;

    // Initialize counter for shape functions
    PetscInt m = 0;
    
    // Define affine coordinates and gradients
    PetscReal Lam[NUM_DIMENSIONS + 1] = {0.0};
    PetscReal DLam[NUM_DIMENSIONS][NUM_DIMENSIONS + 1] = {{0.0}};

    PetscCall(AffineTetrahedron(X, Lam, DLam));
    
    /* Shape functions over edges */
    PetscReal LampE[NUM_EDGES_PER_ELEMENT][2];
    PetscReal DLampE[NUM_EDGES_PER_ELEMENT][NUM_DIMENSIONS][2];
    PetscBool IdecE;
    
    // Compute edges projection
    PetscCall(ProjectTetE(Lam, DLam, LampE, DLampE, &IdecE));

    // Polynomial order and orientation for faces and edges 
    PetscInt NoriF[4];  // Orientation for faces
    PetscInt NoriE[6];  // Orientation for 


    // Extract orientation for faces
    for (PetscInt i = 0; i < NUM_FACES_PER_ELEMENT; ++i){
        NoriF[i] = cellOrientation[i];
    }

    /* Extract orientation for edges */
    for (PetscInt i = 0; i < NUM_EDGES_PER_ELEMENT; ++i){
        NoriE[i] = cellOrientation[i+4]; 
    }


    /* Allocate memory for shape functions on edges */
    PetscReal **EE;       // EE[NUM_DIMENSIONS][nordEdge];
    PetscReal **CurlEE;   // CurlEE[2*NUM_DIMENSIONS-3][nordEdge];

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &EE));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(nord, &EE[i]));
    }

    PetscCall(PetscCalloc1(2*NUM_DIMENSIONS-3, &CurlEE));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(nord, &CurlEE[i]));
    }


    /* Reset matrices */
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < NUM_DIMENSIONS; j++){
            ShapE[i][j] = 0.0;
            CurlE[i][j] = 0.0;
        }
    }

    
    /* Shape functions for edges */
    PetscInt nordEdge = 0;
    PetscInt numDofEdge = 0;    
    for(PetscInt i = 0; i < NUM_EDGES_PER_ELEMENT; i++){
        // Local parameters
        nordEdge = nord;
        numDofEdge = nordEdge;
        if(numDofEdge > 0){
            // Local parameters
            PetscInt maxI = nordEdge - 1;
            // Orient 
            PetscReal GLampE[2] = {0.0};
            PetscReal GDLampE[3][2] = {{0.0}};
            PetscReal S[2] = {0.0} ;
            PetscReal D[NUM_DIMENSIONS][2] = {{0.0}};

            S[0] = LampE[i][0];
            S[1] = LampE[i][1];

            /* Extract the slice into D */
            for (PetscInt j = 0; j < NUM_DIMENSIONS; j++) {
                for (PetscInt k = 0; k < 2; k++) {
                    D[j][k] = DLampE[i][j][k];
                }
            }

            PetscCall(OrientE(S, D, NoriE[i], GLampE, GDLampE));

            /* Construct the shape functions */
            PetscCall(AncEE(GLampE, GDLampE, nordEdge, IdecE, EE, CurlEE));

            for(PetscInt j = minI; j < maxI + 1; j++){
                for(PetscInt k = 0; k < NUM_DIMENSIONS; k++){
                    ShapE[k][m] = EE[k][j];
                    CurlE[k][m] = CurlEE[k][j];
                }
                m += 1;
            }
        }
    }

    // Free memory for shape functions on edges
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(EE[i]));
    }
    PetscCall(PetscFree(EE));

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(CurlEE[i]));
    }
    PetscCall(PetscFree(CurlEE));

    /* Shape functions over faces */ 
    PetscReal LampF[NUM_FACES_PER_ELEMENT][NUM_DIMENSIONS];
    PetscReal DLampF[NUM_FACES_PER_ELEMENT][NUM_DIMENSIONS][NUM_DIMENSIONS];
    PetscBool IdecF;    

    // Compute faces projection
    PetscCall(ProjectTetF(Lam, DLam, LampF, DLampF, &IdecF));

    // Allocate memory for funcions on faces
    PetscReal ***ETri;       // ETri[NUM_DIMENSIONS][nord - 1][nord - 1]
    PetscReal ***CurlETri;   // CurlETri[2*NUM_DIMENSIONS-3][nord-1][nord-1]

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &ETri));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(nord-1, &ETri[i]));
        for (PetscInt j = 0; j < nord-1; j++){
            PetscCall(PetscCalloc1(nord-1, &ETri[i][j]));
        }
    }

    PetscCall(PetscCalloc1(2*NUM_DIMENSIONS-3, &CurlETri));
    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        PetscCall(PetscCalloc1(nord-1, &CurlETri[i]));
        for (PetscInt j = 0; j < nord-1; j++){
            PetscCall(PetscCalloc1(nord-1, &CurlETri[i][j]));
        }
    }    

    // Shape functions for faces
    PetscInt nordFace = 0;
    PetscInt numDofFace = 0;

    for(PetscInt i = 0; i < NUM_FACES_PER_ELEMENT; i++){
        // Local parameters
        nordFace = nord;
        numDofFace = nordFace*(nordFace-1)/2;
        if(numDofFace > 0){
            // Local parameters (again)
            PetscInt maxIJ = nordFace - 1;

            // Orient
            PetscReal GLampF[3];
            PetscReal GDLampF[3][3];
            PetscReal tmpLampF[3];
            PetscReal tempDLampF[3][3];

            // Prepare input matrices
            for (PetscInt j = 0; j<3; j++){
                tmpLampF[j] = LampF[i][j];
                for (PetscInt k = 0; k<3; k++){
                    tempDLampF[j][k] = DLampF[i][j][k];
                }
            }

            PetscCall(OrientTri(tmpLampF, tempDLampF, NoriF[i], GLampF, GDLampF));

            // Loop over families
            PetscInt famctr = m;
            for(PetscInt j = 0; j < 2; j++){
                m = famctr + j - 1;
                PetscInt abc[3];
                for(PetscInt k = 0; k < 3; k++){
                    PetscInt pos = (k - j) % 3;
                    if (pos < 0){
                        pos += 3;
                    } 
                    abc[pos] = k;
                }
        
                PetscReal tempGLampF[3];
                PetscReal tempGDLampF[3][3];
                for(PetscInt k = 0; k < 3; k++){
                    tempGLampF[k] = GLampF[abc[k]];
                    for(PetscInt t = 0; t < NUM_DIMENSIONS; t++){
                        tempGDLampF[t][k] = GDLampF[t][abc[k]];
                    }
                }

                // Construct the shape functions
                PetscCall(AncETri(tempGLampF, tempGDLampF, nordFace, IdecF, ETri, CurlETri));

                 for(PetscInt k = minIJ; k < maxIJ + 1; k++){
                    for(PetscInt r = minI; r < k-minJ+1; r++){
                        PetscInt p = k - r;
                        m += 2;
                        for(PetscInt t = 0; t < NUM_DIMENSIONS; t++){
                            ShapE[t][m-1] = ETri[t][r][p-1];
                            CurlE[t][m-1] = CurlETri[t][r][p-1];
                        }
                    }
                }
            }
        }
    }

    // Free memory
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < nord-1; j++){
            PetscCall(PetscFree(ETri[i][j]));
        }
        PetscCall(PetscFree(ETri[i]));
    }
    PetscCall(PetscFree(ETri));   

    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        for (PetscInt j = 0; j < nord-1; j++){
            PetscCall(PetscFree(CurlETri[i][j]));
        }
        PetscCall(PetscFree(CurlETri[i]));
    }
    PetscCall(PetscFree(CurlETri));

    /* Shape functions over volume */
    PetscInt nordB = nord;
    PetscInt ndofB = nordB*(nordB-1)*(nordB-2)/6;
    PetscInt minbeta = 2*minIJ;
    PetscInt maxIJK = nordB-1;
    PetscInt maxK = maxIJK-minIJ;

    // Allocate memory for shape functions in volume
    PetscReal ***ETriV;         // ETri[NUM_DIMENSIONS][nord-minK-1][nord-minK-1]
    PetscReal ***CurlETriV;     // CurlETri[2*NUM_DIMENSIONS-3][nord-minK-1][nord-minK-1]
    PetscReal **homLbet;        // homLbet[maxK][maxK] 
    PetscReal ***DhomLbet;      // DhomLbet[NUM_DIMENSIONS][maxK][maxK]    

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &ETriV));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(nord-minK-1, &ETriV[i]));
        for (PetscInt j = 0; j < nord-minK-1; j++){
            PetscCall(PetscCalloc1(nord-minK-1, &ETriV[i][j]));
        }
    }

    PetscCall(PetscCalloc1(2*NUM_DIMENSIONS-3, &CurlETriV));
    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        PetscCall(PetscCalloc1(nord-minK-1, &CurlETriV[i]));
        for (PetscInt j = 0; j < nord-minK-1; j++){
            PetscCall(PetscCalloc1(nord-minK-1, &CurlETriV[i][j]));
        }
    }

    PetscCall(PetscCalloc1(maxK, &homLbet));
    for (PetscInt i = 0; i < maxK; i++){
        PetscCall(PetscCalloc1(maxK, &homLbet[i]));
    }

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &DhomLbet));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(maxK, &DhomLbet[i]));
        for (PetscInt j = 0; j < maxK; j++){
            PetscCall(PetscCalloc1(maxK, &DhomLbet[i][j]));
        }
    }    

    // If necessary, create bubbles
    if(ndofB > 0){
        // Local parameters (again)
        IdecB[0] = IdecF;
        IdecB[1] = PETSC_TRUE;
        
        // Loop over families
        PetscInt famctr = m;
        for(PetscInt i = 0; i < 3; i++){
            m = famctr + i - 2;
            PetscInt abcd[4];
            for(PetscInt j = 0; j < 4; j++){
                PetscInt pos = (j - i) % 4;
                if (pos < 0){
                    pos += 4;
                } 
                abcd[pos] = j;
            }
            PetscInt abc[3] = {abcd[0], abcd[1], abcd[2]};
            PetscInt d = abcd[3];

            // Now construct the shape functions (no need to orient)
            PetscReal tempLam[3];
            PetscReal tempDLam[3][3];
            for(PetscInt j = 0; j < 3; j++){
                tempLam[j] = Lam[abc[j]];
                for(PetscInt k = 0; k < NUM_DIMENSIONS; k++){
                    tempDLam[k][j] = DLam[k][abc[j]];
                }
            }
            
            PetscCall(AncETri(tempLam, tempDLam, nordB-minK, IdecB[0], ETriV, CurlETriV));

            PetscReal tmp1[2] = {1-Lam[d], Lam[d]};
            PetscReal tmp2[3][2];

            // Initialize input matrix
            for(PetscInt j = 0; j < NUM_DIMENSIONS; j++){
                tmp2[j][0] = -DLam[j][d]; 
                tmp2[j][1] = DLam[j][d]; 
            }

            PetscCall(HomIJacobi(tmp1, tmp2, maxK, minbeta, IdecB[1], homLbet, DhomLbet));

            for(PetscInt j = minIJK; j < maxIJK+1; j++){
                for(PetscInt k = minIJ; k < j-minK+1; k++){
                    for(PetscInt r = minI; r < k-minJ+1; r++){
                        PetscInt p = k - r;
                        PetscInt q = j - k;
                        m += 3;

                        for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                            ShapE[n][m-1] = ETriV[n][r][p-1]*homLbet[k-1][q-1];
                        }

                        PetscReal DhomLbetxETri[NUM_DIMENSIONS];
                        PetscReal v1[NUM_DIMENSIONS], v2[NUM_DIMENSIONS];
    
                        for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                            v1[n] = DhomLbet[n][k-1][q-1];
                            v2[n] = ETriV[n][r][p-1];
                        }
    
                        PetscCall(crossProduct(v1, v2, DhomLbetxETri));
    
                        for(PetscInt n = 0; n < NUM_DIMENSIONS; n++){
                            CurlE[n][m-1] = homLbet[k-1][q-1]*CurlETriV[n][r][p-1] + DhomLbetxETri[n];
                        }
                    }
                }
            }
        }
    }

    /* Compute reordering vector (PETGEM to PETSc convention) */
    /*PetscInt tmp1[] = {1, 2, 3, 4, 5, 6};
    PetscInt tmp2[] = {13, 14, 15, 16, 19, 20, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    PetscInt tmp3[] = {43, 44, 45, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 37, 38, 39, 40, 41, 42, 31, 32, 33, 34, 35, 36, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    PetscInt tmp4[] = {73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
    PetscInt tmp5[] = {111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    PetscInt tmp6[] = {157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36};*/


    // PetscInt tmp1[] = {0, 1, 2, 3, 4, 5}; 
    // PetscInt tmp2[] = {12, 13, 14, 15, 18, 19, 16, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; // Subtracted values from original array
    // PetscInt tmp3[] = {42, 43, 44, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 36, 37, 38, 39, 40, 41, 30, 31, 32, 33, 34, 35, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}; // Subtracted values from original array
    // PetscInt tmp4[] = {72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23}; // Subtracted values from original array
    // PetscInt tmp5[] = {110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29}; // Subtracted values from original array
    // PetscInt tmp6[] = {156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35}; // Subtracted values from original array

    // // Print the updated arrays
 


    // PetscInt numDofInCell; 
    // numDofInCell = nord*(nord+2)*(nord+3)/2;
    // PetscInt orderPermutation[numDofInCell]; 

    // switch (nord){
    // case 1:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp1[i];
    //     }
    //     break;
    // case 2:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp2[i];
    //     }
    //     break;
    // case 3:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp3[i];
    //     }
    //     break;
    // case 4:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp4[i];
    //     }
    //     break;
    // case 5:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp5[i];
    //     }
    //     break;
    // case 6:
    //     for (PetscInt i=0; i<numDofInCell; i++){
    //         orderPermutation[i] = tmp6[i];
    //     }
    //     break;        
    // }


    // /* Copy basis and curl from PETGEM order convention */
    // PetscReal tmpShapE[NUM_DIMENSIONS][numDofInCell], tmpCurlE[NUM_DIMENSIONS][numDofInCell]; 

    // for (PetscInt i = 0; i<NUM_DIMENSIONS; i++){
    //     for (PetscInt j = 0; j<numDofInCell; j++){
    //         tmpShapE[i][j] = ShapE[i][j];
    //         tmpCurlE[i][j] = CurlE[i][j];        
    //     }    
    // }

    // /* Apply PETSc ordering */
    // for (PetscInt i = 0; i<NUM_DIMENSIONS; i++){
    //     for (PetscInt j = 0; j<numDofInCell; j++){
    //         ShapE[i][j] = tmpShapE[i][orderPermutation[j]];
    //         CurlE[i][j] = tmpCurlE[i][orderPermutation[j]];
    //     }
    // }

    /* Free memory */
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < nordB-minK-1; j++){
            PetscCall(PetscFree(ETriV[i][j]));
        }
        PetscCall(PetscFree(ETriV[i]));
    }
    PetscCall(PetscFree(ETriV));   

    for (PetscInt i = 0; i < 2*NUM_DIMENSIONS-3; i++){
        for (PetscInt j = 0; j < nordB-minK-1; j++){
            PetscCall(PetscFree(CurlETriV[i][j]));
        }
        PetscCall(PetscFree(CurlETriV[i]));
    }
    PetscCall(PetscFree(CurlETriV));
    
    for (PetscInt i = 0; i < maxK; i++){
        PetscCall(PetscFree(homLbet[i]));   
    }
    PetscCall(PetscFree(homLbet));

    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        for (PetscInt j = 0; j < maxK; j++){
            PetscCall(PetscFree(DhomLbet[i][j]));   
        }
        PetscCall(PetscFree(DhomLbet[i]));       
    }
    PetscCall(PetscFree(DhomLbet));           

    PetscFunctionReturn(PETSC_SUCCESS);
}

        
PetscErrorCode computeElementalMatrix(PetscInt nord, PetscInt cellOrientation[10], PetscReal jacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal invJacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscInt numGaussPoints, PetscReal **gaussPoints, PetscReal *weigths, PetscReal *cellResistivity, PetscReal **Me, PetscReal **Ke){
    PetscFunctionBeginUser;

    /* Initial declarations */
    PetscReal e_r[NUM_DIMENSIONS][NUM_DIMENSIONS]  = {{0.0}}; 
    PetscReal mu_r[NUM_DIMENSIONS][NUM_DIMENSIONS] = {{0.0}}; 
    PetscReal iPoint[NUM_DIMENSIONS] = {0.0};
    PetscReal det;
    PetscReal temp1[NUM_DIMENSIONS];
    PetscReal temp2[NUM_DIMENSIONS];
    PetscReal temp3[NUM_DIMENSIONS];
    PetscReal temp4[NUM_DIMENSIONS];
    PetscReal dotResult; 

    /* Tensor for integration (Vertical transverse electric permitivity) */
    e_r[0][0] = cellResistivity[0];
    e_r[1][1] = cellResistivity[1];
    e_r[2][2] = cellResistivity[2];

    /* Tensor for integration (Constant magnetic permittivity) */
    mu_r[0][0] = 1.0;
    mu_r[1][1] = 1.0;
    mu_r[2][2] = 1.0;

    /* Compute the determinant of the Jacobian */ 
    det = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
        - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
        + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);

    /* Compute number of dofs per cell */
    PetscInt numDofInCell; 
    numDofInCell = nord*(nord+2)*(nord+3)/2;    

    /* Allocate matrices for shape functions */    
    PetscReal **ShapE;
    PetscReal **CurlE;
    PetscReal **NiReal;

    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &ShapE));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &CurlE));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &NiReal));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(numDofInCell, &ShapE[i]));
        PetscCall(PetscCalloc1(numDofInCell, &CurlE[i]));
        PetscCall(PetscCalloc1(numDofInCell, &NiReal[i]));
        
    }

    /* Reset elemental matrices */
    for (PetscInt i = 0; i < numDofInCell; ++i){
        for (PetscInt j = 0; j < numDofInCell; ++j){
                Me[i][j] = 0.0;
                Ke[i][j] = 0.0;
        }
    }

    /* Compute elemental matrices (mass and stifness matrix)*/    
    for (PetscInt i = 0; i < numGaussPoints; ++i){
        // Get gauss for i point 
        iPoint[0] = gaussPoints[i][0];
        iPoint[1] = gaussPoints[i][1];
        iPoint[2] = gaussPoints[i][2]; 

        /* Compute basis function for i point */
        PetscCall(shape3DETet(iPoint, nord, cellOrientation, ShapE, CurlE));

        // NiReal = Ni in real element        
        for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
            for (PetscInt k = 0; k < numDofInCell; ++k){
                NiReal[j][k] = 0.0;
            }
        }

        for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
            for (PetscInt k = 0; k < numDofInCell; ++k){
                for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                    NiReal[j][k] += (invJacobian[j][m] * ShapE[m][k]);
                }
            }
        }

        /* Perform mass matrix integration */
        for (PetscInt j = 0; j < numDofInCell; ++j){
            for (PetscInt k = 0; k < numDofInCell; ++k){
                /* Prepare data */
                for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                    /* Extract slide for row j */
                    temp1[m] = NiReal[m][j];
                    // Extract slide for row k
                    temp2[m] = NiReal[m][k];
                }

                // Perform matrix vector multiplication
                PetscCall(matrixVectorProduct(temp1, e_r, temp3));
                 
                // Perform dot product
                PetscCall(dotProduct(temp3, temp2, &dotResult));

                // Integration
                Me[j][k] += weigths[i] * dotResult * det;
            }
        }

        // Transform curl on reference element to real element
        for (PetscInt j = 0; j < numDofInCell; ++j){
            // Prepare data
            for (PetscInt k = 0; k < NUM_DIMENSIONS; ++k){
                // Extract slide for row j
                temp1[k] = CurlE[k][j];            
            }

            // Perform vector matrix multiplication
            PetscCall(vectorMatrixProduct(temp1, jacobian, temp2));

            // Update data
            for (PetscInt k = 0; k < NUM_DIMENSIONS; ++k){
                // Update slide for row j
                CurlE[k][j] = temp2[k] / det;
            }
        }

        // Perform stiffness matrix integration
        for (PetscInt j = 0; j < numDofInCell; ++j){
            for (PetscInt k = 0; k < numDofInCell; ++k){
                // Prepare data
                for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                    // Extract slide for row j
                    temp1[m] = CurlE[m][j];
                    temp2[m] = CurlE[m][k];                                        
                }

                // Perform matrix vector multiplication
                PetscCall(matrixVectorProduct(temp1, mu_r, temp3));
                 
                // Perform point-wise multiplication
                for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                    temp4[m] = temp2[m] * det;                 
                }

                // Perform dot product 
                PetscCall(dotProduct(temp3, temp4, &dotResult));

                // Integration
                Ke[j][k] += weigths[i] * dotResult;
            }
        }
    }

    /* Free memory */ 
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(ShapE[i]));
        PetscCall(PetscFree(CurlE[i]));
        PetscCall(PetscFree(NiReal[i]));
    }
    PetscCall(PetscFree(ShapE));
    PetscCall(PetscFree(CurlE));   
    PetscCall(PetscFree(NiReal));   

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode computeBasisFunctions(PetscInt nord, PetscInt orientation[10], PetscReal jacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal invJacobian[NUM_DIMENSIONS][NUM_DIMENSIONS], PetscReal *point, PetscReal **basisFunctions, PetscReal **curlBasisFunctions){
    PetscFunctionBeginUser;

    /* Initial declarations */
    PetscReal iPoint[NUM_DIMENSIONS]  = {0.0};
    PetscReal temp1[NUM_DIMENSIONS];
    PetscReal temp2[NUM_DIMENSIONS];
    PetscReal det;

    /* Compute number of dofs per cell */
    PetscInt numDofInCell; 
    numDofInCell = nord*(nord+2)*(nord+3)/2;    

    /* Compute the determinant of the Jacobian */ 
    det = jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
        - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
        + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);


    /* Allocate matrices for shape functions */    
    PetscReal **ShapE;
    PetscReal **CurlE;
    
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &ShapE));
    PetscCall(PetscCalloc1(NUM_DIMENSIONS, &CurlE));
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscCalloc1(numDofInCell, &ShapE[i]));
        PetscCall(PetscCalloc1(numDofInCell, &CurlE[i]));
    }

    // Compute basis functions for iPoint
    iPoint[0] = point[0];
    iPoint[1] = point[1];
    iPoint[2] = point[2];

    // Compute basis function for i point
    PetscCall(shape3DETet(iPoint, nord, orientation, ShapE, CurlE));    

    /* Reset basis functions and curl functions */
    for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
        for (PetscInt k = 0; k < numDofInCell; ++k){
            basisFunctions[j][k] = 0.0;
            curlBasisFunctions[j][k] = 0.0;
        }
    }

    /* NiReal = Ni in real element */        
    for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
        for (PetscInt k = 0; k < numDofInCell; ++k){
            for (PetscInt m = 0; m < NUM_DIMENSIONS; ++m){
                basisFunctions[j][k] += invJacobian[j][m] * ShapE[m][k];
            }
        }
    }

    /* Transform curl on reference element to real element */
    for (PetscInt i = 0; i < numDofInCell; ++i){
        // Prepare data
        for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
            // Extract slide for row j
            temp1[j] = CurlE[j][i];            
        }

        // Perform vector matrix multiplication
        PetscCall(vectorMatrixProduct(temp1, jacobian, temp2));

        // Update data
        for (PetscInt j = 0; j < NUM_DIMENSIONS; ++j){
            // Update slide for row j
            curlBasisFunctions[j][i] = temp2[j] / det;
        }
    }

    /* Free memory */ 
    for (PetscInt i = 0; i < NUM_DIMENSIONS; i++){
        PetscCall(PetscFree(ShapE[i]));
        PetscCall(PetscFree(CurlE[i]));
    }
    PetscCall(PetscFree(ShapE));
    PetscCall(PetscFree(CurlE));   

    PetscFunctionReturn(PETSC_SUCCESS);
}

