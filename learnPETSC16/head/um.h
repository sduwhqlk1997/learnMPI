#ifndef UM_H_
#define UM_H_
#include<petsc.h>
//STARTSTRUCT
// location of one node
typedef struct
{
    PetscReal x,y;
} Node;

//data type for an Unstructured Mesh
typedef struct
{
    PetscInt N, // Number of node
             K, // number of elements
             P; // number of Neumann boundary segments; may be 0
    Vec loc; // nodal locations; length N, dof =2 Vec
    IS e,  // element triples; length 3K
           //    values e[3*k+0], e[3*k+1], e[3*k+2]
       bf, // flag for boundary nodes; length N
           // if bf[i] > 0 then node i is on boundary
           // if bf[i] == 2 then node i is Dirichlet
       ns; // Neumann boundary segment pairs; length 2P
           // may be a null ptr; values s[s*p+0],s[2*p+1]
}UM;
//ENDSTRUCT

//methods below are listed in typical call order

//STARTDECLARE
PetscErrorCode UMInitialize(UM *mesh);
PetscErrorCode UMDestory(UM *mesh);

//create Vec and then read node coordinates from file into it
PetscErrorCode UMReadNodes(UM *mesh, char *filename);

// create ISs and then read element triples, Neumann boundary segments, and boundary flags into them; call UMReadNodes() first
PetscErrorCode UMReadISs(UM *mesh, char *filename);
//ENDDECLARE


#endif