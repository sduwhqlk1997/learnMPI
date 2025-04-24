static char help[] = "Poisson Problem in 2d and 3d with finite elements.\n\
We solve the Poisson problem in a rectangular\n\
domain, using a parallel unstructured mesh (DMPLEX) to discretize it.\n\
This example supports automatic convergence estimation\n\
and eventually adaptivity.\n\n\n";

#include <petscdmplex.h>
#include <petscdmceed.h>
#include <petscsnes.h>
#include <petscds.h>
#include <petscconvest.h>

typedef struct{
    /* Domain and mesh definnition*/
    PetscBool spectral; /* Look at the spectrum along planes in the solution */
    PetscBool shear; /* Shear the domain */
    PetscBool adjoint; /* Solve the adjoint problem */
    PetscBool homogeneous; /* Use homogeneous boundary conditions */
    PetscBool viewError; /* Output the solution error */
} AppCtx;

static PetscErrorCode zero(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx){
    *u = 0.0;
    return PETSC_SUCCESS;
}

static PetscErrorCode trig_inhomogeneous_u(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx){
    PetscInt d;
    *u = 0.0;
    for (d = 0; d < dim; ++d){
        *u += PetscSinReal(2.0 * PETSC_PI *x[d]);
    }
    return PETSC_SUCCESS;
}

static PetscErrorCode trig_homogeneous_u(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx){
  PetscInt d;
  *u = 1.0;
  for (d = 0; d < dim; ++d) *u *= PetscSinReal(2.0 * PETSC_PI * x[d]);
  return PETSC_SUCCESS;
}

/* Compute integral of (residual of solution)*(adjoint solution - projection of adjoint solution) */
static void obj_error_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar obj[])
{
  obj[0] = a[aOff[0]] * (u[0] - a[aOff[1]]);
}

static void f0_trig_inhomogeneous_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) f0[0] += -4.0 * PetscSqr(PETSC_PI) * PetscSinReal(2.0 * PETSC_PI * x[d]);
}

static void f0_trig_homogeneous_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) {
    PetscScalar v = 1.;
    for (PetscInt e = 0; e < dim; e++) {
      if (e == d) {
        v *= -4.0 * PetscSqr(PETSC_PI) * PetscSinReal(2.0 * PETSC_PI * x[d]);
      } else {
        v *= PetscSinReal(2.0 * PETSC_PI * x[d]);
      }
    }
    f0[0] += v;
  }
}

static void f0_unity_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = 1.0;
}

static void f0_identityaux_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = a[0];
}

static void f1_u(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) f1[d] = u_x[d];
}

static void g3_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[], PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d * dim + d] = 1.0;
}

PLEXFE_QFUNCTION(Laplace, f0_trig_inhomogeneous_u, f1_u)

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options){
    PetscFunctionBeginUser;
    options->shear       = PETSC_FALSE;
    options->spectral    = PETSC_FALSE;
    options->adjoint     = PETSC_FALSE;
    options->homogeneous = PETSC_FALSE;
    options->viewError   = PETSC_FALSE;

    PetscOptionsBegin(comm, "", "Poisson Problem Options", "DMPLEX");
    PetscCall(PetscOptionsBool("-shear", "Shear the domain", "ex13.c", options->shear, &options->shear, NULL));
    PetscCall(PetscOptionsBool("-spectral", "Look at the spectrum along planes of the solution", "ex13.c", options->spectral, &options->spectral, NULL));
    PetscCall(PetscOptionsBool("-adjoint", "Solve the adjoint problem", "ex13.c", options->adjoint, &options->adjoint, NULL));
    PetscCall(PetscOptionsBool("-homogeneous", "Use homogeneous boundary conditions", "ex13.c", options->homogeneous, &options->homogeneous, NULL));
    PetscCall(PetscOptionsBool("-error_view", "Output the solution error", "ex13.c", options->viewError, &options->viewError, NULL));
    PetscOptionsEnd();
    PetscFunctionReturn(PETSC_SUCCESS);

}

static PetscErrorCode CreateSpectralPlanes(DM dm, PetscInt numPlanes, const PetscInt planeDir[], const PetscReal planeCoord[], AppCtx *user)
{
  PetscSection       coordSection;
  Vec                coordinates;
  const PetscScalar *coords;
  PetscInt           dim, p, vStart, vEnd, v;

  PetscFunctionBeginUser;
  PetscCall(DMGetCoordinateDim(dm, &dim));
  PetscCall(DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd));
  PetscCall(DMGetCoordinatesLocal(dm, &coordinates));
  PetscCall(DMGetCoordinateSection(dm, &coordSection));
  PetscCall(VecGetArrayRead(coordinates, &coords));
  for (p = 0; p < numPlanes; ++p) {
    DMLabel label;
    char    name[PETSC_MAX_PATH_LEN];

    PetscCall(PetscSNPrintf(name, PETSC_MAX_PATH_LEN, "spectral_plane_%" PetscInt_FMT, p));
    PetscCall(DMCreateLabel(dm, name));
    PetscCall(DMGetLabel(dm, name, &label));
    PetscCall(DMLabelAddStratum(label, 1));
    for (v = vStart; v < vEnd; ++v) {
      PetscInt off;

      PetscCall(PetscSectionGetOffset(coordSection, v, &off));
      if (PetscAbsReal(planeCoord[p] - PetscRealPart(coords[off + planeDir[p]])) < PETSC_SMALL) PetscCall(DMLabelSetValue(label, v, 1));
    }
  }
  PetscCall(VecRestoreArrayRead(coordinates, &coords));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode CreateMesh(MPI_Comm comm, AppCtx *user, DM *dm){
    PetscFunctionBeginUser; // 标记用户自定义函数
    PetscCall(DMCreate(comm, dm));
    PetscCall(DMSetType(*dm, DMPLEX));
    PetscCall(DMSetFromOptions(*dm));
    if (user->shear) 
        PetscCall(DMPlexShearGeometry(*dm, DM_X, NULL));
    PetscCall(DMSetApplicationContext(*dm, user));
    PetscCall(DMViewFromOptions(*dm, NULL, "-dm_view"));
    if (user->spectral){
        PetscInt planeDir[2] = {0, 1};
        PetscReal planeCoord[2] = {0., 1.};

        PetscCall(CreateSpectralPlanes(*dm, 2, planeDir, planeCoord, user));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

