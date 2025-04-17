static char help[] =
"Solves structured-grid Poisson problem in 1D, 2D, 3D.  Option prefix fsh_.\n"
"Equation is\n"
"    - cx u_xx - cy u_yy - cz u_zz = f,\n"
"subject to Dirichlet boundary conditions.  Solves three different problems\n"
"where exact solution is known.  Uses DMDA and SNES.  Equation is put in form\n"
"F(u) = - grad^2 u - f.  Call-backs fully-rediscretize for the supplied grid.\n"
"Defaults to 2D, a SNESType of KSPONLY, and a KSPType of CG.\n\n";

# include <petsc.h>
# include "poissonfunctions.h"

// exact solutions  u(x,y),  for boundary condition and error calculation

static PetscReal u_exact_1Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return x*x * (1.0 - x*x);
}

static PetscReal u_exact_2Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return x*x * (1.0 - x*x) * y*y *(y*y - 1.0);
}

static PetscReal u_exact_3Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return x*x * (1.0 - x*x) * y*y * (y*y - 1.0) * z*z * (z*z - 1.0);
}

static PetscReal u_exact_1Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return - PetscExpReal(x);
}

static PetscReal u_exact_2Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return - x * PetscExpReal(y);
}

static PetscReal u_exact_3Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return - x * PetscExpReal(y + z);
}

static PetscReal zero(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return 0.0;
}

// right-hand-side functions  f(x,y) = - laplacian u

static PetscReal f_rhs_1Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx){ // 这里 void型指针 ctx 可以指向任何类型的数据
    PoissonCtx* user = (PoissonCtx*) ctx;
    return user->cx * 12.0 * x*x - 2.0;
}

static PetscReal f_rhs_2Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    PoissonCtx* user = (PoissonCtx*)ctx;
    PetscReal   aa, bb, ddaa, ddbb;
    aa = x*x * (1.0 - x*x);
    bb = y*y * (y*y - 1.0);
    ddaa = 2.0 * (1.0 - 6.0 * x*x);
    ddbb = 2.0 * (6.0 * y*y - 1.0);
    return - (user->cx * ddaa * bb + user->cy * aa * ddbb);
}

static PetscReal f_rhs_3Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    PoissonCtx* user = (PoissonCtx*)ctx;
    PetscReal   aa, bb, cc, ddaa, ddbb, ddcc;
    aa = x*x * (1.0 - x*x);
    bb = y*y * (y*y - 1.0);
    cc = z*z * (z*z - 1.0);
    ddaa = 2.0 * (1.0 - 6.0 * x*x);
    ddbb = 2.0 * (6.0 * y*y - 1.0);
    ddcc = 2.0 * (6.0 * z*z - 1.0);
    return - (user->cx * ddaa * bb * cc + user->cy * aa * ddbb * cc + user->cz * aa * bb * ddcc);
}

static PetscReal f_rhs_1Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return PetscExpReal(x);
}

static PetscReal f_rhs_2Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return x * PetscExpReal(y);  // note  f = - (u_xx + u_yy) = - u
}

static PetscReal f_rhs_3Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return 2.0 * x * PetscExpReal(y + z);  // note  f = - laplacian u = - 2 u
}

// functions simply to put u_exact()=g_bdry() into a grid
// these are irritatingly-dimension-dependent inside ...
extern PetscErrorCode Form1DUExact(DMDALocalInfo*, Vec, PoissonCtx*);
extern PetscErrorCode Form2DUExact(DMDALocalInfo*, Vec, PoissonCtx*);
extern PetscErrorCode Form3DUExact(DMDALocalInfo*, Vec, PoissonCtx*);

//STARTPTRARRAYS
// arrays of pointers to functions
/*从参考代码第92行开始*/

int main(int argc, char* argv[]){
    DM da, da_after;
    SNES snes;
    KSP ksp;
    Vec u_initial, u, u_exact;
    return 0;
}