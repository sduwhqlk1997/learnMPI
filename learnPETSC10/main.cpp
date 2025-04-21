static char help[] =
    "Solves the p-Helmholtz equation in 2D using Q_1 FEM.  Option prefix -ph_.\n"
    "Problem is posed as minimizing this objective functional over W^{1,p}\n"
    "for p>1:\n"
    "    I[u] = int_Omega (1/p) |grad u|^p + (1/2) u^2 - f u.\n"
    "The strong form equation, namely setting the gradient to zero, is a PDE\n"
    "    - div( |grad u|^{p-2} grad u ) + u = f\n"
    "subject to homogeneous Neumann boundary conditions.  Implements objective\n"
    "and gradient (residual) but no Hessian (Jacobian).  Defaults to linear\n"
    "problem (p=2) and quadrature degree 2.  Can be run with only an objective\n"
    "function; use -ph_no_gradient -snes_fd_function.\n\n";

#include <petsc.h>
#include "quadrature.h"

typedef struct
{
    PetscReal p, eps;
    PetscInt quadpts;
    PetscReal (*f)(PetscReal x, PetscReal y, PetscReal p, PetscReal eps);
} PHelmCtx;

static PetscReal f_constant(PetscReal x, PetscReal y, PetscReal p, PetscReal eps)
{
    return 1.0;
}

static PetscReal u_exact_cosines(PetscReal x, PetscReal y, PetscReal p, PetscReal eps)
{
    return PetscCosReal(PETSC_PI * x) * PetscCosReal(PETSC_PI * y);
}

static PetscReal f_cosines(PetscReal x, PetscReal y, PetscReal p, PetscReal eps)
{
    const PetscReal uu = u_exact_cosines(x, y, p, eps),
                    pi2 = PETSC_PI * PETSC_PI,
                    lapu = -2 * pi2 * uu;
    if (p == 2.0)
    {
        return -lapu + uu;
    }
    else
    {
        const PetscReal
            ux = -PETSC_PI * PetscSinReal(PETSC_PI * x) * PetscCosReal(PETSC_PI * y),
            uy = -PETSC_PI * PetscCosReal(PETSC_PI * x) * PetscSinReal(PETSC_PI * y),
            // note regularization changes f(x,y) but not u(x,y):
            w = ux * ux + uy * uy + eps * eps,
            pi3 = pi2 * PETSC_PI,
            wx = pi3 * PetscSinReal(2 * PETSC_PI * x) * PetscCosReal(2 * PETSC_PI * y),
            wy = pi3 * PetscCosReal(2 * PETSC_PI * x) * PetscSinReal(2 * PETSC_PI * y);
        const PetscReal s = (p - 2) / 2; //  -1/2 <= s <= 0
        return -s * PetscPowScalar(w, s - 1) * (wx * ux + wy * uy) - PetscPowScalar(w, s) * lapu + uu;
    }
}

typedef enum
{
    CONSTANT,
    COSINES
} ProblemType;
static const char *ProblemTypes[] = {"constant", "cosines",
                                     "ProblemType", "", NULL};

extern PetscErrorCode GetVecFromFunction(DMDALocalInfo *, Vec,
                                         PetscReal (*)(PetscReal, PetscReal, PetscReal, PetscReal), PHelmCtx *);
extern PetscErrorCode FormObjectiveLocal(DMDALocalInfo *, PetscReal **, PetscReal *, PHelmCtx *);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo *, PetscReal **, PetscReal **, PHelmCtx *);

int main(int argc, char *argv[])
{
    DM da;
    SNES snes;
    Vec u_initial, u, u_exact;
    PHelmCtx user;
    DMDALocalInfo info;
    ProblemType problem = COSINES;
    PetscBool no_objective = PETSC_FALSE,
              no_gradient = PETSC_FALSE,
              exact_init = PETSC_FALSE,
              view_f = PETSC_FALSE;
    PetscReal err;

    PetscCall(PetscInitialize(&argc, &argv, NULL, help));

    user.p = 2.0;
    user.eps = 0.0;
    user.quadpts = 2;
    PetscOptionsBegin(PETSC_COMM_WORLD, "ph_",
                      "p-Helmholtz solver options", "");
    PetscCall(PetscOptionsReal("-eps",
                               "regularization parameter eps",
                               "main.cpp", user.eps, &(user.eps), NULL));
    PetscCall(PetscOptionsBool("-exact_init",
                               "use exact solution to initialize",
                               "main.cpp", exact_init, &(exact_init), NULL));
    PetscCall(PetscOptionsBool("-no_objective",
                               "do not set the objective evaluation function",
                               "main.cpp", no_objective, &(no_objective), NULL));
    PetscCall(PetscOptionsBool("-no_gradient",
                               "do not set the residual evaluation function",
                               "main.cpp", no_gradient, &(no_gradient), NULL));
    PetscCall(PetscOptionsReal("-p",
                               "exponent p > 1",
                               "main.cpp", user.p, &(user.p), NULL));
    if (user.p < 1.0)
    {
        SETERRQ(PETSC_COMM_SELF, 1, "p >= 1 required");
    }
    if (user.p == 1.0)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                              "WARNING: well-posedness only known for p > 1\n"));
    }
    PetscCall(PetscOptionsEnum("-problem",
                               "problem type determines right side f(x,y)",
                               "main.cpp", ProblemTypes, (PetscEnum)problem, (PetscEnum *)&problem,
                               NULL));
    PetscCall(PetscOptionsInt("-quadpts",
                              "number n of quadrature points in each direction (= 1,2,3 only)",
                              "main.cpp", user.quadpts, &(user.quadpts), NULL));
    if ((user.quadpts < 1) || (user.quadpts > 3))
    {
        SETERRQ(PETSC_COMM_SELF, 3, "quadrature points n=1,2,3 only");
    }
    PetscCall(PetscOptionsBool("-view_f",
                               "view right-hand side to STDOUT",
                               "main.cpp", view_f, &(view_f), NULL));
    PetscOptionsEnd();

    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, 
        DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,2,2,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da));
    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da));
    PetscCall(DMSetApplicationContext(da,&user));
    PetscCall(DMDASetUniformCoordinates(da,0.0,1.0,0.0,1.0,-1.0,-1.0));
    PetscCall(DMDAGetLocalInfo(da,&info));

    PetscCall(SNESCreate(PETSC_COMM_WORLD,&snes));
    PetscCall(SNESSetDM(snes,da));
    if (!no_objective) {
        PetscCall(DMDASNESSetObjectiveLocal(da, (DMDASNESObjectiveFn *)FormObjectiveLocal, &user));
    }
    if (no_gradient) {
        PetscCall(PetscOptionsSetValue(NULL,"snes_fd_function_eps","0,0"));
    } else {
        PetscCall(DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunctionFn *)FormFunctionLocal,&user));
    }
    PetscCall(SNESSetFromOptions(snes));

    // set initial iterate and right-hand side
    // (示例代码145行)
    return 0;
}