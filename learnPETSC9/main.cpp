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
static DMDASNESFunctionFn* residual_ptr[3]={
    (DMDASNESFunctionFn*)&Poisson1DFunctionLocal,
    (DMDASNESFunctionFn*)&Poisson2DFunctionLocal,
    (DMDASNESFunctionFn*)&Poisson1DFunctionLocal
};

static DMDASNESJacobianFn* jacobian_ptr[3]={
    (DMDASNESJacobianFn*)&Poisson1DJacobianLocal,
    (DMDASNESJacobianFn*)&Poisson2DJacobianLocal,
    (DMDASNESJacobianFn*)&Poisson3DJacobianLocal
};

typedef PetscErrorCode (*ExactFcnVec)(DMDALocalInfo*, Vec, PoissonCtx*);

static ExactFcnVec getuexact_ptr[3]={&Form1DUExact,&Form2DUExact,&Form3DUExact};
//ENDPTRARRAYS

typedef enum {MANUPOLY, MANUEXP, ZERO} ProblemType;
static const char* ProblemTypes[] = {"manupoly","manuexp","zero",
                                    "ProblemType", "", NULL};
// more arrays of pointers to functions:   ..._ptr[DIMS][PROBLEMS]
typedef PetscReal (*PointwiseFcn)(PetscReal,PetscReal,PetscReal,void*);

static PointwiseFcn g_bdry_ptr[3][3] = {
    {&u_exact_1Dmanupoly, &u_exact_1Dmanuexp, &zero},
    {&u_exact_2Dmanupoly, &u_exact_2Dmanuexp, &zero},
    {&u_exact_3Dmanupoly, &u_exact_3Dmanuexp, &zero}
};

static PointwiseFcn f_rhs_ptr[3][3] = {
    {&f_rhs_1Dmanupoly, &f_rhs_1Dmanuexp, &zero},
    {&f_rhs_2Dmanupoly, &f_rhs_2Dmanuexp, &zero},
    {&f_rhs_3Dmanupoly, &f_rhs_3Dmanuexp, &zero}
};

static const char* InitialTypes[] = {
    "zeros","random","InitialType", "", NULL
};

int main(int argc, char* argv[]){
    DM da, da_after;
    SNES snes;
    KSP ksp;
    Vec u_initial, u, u_exact;
    PoissonCtx user;
    DMDALocalInfo info;
    PetscReal errinf, normconst2h, err2h;
    char           gridstr[99];
    ExactFcnVec    getuexact;

    // fish defaults:
    PetscInt       dim = 2;                  // 2D
    ProblemType    problem = MANUEXP;        // manufactured problem using exp()
    InitialType    initial = ZEROS;          // set u=0 for initial iterate
    PetscBool      gonboundary = PETSC_TRUE; // initial iterate has u=g on boundary

    PetscCall(PetscInitialize(&argc,&argv,NULL,help));

    // get options and configure context
    user.Lx = 1.0;
    user.Ly = 1.0;
    user.Lz = 1.0;
    user.cx = 1.0;
    user.cy = 1.0;
    user.cz = 1.0;
    PetscOptionsBegin(PETSC_COMM_WORLD,"fsh_","options for main.cpp","");
    PetscCall(PetscOptionsReal("-cx",
        "set coefficient of x term u_xx in equation",
        "main.cpp",user.cx,&user.cx,NULL));
   PetscCall(PetscOptionsReal("-cy",
        "set coefficient of y term u_yy in equation",
        "main.cpp",user.cy,&user.cy,NULL));
   PetscCall(PetscOptionsReal("-cz",
        "set coefficient of z term u_zz in equation",
        "main.cpp",user.cz,&user.cz,NULL));
    PetscCall(PetscOptionsInt("-dim",
        "dimension of problem (=1,2,3 only)",
        "main.cpp",dim,&dim,NULL));
    PetscCall(PetscOptionsBool("-initial_gonboundary",
        "set initial iterate to have correct boundary values",
        "main.cpp",gonboundary,&gonboundary,NULL));
    PetscCall(PetscOptionsEnum("-initial_type",
        "type of initial iterate",
        "main.cpp",InitialTypes,(PetscEnum)initial,(PetscEnum*)&initial,NULL));
    PetscCall(PetscOptionsReal("-Lx",
        "set Lx in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
        "main.cpp",user.Lx,&user.Lx,NULL));
    PetscCall(PetscOptionsReal("-Ly",
        "set Ly in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
        "main.cpp",user.Ly,&user.Ly,NULL));
    PetscCall(PetscOptionsReal("-Lz",
        "set Ly in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
        "main.cpp",user.Lz,&user.Lz,NULL));
    PetscCall(PetscOptionsEnum("-problem",
        "problem type; determines exact solution and RHS",
        "main.cpp",ProblemTypes,(PetscEnum)problem,(PetscEnum*)&problem,NULL));
    PetscOptionsEnd();
    user.g_bdry = g_bdry_ptr[dim-1][problem];
    user.f_rhs = f_rhs_ptr[dim-1][problem];

    if ( user.cx <= 0.0 || user.cy <= 0.0 || user.cz <= 0.0 ) {
        SETERRQ(PETSC_COMM_SELF,2,"positivity required for coefficients cx,cy,cz\n"); // 终止程序运行，PETSC_COMM_SELF表示当前进程
    }
    if ((problem == MANUEXP) && ( user.cx != 1.0 || user.cy != 1.0 || user.cz != 1.0)) {
        SETERRQ(PETSC_COMM_SELF,3,"cx=cy=cz=1 required for problem MANUEXP\n");
    }
    //STARTCREATE
    //create DMDA in chosen dimension
    switch (dim){
        case 1:
            PetscCall(DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,3,1,1,NULL,&da));
            break;
        case 2:
            PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,3,3,
            PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da));
            break;
        case 3:
            PetscCall(DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
            DMDA_STENCIL_STAR,3,3,3,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
            1,1,NULL,NULL,NULL,&da));
            break;
        default:
            SETERRQ(PETSC_COMM_SELF,1,"invalid dim for DMDA creation\n");
    }
    PetscCall(DMSetApplicationContext(da,&user));
    PetscCall(DMSetFromOptions(da));
    PetscCall(DMSetUp(da)); // call BEFORE SetUniformCoordinates
    PetscCall(DMDASetUniformCoordinates(da,0.0,user.Lx,0.0,user.Ly,0.0,user.Lz)); // 均匀分布的网格

    // set SNES call-backs
    PetscCall(SNESCreate(PETSC_COMM_WORLD,&snes));
    PetscCall(SNESSetDM(snes,da));
    PetscCall(DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunctionFn *)(residual_ptr[dim-1]),&user));
    PetscCall(DMDASNESSetJacobianLocal(da,(DMDASNESJacobianFn *)(jacobian_ptr[dim-1]),&user));

    // default to KSPONLY+CG because problem is linear and SPD
    PetscCall(SNESSetType(snes,SNESKSPONLY)); // 表示仅使用 KSP（Krylov Subspace Solver）求解线性系统。
    PetscCall(SNESGetKSP(snes,&ksp));
    PetscCall(KSPSetType(ksp,KSPCG));
    PetscCall(SNESSetFromOptions(snes));

    // set initial iterate and then solve
    PetscCall(DMGetGlobalVector(da,&u_initial));
    PetscCall(InitialState(da, initial, gonboundary, u_initial, &user));
    PetscCall(SNESSolve(snes,NULL,u_initial));
    //ENDCREATE

    //STARTGETSOLUTION (242 行)
    return 0;
}