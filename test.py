    c-----------------------------------------------------------------------
    c\BeginDoc
    c
    c\Name: dsaupd 
    c
    c\Description: 
    ...
    c  dsaupd  is usually called iteratively to solve one of the 
    c  following problems:
    c
    c  Mode 1:  A*x = lambda*x, A symmetric 
    c           ===> OP = A  and  B = I.
    ...
    c\Usage:
    c  call dsaupd  
    c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
    c       IPNTR, WORKD, WORKL, LWORKL, INFO )
    c
    c\Arguments
    c  IDO     Integer.  (INPUT/OUTPUT)
    ...             
    c  BMAT    Character*1.  (INPUT)
    ...
    c  N       Integer.  (INPUT)
    ...
    c  WHICH   Character*2.  (INPUT)
    ...
    c  NEV     Integer.  (INPUT)
    ...
    c  TOL     Double precision  scalar.  (INPUT)
    ...
    c  RESID   Double precision  array of length N.  (INPUT/OUTPUT)
    ...
    c  NCV     Integer.  (INPUT)
    ...
    c  V       Double precision  N by NCV array.  (OUTPUT)
    ...
    c  LDV     Integer.  (INPUT)
    ...
    c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
    ...
    c  IPNTR   Integer array of length 11.  (OUTPUT)         -------------------------------------------------------------
    ...
    c  WORKD   Double precision  work array of length 3*N.  (REVERSE COMMUNICATION)
    ...
    c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
    ...
    c  LWORKL  Integer.  (INPUT)
    ...
    c  INFO    Integer.  (INPUT/OUTPUT)
    ...
    c-----------------------------------------------------------------------
    c
          subroutine dsaupd 
         &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
         &     ipntr, workd, workl, lworkl, info )
    c
    c     %----------------------------------------------------%
    c     | Include files for debugging and timing information |
    c     %----------------------------------------------------%
    c
          include   'debug.h'
          include   'stat.h'
    c
    c     %------------------%
    c     | Scalar Arguments |
    c     %------------------%
    c
          character  bmat*1, which*2
          integer    ido, info, ldv, lworkl, n, ncv, nev
          Double precision 
         &           tol
    c
    c     %-----------------%
    c     | Array Arguments |
    c     %-----------------%
    c
          integer    iparam(11), ipntr(11)
          Double precision 
         &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
    c
    c     %------------%
    c     | Parameters |
    c     %------------%
    c
          Double precision 
         &           one, zero
          parameter (one = 1.0D+0 , zero = 0.0D+0 )
    c
    c     %---------------%
    c     | Local Scalars |
    c     %---------------%
    c
          integer    bounds, ierr, ih, iq, ishift, iupd, iw, 
         &           ldh, ldq, msglvl, mxiter, mode, nb,
         &           nev0, next, np, ritz, j
          save       bounds, ierr, ih, iq, ishift, iupd, iw,
         &           ldh, ldq, msglvl, mxiter, mode, nb,
         &           nev0, next, np, ritz
    


    #include <iostream>
    #include <bitset>
    #include <string>
    #include <vector>
    using namespace std;

    ...

    extern"C" {
    void dsaupd_(int IDO, char BMAT, int N, char WHICH[2], int NEV,
            double TOL, double * RESID, int NCV, double * V, int LDV,
            int IPARAM[11], int IPNTR[11], double * WORKD,
            double * WORKL, int LWORKL, int INFO );
    }

    int main(int argc, char *argv[]) {
        ...

        // variables to try whether dsaupd can be called
        int i;
        double d;
        char c;
        int * ip;
        double * dp;
        char * cp;
        dsaupd_(i,c,i,cp,i,d,dp,i,dp,i,ip,ip,dp,dp,i,i);

        ...
    }

    ...
















    c-----------------------------------------------------------------------
    c\BeginDoc
    c
    c\Name: dsaupd 
    c
    c\Description: 
    ...
    c  dsaupd  is usually called iteratively to solve one of the 
    c  following problems:
    c
    c  Mode 1:  A*x = lambda*x, A symmetric 
    c           ===> OP = A  and  B = I.
    ...
    c\Usage:
    c  call dsaupd  
    c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
    c       IPNTR, WORKD, WORKL, LWORKL, INFO )
    c
    c\Arguments
    c  IDO     Integer.  (INPUT/OUTPUT)
    ...             
    c  BMAT    Character*1.  (INPUT)
    ...
    c  N       Integer.  (INPUT)
    ...
    c  WHICH   Character*2.  (INPUT)
    ...
    c  NEV     Integer.  (INPUT)
    ...
    c  TOL     Double precision  scalar.  (INPUT)
    ...
    c  RESID   Double precision  array of length N.  (INPUT/OUTPUT)
    ...
    c  NCV     Integer.  (INPUT)
    ...
    c  V       Double precision  N by NCV array.  (OUTPUT)
    ...
    c  LDV     Integer.  (INPUT)
    ...
    c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
    ...
    c  IPNTR   Integer array of length 11.  (OUTPUT)         -------------------------------------------------------------
    ...
    c  WORKD   Double precision  work array of length 3*N.  (REVERSE COMMUNICATION)
    ...
    c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
    ...
    c  LWORKL  Integer.  (INPUT)
    ...
    c  INFO    Integer.  (INPUT/OUTPUT)
    ...
    c-----------------------------------------------------------------------
    c
          subroutine dsaupd 
         &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
         &     ipntr, workd, workl, lworkl, info )
    c
    c     %----------------------------------------------------%
    c     | Include files for debugging and timing information |
    c     %----------------------------------------------------%
    c
          include   'debug.h'
          include   'stat.h'
    c
    c     %------------------%
    c     | Scalar Arguments |
    c     %------------------%
    c
          character  bmat*1, which*2
          integer    ido, info, ldv, lworkl, n, ncv, nev
          Double precision 
         &           tol
    c
    c     %-----------------%
    c     | Array Arguments |
    c     %-----------------%
    c
          integer    iparam(11), ipntr(11)
          Double precision 
         &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
    c
    c     %------------%
    c     | Parameters |
    c     %------------%
    c
          Double precision 
         &           one, zero
          parameter (one = 1.0D+0 , zero = 0.0D+0 )
    c
    c     %---------------%
    c     | Local Scalars |
    c     %---------------%
    c
          integer    bounds, ierr, ih, iq, ishift, iupd, iw, 
         &           ldh, ldq, msglvl, mxiter, mode, nb,
         &           nev0, next, np, ritz, j
          save       bounds, ierr, ih, iq, ishift, iupd, iw,
         &           ldh, ldq, msglvl, mxiter, mode, nb,
         &           nev0, next, np, ritz





        int ido = 0;           // Reverse communication flag. (in/out-)
        char bmat = 'I';        // Normal eigenvalue equation. (in)
        int n = 1<<L;           // Dimension of the eigenproblem. (in)
        char which[3] = "SM";   // Compute smallest (in magnitude) eigenvalues. (in)
        int nev = 4;            // Number of eigenvalues to be computed. (in)
        double tol = 0;         // Tolerated error on the eigenvalues. (in)

        double resid[n];        // Will contain the final residual vector. (in/out)
        int ncv = min(n, max(2*nev + 1, 20));   // Number of Lanczos vectors. (in)
        double v[n*ncv];        // Will contain the Lanczos vectors. (out)
        int ldv = n;            //???? Leading dimension of v.??? (in)
        int iparam[11];         // Itaration parameters: (in)
            iparam[0] = 1;          // Exact shifts are used.
            iparam[2] = 300;        // Maximum number of Arnoldi iterations allowed.
            iparam[6] = 1;          // Mode 1 of dsaupd (normal eigenvalue equation).
        int ipntr[11];          // Pointers inside working spaces. (work)
        double workd[3*n];      // Workspace that will contain residuals. (out/work)
        int lworkl = ncv*(ncv+8);   // Size of private workspace. (in/work)
        double workl[lworkl];       // Private workspace. (work)
        int info = 0;           // Randomly initial residual vector is used. (in/out)

        dsaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                iparam, ipntr, workd, workl, lworkl, info );


        int ido = 0;            // Reverse communication flag. (in/out-)
        char bmat = 'I';        // Normal eigenvalue equation. (in)
        int n = 1<<L;           // Dimension of the eigenproblem. (in)
        char which[3] = "SM";   // Compute smallest (in magnitude) eigenvalues. (in)
        int nev = 4;            // Number of eigenvalues to be computed. (in)
        double tol = 0;         // Tolerated error on the eigenvalues. (in)

        double resid[n];        // Will contain the final residual vector. (in/out)
        int ncv = min(n, max(2*nev + 1, 20));   // Number of Lanczos vectors. (in)
        double v[n*ncv];        // Will contain the Lanczos vectors. (out)
        int ldv = n;            //???? Leading dimension of v.??? (in)
        int iparam[11];         // Itaration parameters: (in)
            iparam[0] = 1;          // Exact shifts are used.
            iparam[2] = 300;        // Maximum number of Arnoldi iterations allowed.
            iparam[6] = 1;          // Mode 1 of dsaupd (normal eigenvalue equation).
        int ipntr[11];          // Pointers inside working spaces. (work)
        double workd[3*n];      // Workspace that will contain residuals. (out/work)
        int lworkl = ncv*(ncv+8);   // Size of private workspace. (in/work)
        double workl[lworkl];       // Private workspace. (work)
        int info = 0;           // Randomly initial residual vector is used. (in/out)

        dsaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                iparam, ipntr, workd, workl, lworkl, info );