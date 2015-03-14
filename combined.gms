
$Title Mix and Match General Equilibrium with Iceberg Costs

*Edward J. Balistreri, Colorado School of Mines (ebalistr@mines.edu)
*Thomas F. Rutherford, ETH Z\"{u}rich (tom@mpsge.org).
*March 2011
*        This formulation allows for Armington, Krugman, or Melitz
*        trade depending on the user defined subset j(i),k(i),or h(i).

$onempty
Set
        r        countries or regions      /R1*R3/
        f        factors of production     /L1*L3/,
        i        goods                /G1, G2, G3/
         j(i)     Armingtion goods     /G1/,
        k(i)     Krugman goods            /G2/,
        h(i)     Melitz goods                 /G3/;

;
display i,j,k,h;
Alias (r,s),(f,g);

Parameters
        alpha        top level elasticity of substitution      /2.0/,
        a            Pareto shape parameter                    /4.6/,
        b            Pareto lower support                      /0.5/,
        sig_j        industry elasticity of substitution (a+1) /5.6/,
        sig_k        industry elasticity of substitution (a+1) /5.6/,
        sig_h        industry elasticity of substitution       /3.8/,
         esub          elasticity of substitution in production  /1.0/
        vx0(i,r,s)   arbitrary benchmark export values
        beta(i,r)    expenditure weights
        gamma(i,f,r) primary factor value shares
        lbar(f,r)    primary factor supply
        ra0(r)       benchmark income
        c0(i,r)      benchmark input cost,
         w0(f,r)       benchmark wage,
        y0(i,r)      benchmark input supply,
        q0(i,r)      benchmark aggregate quantity,
        p0(i,r)      benchmark price index,
        m0(*,r)      benchmark number of entered firms,
        n0(*,r,s)    benchmark number of operating firms,
        qf0(i,r,s)   benchmark avg firm-level quantity,
        pf0(i,r,s)   benchmark avg firm-level pricing (gross),
        phi0(i,r,s)  benchmark avg productivity
        fc(*,r,s)    bilateral fixed costs,
        delt_fs(*,r) annualized sunk cost,
        nk0(k,r)     benchmark number of (krugman) firms
        fck(k,r)     krugman fixed costs,
         t(i,r,s)     bilateral tariffs
        tau(i,r,s)   iceberg transport cost factor,
        lambda(i,r,s)  preference weight parameter
        ;

*        Finite variance restriction
Abort$(a le (sig_h - 1))
     "Firm size distribution must have a finite variance. a > sig-1";

*        Setup the benchmark with arbitrary data
t(i,r,s)   = 0;
vx0(i,r,s) = 1;
vx0(i,r,r) = 3;
*        Unit choice
c0(i,r)    = 1;
p0(i,r)    = 1;
w0(f,r)     = 1;
y0(i,r)    = sum(s, vx0(i,r,s))/c0(i,r);
q0(i,r)    = sum(s, vx0(i,s,r))/p0(i,r);
ra0(r)     = sum(i,q0(i,r)*p0(i,r));
beta(i,r)  = p0(i,r) * (q0(i,r)/RA0(r))**(1/alpha);

*        randomly distribute the factor shares
gamma(i,f,r)=0;
loop(f$(ord(f) ne card(f)),
        gamma(i,f,r) = uniform(0,(1-sum(g,gamma(i,g,r))));
);
gamma(i,f,r)$(ord(f) eq card(f)) = 1-sum(g,gamma(i,g,r));
display gamma;

lbar(f,r)  = sum(i,gamma(i,f,r)*y0(i,r)*c0(i,r));

*---- Melitz model Calibration ----*
M0(i,r)    = 10;
N0(i,r,r)  = 9;
*     Assume that foreign operation falls off using the following
N0(i,r,s)  = (vx0(i,r,s)/vx0(i,r,r))**2 * N0(i,r,r);
*     Calibrate the sunk cost based on free entry
delt_fs(i,r)= y0(i,r)/M0(i,r) * (sig_h-1)/(a*sig_h);
*     Calibrate the fixed cost based on zero cutoff profit
fc(i,r,s)   = vx0(i,r,s)/(N0(i,r,s)*c0(i,r)) * (a + 1 - sig_h)/(a*sig_h);
pf0(i,r,s) = (vx0(i,r,s)/(N0(i,r,s)*q0(i,s)*p0(i,s)))**(1/(1-sig_h));
qf0(i,r,s) = (p0(i,s)*q0(i,s))*pf0(i,r,s)**(-sig_h);
lambda(i,r,s)= p0(i,s)**(1-sig_h);
phi0(i,r,s)= b * (a/(a+1-sig_h))**(1/(sig_h-1)) *
            (N0(i,r,s)/M0(i,r))**(-1/a);
*     Calibrated tau (used for all sectors to give us a
*     consistent simulation instrument):
tau(i,r,s) = (1-1/sig_h)*pf0(i,r,s)*phi0(i,r,s)/c0(i,r);
*-----------------------------------*

*---- Krugman model Calibration-----*
* Assume tau(k,r,s) is given by the melitz calibration
NK0(k,r)   = 10;
fcK(k,r)   = sum(s,vx0(k,r,s))/(sig_k*NK0(k,r)*c0(k,r));
pf0(k,r,s) = tau(k,r,s)*c0(k,r)/(1-1/sig_k);
qf0(k,r,s) = vx0(k,r,s)/(NK0(k,r)*pf0(k,r,s));
lambda(k,r,s)= qf0(k,r,s)/q0(k,s) * (pf0(k,r,s)/p0(k,s))**sig_k;

$ontext
*Note: We could have used the following code to calibrate tau
NK0(k,r)   = 10;
fcK(k,r)   = sum(s,vx0(k,r,s))/(sig_k*NK0(k,r)*c0(k,r));
pf0(k,r,s) = (vx0(k,r,s)/(NK0(k,r)*q0(k,s)*p0(k,s)))**(1/(1-sig_k));
qf0(k,r,s) = (p0(k,s)*q0(k,s))*pf0(k,r,s)**(-sig_k);
lambda(k,r,s)= p0(k,s)**(1-sig_k);
tau(k,r,s) = (1-1/sig_k)*pf0(k,r,s)/c0(k,r);
$offtext
*-----------------------------------*

*---- Armington model Calibration---*
* Assume tau(j,r,s) is given by the melitz calibration
lambda(j,r,s)=(vx0(j,r,s)/(c0(j,r)*q0(j,s)))**(1/sig_j) *
            (p0(j,s)/(c0(j,r)))**(-1)  *
             tau(j,r,s)**((sig_j-1)/sig_j);
*-----------------------------------*

display tau;

Positive Variables
        U(r)          Welfare
        E(r)          True cost of living index
        Q(i,r)        Composite Quantity,
        P(i,r)        Composite price index,
        M(h,r)        Number of Entered firms
        N(h,r,s)      Number of Operating firms (varieties)
        NK(k,r)       Number of Krugman firms (varieties)
        QF(*,r,s)     Avg Firm output in s-market
        PF(*,r,s)     Avg Firm (gross) pricing in s-market
        PHI(h,r,s)    Avg Firm productivity
        c(i,r)        Composite input price (marginal cost),
        Y(i,r)        Composite input supply (output)
        w(f,r)        Primary factor price
        RA(r)         Income
         pi(h,r,s)     Capacity rents;

Equations
        EXPFUN(r)     Unit expenditure function,
        DEM(i,r)      Aggregate demand,
        PRC_h(h,r)    Price index,
        PRC_k(k,r)    Price index,
        PRC_j(j,r)    Price index,
        FE(h,r)       Free entry,
        FEK(k,r)      Free entry,
        ZCP(h,r,s)    Zero cutoff profits
        DEMF(h,r,s)   Firm demand,
        DEMFK(k,r,s)  Firm demand,
        MKUP(h,r,s)   Optimal firm pricing,
        MKUPK(k,r,s)  Optimal firm pricing,
        PAR(h,r,s)    Pareto Productivity
        MKT_h(h,r)    Input market clearance,
        MKT_k(k,r)    Input market clearance,
        MKT_j(j,r)    Input market clearance,
        COST(i,r)     Unit cost functions
        LMKT(f,r)     Primary factor markets
        FINAL(r)      Final demand
        BC(r)         Budget constraint
         Capacity(h,r,s) Capacity constraint M gt N
         ;

EXPFUN(r)..
    (sum(i,beta(i,r)**alpha *P(i,r)**(1-alpha))**(1/(1-alpha)))$(alpha ne 1)
    + prod(i,(P(i,r)/p0(i,r))**beta(i,r))$(alpha eq 1)
    - E(r) =g= 0;

DEM(i,r)..
    Q(i,r) - RA0(r)*U(r)*(beta(i,r)*E(r)/P(i,r))**alpha =g= 0;

PRC_h(h,s)..
    sum(r,lambda(h,r,s)*N(h,r,s)*PF(h,r,s)**(1-sig_h))**(1/(1-sig_h)) -
    P(h,s) =g= 0;

PRC_k(k,s)..
    sum(r,lambda(k,r,s)*NK(k,r)*PF(k,r,s)**(1-sig_k))**(1/(1-sig_k)) -
    P(k,s) =g= 0;

PRC_j(j,s)..
    sum(r,lambda(j,r,s)**(sig_j) *((1+t(j,r,s))*tau(j,r,s)*c(j,r))**(1-sig_j)
       )**(1/(1-sig_j)) - P(j,s) =g= 0;

DEMF(h,r,s)..
    QF(h,r,s) -  lambda(h,r,s)*Q(h,s)*(P(h,s)/PF(h,r,s))**sig_h =g= 0;

DEMFK(k,r,s)..
    QF(k,r,s) -  lambda(k,r,s)*Q(k,s)*(P(k,s)/PF(k,r,s))**sig_k =g= 0;

MKUP(h,r,s)..
    (1+t(h,r,s))*tau(h,r,s)*c(h,r)/PHI(h,r,s) - (1 - 1/sig_h)*PF(h,r,s) =g= 0;

MKUPK(k,r,s)..
    (1+t(k,r,s))*tau(k,r,s)*c(k,r)            - (1 - 1/sig_k)*PF(k,r,s) =g= 0;

FE(h,r)..
    c(h,r)*delt_fs(h,r) -
    sum(s,(N(h,r,s)/M(h,r))*PF(h,r,s)*QF(h,r,s)*(sig_h-1)/((1+t(h,r,s))*a*sig_h)
          +pi(h,r,s))
    =g= 0;

FEK(k,r)..
    c(k,r)*fcK(k,r) -
    sum(s,PF(k,r,s)*QF(k,r,s)/((1+t(k,r,s))*sig_k))
    =g= 0;

ZCP(h,r,s)..
    c(h,r)*fc(h,r,s) + pi(h,r,s) -
    (PF(h,r,s)*QF(h,r,s)*(a+1-sig_h))/((1+t(h,r,s))*a*sig_h)
    =g= 0;

PAR(h,r,s)..
    PHI(h,r,s) * (N(h,r,s)/M(h,r))**(1/a) -
    b * (a/(a+1-sig_h))**(1/(sig_h-1))
    =g= 0;

MKT_j(j,r)..
    Y(j,r) -
    sum(s,tau(j,r,s)*Q(j,s)*
         (lambda(j,r,s)*P(j,s)/((1+t(j,r,s))*tau(j,r,s)*c(j,r)))**(sig_j)
       ) =g= 0;

MKT_k(k,r)..
    Y(k,r) -
    NK(k,r)*(fcK(k,r) + sum(s,tau(k,r,s)*QF(k,r,s))) =g= 0;

MKT_h(h,r)..
    Y(h,r) - (delt_fs(h,r)*M(h,r) +
    sum(s,N(h,r,s)*(fc(h,r,s) + tau(h,r,s)*QF(h,r,s)/PHI(h,r,s)))
             ) =g= 0;

COST(i,r)..
    c(i,r) - (
    (sum(f,gamma(i,f,r)*(w(f,r)/w0(f,r))**(1-esub))**(1/(1-esub)))$(esub ne 1) +
    prod(f,w(f,r)**gamma(i,f,r))$(esub eq 1)
    )=g=0;

LMKT(f,r)..
    lbar(f,r) -
    sum(i,gamma(i,f,r)*Y(i,r)*(c(i,r)*w0(f,r)/(c0(i,r)*w(f,r)))**esub) =g= 0;

FINAL(r)..  RA0(r)*U(r)*E(r) - RA(r) =g= 0;

BC(s)..     RA(s) =e=
*   Factor income
    sum(f,w(f,s)*lbar(f,s)) +
*   plus tariff revenue
    sum(r,
       sum(j,t(j,r,s)*c(j,r)*tau(j,r,s)*Q(j,s)*
             (lambda(j,r,s)*P(j,s)/((1+t(j,r,s))*tau(j,r,s)*c(j,r)))**(sig_j))
     + sum(k,t(k,r,s)*pf(k,r,s)*qf(k,r,s)* NK(k,r)/(1+t(k,r,s)))
     + sum(h,t(h,r,s)*pf(h,r,s)*qf(h,r,s)*N(h,r,s)/(1+t(h,r,s)))
       );

Capacity(h,r,s).. M(h,r) =g= N(h,r,s);

Model b_1  /
*          GE          Armington    Krugman      Melitz
*          ---------   ---------    ---------    --------
           expfun.U,
           DEM.P,
                       PRC_j.Q,     PRC_k.Q,     PRC_h.Q,
                                    FEK.NK,      FE.M,
                                                 ZCP.N,
                                    DEMFK.PF,    DEMF.PF,
                                    MKUPK.QF,    MKUP.QF,
                                                 PAR.PHI,
                       MKT_j.c,     MKT_k.c,     MKT_h.c,
                                                       capacity.pi,
           COST.Y,
           LMKT.w
           FINAL.E,
           BC.RA/;

*Set the level values and check for benchmark consistency
Q.l(i,r)    = q0(i,r)   ;
P.l(i,r)    = p0(i,r)   ;
M.l(h,r)    = M0(h,r)   ;
N.l(h,r,s)  = N0(h,r,s) ;
NK.l(k,r)   = NK0(k,r)  ;
QF.l(h,r,s) = QF0(h,r,s);
PF.l(h,r,s) = PF0(h,r,s);
QF.l(k,r,s) = QF0(k,r,s);
PF.l(k,r,s) = PF0(k,r,s);
PHI.l(h,r,s)= PHI0(h,r,s);
c.l(i,r)    = c0(i,r)   ;
Y.l(i,r)    = y0(i,r)   ;
w.l(f,r)    = 1         ;
U.l(r)      = 1         ;
E.l(r)      = 1         ;
RA.l(r)     = RA0(r)    ;

*Set lower bounds to avoid bad function calls:
Q.lo(i,r)    =1e-6;
P.lo(i,r)    =1e-6;
M.lo(h,r)    =1e-6;
N.lo(h,r,s)  =1e-6;
NK.lo(k,r)   =1e-6;
QF.lo(h,r,s) =1e-6;
PF.lo(h,r,s) =1e-6;
QF.lo(k,r,s) =1e-6;
PF.lo(k,r,s) =1e-6;
PHI.lo(h,r,s)=1e-6;
c.lo(i,r)    =1e-6;


b_1.iterlim = 0;
Solve b_1 using MCP;
Abort$(b_1.objval > 1e-6) "Benchmark Replication Failed";


