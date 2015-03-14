
$Title Melitz Trade Equilibrium with Iceberg Costs

*Edward J. Balistreri, Colorado School of Mines (ebalistr@mines.edu)
*Thomas F. Rutherford, ETH Z\"{u}rich (tom@mpsge.org).
*March 2011

*Warning: This partial equilibrium illustrative model does not
*         include a logical constraint to ensure M > N.

Set
        r        countries or regions /R1,R2,R3/
        h        goods                      /G1/;
Alias (r,s);

Parameters
        sig          elasticity of substitution       /3.8/,
        eta          demand elasticity                /2/,
        mu           supply elasticity                /0.5/,
         a            Pareto shape parameter           /4.6/,
         b            Pareto lower support             /0.5/,
        Q0(h,r)      benchmark aggregate quantity,
        P0(h,r)      benchmark price index,
        M0(h,r)      benchmark number of entered firms,
        N0(h,r,s)    benchmark number of operating firms,
        qf0(h,r,s)   benchmark avg firm-level quantity,
        pf0(h,r,s)   benchmark avg firm-level pricing (gross),
         phi0(h,r,s)  benchmark avg productivity
        c0(h,r)      benchmark input cost,
        Y0(h,r)      benchmark input supply,
        fc(h,r,s)    bilateral fixed costs,
        delt_fs(h,r) annualized sunk cost,
        tau(h,r,s)   iceberg transport cost factor,
        vx0(h,r,s)   arbitrary benchmark export values
        ;


c0(h,r)    = 1;
vx0(h,r,s) = 1;
vx0(h,r,r) = 3;
Y0(h,r)    = sum(s, vx0(h,r,s))/c0(h,r);
M0(h,r)    = 10;
N0(h,r,r)  = 9;
N0(h,r,s)  = (vx0(h,r,s)/vx0(h,r,r))**2 * N0(h,r,r);
*     Calibrate the sunk cost based on free entry
delt_fs(h,r)= Y0(h,r)/M0(h,r) * (sig-1)/(a*sig);
*     Calibrate the fixed cost based on zero cutoff profit
fc(h,r,s)   = vx0(h,r,s)/(N0(h,r,s)*c0(h,r)) * (a + 1 - sig)/(a*sig);
P0(h,r)    = 1;
Q0(h,r)    = sum(s, vx0(h,s,r))/P0(h,r);
pf0(h,r,s) = (vx0(h,r,s)/(N0(h,r,s)*Q0(h,s)))**(1/(1-sig));
qf0(h,r,s) = Q0(h,s)*pf0(h,r,s)**(-sig);
phi0(h,r,s)= b * (a/(a+1-sig))**(1/(sig-1)) *
            (N0(h,r,s)/M0(h,r))**(-1/a);
tau(h,r,s) = (1-1/sig)*pf0(h,r,s)*phi0(h,r,s)/c0(h,r);
display N0,tau;

Positive Variables
        Q(h,r)        Composite Quantity,
        P(h,r)        Composite price index,
        M(h,r)        Number of Entered firms
        N(h,r,s)      Number of Operating firms (varieties)
        QF(h,r,s)     Avg Firm output in s-market
        PF(h,r,s)     Avg Firm (gross) pricing in s-market
         PHI(h,r,s)    Avg Firm productivity
        c(h,r)        Composite input price (marginal cost),
        Y(h,r)        Composite input supply (output);

Equations
        DEM(h,r)      Aggregate demand,
        DS(h,r)       Dixit-Stiglitz price index,
        FE(h,r)       Free entry,
         ZCP(h,r,s)    Zero cutoff profits
        DEMF(h,r,s)   Firm demand,
        MKUP(h,r,s)   Optimal firm pricing,
         PAR(h,r,s)    Pareto Productivity
        MKT(h,r)      Input market clearance,
        SUP(h,r)      Input supply (output);

DEM(h,r)..  Q(h,r) - Q0(h,r)*(P0(h,r)/P(h,r))**eta =g= 0;

DS(h,s).. sum(r,N(h,r,s)*PF(h,r,s)**(1-sig))**(1/(1-sig)) -
           P(h,s) =g= 0;

DEMF(h,r,s).. QF(h,r,s) -  Q(h,s)*(P(h,s)/PF(h,r,s))**sig =g= 0;

MKUP(h,r,s).. tau(h,r,s)*c(h,r)/PHI(h,r,s) -
             (1 - 1/sig)*PF(h,r,s) =g= 0;

FE(h,r)..  c(h,r)*delt_fs(h,r) -
           sum(s,(N(h,r,s)/M(h,r))*PF(h,r,s)*QF(h,r,s)*(sig-1)/(a*sig))
            =g= 0;

ZCP(h,r,s).. c(h,r)*fc(h,r,s) -
            (PF(h,r,s)*QF(h,r,s)*(a+1-sig))/(a*sig) =g= 0;

PAR(h,r,s).. PHI(h,r,s) -
             b * (a/(a+1-sig))**(1/(sig-1)) * (N(h,r,s)/M(h,r))**(-1/a)
              =g= 0;

MKT(h,r).. Y(h,r) - (
            delt_fs(h,r)*M(h,r) +
            sum(s,N(h,r,s)*(fc(h,r,s) + tau(h,r,s)*QF(h,r,s)/PHI(h,r,s)))
                     ) =g= 0;

SUP(h,r).. Y0(h,r)*(c(h,r)/c0(h,r))**mu - Y(h,r) =g= 0;

model A_3 /DEM.P,DS.Q,FE.M,ZCP.N,DEMF.PF,MKUP.QF,PAR.PHI,MKT.c,SUP.Y/;

*Set the level values and check for benchmark consistency
Q.l(h,r)   =Q0(h,r)   ;
P.l(h,r)   =P0(h,r)   ;
M.l(h,r)   =M0(h,r)   ;
N.l(h,r,s) =N0(h,r,s) ;
QF.l(h,r,s)=QF0(h,r,s);
PF.l(h,r,s)=PF0(h,r,s);
PHI.l(h,r,s)=PHI0(h,r,s);
c.l(h,r)   =c0(h,r)   ;
Y.l(h,r)   =Y0(h,r)   ;

A_3.iterlim = 0;
Solve A_3 using MCP;
Abort$(A_3.objval > 1e-6) "Benchmark Replication Failed";

