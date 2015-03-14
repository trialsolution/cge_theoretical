
$Title Krugman Trade Equilibrium with Iceberg Costs

*Edward J. Balistreri, Colorado School of Mines (ebalistr@mines.edu)
*Thomas F. Rutherford, ETH Z\"{u}rich (tom@mpsge.org).
*March 2011

Set
        r        countries or regions /R1,R2,R3/
        k        goods                      /G1/;
Alias (r,s);

Parameters
        sig          elasticity of substitution       /5.6/,
        eta          demand elasticity                /2/,
        mu           supply elasticity                /0.5/,
        Q0(k,r)      benchmark aggregate quantity,
        P0(k,r)      benchmark price index,
        N0(k,r)      benchmark number of firms,
        qf0(k,r,s)   benchmark firm-level quantity,
        pf0(k,r,s)   benchmark firm-level pricing (gross of tau),
        c0(k,r)      benchmark input cost,
        Y0(k,r)      benchmark input supply,
        fc(k,r)      fixed costs,
        tau(k,r,s)   iceberg transport cost factor,
        vx0(k,r,s)   arbitrary benchmark export values
        ;
c0(k,r)    = 1;
vx0(k,r,s) = 1;
vx0(k,r,r) = 3;
Y0(k,r)    = sum(s, vx0(k,r,s))/c0(k,r);
N0(k,r)    = 10;
*        Calibrate the fixed cost based on zero profit
fc(k,r)    = sum(s,vx0(k,r,s))/(sig*N0(k,r)*c0(k,r));
P0(k,r)    = 1;
Q0(k,r)    = sum(s, vx0(k,s,r))/P0(k,r);
pf0(k,r,s) = (vx0(k,r,s)/(N0(k,r)*Q0(k,s)))**(1/(1-sig));
qf0(k,r,s) = Q0(k,s)*pf0(k,r,s)**(-sig);
tau(k,r,s) = (1-1/sig)*pf0(k,r,s)/c0(k,r);
display tau;

Positive Variables
        Q(k,r)        Composite Quantity,
        P(k,r)        Composite price index,
        N(k,r)        Number of firms (varieties)
        QF(k,r,s)     Firm-level output in s-market
        PF(k,r,s)     Firm-level (gross) pricing in s-market
        c(k,r)        Composite input price (marginal cost),
        Y(k,r)        Composite input supply (output);

Equations
        DEM(k,r)      Aggregate demand,
        DS(k,r)       Dixit-Stiglitz price index,
        FE(k,r)       Free entry,
        DEMF(k,r,s)   Firm demand,
        MKUP(k,r,s)   Optimal firm pricing,
        MKT(k,r)      Input market clearance,
        SUP(k,r)      Input supply (output);

DEM(k,r)..  Q(k,r) - Q0(k,r)*(P0(k,r)/P(k,r))**eta =g= 0;

DS(k,s).. sum(r,N(k,r)*PF(k,r,s)**(1-sig))**(1/(1-sig)) -
           P(k,s) =g= 0;

FE(k,r)..     c(k,r)*fc(k,r) - sum(s,PF(k,r,s)*QF(k,r,s)/sig) =g= 0;

DEMF(k,r,s).. QF(k,r,s) -  Q(k,s)*(P(k,s)/PF(k,r,s))**sig =g= 0;

MKUP(k,r,s).. tau(k,r,s)*c(k,r) - (1 - 1/sig)*PF(k,r,s) =g= 0;

MKT(k,r).. Y(k,r) -
           N(k,r)*(fc(k,r) + sum(s,tau(k,r,s)*QF(k,r,s)))
           =g= 0;

SUP(k,r).. Y0(k,r)*(c(k,r)/c0(k,r))**mu - Y(k,r) =g= 0;

model A_2 /DEM.P,DS.Q,FE.N,DEMF.PF,MKUP.QF,MKT.c,SUP.Y/;

*Set the level values and check for benchmark consistency
Q.l(k,r)   =Q0(k,r)   ;
P.l(k,r)   =P0(k,r)   ;
N.l(k,r)   =N0(k,r)   ;
QF.l(k,r,s)=QF0(k,r,s);
PF.l(k,r,s)=PF0(k,r,s);
c.l(k,r)   =c0(k,r)   ;
Y.l(k,r)   =Y0(k,r)   ;

A_2.iterlim = 0;
Solve A_2 using MCP;
Abort$(A_2.objval > 1e-6) "Benchmark Replication Failed";

