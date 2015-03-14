
$Title Armington Trade Equilibrium with Iceberg Costs

*Edward J. Balistreri, Colorado School of Mines (ebalistr@mines.edu)
*Thomas F. Rutherford, ETH Z\"{u}rich (tom@mpsge.org).
*March 2011

Set
        r        countries or regions /R1,R2,R3/
        j        goods                     /G1/;
Alias (r,s);

Parameters
        sig         elasticity of substitution       /5.6/,
        eta         demand elasticity                /1/,
        mu          supply elasticity                /0/,
        Q0(j,r)     benchmark aggregate quantity,
        P0(j,r)     benchmark price index,
        c0(j,r)     benchmark input cost,
        Y0(j,r)     benchmark input supply,
        tau(j,r,s)  iceberg transport cost factor,
        vx0(j,r,s)  arbitrary benchmark export values,
        lambda(j,r,s) bilateral preference weights
        ;
P0(j,r)    = 1;
c0(j,r)    = 1;
vx0(j,r,s) = 1;
vx0(j,r,r) = 3;
Q0(j,r)    = sum(s, vx0(j,s,r))/P0(j,r);
Y0(j,r)    = sum(s, vx0(j,r,s))/c0(j,r);

*        Assume neutral preference weights and calibrate tau
lambda(j,r,s)= 1;
tau(j,r,s) = (vx0(j,r,s)/(c0(j,r)*Q0(j,s)))**(1/(1-sig)) *
             (lambda(j,r,s)*P0(j,s)/(c0(j,r)))**(sig/(sig-1));

*        Alternatively we could specify tau and calibrate lambda
*tau(j,r,s)  = 1;
*lambda(j,r,s) = c0(j,r)/P0(j,s) *
*             (vx0(j,r,s)/(c0(j,r)*Q0(j,s)))**(1/sig) *
*              tau(j,r,s)**((sig-1)/sig);
Display lambda,tau;

Positive Variables
        Q(j,r)        Composite Quantity,
        P(j,r)        Composite price index,
        c(j,r)        Composite input price (marginal cost),
        Y(j,r)        Composite input supply (output);

Equations
        DEM(j,r)        Aggregate demand,
        ARM(j,r)        Armington unit cost function,
        MKT(j,r)        Market clearance,
        SUP(j,r)        Input supply (output);

DEM(j,r)..  Q(j,r) - Q0(j,r)*(P0(j,r)/P(j,r))**eta =g= 0;

ARM(j,s).. sum(r,lambda(j,r,s)**(sig) *(tau(j,r,s)*c(j,r))**(1-sig)
              )**(1/(1-sig)) -
           P(j,s) =g= 0;

MKT(j,r).. Y(j,r) -
           sum(s,tau(j,r,s)*Q(j,s)*
                 (lambda(j,r,s)*P(j,s)/(tau(j,r,s)*c(j,r)))**(sig)
              ) =g= 0;

SUP(j,r).. Y0(j,r)*(c(j,r)/c0(j,r))**mu - Y(j,r) =g= 0;

model A_1 /DEM.P,ARM.Q,MKT.c,SUP.Y/;

*Set the level values and check for benchmark consistency
Q.l(j,r) = Q0(j,r);
P.l(j,r) = P0(j,r);
c.l(j,r) = c0(j,r);
Y.l(j,r) = Y0(j,r);

A_1.iterlim = 0;
Solve A_1 using MCP;
Abort$(A_1.objval > 1e-6) "Benchmark Replication Failed";


