set VV;
set VS;

#processing power
param cv {v in VV};
param cs {s in VS};

#bandwidth
param bv {k in VV, l in VV};
param bs {n in VS, m in VS};

var M {v in VV, s in VS} binary;
var D {k in VV, l in VV, m in VS, n in VS} binary;

minimize Cost: sum {k in VV, l in VV, m in VS, n in VS} D[k, l, m, n] * bv[k,l];

subject to ValidPaths  {k in VV, l in VV, m in VS, n in VS}: 
                        D[k,l,m,n] <= bs[m,n];
subject to CPU {s in VS}: sum {v in VV} M[v, s] * cv[v] <= cs[s];
subject to VirtualOne {v in VV}: sum {s in VS} M[v,s] = 1;
subject to SubstratelOne {s in VS}: sum {v in VV} M[v,s] <= 1;

subject to Bandwidth {m in VS, n in VS}: sum {k in VV, l in VV} 
              (D[k,l,m,n] * bv[k,l] + D[k,l,n,m] * bv[k,l]) <= bs[m,n];

subject to Path {k in VV, l in VV, s in VS}: 
( sum {j in VS} D[k,l,s,j] ) - ( sum {j in VS} D[k,l,j,s] ) = M[k,s] - M[l,s];
