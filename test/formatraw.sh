egrep "(time)|(cp)|^@|^cplex" | tformat --sum @xfrac @ffrac @cguFeasible @cguInfeasible cplex@cpFeasible @cplex@cpInfeasible --mean cplex@cpedges | sort -n
