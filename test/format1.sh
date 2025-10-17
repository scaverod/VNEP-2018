egrep "(cpCost)|(cguCost)|(compactTime)|(cgtime)" | tformat --latex --meanstddev cplex@compactTime @cgtime | sort -n
