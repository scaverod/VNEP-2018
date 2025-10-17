model embed.mod;
data ../data/ampl/small.data;
option solver cplexamp;
option send_statuses 0;
solve;
display solve_result;
display Cost;
