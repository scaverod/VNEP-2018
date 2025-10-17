model embed.mod;
data d.data
option solver cplexamp;
solve;
#display cv;
#display cs;
#display M;
#display D;
display solve_result;
display Cost;
