
  in"bin/../utile/sum_22.red";

  in"results/symb1.red";
  sum1:=sum;
  in"results_/symb1.red";

  in"sub.red";
  sum1:=sum1;
  sum:=sum;

  diff:=sum1-sum;
  if( not(diff =0) ) then  diff:= (diff where  
                            tp(~m)=>1/(t-m^2),
                            sp(~m)=>1/(s-m^2),
                            up(~m)=>1/(-s-t+p1.p1+p2.p2+p3.p3+p4.p4 -m^2));

   
  out "message"$
  if( diff=0 )     then  write "  OK"  
                   else  write "  ERROR";
  shut "message";
  out T;

end$
