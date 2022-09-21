function [RHS] = applyimpaint(RHS,Value,start_y,end_y,start_x,end_x)
Indy = start_y:end_y;
Indx = start_x:end_x;
if (size(Value,1)==range(Indy)+1 && size(Value,2)==range(Indx)+1)
RHS(Indy,Indx) = RHS(Indy,Indx) + Value;
else
   error('Dimensionn is not matched!'); 
end
end