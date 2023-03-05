
function c=testarg1(a,b)
 
    if(nargin==1)
 
        c=a.^2;
 
    elseif(nargin==2)
 
        c=a*b;
 
    end
end
