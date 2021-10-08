function constr = constraint(x,y,z,R,rE,N)
%---------------------------------------------------------------------------------------------

%This function checks if any catalysts in the arrangement x,y with catalyst radius rE are overlapping
%it returns 1 if there is an overlap and 0 if there is no overlap

%---------------------------------------------------------------------------------------------

flag1 = 0;
flag2 = 0;

for i = 1:N
    
    c = sqrt(x(i)^2+y(i)^2+z(i)^2)-(R-rE);
   
    if c>=0
       
        flag1 = 1;
        
        break
        
    end
    
end

if flag1 == 0

for h = 1:N-1
    
    for k = h+1:N
    
        c = 2*rE - sqrt(((x(h)-x(k))^2)+((y(h)-y(k))^2)++((z(h)-z(k))^2));
    
        if c>=0
            
            flag2 = 1;
            
            break
            
        end
        
    end
    
    if flag2 == 1
        
        break
        
    end
    
end

end

if flag1 == 1 || flag2 == 1
    
    constr = 1;
    
else
    
    constr = 0;
    
end
    
    
    
end