function current_section = get_section_property(W,gh2)
global W


torf=zeros(1,max(size(gh2)+2));
for i = 1:max(size(W))
   
    if(size(W{i,1})==size(gh2))
       torf = W{i,1} == gh2;
       df= all(torf==1);
             if df==1 
             current_section=i; break; end
     end
end

current_section;
