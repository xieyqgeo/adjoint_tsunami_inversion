function dk1=conjugate_direction(it,gk1,gk,dk)
if it==1
    dk1=gk1;
else
    dk1=gk1+sum(gk1(:).*gk1(:))/sum(gk(:).*gk(:))*dk;
end
end