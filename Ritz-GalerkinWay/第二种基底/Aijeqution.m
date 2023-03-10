function Aij = Aijeqution(i,j)
    Aij=i*j/(i+j-1)-2*i*j/(i+j)-1+(i*j+i+j)/(i+j+1)+2.0/(i+j+2)-1.0/(i+j+3);
end