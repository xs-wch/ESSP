function amariindex = amari(mix,estmix,nc)% calculate Amari index
B = ((mix'*mix)^(-1))*mix'*estmix;    

tempB = abs(B);
amariindex = (sum(sum(tempB,1)./max(tempB,[],1))+sum(sum(tempB,2)./max(tempB,[],2))-2*nc)/2/nc;
return