function [DATA,C] = test_RealSamples(data0,data1,Cants)

%Generate a random set of samples. Mixes real data, number from each cluster provided by Cants
CantClas=2;
Inic = 0;
[~,DimeData] = size(data0);

for k=1:CantClas
   C1(1,Inic+1:Inic+Cants(k))=ones(1,Cants(k))*k;
   Inic = Inic+Cants(k);
   if k==1
    samp0 = datasample(data0,Cants(k),'Replace',false);
   elseif k==2
    samp1 = datasample(data1,Cants(k),'Replace',false);   
   end
end

TotalCant = Inic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     mixes the classes     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:TotalCant
   Pos = unidrnd(TotalCant-k+1,1,1);
   C(k)=C1(Pos);
   C1(Pos)=C1(TotalCant-k+1);
   C1(TotalCant-k+1)=C(k);
end

DATA = zeros(TotalCant,DimeData);
zero_counter=1;
one_counter=1;
for k=1:TotalCant

    if C(k) == 1
        DATA(k,:) = samp0(zero_counter,:);
        zero_counter = zero_counter + 1;
    elseif C(k)==2
        DATA(k,:) = samp1(one_counter,:);
        one_counter = one_counter + 1;        
    end
end

end

