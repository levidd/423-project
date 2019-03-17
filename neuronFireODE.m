function dydt = neuronFireODE(t,y, stimFun,j1,j2,j3,j4,j5,ta1,ta2,ta3,ta4,ta5)
stim = stimFun(t);
p = y(1);
actin = y(2);
rna = y(3);

dydt = [(j1.*rna.*(1-p) - p)./ta1;
    ((j2+j3.*p).*(1-actin) - actin)./ta2;
    (j4.*actin.*(p+stim).*(1-rna)-rna)./ta3];
end
