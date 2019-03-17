function dydt = neuronFireODENewTerm(t,y, model, stimFun,stimFun2,j1,j2,j3,j4,j5,j6,ta1,...
    ta2,ta3,ta4,ta5)
type = 0;
if model ~= 1
    type = 1;
end
stim = stimFun(t);
stim2 = stimFun2(t);
hs = y(4);
p = y(1);
actin = y(2);
rna = y(3);
actin = y(2);
% dydt = [dydt; 0] + [0%-j6/(1+hs)./ta1;
%             0;
%             j6/(1+hs)/ta3;
%             (j5.*actin.*stim*(1 - hs)-hs)./ta5];
% dydt(3) = ((j4.*actin.*(p+stim).*(1-rna))./(1+hs)-rna)./ta3;

if type
    % inhibition model
    dydt = [((j1*rna*(1-p)) - p)/ta1;
        ((j2+j3*p)*(1-actin)-actin)/ta2;
        ((j4*actin*(p+stim)*(1-rna))/(1+j5*hs)-rna)/ta3;
        (j5*actin*(1-hs)-hs)/ta5];
else
    % toggle-Off model
    dydt = [(j1.*rna.*(1-p) - p)./ta1;
        ((j2+j3.*p).*(1-actin) - actin)./ta2;
        (j4.*actin.*(p+stim).*(1-rna)/(1+stim2)-rna)./ta3];
end

end