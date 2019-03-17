% takes in paramter values and computes the steady states and the stability
% of the found steady states. Returns everything via cells. Model = 1 is
% toggle switch model, and anything else is the inhibition model
function [pkm, actin, rna, hs, stability, indexes] = computeSSNew(model, j1,j2,j3,j4,j5,j6)

% define returning vectors of length max of j1 and j2 (for looping through
% different parameter values. Only set up for changing j1 and j2
total = max(max(length(j1),length(j2)), length(j5));
pkm = cell(total,1);
actin = cell(total,1);
rna = cell(total,1);
hs = cell(total,1);
stability = cell(total,1);
indexes = cell(total, 1);
type = 0;
if model ~= 1
    type = 1;
end

for j = 1:length(j1)
    for q = 1:length(j2)
        for k = 1:length(j5)
            % find index for ease of indexing
            index = max(max(j,q),k);
            %preload empty arrays for concatenating error bounds
            pkm{index} = [];
            actin{index} = [];
            rna{index} = [];
            hs{index} = [];
            stability{index} = [];

            % start the computations!
            if type
                %inhibition Model
                ssActin = @(pkm) (j2(q)+ j3.*pkm)./(1+j2(q)+j3.*pkm);
                ssHs = @(pkm) (j5(k).*ssActin(pkm))./(1+pkm+j5(k).*ssActin(pkm));
                ssRNA = @(pkm) (j4.*ssActin(pkm).*pkm)/(1+j5(k).*ssHs(pkm)+j4.*ssActin(pkm).*pkm);
                ssPkm =@(pkm) j1(j)*ssRNA(pkm)/(1+ssHs(pkm)+j1(j)*ssRNA(pkm));
                % jacobian entries
                jacMat = @(pkm, actin, rna, hs) ...
                    [-j1(j).*rna-1, 0, j1(j)-j1(j).*pkm, 0;
                    j3*(1-actin), -j2(q)-j3*pkm-1, 0, 0;
                    (j4*actin*(1-rna))/(1+j5(k)*hs), (j4*pkm*(1-rna))/(1+j5(k)*hs),...
                        -(j4*actin*pkm)/(1+j5(k)*hs)-1, -j4*actin*pkm*(1-rna)/(1+j5(k)*hs)^2;
                    0, j6*(1-hs), 0, -j6*actin - 1];
            else
            
                % off-Toggle Model
                ssActin = @(pkm) (j2(q)+ j3.*pkm)./(1+j2(q)+j3.*pkm);
                ssHs = @(pkm) 0;
                ssRNA = @(pkm) (j4.*ssActin(pkm).*pkm)./(1+j4.*ssActin(pkm).*pkm);
                ssPkm =@(pkm) j1(j)*ssRNA(pkm)/(1+ssHs(pkm)+j1(j)*ssRNA(pkm));

                %Jacobian entries
                jacMat = @(pkm, actin, rna, hs) ...
                        [-j1(j).*rna-1, 0, j1(j)-j1(j).*pkm, -1;
                        j3*(1-actin), -j2(q)-j3.*pkm-1, 0, 0;
                        j4.*actin-j4.*actin.*rna, j4.*pkm-j4.*pkm.*rna, -j4...
                            .*actin.*pkm - 1, 0;
                        0, j5(k)-j5(k).*hs, 0, -j5(k).*actin-1];
            end

            % symbolic solving for roots
            syms x;
            ssAns = double(vpasolve(ssPkm(x) - x ==0, x));

            for i = 1:length(ssAns)
                p = ssAns(i);
                if isreal(p) && p>=0
                    a = ssActin(p);
                    r = ssRNA(p);
                    h = ssHs(p);
                    pkm{index} = [pkm{index}, p];
                    actin{index} = [actin{index}, a];
                    rna{index} = [rna{index}, r];
                    hs{index}= [hs{index}, h];
                    jac = jacMat(p, a, r, h);
                    [~,~,isStable] = find(sign(eigs(jac)) == 1, 1); % returns -empty if no postive eigenvalues
                    if isempty(isStable)
                        stability{index} = [stability{index}, 1];
                    else
                        stability{index} = [stability{index}, -1];
                    end
                    indexes{index} = [indexes{index}, index];
                end
            end
        end
    end
end
end
