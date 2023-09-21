function Population = EnvironmentalSelection(Population,V,theta)
    % The environmental selection of LSTPA
    
        PopObj = Population.objs;
        [N,M]  = size(PopObj);
        NV     = size(V,1);
        
        %% Translate the population
        PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
        
        %% Calculate the degree of violation of each solution
        CV = sum(max(0,Population.cons),2);
        
        %% Calculate the smallest angle value between each vector and others
        cosine = 1 - pdist2(V,V,'cosine');
        cosine(logical(eye(length(cosine)))) = 0;
        gamma  = min(acos(cosine),[],2);
    
        %% Associate each solution to a reference vector
        Angle = acos(1-pdist2(PopObj,V,'cosine'));
        [~,associate] = min(Angle,[],2);
        
        [tFrontNo,tFitness] = tNDSortFitness(PopObj,V,theta);
        %% Select one solution for each reference vector
        Next = zeros(1,NV);
        for i = unique(associate)'
            current1 = find(associate==i & CV==0);
            current2 = find(associate==i & CV~=0);
            if ~isempty(current1)
                % Select the one with the minimum APBI value
                APBI = tFitness(current1);
                [~,best] = min(APBI);
                Next(i)  = current1(best);
            elseif ~isempty(current2)
                % Select the one with the minimum CV value
                [~,best] = min(CV(current2));
                Next(i)  = current2(best);
            end
        end
        % Population for next generation
        Population = Population(Next(Next~=0));
    end
    
    
    function [tFrontNo,Fitness] = tNDSortFitness(PopObj,W,thetaGen)
    % Do theta-non-dominated sorting
    
        N  = size(PopObj,1);
        NW = size(W,1);
    
        %% Calculate the d1 and d2 values for each solution to each weight
        normP  = sqrt(sum(PopObj.^2,2));
        Cosine = 1 - pdist2(PopObj,W,'cosine');
        d1     = repmat(normP,1,size(W,1)).*Cosine;
        d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2);
        
        %% Clustering
        [~,class] = min(d2,[],2);
        
        %% Sort
        theta = zeros(1,NW) + 5;
        theta(sum(W>1e-4,2)==1) = 1e6;
        tFrontNo = zeros(1,N);
        Fitness = zeros(1,N);
    
        tmp_min = tan(pi/NW/4);         % (14)
        for i = 1 : NW
            C = find(class==i);
            [~,rank] = sort(d1(C,i)+theta(i)*d2(C,i));
            tFrontNo(C(rank)) = 1 : length(C);
        
            Fitness(C)=d1(C,i)+theta(i)*d2(C,i)*(thetaGen+tmp_min); % (9)
    
        end
    %     Fitness = min(Dis,[],2);
    end
    