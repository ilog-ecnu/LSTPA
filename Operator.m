function Offspring = Operator(Loser,Winner,FitnessDiff)
% The competitive swarm optimizer of LSTPA

    %% Parameter setting
    LoserDec  = Loser.decs;
    WinnerDec = Winner.decs;
    [N,D]     = size(LoserDec);
	LoserVel  = Loser.adds(zeros(N,D));
    WinnerVel = Winner.adds(zeros(N,D));
    Problem   = PROBLEM.Current();

    %% Competitive swarm optimizer

    meanr2 = mean(FitnessDiff);
    coeff = mean(power(10,FitnessDiff));
 
    FitnessDiff = repmat(FitnessDiff,1,D);

    r1     = repmat(rand(N,1),1,D);
    r2     = repmat(rand(N,1),1,D);
    r3     = repmat(rand(N,1),1,D);
    
    tmp = FitnessDiff.*(WinnerDec - LoserDec);      % (6)
    OffVel = r1.*LoserVel + (1-r1).*tmp;            % (4)
    OffDec = LoserDec + r2.*(WinnerDec - LoserDec) + r3.*(OffVel).*(FitnessDiff+1);    %(5)

    
    %% Add the winners
    OffDec = [OffDec;WinnerDec];
    OffVel = [OffVel;WinnerVel];
 
    %% Polynomial mutation
    Lower  = repmat(Problem.lower,2*N,1);
    Upper  = repmat(Problem.upper,2*N,1);
    disM   = 20;
    Site   = rand(2*N,D) < 1/D;
    mu     = rand(2*N,D);
    temp   = Site & mu<=0.5;
    OffDec       = max(min(OffDec,Upper),Lower);
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                   (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp  = Site & mu>0.5; 
    OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                   (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
	Offspring = SOLUTION(OffDec,OffVel);
end