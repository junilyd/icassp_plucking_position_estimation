function [pluckCmFromBridge, Cn,regCoeff] = icasssp19_plucking_position_estimator_LSD(amplitudes,L)

%clear all 
clf
delta=1;
M = min(length(amplitudes),46);
m = (1:M)';
amplitudes=amplitudes(1:M)';
ii=0;
D_LS=[];
pCm = 6:0.0001:ceil(L/2); % search grid in cm
Ppickup = 4.4423/L;
for cnt=pCm
    ii=ii+1;
    P=cnt/L;

    S(:,1) = ones(M,1);
    S(:,2) = 2*delta./(m.^2.*pi^2*P*(1-P)).*abs(sin((m.*pi*P))); 
    %S(:,2) = S(:,2)./max(S(:,2));
  %  S(:,3) = abs(sin((n.*pi*Ppickup))); 

    alphaHat = amplitudes(1:M); %.*(m.^2);
    alphaHat = alphaHat(:);%./max(alphaHat);

    regCoeff(:,ii) = S\alphaHat;
%
    alphaModel(:,ii) =   regCoeff(1,ii)+ regCoeff(2,ii) .*S(:,2);% + regCoeff(3,ii).*S(:,3);
%    D_LS(ii) = sqrt( 1/numAmps * sum( (10*log10(alphaHat./alphaModel(:,ii))).^2) ) ;

    D_LS(ii) = sqrt( 1/M * sum( ( 10*log10(alphaHat./( S(:,2) )) ).^2) ) ;
    
end


[~,argmin] = min(D_LS);
Cn = alphaModel(:,argmin);
regCoeff = regCoeff(:,argmin);
pluckCmFromBridge = pCm(argmin);

%figure(101)
%plot(pCm, D_LS)
% figure(100)
% stem(n,log(y),'b');
% hold on
% stem(n,log(yhat(:,argmin)),'r','linewidth',2);

end