%NCFSYN  H-infinity normalized coprime factor controller synthesis.
%
% [K,CL,GAM,INFO]=NCFSYN(G) 
% [K,CL,GAM,INFO]=NCFSYN(G,W1) 
% [K,CL,GAM,INFO]=NCFSYN(G,W1,W2) computes the Glover-McFarlane 
% 	H-infinity normalized coprime factor (NCF) loop-shaping 
% 	controller K for LTI plant G. with LTI weights W1 and W2. 
%        +         ____      ____      ____     ____
%    r ---> + --->| W2 |--->| Ks |--->| W1 |---| G  |-------> y
%           ^      ----      ----      ----     ----     |
%         + |                                            |
%           |____________________________________________|
% 
% 	The positive feedback controller Ks robustly stabilizes a family 
% 	of systems given by a ball of uncertainty in the normalized coprime
% 	factors of the shaped plant Gs=W2*G*W1.
%  
%     Input:
%     -------
%      G  -   LTI plant to be controlled, which must be square (NxN) 
%
%     Optional Input:
%     -------
%      W1,W2  - stable min-phase LTI weights, either SISO or square (NxN).
%               Default is W1=W2=eye(N)
%     Output:
%     -------
%      K      -  K=W1*Ks*W2, H-infinity optimal loopshaping controller 
%      CL     -  CL = [I; Ks]*inv(I-Gs*Ks)*[I, Gs] , closed-loop 
%      GAM    -  the H-infinity optimal cost 
%      INFO   -  struct array with various information fields, such as
%                INFO.emax   -  the nugap robustness, EMAX=1/GAM 
%                INFO.Gs     -  the 'shaped plant' Gs=W2*G*W1
%                INFO.Ks     -  Ks = NCFSYN(Gs) = NCFSYN(Gs,I,I)
% 
% 	THEORY:  Gs = W2*G*W1 is called the 'shaped plant'. The controller
% 	  Ks simultaneously solves both of the following two H-infinity 
%     optimal control problems for the shaped plant Gs:
%          1).    min  || [I; Ks]*inv(I-Gs*Ks)*[I, Gs] ||
%                  K                                     Hinf
% 	
%          2).    min  || [I; Gs]*inv(I-Ks*Gs)*[I, Ks] ||
%                  Ks                                    Hinf
% 	and GAM is the minimal H-infinity COST for both problems.   
% 	
%   See also: GAPMETRIC, HINFSYN, LOOPSYN, NCFMARGIN

% OLD HELP
% function  [sysk,emax,sysobs]=ncfsyn(sysgw,factor,opt)
%
%   Synthesizes a controller to robustly stabilize
%   a family of systems given by a ball of uncertainty
%   in the normalized coprime factors of the system description.
%
%   Inputs:
%   -------
%    SYSGW  - the weighted system to be controlled.
%    FACTOR - =1 implies that an optimal controller is required.
%             >1 implies that a suboptimal controller is required
%                achieving a performance FACTOR less than optimal.
%    OPT    - 'ref'  includes a reference input (optional).
%
%
%   Outputs:
%   -------
%    SYSK   -  H-infinity loopshaping controller
%    EMAX   -  stability margin as an indication robustness to
%              unstructured perturbations. EMAX is always less
%              than 1 and values of EMAX below 0.3 generally indicate
%              good robustness margins.
%    SYSOBS -  H-infinity loopshaping observer controller. This variable
%              is created only if FACTOR>1 and OPT = 'ref'
%
%   See also: HINFNORM, HINFSYN, and NUGAP

%   Copyright 1991-2004 MUSYN Inc. and The MathWorks, Inc.

 function    [sysk,emax,sysobs]=ncfsyn(sysgw,factor,opt)


if nargin<2,
    disp('usage: [k,cl]=ncfsyn(sysgw,factor,opt)');
    return;
 end;

[type,p,m,n]=minfo(sysgw);

if type~='syst',
    disp('sysgw needs to be of type SYSTEM');
    return;
 end;

if nargin<3, opt='noref'; end


if strcmp(opt,'noref') & (factor==1),
     [sysnlcf,sig] = sncfbal(sysgw);
   %  [tmp1,tmp2,tmp3,n]=minfo(sysnlcf);
     emax = sqrt(1-sig(1)^2);
     sysUV=hankmr(sysnlcf,sig,0,'a');
     sysUV=cjt(sysUV);
     sysk=mscl(cf2sys(sysUV),-1);

 end; %if strcmp(opt,'noref') & (factor==1),

if strcmp(opt,'noref') & (factor>1),
     [a,b,c,d]=unpck(sysgw);
     R=eye(p)+d*d';
     S=eye(m)+d'*d;
     Hamx=[a-b*(S\d')*c, -(b/S)*b'; -c'*(R\c), -(a-b*(S\d')*c)'];
     Hamz=[(a-b*(d'/R)*c)', -c'*(R\c) ; -b*(S\b')  -(a-b*(d'/R)*c)];
     [x1,x2,fail]=ric_schr(Hamx);
     if fail == 0
	X=x2/x1;
     else
	sysk = [];
        disp('Error in solving the X Hamiltonian');
	return
     end
     [z1,z2,fail]=ric_schr(Hamz);
     if fail == 0
	Z=z2/z1;
     else
	sysk = [];
        disp('Error in solving the Z Hamiltonian');
	return
     end
     F=-S\(d'*c+b'*X);
     H=-(b*d'+Z*c')/R;
     sig=eig(Z*X);
     emax=1/sqrt(1+max(real(sig)));
     gamma=factor/emax;
     W2=eye(n)-gamma^(-2)*(eye(n)+Z*X);
     sysk=pck(a+b*F-W2\(Z*c'*(c+d*F)),W2\(Z*c'),-b'*X,-d');
 end; %if strcmp(opt,'noref') & (factor>1)

if (nargin>2) & strcmp(opt,'ref') & factor>1,
     [a,b,c,d]=unpck(sysgw);
     R=eye(p)+d*d';
     S=eye(m)+d'*d;
     Hamx=[a-b*(S\d')*c, -(b/S)*b'; -c'*(R\c), -(a-b*(S\d')*c)'];
     Hamz=[(a-b*(d'/R)*c)', -c'*(R\c) ; -b*(S\b')  -(a-b*(d'/R)*c)];
     [x1,x2,fail]=ric_schr(Hamx);
     if fail == 0
	X=x2/x1;
     else
	sysk = [];
        disp('Error in solving the X Hamiltonian');
	return
     end
     [z1,z2,fail]=ric_schr(Hamz);
     if fail == 0
	Z=z2/z1;
     else
	sysk = [];
        disp('Error in solving the Z Hamiltonian');
	return
     end
     F=-S\(d'*c+b'*X);
     H=-(b*d'+Z*c')/R;
     sig=eig(Z*X);
     emax=1/sqrt(1+max(real(sig)));
     gamma=factor/emax;
     W2=eye(n)-gamma^(-2)*(eye(n)+Z*X);
     Hi=-(b*d'+W2\(Z*c'))/R;
     sqrtS=sqrtm(S); bm=(b-W2\(Z*c'*d))/S;
     sysobs=pck(a+Hi*c,[-Hi, bm, zeros(n,m)],...
              -b'*X,[-d', zeros(m,m), -sqrtS]);
     sysk=pck(a+Hi*c-bm*b'*X,[-Hi-bm*d',-bm*sqrtS],-b'*X,[-d', -sqrtS]);
end; %  if (nargin>2) & (opt=='ref') & factor>1

if  strcmp(opt,'ref') & (factor==1),
     [sysnrcf,sig] = sncfbal(transp(sysgw));
     sysnrcf=transp(sysnrcf);
     emax = sqrt(1-sig(1)^2);
     disp(['emax = ',num2str(emax)]);
     sysUV=hankmr(sysnrcf,sig,0,'a');
     sysUV=cjt(sysUV);
     sysUV=mmult(sysUV,[-eye(p) zeros(p,m);zeros(m,p) eye(m)]);
     sysL=mmult(sysUV,[-eye(p) zeros(p,m); zeros(m,p) eye(m)],sysnrcf);
     sysUVn=mmult(minv(sysL),sysUV);
     sysUVn=sysbal(sysUVn);
     [typeUV,pUV,mUV,nUV]=minfo(sysUV);
     [typeUVn,pUVn,mUVn,nUVn]=minfo(sysUVn);
     sysUVn=strunc(sysUVn,min(nUV,nUVn)); % remove nonminimal states
     sysUVn=madd(mmult(sysUVn,[eye(p) zeros(p,2*m); zeros(m,m+p) eye(m)]),...
                  [zeros(m,p) eye(m) zeros(m,m)]);
     sysk=cf2sys(sysUVn);
  end; %if (nargin>2) & (opt=='ref') & (factor==1)
%
%