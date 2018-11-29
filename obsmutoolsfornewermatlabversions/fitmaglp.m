% function [sys] = fitmaglp(magdata,wt,heading,osysl_g,...
%                                dmdi,upbd,blk,blknum)
%
%  FITMAGLP fits a stable, minimum phase transfer function
%   to magnitude data, MAGDATA, with a supplied frequency
%   domain weighting function, WT. Both these should
%   be VARYING matrices, with identical INDEPENDENT
%   VARIABLE values.
%
%   HEADING (optional) is the title displayed and OSYSL_G
%   (optional) is the FRSP of the old fit.  If given, it will be
%   displayed along with the data.
%
%   The (optional) variable DMDI is the matrix that produced the
%   MAGDATA, WT and UPBD data. Upon fitting the D, the stable,
%   minimum phase system matrix SYS, is absorbed into the original
%   matrix DMDI and plotted along with the UPBD data. BLK is the
%   block structure using in the MU command and BLKNUM defines
%   the current D scaling to be fit.
%
%   A comparison of these two plots provides insight into how
%   well the rational system matrix SYS fits the D data. The newly
%   scaled matrix DMDI, which is CLPG with SYS wrapped is then
%   return. If DMDI, UPBD, BLK and BLKNUM are not provided, the
%   default is to plot the WT variable instead. These options are
%   used with MUSYNFIT.
%
%   Identical to FITMAG except that MAGFIT is used instead
%   of GENPHASE and FITSYS.
%
%   See also: FITSYS, INVFREQS, FITMAG, MUSYNFIT,
%             MUFTBTCH, and MUSYNFLP.

%   Copyright 1991-2004 MUSYN Inc. and The MathWorks, Inc.

function [sys] = fitmaglp(d,wt,heading,osysl_g,dmdi,upbd,blk,blknum)

 if nargin == 0
   disp(['usage: sys = fitmaglp(magdata,weight)']);
   return
 end
 if nargin == 1
   wt = vpck(ones(length(getiv(d)),1),getiv(d));
   heading = 'CURVE FITTING';
 end
 if nargin == 2
   heading = 'CURVE FITTING';
 end
 if nargin == 5 | nargin == 6 | nargin == 7 | nargin >8
   error(['Incorrect input argument']);
   return
 end

 [dtype,drows,dcols,dnum] = minfo(d);
 if dtype ~= 'vary'
   error(['magdata should be a VARYING matrix']);
   return
 end

 [wtype,wrows,wcols,wnum] = minfo(wt);
 if wtype ~= 'vary'
   error(['weight should be a VARYING matrix']);
   return
 end
 code = indvcmp(d,wt);
 if code ~= 1
    error('inconsistent VARYING data in magdata and weight')
    return
 end

 axisvar = 'liv,m';
 laxisvar = 'liv,lm';
 omega = getiv(d);
 mlinsp = max(diff(omega)) - min(diff(omega));
 if abs(mlinsp) < 1e-5
   axisvar = 'iv,m';
   laxisvar = 'iv,lm';
 end

%FIRST: SIMPLE CURVEFITTING
 if nargin <= 3
   clf;subplot(211);
   vplot(laxisvar,d);
   xlabel('mag data ');
   if ~isempty(heading)
     title([heading]);
   end
   subplot(212);
   vplot(laxisvar,wt);  %%
   title('SENSITIVITY WEIGHT FOR FIT');
   xlabel(['NOTE APPROXIMATE ORDER NECESSARY FOR FIT.....']);
   ord = input('ENTER ORDER OF CURVE FIT or ''drawmag''      ');

   while  looptst(ord,'cn')
      ord = input('try again - a nonnegative integer or ''drawmag''...    ');
   end

   go = 1;
   while go == 1
     if ~strcmp('drawmag',ord)
       sysh = magfit(d,[.26,.1,ord,ord]);
       hhead = [',   W/ORDER = ' int2str(ord)];
     else
       [sysh] = drawmag(d);
       [mattype,rodw,cold,num] = minfo(sysh);
       hhead = [',   W/ORDER = ' int2str(num) ' (drawmag)' ];
     end
     sysh_g = frsp(sysh,d);
     clf;subplot(211);
     vplot(laxisvar,d,sysh_g);
     xlabel('  1) mag data    2) newfit ');
     title([heading hhead]);
     subplot(212);
     vplot(laxisvar,wt);  %%
     title('SENSITIVITY WEIGHT FOR FIT');
     xlabel('ENTER NEW ORDER, ''drawmag'', or NEGATIVE NUMBER TO STOP ');
     pause

     ord = input('ENTER NEW ORDER, ''drawmag'', or NEGATIVE TO STOP    ');
     while looptst(ord)
       ord = input('try again - an integer, ''drawmag'', or negative to stop... ');
     end %while any
     go = 0;
     if isstr(ord)
       if strcmp('drawmag',ord)
        go = 1;
       end
     else
       if ord >= 0
         go = 1;
       end
     end
   end %while ord
   sys = sysh;
   %subplot(111);
%  ----------  END OF ONLY_CUREVFITTING  ----------
 elseif nargin == 4
   clf;subplot(211);
   vplot(laxisvar,d,osysl_g);
   xlabel('  data and old fit ');
   if ~isempty(heading)
     title([heading]);
   end
   subplot(212);
   vplot(laxisvar,wt);
   title('WT FOR FIT');
   xlabel(['NOTE APPROXIMATE ORDER NECESSARY FOR FIT.....']);
   ord = input('ENTER ORDER OF CURVE FIT or ''drawmag''      ');

   while  looptst(ord,'cn')
     ord = input('try again - a nonnegative integer or ''drawmag''...    ');
   end

   go = 1;
   while go == 1
     if ~strcmp('drawmag',ord)
       sysh = magfit(d,[.26,.1,ord,ord]);
       hhead = [',   W/ORDER = ' int2str(ord)];
     else
       [sysh] = drawmag(d);
       [mattype,rodw,cold,num] = minfo(sysh);
       hhead = [',   W/ORDER = ' int2str(num) ' (drawmag)' ];
     end
     sysh_g = frsp(sysh,d);
     clf;subplot(211);
     vplot(laxisvar,d,sysh_g,osysl_g);
     xlabel('  1) data    2) newfit   3) oldfit ');
     title([heading hhead]);
     subplot(212);
     vplot(laxisvar,wt);
     title('WEIGHT FOR FIT');

     xlabel('ENTER NEW ORDER, ''drawmag'', or NEGATIVE NUMBER TO STOP ');
     pause;
     ord = input('ENTER NEW ORDER, ''drawmag'', or NEGATIVE TO STOP    ');
     while looptst(ord)
       ord = input('try again - an integer, ''drawmag'', or negative to stop... ');
     end %while any
     go = 0;
     if isstr(ord)
       if strcmp('drawmag',ord)
         go = 1;
       end
     else
       if ord >= 0
         go = 1;
       end
     end
   end %while ord
   sys = sysh;
   %subplot(111);
%  ----------  END OF OLD FITMAGLP  ----------
 else
   clf;subplot(211);
   vplot(laxisvar,d,osysl_g);
   xlabel('mag data and previous D-K rational scaling');
   if ~isempty(heading)
	  title([heading]);
   end
   subplot(212);
   vplot(laxisvar,wt);
   title('SENSITIVITY WEIGHT FOR FIT');
   xlabel(['NOTE APPROXIMATE ORDER NECESSARY FOR FIT.....']);
   ord = input('ENTER ORDER OF CURVE FIT or ''drawmag''      ');

   while  looptst(ord,'cn')
     ord = input('try again - a nonnegative integer or ''drawmag''...    ');
   end

   blksum = cumsum([1 1;blk]);
   inputs = [blksum(blknum,1):blksum(blknum+1,1)-1]';
   outputs = [blksum(blknum,2):blksum(blknum+1,2)-1]';

   go = 1;
   while go == 1
     if ~strcmp('drawmag',ord)
       sysh = magfit(d,[.26,.1,ord,ord]);
       hhead = [',   W/ORDER = ' int2str(ord)];
     else
       sysh = drawmag(d);
       [mattype,rodw,cold,num] = minfo(sysh);
       hhead = [',   W/ORDER = ' int2str(num) ' (drawmag)' ];
     end
     sysh_g = frsp(sysh,d);
     dratodf = vrdiv(sysh_g,d);
     tmpdmdi = sclin(sclout(dmdi,outputs,vinv(dratodf)),inputs,dratodf);
     clf;subplot(211);
     vplot(laxisvar,d,sysh_g,osysl_g);
     xlabel('  1) mag data    2) newfit   3) previous D-K');

     title([heading hhead]);
     subplot(212);
     vplot(axisvar,upbd,vnorm(tmpdmdi));
     title('SCALED TRANSFER FUNCTION: OPTIMAL & RATIONAL');
     xlabel(' 1) MU UPPER BND  2) UPPER BND WITH RATIONAL FIT ');
     pause;
     ord = input('ENTER NEW ORDER, ''drawmag'', or NEGATIVE TO STOP    ');
     while looptst(ord)
       ord = input('try again - an integer, ''drawmag'', or negative to stop... ');
     end %while any
     go = 0;
     if isstr(ord)
       if strcmp('drawmag',ord)
         go = 1;
       end
     else
       if ord >= 0
         go = 1;
       end
     end
   end %while ord
   sys = sysh;
   %subplot(111);
 end
%
%