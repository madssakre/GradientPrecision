%% main
function opt = GRADPREC(spc,khr,opt)
% function opt = GRADPREC(spc,khr,opt)
%
%   GRADPREC described in:
% 
%       Optimal control gradient precision trade-offs: Application to fast generation of DeepControl libraries for MRI
%       Vinding, M.S., et al.
%       J Magn Reson.  2021, 333, 107094
%       https://doi.org/10.1016/j.jmr.2021.107094
%
%       https://github.com/madssakre/GradientPrecision
%
%   is an extension to the blOCh optimization framework:
%
%       https://github.com/madssakre/blOCh 
%
%   which has been used in:
%
%       Application of the limited-memory quasi-Newton algorithm for multi-dimensional, large flip-angle RF pulses at 7T.
%       Vinding, M.S., et al.,
%       Magn Reson Mater Phy. 2017, 30, 29-39.
%       doi:10.1007/s10334-016-0580-1
% 
%       Local SAR, global SAR, and power-constrained large-flip-angle pulses with optimal control and virtual observation points.
%       Vinding, M.S., et al.,
%       Magn Reson Med. 2017, 77, 1, 374-384.
%       doi: 10.1002/mrm.26086
% 
%       Real-time 2D spatially selective MRI experiments: Comparative analysis of optimal control design methods.
%       Maximov, I. et al.,
%       J Magn Reson. 2015, 254, 110-20.
%       doi: 10.1016/j.jmr.2015.03.003
% 
%       Fast numerical design of spatial-selective rf pulses in MRI using Krotov and quasi-Newton based optimal control methods.
%       Vinding, M.S., et al.,
%       J Chem Phys. 2012, 137, 054203.
%       doi: 10.1063/1.4739755
% 
%       and more...
%
%   From the adjacent "main.m" you can run GRADPREC with several
%   different control gradients with different accuracies. You can find the
%   gradient implementations in the "GRAD*" functions below.
%
%   Please cite the appropriate papers if you use GRADPREC or blOCh. And it
%   comes with no warranty etc.
%     
%
%
%     Copyright (C) 2021  Mads Sloth Vinding
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
% 


global history opt2 dt
tic_Opt_Wall_Clock_Time = tic;

% adopt variables from other structs to avoid passing entire structs around
opt.R = spc.R;
opt.Rv = spc.Rv;
opt.D = spc.D;
opt.Dv = spc.Dv;
opt.NS = spc.NS;
opt.Dim = spc.Dim;
opt.pTx = spc.pTx;
opt.P = spc.P;
opt.Md = spc.Md;
opt.idxtot = spc.idxtot;
opt.gamma = khr.gamma;
opt.B1_nom_amp_val = spc.B1_nom_amp_val;
opt.B1_nom_amp_unit = spc.B1_nom_amp_unit;
opt.f_B1_val = spc.f_B1_val;

% default parameters. 
Par.Grad= 'midpoint';
Par.StoreGrad = 'off';

% Masks in this context refers to indices in the control waveforms that we
% don't optimize but keep fixed. % this is not applicable in all gradient
% types
Par.Mask.RF = 'Off';
Par.Mask.GRA  = 'Off';

% This specifies which controls we optimize uv means real and imag RF. s
% means shims, x, y, z means Gx, Gy, Gz. Respectively. Any combination is
% possible e.g. usz.
Par.Controls = 'uv'; % uvxyzs

% for average power and SAR we need to know the pulse duty cycle
Par.Constr.Dutycycle = 0.1;

% about peak RF. This is the one we  pay attention to.
Par.Constr.peak_uv.Type = 'bnd';
Par.Constr.peak_uv.Lim = 1000;
Par.Constr.peak_uv.Unit = 'Hz';
Par.Constr.peak_uv.FreeIterations = 10;
Par.Constr.peak_uv.Cope = 'Constr';

% we can constrain the first and final time slices to zero. Ignore for now
Par.Constr.edges_uv.Type = 'bnd';
Par.Constr.edges_uv.Lim  = 0;
Par.Constr.edges_uv.Unit = 'W';
Par.Constr.edges_uv.FreeIterations = 0;
Par.Constr.edges_uv.Cope = 'Ignore';

% about average power. Ignore for now
Par.Constr.ave_uv.Type = 'nlc';
Par.Constr.ave_uv.Lim = 1e16;
Par.Constr.ave_uv.Unit = '(rad/s)^2';
Par.Constr.ave_uv.FreeIterations = 10;
Par.Constr.ave_uv.Cope = 'Ignore';


% about jaggedness, not really implemented yet. Ignore
Par.Constr.jag_uv.Type = 'nlc';
Par.Constr.jag_uv.Lim = 1;
Par.Constr.jag_uv.Unit = '(rad/s)^2/s';
Par.Constr.jag_uv.FreeIterations = 10;
Par.Constr.jag_uv.Cope = 'Ignore';


% about peak Grad. This is the one we  pay attention to.
Par.Constr.peak_xyz.Type = 'bnd';
Par.Constr.peak_xyz.Lim = 40e-3;
Par.Constr.peak_xyz.Unit = 'T/m';
Par.Constr.peak_xyz.FreeIterations = 10;
Par.Constr.peak_xyz.Cope = 'Constr';

% we can constrain the first and final time slices to zero. Ignore
Par.Constr.edges_xyz.Type = 'bnd';
Par.Constr.edges_xyz.Lim  = 0;
Par.Constr.edges_xyz.Unit = 'T/m';
Par.Constr.edges_xyz.FreeIterations = 0;
Par.Constr.edges_xyz.Cope = 'Ignore';


% about peak dGrad/dt. This is the one we  pay attention to.
Par.Constr.slew_xyz.Type = 'nlc';
Par.Constr.slew_xyz.Lim = 200;
Par.Constr.slew_xyz.Unit = 'T/m/s';
Par.Constr.slew_xyz.FreeIterations = 10;
Par.Constr.slew_xyz.Cope = 'Constr';

% about RF robustness
Par.Robust.Num = 1; % not che
Par.Robust.Min = 0.9;
Par.Robust.Max = 1.1;
Par.Robust.FWHM = 0.1;
Par.Robust.Distribution = 'lorentzian'; % 'evenly'
Par.Robust.Show = 1;
Par.Robust.FixGradient = 1;


% About right or left handedness. % dont try this at home
Par.Turn = 1;

% Here I exhange default parameters with user specified parameters
opt.Par = NEWPAR(opt.Par,Par);

% about sanity checking controls to optimize
[opt.Par.Controls,opt.Par.PosControls,opt.Par.FixControls] = CONTROLS(opt.Par.Controls);

% about masks in time domain.
if strcmp(opt.Par.Mask.RF,'Off') && strcmp(opt.Par.Mask.GRA,'Off')
    
    opt.mon = [1:opt.N]; %#ok<*NBRAK>
    opt.moff = [];
    opt.mon_u = [1:opt.N];
    opt.mon_v = [1:opt.N];
    opt.mon_x =  [1:opt.N];
    opt.mon_y =  [1:opt.N];
    opt.mon_z =  [1:opt.N];
    opt.mon_s =  [1:opt.N];
    
    opt.m_u = ones(opt.pTx,opt.N);
    opt.m_v = ones(opt.pTx,opt.N);
    opt.m_x = ones(1,opt.N);
    opt.m_y = ones(1,opt.N);
    opt.m_z = ones(1,opt.N);
    opt.m_s = ones(opt.NS,opt.N);
    
elseif strcmp(opt.Par.Mask.RF,'Off') && strcmp(opt.Par.Mask.GRA,'On')
    
    opt.mon = find(opt.m);
    opt.moff = find(~opt.m);
    opt.mon_u = [1:opt.N];
    opt.mon_v = [1:opt.N];
    opt.mon_x =  opt.mon;
    opt.mon_y =  opt.mon;
    opt.mon_z =  opt.mon;
    opt.mon_s =  [1:opt.N];
    
    opt.m_u = ones(opt.pTx,opt.N);
    opt.m_v = ones(opt.pTx,opt.N);
    opt.m_x = opt.m;
    opt.m_y = opt.m;
    opt.m_z = opt.m;
    opt.m_s = ones(opt.NS,opt.N);
    
elseif strcmp(opt.Par.Mask.RF,'On') && strcmp(opt.Par.Mask.GRA,'Off')
    
    opt.mon = find(opt.m);
    opt.moff = find(~opt.m);
    opt.mon_u = opt.mon;
    opt.mon_v = opt.mon;
    opt.mon_x =  [1:opt.N];
    opt.mon_y =  [1:opt.N];
    opt.mon_z =  [1:opt.N];
    opt.mon_s =  [1:opt.N];
    
    opt.m_u = repmat(opt.m,[opt.pTx,1]);
    opt.m_v = repmat(opt.m,[opt.pTx,1]);
    opt.m_x = ones(1,opt.N);
    opt.m_y = ones(1,opt.N);
    opt.m_z = ones(1,opt.N);
    opt.m_s = ones(opt.NS,opt.N);
    
elseif strcmp(opt.Par.Mask.RF,'On') && strcmp(opt.Par.Mask.GRA,'On')
    
    opt.mon = find(opt.m);
    opt.moff = find(~opt.m);
    opt.mon_u = opt.mon;
    opt.mon_v = opt.mon;
    opt.mon_x =  opt.mon;
    opt.mon_y =  opt.mon;
    opt.mon_z =  opt.mon;
    opt.mon_s =  [1:opt.N];
    
    opt.m_u = repmat(opt.m,[opt.pTx,1]);
    opt.m_v = repmat(opt.m,[opt.pTx,1]);
    opt.m_x = opt.m;
    opt.m_y = opt.m;
    opt.m_z = opt.m;
    opt.m_s = ones(opt.NS,opt.N);
    
end

opt = CONSTRAINTSANITY(opt);

[opt.Par.Robust.Values,opt.Par.Robust.Weightings] = RFROBUSTNESS(opt.Par.Robust.Num,opt.Par.Robust.Min,opt.Par.Robust.Max,opt.Par.Robust.Distribution,opt.Par.Robust.FWHM);


opt.k = 0;

opt = CONSTRAINHOW(opt);

opt = ALLOCATE(spc,khr,opt);

opt.s = zeros(opt.NS,opt.N);

arraymask = [];
if contains(opt.Par.Controls,'v')
    arraymask = [arraymask,REARRANGE('u','array',opt.m_u)];
end
if contains(opt.Par.Controls,'v')
    arraymask = [arraymask,REARRANGE('v','array',opt.m_v)];
end

if contains(opt.Par.Controls,'x')
    arraymask = [arraymask,REARRANGE('x','array',opt.m_x)];
end
if contains(opt.Par.Controls,'y')
    arraymask = [arraymask,REARRANGE('y','array',opt.m_y)];
end
if contains(opt.Par.Controls,'z')
    arraymask = [arraymask,REARRANGE('z','array',opt.m_z)];
end
if contains(opt.Par.Controls,'s')
    arraymask = [arraymask,REARRANGE('s','array',opt.m_s)];
end

opt.arraymask = arraymask;


opt.gx = opt.g(1,:);
opt.gy = opt.g(2,:);
opt.gz = opt.g(3,:);

opt.gxo = zeros(1,opt.N,opt.MaxIter+1);
opt.gyo = zeros(1,opt.N,opt.MaxIter+1);
opt.gzo = zeros(1,opt.N,opt.MaxIter+1);

opt.gxo(1,:,1) = opt.gx;
opt.gyo(1,:,1) = opt.gy;
opt.gzo(1,:,1) = opt.gz;


opt.array0 = ARRAY(opt.Par.Controls,opt.u,opt.v,opt.gx,opt.gy,opt.gz,opt.s);

opt.array0(isinf(abs(opt.array0))) = 0;
opt.array0(isnan(opt.array0)) = 0;

opt.max_uv = opt.w1m;

opt.max_g = opt.Par.Constr.peak_xyz.Lim;

[opt.arrayub,opt.arraylb,opt.arraysc] = UB_LB_SC(opt);

opt.array0 = ARRAY2NORM(opt.array0,opt.arraysc);

options=optimset('TolX',opt.TolX,'TolFun',opt.TolFun,'Display','off','MaxIter',opt.MaxIter,'SubproblemAlgorithm','cg',...
    'MaxFunEvals',opt.MaxFunEvals,'Algorithm','interior-point','Hessian',{'lbfgs',10},'ScaleProblem','none',...
    'GradObj','on','DerivativeCheck',opt.deriv_check,'FinDiffType','central','LargeScale','off','GradConstr','on','OutputFcn',@OUTPUT);

dt = opt.dt;
opt.Go = true;
opt2 = opt;
[opt.array,fval,opt.exitflag,opt.output] =fmincon(@FUN,opt.array0,[],[],[],[],opt.arraylb,opt.arrayub,@CONSTR,options,opt);

fprintf(opt.output.message);

opt.Go = false;
opt.k = 1;
opt.ksafe = history.safeiter;

opt.Par.Constr.sar_l.Viol = history.Constr.sar_l.Viol;
opt.Par.Constr.sar_g.Viol = history.Constr.sar_g.Viol;

opt.Par.Constr.peak_uv.Violu = history.Constr.peak_uv.Violu;
opt.Par.Constr.peak_uv.Violv = history.Constr.peak_uv.Violv;

opt.Par.Constr.ave_uv.Violuv = history.Constr.ave_uv.Violuv;

opt.Par.Constr.jag_uv.Viol = history.Constr.jag_uv.Viol;

opt.Par.Constr.peak_xyz.Violx = history.Constr.peak_xyz.Violx;
opt.Par.Constr.peak_xyz.Violy = history.Constr.peak_xyz.Violy;
opt.Par.Constr.peak_xyz.Violz = history.Constr.peak_xyz.Violz;

opt.Par.Constr.slew_xyz.Violx = history.Constr.slew_xyz.Violx;
opt.Par.Constr.slew_xyz.Violy = history.Constr.slew_xyz.Violy;
opt.Par.Constr.slew_xyz.Violz = history.Constr.slew_xyz.Violz;

Copes = {opt.Par.Constr.peak_uv.Cope,opt.Par.Constr.peak_uv.Cope,opt.Par.Constr.ave_uv.Cope,opt.Par.Constr.jag_uv.Cope,opt.Par.Constr.peak_xyz.Cope,opt.Par.Constr.peak_xyz.Cope,opt.Par.Constr.peak_xyz.Cope,opt.Par.Constr.slew_xyz.Cope,opt.Par.Constr.slew_xyz.Cope,opt.Par.Constr.slew_xyz.Cope}; %,opt.Par.NRMSE.Cope

Viols = [opt.Par.Constr.peak_uv.Violu,opt.Par.Constr.peak_uv.Violv,opt.Par.Constr.ave_uv.Violuv,opt.Par.Constr.jag_uv.Viol,opt.Par.Constr.peak_xyz.Violx,opt.Par.Constr.peak_xyz.Violy,opt.Par.Constr.peak_xyz.Violz,opt.Par.Constr.slew_xyz.Violx,opt.Par.Constr.slew_xyz.Violy,opt.Par.Constr.slew_xyz.Violz]; % ,~opt.Par.NRMSE.Viol

opt.ksafe = KSAFE(opt.ksafe,Copes,Viols);

opt.StopCriteria = history.StopCriteria;
opt.Fun = history.Fun(1:opt.ksafe);
opt.dFun = history.dFun(1:opt.ksafe);

opt.con.peak_u   = zeros(1,opt.ksafe);
opt.con.peak_v   = zeros(1,opt.ksafe);
opt.con.ave_uv   = zeros(1,opt.ksafe);
opt.con.jag_uv   = zeros(1,opt.ksafe);
opt.con.peak_x   = zeros(1,opt.ksafe);
opt.con.peak_y   = zeros(1,opt.ksafe);
opt.con.peak_z   = zeros(1,opt.ksafe);
opt.con.slew_x   = zeros(1,opt.ksafe);
opt.con.slew_y   = zeros(1,opt.ksafe);
opt.con.slew_z   = zeros(1,opt.ksafe);


if isfield(history,'Durations')
    opt.Durations = history.Durations;
else
    opt.Durations = -1;
end

opt.uo = zeros(opt.pTx,opt.N,opt.ksafe);

opt.vo = zeros(opt.pTx,opt.N,opt.ksafe);

opt.gxo = zeros(1,opt.N,opt.ksafe);
opt.gyo = zeros(1,opt.N,opt.ksafe);
opt.gzo = zeros(1,opt.N,opt.ksafe);

for p = 1:opt.ksafe
    
    opt.k = p;
    array_scaled = ARRAY2PHYS(history.array(p,:),opt.arraysc);
    
    [u_,v_,gx_,gy_,gz_,s_] = iARRAY(opt.Par.Controls,array_scaled,opt.N,opt.pTx,opt.NS,opt.u,opt.v,opt.gx,opt.gy,opt.gz,opt.s);
    opt.uo(:,:,p) = u_;
    opt.vo(:,:,p) = v_;
    
    opt.gxo(1,:,p) = gx_;
    opt.gyo(1,:,p) = gy_;
    opt.gzo(1,:,p) = gz_;
    
    opt.con.peak_u(p)   = history.con(p).peak_u;
    opt.con.peak_v(p)   = history.con(p).peak_v;
    opt.con.ave_uv(p)   = history.con(p).ave_uv;
    opt.con.jag_uv(p)   = history.con(p).jag_uv;
    opt.con.peak_x(p)   = history.con(p).peak_x;
    opt.con.peak_y(p)   = history.con(p).peak_y;
    opt.con.peak_z(p)   = history.con(p).peak_z;
    opt.con.slew_x(p)   = history.con(p).slew_x;
    opt.con.slew_y(p)   = history.con(p).slew_y;
    opt.con.slew_z(p)   = history.con(p).slew_z;
    
    
end

% for saving

opt.Path2Package = strrep(mfilename('fullpath'),mfilename,'');

opt = rmfield(opt,{'u','v','g','Md','mon','mon_u','mon_v','mon_x','mon_y','mon_z','mon_s','m_u','m_v','m_x','m_y','m_z','m_s','sr','si','w0','yX','yY','yZ','yS','Gm','Sm','L_t','M_t','s','arraymask','gx','gy','gz','array0','arrayub','arraylb','arraysc','array','ksaveintermediate'});



opt.Durations(opt.ksafe+1:end) = [];
clearvars -global history
opt.Opt_Wall_Clock_Time = toc(tic_Opt_Wall_Clock_Time);
opt.Go = false;
end

%% Optimization and physics functions

function [Fun,arraygrad] = FUN(array,opt)
global M_T Duration

array_scaled = ARRAY2PHYS(array,opt.arraysc);

[u,v,gx,gy,gz,s] = iARRAY(opt.Par.Controls,array_scaled,opt.N,opt.pTx,opt.NS,opt.u,opt.v,opt.gx,opt.gy,opt.gz,opt.s);

[M_t,L_t] = PROPAGATE(opt.N,opt.pTx,opt.dt,u,v,gx,gy,gz,s,opt.yX,opt.yY,opt.yZ,opt.yS,opt.sr,opt.si,opt.w0,opt.M_t,opt.L_t,opt.Par.Turn);

if nargout >= 2
    
    tic_One_Iteration = tic;
    if strcmp(opt.Par.Grad,'standard_0')
        [Gu,Gv] = GRAD_STANDARD_0(opt.Par.Controls,opt.N,opt.P,opt.pTx,opt.dt,opt.sr,opt.si,M_t,L_t);
    elseif strcmp(opt.Par.Grad,'standard_1')
        [Gu,Gv] = GRAD_STANDARD_1(opt.Par.Controls,opt.N,opt.P,opt.pTx,opt.dt,opt.sr,opt.si,M_t,L_t);
    elseif contains(opt.Par.Grad,'midpoint')
        [Gu,Gv,Gx,Gy,Gz,Gs] = GRAD_MIDPOINT(opt.Par.Controls,opt.N,opt.P,opt.pTx,opt.NS,opt.dt,opt.yX,opt.yY,opt.yZ,opt.yS,opt.sr,opt.si,M_t,L_t,opt.Par.Turn);
    elseif contains(opt.Par.Grad,'exact')
        [Gu,Gv,Gx,Gy,Gz,Gs] = GRAD_EXACT(opt.Par.Controls,opt.par_Ncores,opt.N,opt.P,opt.pTx,opt.NS,opt.dt,opt.yX,opt.yY,opt.yZ,opt.yS,opt.sr,opt.si,opt.w0,M_t,L_t,u,v,gx,gy,gz,s,opt.mon_u,opt.mon_v,opt.mon_x,opt.mon_y,opt.mon_z,opt.mon_s,opt.Par.Turn,opt.gamma);
    elseif strcmp(opt.Par.Grad,'standard_0_loop')
        Gu = GRAD_STANDARD_0_forloop('u',opt.N,opt.pTx,opt.P,opt.dt,opt.sr,opt.si,opt.mon_u,M_t,L_t);
        Gv = GRAD_STANDARD_0_forloop('v',opt.N,opt.pTx,opt.P,opt.dt,opt.sr,opt.si,opt.mon_v,M_t,L_t);
    elseif strcmp(opt.Par.Grad,'standard_1_loop')
        Gu = GRAD_STANDARD_1_forloop('u',opt.N,opt.pTx,opt.P,opt.dt,opt.sr,opt.si,opt.mon_u,M_t,L_t);
        Gv = GRAD_STANDARD_1_forloop('v',opt.N,opt.pTx,opt.P,opt.dt,opt.sr,opt.si,opt.mon_v,M_t,L_t);
    end
    Duration = toc(tic_One_Iteration);
    arraygrad = [];
    if contains(opt.Par.Controls,'u')
        arraygrad = [arraygrad,REARRANGE('u','array',Gu)];
    end
    if contains(opt.Par.Controls,'v')
        arraygrad = [arraygrad,REARRANGE('v','array',Gv)];
    end
    
    if contains(opt.Par.Controls,'x')
        arraygrad = [arraygrad,REARRANGE('x','array',Gx)];
    end
    if contains(opt.Par.Controls,'y')
        arraygrad = [arraygrad,REARRANGE('y','array',Gy)];
    end
    if contains(opt.Par.Controls,'z')
        arraygrad = [arraygrad,REARRANGE('z','array',Gz)];
    end
    
    
    arraygrad = arraygrad.*opt.arraymask;
    arraygrad = ARRAY2PHYS(arraygrad,opt.arraysc);
    
    
    if strcmp(opt.Par.StoreGrad,'on')
        
        GradFileID = fopen('GradFile1.txt','a');
        
        
        for h = 1:length(arraygrad)
            fprintf(GradFileID,'%02.6e ',arraygrad(h));
        end
        fprintf(GradFileID,'\n');
        fclose(GradFileID);
        
    end
    
    
end

M_T = M_t(:,end);


Fun = -EFFICIENCY(opt.Md,M_T,opt.P);



end

function [M_t,L_t] = PROPAGATE(N,pTx,dt,u,v,gx,gy,gz,s,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,Turn)


R11 = cell(1,N);R12 = cell(1,N);R13 = cell(1,N);
R21 = cell(1,N);R22 = cell(1,N);R23 = cell(1,N);
R31 = cell(1,N);R32 = cell(1,N);R33 = cell(1,N);


for n = 1:N
    [R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n}] = ROTATOR(u,v,gx,gy,gz,s,w0,yX,yY,yZ,yS,pTx,sr,si,dt,n,Turn);
end



for n = 1:N
    M_t(:,n+1) = ROTATE(M_t(:,n),R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n},'Forward');
end


for n = N:-1:1
    L_t(:,n) = ROTATE(L_t(:,n+1),R11{n},R12{n},R13{n},R21{n},R22{n},R23{n},R31{n},R32{n},R33{n},'Backward');
end


end

function [R11,R12,R13,R21,R22,R23,R31,R32,R33] = ROTATOR(u,v,gx,gy,gz,s,w0,yX,yY,yZ,yS,pTx,sr,si,dt,n,varargin)

P = length(w0);
wx = zeros(size(w0));
wy = zeros(size(w0));
for ns = 1:pTx
    wx = wx+(sr(:,ns).*u(ns,n)-si(:,ns).*v(ns,n)).*dt;
    wy = wy+(si(:,ns).*u(ns,n)+sr(:,ns).*v(ns,n)).*dt;
end


wz = (yX*gx(n) + yY*gy(n) + yZ*gz(n) + w0).*dt;

if nargin == 16
    Turn = -1;
elseif nargin >= 16
    Turn = varargin{1};
end
switch Turn
    case -1
        wx = -wx+eps;
        wy = -wy+eps;
        wz = wz+eps;
    case 1
        wx = -wx+eps;
        wy = -wy+eps;
        wz = -wz+eps;
end

theta = sqrt(wx.^2+wy.^2+wz.^2);

x = wx./theta;
y = wy./theta;
z = wz./theta;


c = cos(theta);
s = sin(theta);

C =1-c;

zs = z.*s;
xs = x.*s;
ys = y.*s;

xyC = x.*y.*C;
yzC = y.*z.*C;
xzC = x.*z.*C;

R11 = x.*x.*C+c;
R12 = xyC-zs;
R13 = xzC+ys;

R21 = xyC+zs;
R22 = y.*y.*C+c;
R23 = yzC-xs;

R31 = xzC-ys;
R32 = yzC+xs;
R33 = z.*z.*C+c;

end

function ML_out = ROTATE(ML_in,R11,R12,R13,R21,R22,R23,R31,R32,R33,Type)



[P,N] = size(R11);

ML = ML_in(:)';

ML = reshape(ML,3,P);

ML = permute(ML,[2 1]);

MLx = ML(:,1);
MLy = ML(:,2);
MLz = ML(:,3);


if strcmp(Type,'Forward')
    
    if N > 1
        
        for n = 1:N
            
            MLx_ = R11(:,n).*MLx + R12(:,n).*MLy + R13(:,n).*MLz;
            MLy_ = R21(:,n).*MLx + R22(:,n).*MLy + R23(:,n).*MLz;
            MLz_ = R31(:,n).*MLx + R32(:,n).*MLy + R33(:,n).*MLz;
            
            MLx = MLx_;
            MLy = MLy_;
            MLz = MLz_;
            
            
            
        end
        
        
    else
        
        MLx_ = R11.*MLx + R12.*MLy + R13.*MLz;
        MLy_ = R21.*MLx + R22.*MLy + R23.*MLz;
        MLz_ = R31.*MLx + R32.*MLy + R33.*MLz;
        
        MLx = MLx_;
        MLy = MLy_;
        MLz = MLz_;
        
    end
    
    
    
elseif strcmp(Type,'Backward')
    
    if N > 1
        
        for n = N:-1:1
            
            MLx_ = R11(:,n).*MLx + R21(:,n).*MLy + R31(:,n).*MLz;
            MLy_ = R12(:,n).*MLx + R22(:,n).*MLy + R32(:,n).*MLz;
            MLz_ = R13(:,n).*MLx + R23(:,n).*MLy + R33(:,n).*MLz;
            
            MLx = MLx_;
            MLy = MLy_;
            MLz = MLz_;
            
            
            
        end
        
        
    else
        
        MLx_ = R11.*MLx + R21.*MLy + R31.*MLz;
        MLy_ = R12.*MLx + R22.*MLy + R32.*MLz;
        MLz_ = R13.*MLx + R23.*MLy + R33.*MLz;
        
        MLx = MLx_;
        MLy = MLy_;
        MLz = MLz_;
        
        
    end
    
    
    
end

ML(:,1) = MLx;
ML(:,2) = MLy;
ML(:,3) = MLz;

ML = permute(ML,[2 1]);

ML = reshape(ML,3*P,1);

ML_out = ML(:);
end

function [Eff] = EFFICIENCY(Md,M_T,P,varargin)


A_x = Md(1:3:end);
A_y = Md(2:3:end);
A_z = Md(3:3:end);

B_x = M_T(1:3:end);
B_y = M_T(2:3:end);
B_z = M_T(3:3:end);
Phi = sum(A_x.*B_x+A_y.*B_y+A_z.*B_z);


Eff = Phi/P;



end

function opt = PROGRESS(opt)


Stop = zeros(14,1);

if strcmp(opt.Par.Constr.peak_uv.Cope,'Stop') && opt.Par.Constr.peak_uv.Violu && opt.k > opt.Par.Constr.peak_uv.FreeIterations
    Stop(3) = 1;
    fprintf('%s: Will stop because rf u peak is violated\n',mfilename);
end
if strcmp(opt.Par.Constr.peak_uv.Cope,'Stop') && opt.Par.Constr.peak_uv.Violv && opt.k > opt.Par.Constr.peak_uv.FreeIterations
    Stop(4) = 1;
    fprintf('%s: Will stop because rf v peak is violated\n',mfilename);
end
if strcmp(opt.Par.Constr.ave_uv.Cope,'Stop') && opt.Par.Constr.ave_uv.Viol && opt.k > opt.Par.Constr.ave_uv.FreeIterations
    Stop(5) = 1;
    fprintf('%s: Will stop because uv average is violated\n',mfilename);
end
if strcmp(opt.Par.Constr.peak_xyz.Cope,'Stop') && opt.Par.Constr.peak_xyz.Violx && opt.k > opt.Par.Constr.peak_xyz.FreeIterations
    Stop(6) = 1;
    fprintf('%s: Will stop because gradient x peak is violated\n',mfilename);
end
if strcmp(opt.Par.Constr.peak_xyz.Cope,'Stop') && opt.Par.Constr.peak_xyz.Violy && opt.k > opt.Par.Constr.peak_xyz.FreeIterations
    Stop(7) = 1;
    fprintf('%s: Will stop because gradient y peak is violated\n',mfilename);
end
if strcmp(opt.Par.Constr.peak_xyz.Cope,'Stop') && opt.Par.Constr.peak_xyz.Violz && opt.k > opt.Par.Constr.peak_xyz.FreeIterations
    Stop(8) = 1;
    fprintf('%s: Will stop because gradient z peak is violated\n',mfilename);
end
if strcmp(opt.Par.Constr.slew_xyz.Cope,'Stop') && opt.Par.Constr.slew_xyz.Violx && opt.k > opt.Par.Constr.slew_xyz.FreeIterations
    Stop(9) = 1;
    fprintf('%s: Will stop because gradient x slew is violated\n',mfilename);
end
if strcmp(opt.Par.Constr.slew_xyz.Cope,'Stop') && opt.Par.Constr.slew_xyz.Violy && opt.k > opt.Par.Constr.slew_xyz.FreeIterations
    Stop(10) = 1;
    fprintf('%s: Will stop because gradient y slew is violated\n',mfilename);
end
if strcmp(opt.Par.Constr.slew_xyz.Cope,'Stop') && opt.Par.Constr.slew_xyz.Violz && opt.k > opt.Par.Constr.slew_xyz.FreeIterations
    Stop(11) = 1;
    fprintf('%s: Will stop because gradient z slew is violated\n',mfilename);
end




% save intermediate
if opt.k-1 > 0
    if mod(opt.k-1,opt.ksaveintermediate) == 0
        
        %         opt.arrayintermed = Rearrange_controls(opt.u,opt.v);
        tempsave = opt.Save; % take backup
        
        opt.Save = struct('Bundle','reduced','Data',tempsave.Data,'Controls',tempsave.Controls,'Figures','none','Scripts','none');
        try
            b6_Save_Job([],[],opt,[]);
            
        catch me
            fprintf('%s: %s',mfilename,me.message);
        end
        
        opt.Save = tempsave;
        %         opt.k = tempk;
        
    end
end

if sum(Stop)>0
    opt.Go = false;
    opt.StopCriteria = Stop;
    opt.ksafe = opt.k;
    opt.k = opt.k+1;
    return
else
    opt.StopCriteria = Stop;
    opt.Go = true;
    opt.ksafe = opt.k;
    opt.k = opt.k+1;
end




end

function [stop] = OUTPUT(array,optimValues,state,opt)
global M_T history Duration opt2

stop = false;

switch state
    case 'init'

        history.array = zeros(opt2.MaxIter,length(array)); 
        
        history.array(1,:) = array;
        
        history.fval(1)= optimValues.fval;
        
        history.safeiter  = optimValues.iteration+1;
        
        array_scaled = ARRAY2PHYS(array,opt.arraysc);
        opt2.arrayiintermed = array;
        
        [opt2.u,opt2.v,opt2.gx,opt2.gy,opt2.gz,opt2.s] = iARRAY(opt.Par.Controls,array_scaled,opt.N,opt.pTx,opt.NS,opt.u,opt.v,opt.gx,opt.gy,opt.gz,opt.s);
        
        opt2.M_T = M_T;
        history.M_T{1} = M_T;
        
        [Fun] = EFFICIENCY(opt.Md,M_T,opt.P);
        
        [con,opt2] = CONSTR_dummyfunc(opt2);
        
        history.con(1) = con;
        
        opt2.con = history.con;
        history.Fun(1) = Fun;
        
        opt2.dFun = [];
        opt2 = PROGRESS(opt2);
        PRINT(opt2,'all',Fun,con);
        opt2.dFun = 0;
        history.dFun(1) = opt2.dFun;
        
        history.Constr.sar_l.Viol = opt2.Par.Constr.sar_l.Viol;
        history.Constr.sar_g.Viol = opt2.Par.Constr.sar_g.Viol;
        history.Constr.peak_uv.Violu = opt2.Par.Constr.peak_uv.Violu;
        history.Constr.peak_uv.Violv = opt2.Par.Constr.peak_uv.Violv;
        history.Constr.ave_uv.Violuv = opt2.Par.Constr.ave_uv.Violuv;
        history.Constr.jag_uv.Viol = opt2.Par.Constr.jag_uv.Viol;
        history.Constr.peak_xyz.Violx = opt2.Par.Constr.peak_xyz.Violx;
        history.Constr.peak_xyz.Violy= opt2.Par.Constr.peak_xyz.Violy;
        history.Constr.peak_xyz.Violz = opt2.Par.Constr.peak_xyz.Violz;
        history.Constr.slew_xyz.Violx = opt2.Par.Constr.slew_xyz.Violx;
        history.Constr.slew_xyz.Violy = opt2.Par.Constr.slew_xyz.Violy;
        history.Constr.slew_xyz.Violz = opt2.Par.Constr.slew_xyz.Violz;
        
        history.StopCriteria = opt2.StopCriteria;
        
    case 'iter'

        if optimValues.iteration > 0
            
            iter = optimValues.iteration+1;
            
            history.array(iter,:) = array;
            
            history.fval(iter)= optimValues.fval;
            history.M_T{iter} = M_T;
            
            history.safeiter  = iter;
            
            array_scaled = ARRAY2PHYS(array,opt.arraysc);
            
            opt2.arrayiintermed = array;
            
            [opt2.u,opt2.v,opt2.gx,opt2.gy,opt2.gz,opt2.s] = iARRAY(opt.Par.Controls,array_scaled,opt.N,opt.pTx,opt.NS,opt.u,opt.v,opt.gx,opt.gy,opt.gz,opt.s);
            
            opt2.M_T = M_T;
            
            [Fun] = EFFICIENCY(opt.Md,M_T,opt.P);
            
            [con,opt2] = CONSTR_dummyfunc(opt2);
            
            
            history.Fun(iter) = Fun;
            
            history.con(iter) = con;
            opt2.con = history.con;
            opt2.Fun = history.Fun;
            
            opt2.dFun = history.Fun(iter)-history.Fun(iter-1);
            opt2.arrayintermed = array;
            opt2 = PROGRESS(opt2);
            PRINT(opt2,'all',Fun,con);
            history.dFun(iter) = opt2.dFun;
            
            
            history.Constr.sar_l.Viol = opt2.Par.Constr.sar_l.Viol;
            history.Constr.sar_g.Viol = opt2.Par.Constr.sar_g.Viol;
            history.Constr.peak_uv.Violu = opt2.Par.Constr.peak_uv.Violu;
            history.Constr.peak_uv.Violv = opt2.Par.Constr.peak_uv.Violv;
            history.Constr.ave_uv.Violuv = opt2.Par.Constr.ave_uv.Violuv;
            history.Constr.jag_uv.Viol = opt2.Par.Constr.jag_uv.Viol;
            history.Constr.peak_xyz.Violx = opt2.Par.Constr.peak_xyz.Violx;
            history.Constr.peak_xyz.Violy= opt2.Par.Constr.peak_xyz.Violy;
            history.Constr.peak_xyz.Violz = opt2.Par.Constr.peak_xyz.Violz;
            history.Constr.slew_xyz.Violx = opt2.Par.Constr.slew_xyz.Violx;
            history.Constr.slew_xyz.Violy = opt2.Par.Constr.slew_xyz.Violy;
            history.Constr.slew_xyz.Violz = opt2.Par.Constr.slew_xyz.Violz;
            
            history.Durations(iter) = Duration;
            history.StopCriteria = opt2.StopCriteria;
            if ~opt2.Go
                stop = true;
                
            end
            
        end
    case 'done'

    otherwise
end
end
%% midpoint gradient functions

function [Gu,Gv,Gx,Gy,Gz,Gs] = GRAD_MIDPOINT(Controls,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,M_t,L_t,Turn)

if contains(Controls,'u') || contains(Controls,'v')
    LzMx = L_t(3:3:end,:).*M_t(1:3:end,:);
    LzMy = L_t(3:3:end,:).*M_t(2:3:end,:);
    LxMz = L_t(1:3:end,:).*M_t(3:3:end,:);
    LyMz = L_t(2:3:end,:).*M_t(3:3:end,:);
end

switch Turn
    case -1
        
        if contains(Controls,'x') || contains(Controls,'y') || contains(Controls,'z')
            LxMy = -L_t(1:3:end,:).*M_t(2:3:end,:);
            LyMx = -L_t(2:3:end,:).*M_t(1:3:end,:);
        end
    case 1
        
        if contains(Controls,'x') || contains(Controls,'y') || contains(Controls,'z')
            LxMy = L_t(1:3:end,:).*M_t(2:3:end,:);
            LyMx = L_t(2:3:end,:).*M_t(1:3:end,:);
        end
end


if contains(Controls,'u')
    Gu = GRAD_u_MIDPOINT(pTx,dt,N,P,si,sr,LzMx,LzMy,LxMz,LyMz);
else
    Gu = [];
end
if contains(Controls,'v')
    Gv = GRAD_v_MIDPOINT(pTx,dt,N,P,si,sr,LzMx,LzMy,LxMz,LyMz);
else
    Gv = [];
end

if contains(Controls,'x')
    Gx = GRAD_xyz_MIDPOINT(dt,P,yX,LxMy,LyMx);
else
    Gx = [];
end
if contains(Controls,'y')
    Gy = GRAD_xyz_MIDPOINT(dt,P,yY,LxMy,LyMx);
else
    Gy = [];
end
if contains(Controls,'z')
    Gz = GRAD_xyz_MIDPOINT(dt,P,yZ,LxMy,LyMx);
else
    Gz = [];
end

Gs = [];




end

function Ga = GRAD_xyz_MIDPOINT(dt,P,yA,LxMy,LyMx)

temp1a=yA.*LxMy;
temp2a = sum(temp1a);
temp1b=yA.*LyMx;
temp2b = sum(temp1b);


Ga = temp2a-temp2b;

Ga = -(Ga(:,1:end-1)+Ga(:,2:end))./2*dt/P;

end

function Gu = GRAD_u_MIDPOINT(pTx,dt,N,P,si,sr,LzMx,LzMy,LxMz,LyMz)
Gu = zeros(pTx,N+1);
for s = 1:pTx
    temp1a=si(:,s).*LzMx;
    temp2a = sum(temp1a);
    temp1b=sr(:,s).*LzMy;
    temp2b = sum(temp1b);
    temp1c=si(:,s).*LxMz;
    temp2c = sum(temp1c);
    temp1d=sr(:,s).*LyMz;
    temp2d = sum(temp1d);
    
    
    Gu(s,:)  = temp2a-temp2b-temp2c+temp2d;
end
Gu = -(Gu(:,1:end-1)+Gu(:,2:end))./2*dt/P;
end

function Gv = GRAD_v_MIDPOINT(pTx,dt,N,P,si,sr,LzMx,LzMy,LxMz,LyMz)
Gv = zeros(pTx,N+1);
for s = 1:pTx
    temp1a=si(:,s).*LzMy;
    temp2a = sum(temp1a);
    temp1b=sr(:,s).*LzMx;
    temp2b = sum(temp1b);
    temp1c=si(:,s).*LyMz;
    temp2c = sum(temp1c);
    temp1d=sr(:,s).*LxMz;
    temp2d = sum(temp1d);
    
    Gv(s,:)  = temp2a+temp2b-temp2c-temp2d;
end
Gv = -(Gv(:,1:end-1)+Gv(:,2:end))./2*dt/P;
end


%% standard gradients (vectorized)

function [Gu,Gv] = GRAD_STANDARD_0(Controls,N,P,pTx,dt,sr,si,M_t,L_t)

LzMx = L_t(3:3:end,2:end).*M_t(1:3:end,1:end-1);
LzMy = L_t(3:3:end,2:end).*M_t(2:3:end,1:end-1);
LxMz = L_t(1:3:end,2:end).*M_t(3:3:end,1:end-1);
LyMz = L_t(2:3:end,2:end).*M_t(3:3:end,1:end-1);

if contains(Controls,'u')
    Gu = GRAD_u_STANDARD(N,P,pTx,dt,sr,si,LzMx,LzMy,LxMz,LyMz);
else
    Gu = [];
end
if contains(Controls,'v')
    Gv = GRAD_v_STANDARD(N,P,pTx,dt,sr,si,LzMx,LzMy,LxMz,LyMz);
else
    Gv = [];
end

end

function [Gu,Gv] = GRAD_STANDARD_1(Controls,N,P,pTx,dt,sr,si,M_t,L_t)

LzMx = L_t(3:3:end,2:end).*M_t(1:3:end,2:end);
LzMy = L_t(3:3:end,2:end).*M_t(2:3:end,2:end);
LxMz = L_t(1:3:end,2:end).*M_t(3:3:end,2:end);
LyMz = L_t(2:3:end,2:end).*M_t(3:3:end,2:end);

if contains(Controls,'u')
    Gu = GRAD_u_STANDARD(N,P,pTx,dt,sr,si,LzMx,LzMy,LxMz,LyMz);
else
    Gu = [];
end
if contains(Controls,'v')
    Gv = GRAD_v_STANDARD(N,P,pTx,dt,sr,si,LzMx,LzMy,LxMz,LyMz);
else
    Gv = [];
end

end

function Gu = GRAD_u_STANDARD(N,P,pTx,dt,sr,si,LzMx,LzMy,LxMz,LyMz)
Gu = zeros(pTx,N);

for s = 1:pTx
    temp1a=si(:,s).*LzMx;
    temp2a = sum(temp1a);
    temp1b=sr(:,s).*LzMy;
    temp2b = sum(temp1b);
    temp1c=si(:,s).*LxMz;
    temp2c = sum(temp1c);
    temp1d=sr(:,s).*LyMz;
    temp2d = sum(temp1d);
    
    
    Gu(s,:)  = temp2a-temp2b-temp2c+temp2d;
end

Gu = -Gu*dt/P;
end

function Gv = GRAD_v_STANDARD(N,P,pTx,dt,sr,si,LzMx,LzMy,LxMz,LyMz)
Gv = zeros(pTx,N);
for s = 1:pTx
    temp1a=si(:,s).*LzMy;
    temp2a = sum(temp1a);
    temp1b=sr(:,s).*LzMx;
    temp2b = sum(temp1b);
    temp1c=si(:,s).*LyMz;
    temp2c = sum(temp1c);
    temp1d=sr(:,s).*LxMz;
    temp2d = sum(temp1d);
    
    Gv(s,:)  = temp2a+temp2b-temp2c-temp2d;
end
Gv = -Gv*dt/P;
end

%% standard gradients (for loop)

function Guv = GRAD_STANDARD_0_forloop(uorv,N,pTx,P,dt,sr,si,mon,M_t,L_t)
Guv = zeros(pTx,N);
for n = 1:N
    if ismember(n,mon)
        switch uorv
            case 'u'
                for s = 1:pTx
                    Guv(s,n) = sum(si(:,s).*L_t(3:3:end,n+1).*M_t(1:3:end,n))-sum(sr(:,s).*L_t(3:3:end,n+1).*M_t(2:3:end,n))-sum(si(:,s).*L_t(1:3:end,n+1).*M_t(3:3:end,n))+sum(sr(:,s).*L_t(2:3:end,n+1).*M_t(3:3:end,n));
                end
            case 'v'
                for s = 1:pTx
                    Guv(s,n) = sum(si(:,s).*L_t(3:3:end,n+1).*M_t(2:3:end,n))+sum(sr(:,s).*L_t(3:3:end,n+1).*M_t(1:3:end,n))-sum(si(:,s).*L_t(2:3:end,n+1).*M_t(3:3:end,n))-sum(sr(:,s).*L_t(1:3:end,n+1).*M_t(3:3:end,n));
                end
        end
    end
end
Guv = -Guv/P*dt;
end

function Guv = GRAD_STANDARD_1_forloop(uorv,N,pTx,P,dt,sr,si,mon,M_t,L_t)
Guv = zeros(pTx,N);
for n = 1:N
    if ismember(n,mon)
        switch uorv
            case 'u'
                for s = 1:pTx
                    Guv(s,n) = sum(si(:,s).*L_t(3:3:end,n+1).*M_t(1:3:end,n+1))-sum(sr(:,s).*L_t(3:3:end,n+1).*M_t(2:3:end,n+1))-sum(si(:,s).*L_t(1:3:end,n+1).*M_t(3:3:end,n+1))+sum(sr(:,s).*L_t(2:3:end,n+1).*M_t(3:3:end,n+1));
                end
            case 'v'
                for s = 1:pTx
                    Guv(s,n) = sum(si(:,s).*L_t(3:3:end,n+1).*M_t(2:3:end,n+1))+sum(sr(:,s).*L_t(3:3:end,n+1).*M_t(1:3:end,n+1))-sum(si(:,s).*L_t(2:3:end,n+1).*M_t(3:3:end,n+1))-sum(sr(:,s).*L_t(1:3:end,n+1).*M_t(3:3:end,n+1));
                end
        end
    end
end
Guv = -Guv/P*dt;
end

%% exact gradients functions

function [Gu,Gv,Gx,Gy,Gz,Gs] = GRAD_EXACT(Controls,par_Ncores,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon_u,mon_v,mon_x,mon_y,mon_z,mon_s,Turn,gamma)

if contains(Controls,'u')
    Gu = GRAD_uv_EXACT('u',par_Ncores,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon_u,Turn);
else
    Gu = [];
end
if contains(Controls,'v')
    Gv = GRAD_uv_EXACT('v',par_Ncores,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon_v,Turn);
else
    Gv = [];
end
if contains(Controls,'x')
    Gx = GRAD_xyz_EXACT('x',par_Ncores,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon_x,Turn,gamma);
else
    Gx = [];
end
if contains(Controls,'y')
    Gy = GRAD_xyz_EXACT('y',par_Ncores,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon_y,Turn,gamma);
else
    Gy = [];
end
if contains(Controls,'z')
    Gz = GRAD_xyz_EXACT('z',par_Ncores,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon_z,Turn,gamma);
else
    Gz = [];
end

Gs = [];

end

function Guv = GRAD_uv_EXACT(uorv,par_Ncores,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon,Turn)

Guv = zeros(pTx,N);
if par_Ncores > 1
    parfor (n = 1:N,par_Ncores)
        if ismember(n,mon)
            
            Guv(:,n) = GRAD_subfunc_uv_EXACT(uorv,n,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon,Turn)
        end
    end
else
    for n = 1:N
        if ismember(n,mon)
            Guv(:,n) = GRAD_subfunc_uv_EXACT(uorv,n,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon,Turn);
        end
    end
end
Guv = -Guv/P;
end

function Guv = GRAD_subfunc_uv_EXACT(uorv,n,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon,Turn) %#ok<*INUSL>


z = yX*gx(n) + yY*gy(n)+yZ*gz(n)+w0;

x = sum((sr.*repmat(u(:,n).',[P,1])-si.*repmat(v(:,n).',[P,1])),2);
y = sum((si.*repmat(u(:,n).',[P,1])+sr.*repmat(v(:,n).',[P,1])),2);
switch uorv
    
    case 'u'
        
        a = si; b = -sr;
    case 'v'
        a = sr; b = si;
        
end
switch Turn
    case -1
        z = -z;
    case 1
end

Guv = zeros(pTx,1);

for s = 1:pTx
    %     tic
    Dm2 = zeros(P,3);
    Dm2(:,1) = y+1i.*a(:,s);
    
    Dm1 = zeros(P,3);
    Dm1(:,1) = -z;
    Dm1(:,2) = -x+1i.*b(:,s);
    
    Dp1 = zeros(P,3);
    
    Dp1(:,1) = z;
    Dp1(:,2) = x-1i.*b(:,s);
    
    Dp2 = zeros(P,3);
    Dp2(:,1) = -y-1i.*a(:,s);
    
    Dm2_ = Dm2.';Dm2_ = Dm2_(:);
    Dm1_ = Dm1.';Dm1_ = Dm1_(:);
    Dp1_ = Dp1.';Dp1_ = Dp1_(:);
    Dp2_ = Dp2.';Dp2_ = Dp2_(:);
    
    Dm2_ = [Dm2_;0;0]; %#ok<*AGROW>
    Dm1_ = [Dm1_;0;0];
    Dp1_ = [0;Dp1_;0];
    Dp2_ = [0;0;Dp2_];
    
    matrix = SPDIAGS_old([Dm2_,Dm1_,Dp1_,Dp2_],[-2,-1,1,2],P*3,P*3);
    
    vector_prop2=EXPV_IK(matrix,M_t(:,n),dt);
    
    Guv(s) = L_t(:,n+1)'*imag(vector_prop2);
    
end
end

function Gg = GRAD_xyz_EXACT(gxgyorgz,par_Ncores,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon,Turn,gamma)

Gg = zeros(1,N);
if par_Ncores > 1
    parfor (n = 1:N,par_Ncores)
        if ismember(n,mon)
            Gg(n) = GRAD_subfunc_xyz_EXACT(gxgyorgz,n,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon,Turn,gamma);
        end
    end
    
else
    for n = 1:N
        if ismember(n,mon)
            Gg(n) = GRAD_subfunc_xyz_EXACT(gxgyorgz,n,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon,Turn,gamma);
            
        end
    end
    
end
Gg = -Gg./P;
end

function [Gg] = GRAD_subfunc_xyz_EXACT(xyz,n,N,P,pTx,NS,dt,yX,yY,yZ,yS,sr,si,w0,M_t,L_t,u,v,gx,gy,gz,s,mon,Turn,gamma)


z = yX*gx(n) + yY*gy(n)+yZ*gz(n)+w0;

x = sum((sr.*repmat(u(:,n).',[P,1])-si.*repmat(v(:,n).',[P,1])),2);
y = sum((si.*repmat(u(:,n).',[P,1])+sr.*repmat(v(:,n).',[P,1])),2);


switch xyz
    case 'x'
        a = yX./gamma;
    case 'y'
        a = yY./gamma;
    case 'z'
        a = yZ./gamma;
end

switch Turn
    case -1
    case 1
        %         a = -a;
end

Dm2 = zeros(P,3);
Dm2(:,1) = y;

Dm1 = zeros(P,3);
Dm1(:,1) = -z-1i.*a;
Dm1(:,2) = -x;

Dp1 = zeros(P,3);

Dp1(:,1) = z+1i.*a;
Dp1(:,2) = x;

Dp2 = zeros(P,3);
Dp2(:,1) = -y;

Dm2_ = Dm2.';Dm2_ = Dm2_(:);
Dm1_ = Dm1.';Dm1_ = Dm1_(:);
Dp1_ = Dp1.';Dp1_ = Dp1_(:);
Dp2_ = Dp2.';Dp2_ = Dp2_(:);

Dm2_ = [Dm2_;0;0]; %#ok<*AGROW>
Dm1_ = [Dm1_;0;0];
Dp1_ = [0;Dp1_;0];
Dp2_ = [0;0;Dp2_];


matrix = SPDIAGS_old([Dm2_,Dm1_,Dp1_,Dp2_],[-2,-1,1,2],P*3,P*3);

vector_prop2=EXPV_IK(matrix,M_t(:,n),dt);

Gg = L_t(:,n+1)'*imag(vector_prop2).*gamma;

end

%% hard boundaries and scaling
function [arrayub,arraylb,arraysc] = UB_LB_SC(opt)
if contains(opt.Par.Controls,'u')|| contains(opt.Par.Controls,'v')
    
    if ~strcmp(opt.Par.Constr.peak_uv.Cope,'Constr')
        uub =  ones(size(opt.u)).*10000;
        ulb = -ones(size(opt.u)).*10000;
        
        usc = ones(size(opt.u)).*opt.max_uv;
        
        vub =  ones(size(opt.v)).*10000;
        vlb = -ones(size(opt.v)).*10000;
        
        vsc = ones(size(opt.v)).*opt.max_uv;
        
    else
        
        if strcmp(opt.Par.Constr.peak_uv.Type,'bnd')
            uub =  ones(size(opt.u));
            ulb = -ones(size(opt.u));
            
            usc = ones(size(opt.u)).*opt.max_uv;
            
            vub =  ones(size(opt.v));
            vlb = -ones(size(opt.v));
            
            vsc = ones(size(opt.v)).*opt.max_uv;
            
        end
        
    end
    
    if strcmp(opt.Par.Constr.edges_uv.Cope,'Constr') && strcmp(opt.Par.Constr.edges_uv.Type,'bnd')
        
        if opt.pTx > 1
            uub(:,1) = 0;
            uub(:,end) = 0;
            
            ulb(:,1) = 0;
            ulb(:,end) = 0;
            
            vub(:,1) = 0;
            vub(:,end) = 0;
            
            vlb(:,1) = 0;
            vlb(:,end) = 0;
        else
            uub(1) = 0;
            uub(end) = 0;
            
            ulb(1) = 0;
            ulb(end) = 0;
            
            vub(1) = 0;
            vub(end) = 0;
            
            vlb(1) = 0;
            vlb(end) = 0;
        end
        
        
        
    end
    
    
    
    
end


if contains(opt.Par.Controls,'x')
    
    if ~strcmp(opt.Par.Constr.peak_xyz.Cope,'Constr')
        
        xub =  ones(size(opt.gx)).*100;
        xlb = -ones(size(opt.gx)).*100;
        
        xsc = ones(size(opt.gx)).*opt.max_g;
        
    else
        if strcmp(opt.Par.Constr.peak_xyz.Type,'bnd')
            xub =  ones(size(opt.gx));
            xlb = -ones(size(opt.gx));
            
            xsc = ones(size(opt.gx)).*opt.max_g;
            
        end
        
    end
    
    if strcmp(opt.Par.Constr.edges_xyz.Cope,'Constr') && strcmp(opt.Par.Constr.edges_xyz.Type,'bnd')
        
        
        
        xub(1) = 0;
        xub(end) = 0;
        
        xlb(1) = 0;
        xlb(end) = 0;
        
        
    end
    
    
    
    
end

if contains(opt.Par.Controls,'y')
    
    if ~strcmp(opt.Par.Constr.peak_xyz.Cope,'Constr')
        
        yub =  ones(size(opt.gy)).*100;
        ylb = -ones(size(opt.gy)).*100;
        
        ysc = ones(size(opt.gy)).*opt.max_g;
        
    else
        if strcmp(opt.Par.Constr.peak_xyz.Type,'bnd')
            yub =  ones(size(opt.gy));
            ylb = -ones(size(opt.gy));
            
            ysc = ones(size(opt.gy)).*opt.max_g;
            
        end
        
    end
    
    
    if strcmp(opt.Par.Constr.edges_xyz.Cope,'Constr') && strcmp(opt.Par.Constr.edges_xyz.Type,'bnd')
        
        yub(1) = 0;
        yub(end) = 0;
        
        ylb(1) = 0;
        ylb(end) = 0;
        
        
    end
    
    
end
if contains(opt.Par.Controls,'z')
    
    if ~strcmp(opt.Par.Constr.peak_xyz.Cope,'Constr')
        
        zub =  ones(size(opt.gz)).*100;
        zlb = -ones(size(opt.gz)).*100;
        
        zsc = ones(size(opt.gz)).*opt.max_g;
        
    else
        if strcmp(opt.Par.Constr.peak_xyz.Type,'bnd')
            zub =  ones(size(opt.gz));
            zlb = -ones(size(opt.gz));
            
            zsc = ones(size(opt.gz)).*opt.max_g;
            
        end
        
    end
    
    if strcmp(opt.Par.Constr.edges_xyz.Cope,'Constr') && strcmp(opt.Par.Constr.edges_xyz.Type,'bnd')
        
        zub(1) = 0;
        zub(end) = 0;
        
        zlb(1) = 0;
        zlb(end) = 0;
        
        
    end
    
    
end



arrayub = [];
arraylb = [];
arraysc = [];

if exist('uub','var')
    arrayub = [arrayub,REARRANGE('u','array',uub)];
end
if exist('vub','var')
    arrayub = [arrayub,REARRANGE('v','array',vub)];
end


if exist('xub','var')
    arrayub = [arrayub,REARRANGE('x','array',xub)];
end
if exist('yub','var')
    arrayub = [arrayub,REARRANGE('y','array',yub)];
end
if exist('zub','var')
    arrayub = [arrayub,REARRANGE('z','array',zub)];
end

if exist('sub','var')
    arrayub = [arrayub,REARRANGE('s','array',sub)];
end



if exist('ulb','var')
    arraylb = [arraylb,REARRANGE('u','array',ulb)];
end
if exist('vlb','var')
    arraylb = [arraylb,REARRANGE('v','array',vlb)];
end

if exist('xlb','var')
    arraylb = [arraylb,REARRANGE('x','array',xlb)];
end
if exist('ylb','var')
    arraylb = [arraylb,REARRANGE('y','array',ylb)];
end
if exist('zlb','var')
    arraylb = [arraylb,REARRANGE('z','array',zlb)];
end

if exist('slb','var')
    arraylb = [arraylb,REARRANGE('s','array',slb)];
end


if exist('usc','var')
    arraysc = [arraysc,REARRANGE('u','array',usc)];
end
if exist('vsc','var')
    arraysc = [arraysc,REARRANGE('v','array',vsc)];
end

if exist('xsc','var')
    arraysc = [arraysc,REARRANGE('x','array',xsc)];
end
if exist('ysc','var')
    arraysc = [arraysc,REARRANGE('y','array',ysc)];
end
if exist('zsc','var')
    arraysc = [arraysc,REARRANGE('z','array',zsc)];
end

if exist('ssc','var')
    arraysc = [arraysc,REARRANGE('s','array',ssc)];
end
end

%% Utility functions
function varargout = REARRANGE(I,O,varargin)

% non exhaustive sanity check
if ~ischar(I) && ~ischar(O)
    fprintf(2,'REARRANGE: first two arguments must be strings.\n' )
    return
end


if strcmp(O,'array') && nargout ~= 1
    fprintf(2,'REARRANGE: you are requesting a wrong number of outputs\n' )
    return
end

if strcmp(I,'array') && nargin ~= 6
    fprintf(2,'REARRANGE: you are supplying a wrong number of inputs\n' )
    return
end
%
Type = strcat(I,'2',O);

switch Type
    % without s
    case  {'uvxyz2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        i3 = varargin{3};
        i4 = varargin{4};
        i5 = varargin{5};
        varargout{1} = [i1(:).',i2(:).',i3(:).',i4(:).',i5(:).'];
        
        
    case  {'uvyz2array','uvxz2array','uvxy2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        i3 = varargin{3};
        i4 = varargin{4};
        varargout{1} = [i1(:).',i2(:).',i3(:).',i4(:).'];
        
        
        
    case  {'uvx2array','uvy2array', 'uvz2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        i3 = varargin{3};
        varargout{1} = [i1(:).',i2(:).',i3(:).'];
        
    case {'u2array','v2array','s2array'}
        
        i1 = permute(varargin{1},[2,1]);
        varargout{1} = i1(:).';
        
    case  {'uv2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        varargout{1} = [i1(:).',i2(:).'];
        
        
    case  {'uxyz2array','vxyz2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i3 = varargin{2};
        i4 = varargin{3};
        i5 = varargin{4};
        varargout{1} = [i1(:).',i3(:).',i4(:).',i5(:).'];
        
        
    case  {'uy2array','ux2array', 'uz2array','vy2array', 'vx2array','vz2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i3 = varargin{2};
        varargout{1} = [i1(:).',i3(:).'];
        
        
    case  {'uyz2array','uxz2array','uxy2array','vxy2array','vyz2array','vxz2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i3 = varargin{3};
        i4 = varargin{4};
        varargout{1} = [i1(:).',i3(:).',i4(:).'];
        
    case  {'x2array','y2array','z2array'}
        
        
        i3 = varargin{1};
        varargout{1} = i3(:).';
        
    case  {'yz2array','xz2array', 'xy2array'}
        
        i3 = varargin{1};
        i4 = varargin{2};
        varargout{1} = [i3(:).',i4(:).'];
        
    case  'xyz2array'
        
        i3 = varargin{1};
        i4 = varargin{2};
        i5 = varargin{3};
        varargout{1} = [i3(:).',i4(:).',i5(:).'];
        % ----------------------
        
    case  {'array2u','array2v'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        
        varargout{1} = permute(reshape(array,N,pTx),[2,1]);
        
    case  {'array2s'}
        
        array = varargin{1};
        N = varargin{2};
        NS = varargin{4};
        
        varargout{1} = permute(reshape(array,N,NS),[2,1]);
        
        
    case  {'array2uv'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:2*N*pTx),N,pTx),[2,1]);
        
    case  {'array2uvx','array2uvy','array2uvz'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:2*N*pTx),N,pTx),[2,1]);
        varargout{3} = array(2*N*pTx+1:2*N*pTx+N);
        
        
    case  {'array2uvyz','array2uvxy','array2uvxz'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:2*N*pTx),N,pTx),[2,1]);
        varargout{3} = array(2*N*pTx+1:2*N*pTx+N);
        varargout{4} = array(2*N*pTx+N+1:2*N*pTx+2*N);
        
    case  {'array2uvxyz'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:2*N*pTx),N,pTx),[2,1]);
        varargout{3} = array(2*N*pTx+1:2*N*pTx+N);
        varargout{4} = array(2*N*pTx+N+1:2*N*pTx+2*N);
        varargout{5} = array(2*N*pTx+2*N+1:2*N*pTx+3*N);
        
    case  {'array2ux', 'array2uy','array2uz','array2vx','array2vy', 'array2vz'}
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = array(N*pTx+1:N*pTx+N);
        
    case  {'array2uxz','array2uxy','array2uyz','array2vxz','array2vxy','array2vyz'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = array(N*pTx+1:N*pTx+N);
        varargout{3} = array(N*pTx+N+1:N*pTx+N*2);
        
        
    case  {'array2uxyz','array2vxyz'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = array(N*pTx+1:N*pTx+N);
        varargout{3} = array(N*pTx+N+1:N*pTx+N*2);
        varargout{4} = array(N*pTx+N*2+1:end);
        
    case  {'array2y','array2x','array2z'}
        
        array = varargin{1};
        varargout{1} = array;
        
    case  {'array2yz','array2xz','array2xy'}
        
        array = varargin{1};
        N = varargin{2};
        varargout{1} = array(1:N);
        varargout{2} = array(N+1:end);
        
    case  'array2xyz'
        
        array = varargin{1};
        N = varargin{2};
        varargout{1} = array(1:N);
        varargout{2} = array(N+1:2*N);
        varargout{3} = array(2*N+1:end);
        
        %         with s
        
    case  {'uvxyzs2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        i3 = varargin{3};
        i4 = varargin{4};
        i5 = varargin{5};
        i6 = permute(varargin{6},[2,1]);
        
        varargout{1} = [i1(:).',i2(:).',i3(:).',i4(:).',i5(:).',i6(:).'];
        
        
    case  {'uvyzs2array','uvxzs2array','uvxys2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        i3 = varargin{3};
        i4 = varargin{4};
        i5 = permute(varargin{5},[2,1]);
        varargout{1} = [i1(:).',i2(:).',i3(:).',i4(:).',i5(:).'];
        
        
        
    case  {'uvxs2array','uvys2array', 'uvzs2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        i3 = varargin{3};
        i4 = permute(varargin{4},[2,1]);
        varargout{1} = [i1(:).',i2(:).',i3(:).',i4(:).'];
        
    case {'us2array','vs2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        varargout{1} = [i1(:).',i2(:).'];
        
    case  {'uvs2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i2 = permute(varargin{2},[2,1]);
        i3 = permute(varargin{3},[2,1]);
        varargout{1} = [i1(:).',i2(:).',i3(:).'];
        
        
    case  {'uxyzs2array','vxyzs2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i3 = varargin{2};
        i4 = varargin{3};
        i5 = varargin{4};
        i6 = permute(varargin{5},[2,1]);
        varargout{1} = [i1(:).',i3(:).',i4(:).',i5(:).',i6(:).'];
        
        
    case  {'uys2array','uxs2array', 'uzs2array','vys2array', 'vxs2array','vzs2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i3 = varargin{2};
        i6 = permute(varargin{3},[2,1]);
        varargout{1} = [i1(:).',i3(:).',i6(:).'];
        
        
    case  {'uyzs2array','uxzs2array','uxys2array','vxys2array','vyzs2array','vxzs2array'}
        
        i1 = permute(varargin{1},[2,1]);
        i3 = varargin{2};
        i4 = varargin{3};
        i6 = permute(varargin{4},[2,1]);
        varargout{1} = [i1(:).',i3(:).',i4(:).',i6(:).'];
        
    case  {'xs2array','ys2array','zs2array'}
        
        
        i3 = varargin{1};
        i6 = permute(varargin{2},[2,1]);
        
        varargout{1} = [i3(:).',i6(:).'];
        
    case  {'yzs2array','xzs2array', 'xys2array'}
        
        i3 = varargin{1};
        i4 = varargin{2};
        i6 = permute(varargin{3},[2,1]);
        varargout{1} = [i3(:).',i4(:).',i6(:).'];
        
    case  'xyzs2array'
        
        i3 = varargin{1};
        i4 = varargin{2};
        i5 = varargin{3};
        i6 = permute(varargin{4},[2,1]);
        varargout{1} = [i3(:).',i4(:).',i5(:).',i6(:).'];
        % ----------------------
        
    case  {'array2us','array2vs'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        NS = varargin{4};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:end),N,NS),[2,1]);
        
    case  {'array2uvs'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        NS = varargin{4};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:2*N*pTx),N,pTx),[2,1]);
        varargout{3} = permute(reshape(array(2*N*pTx+1:end),N,NS),[2,1]);
    case  {'array2uvxs','array2uvys','array2uvzs'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        NS = varargin{4};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:2*N*pTx),N,pTx),[2,1]);
        varargout{3} = array(2*N*pTx+1:2*N*pTx+N);
        varargout{4} = permute(reshape(array(2*N*pTx+N+1:end),N,NS),[2,1]);
        
    case  {'array2uvyzs','array2uvxys','array2uvxzs'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        NS = varargin{4};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:2*N*pTx),N,pTx),[2,1]);
        varargout{3} = array(2*N*pTx+1:2*N*pTx+N);
        varargout{4} = array(2*N*pTx+N+1:2*N*pTx+2*N);
        varargout{5} = permute(reshape(array(2*N*pTx+2*N+1:end),N,NS),[2,1]);
        
    case  {'array2uvxyzs'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        NS = varargin{4};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = permute(reshape(array(N*pTx+1:2*N*pTx),N,pTx),[2,1]);
        varargout{3} = array(2*N*pTx+1:2*N*pTx+N);
        varargout{4} = array(2*N*pTx+N+1:2*N*pTx+2*N);
        varargout{5} = array(2*N*pTx+2*N+1:2*N*pTx+3*N);
        varargout{6} = permute(reshape(array(2*N*pTx+3*N+1:end),N,NS),[2,1]);
        
    case  {'array2uxs', 'array2uys','array2uzs','array2vxs','array2vys', 'array2vzs'}
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        NS = varargin{4};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = array(N*pTx+1:N*pTx+N);
        varargout{3} = permute(reshape(array(N*pTx+N+1:end),N,NS),[2,1]);
        
    case  {'array2uxzs','array2uxys','array2uyzs','array2vxzs','array2vxys','array2vyzs'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        NS = varargin{4};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = array(N*pTx+1:N*pTx+N);
        varargout{3} = array(N*pTx+N+1:N*pTx+N*2);
        varargout{4} = permute(reshape(array(N*pTx+N*2+1:end),N,NS),[2,1]);
        
    case  {'array2uxyzs','array2vxyzs'}
        
        array = varargin{1};
        N = varargin{2};
        pTx = varargin{3};
        NS = varargin{4};
        varargout{1} = permute(reshape(array(1:N*pTx),N,pTx),[2,1]);
        varargout{2} = array(N*pTx+1:N*pTx+N);
        varargout{3} = array(N*pTx+N+1:N*pTx+N*2);
        varargout{4} = array(N*pTx+N*2+1:N*pTx+N*2+N);
        varargout{5} = permute(reshape(array(N*pTx+N*2+N+1:end),N,NS),[2,1]);
        
    case  {'array2ys','array2xs','array2zs'}
        
        array = varargin{1};
        N = varargin{2};
        NS = varargin{4};
        varargout{1} = array(1:N);
        varargout{2} = permute(reshape(array(N+1:end),N,NS),[2,1]);
        
    case  {'array2yzs','array2xzs','array2xys'}
        
        array = varargin{1};
        N = varargin{2};
        NS = varargin{4};
        varargout{1} = array(1:N);
        varargout{2} = array(N+1:2*N);
        varargout{3} = permute(reshape(array(2*N+1:end),N,NS),[2,1]);
        
    case  'array2xyzs'
        
        array = varargin{1};
        N = varargin{2};
        NS = varargin{4};
        varargout{1} = array(1:N);
        varargout{2} = array(N+1:2*N);
        varargout{3} = array(2*N+1:3*N);
        varargout{4} = permute(reshape(array(3*N+1:end),N,NS),[2,1]);
        
    otherwise
        
        fprintf(2,'REARRANGE: can''t rearrange %s to %s\n',I,O)
        
        
        
end
end

function [array]= ARRAY(Controls,u,v,gx,gy,gz,s)

if ~isempty(Controls)
    switch Controls
        case 'uvxyz'
            array =REARRANGE(Controls,'array',u,v,gx,gy,gz);
        case 'uvyz'
            array =REARRANGE(Controls,'array',u,v,gy,gz);
        case 'uvxz'
            array =REARRANGE(Controls,'array',u,v,gx,gz);
        case 'uvxy'
            array =REARRANGE(Controls,'array',u,v,gx,gy);
        case 'uvy'
            array =REARRANGE(Controls,'array',u,v,gy);
        case 'uvx'
            array =REARRANGE(Controls,'array',u,v,gx);
        case 'uvz'
            array =REARRANGE(Controls,'array',u,v,gz);
        case 'uv'
            array =REARRANGE(Controls,'array',u,v);
        case 'uxyz'
            array =REARRANGE(Controls,'array',u,gx,gy,gz);
        case 'uyz'
            array =REARRANGE(Controls,'array',u,gy,gz);
        case 'uxz'
            array =REARRANGE(Controls,'array',u,gx,gz);
        case 'uxy'
            array =REARRANGE(Controls,'array',u,gx,gy);
        case 'uy'
            array =REARRANGE(Controls,'array',u,gy);
        case 'ux'
            array =REARRANGE(Controls,'array',u,gx);
        case 'uz'
            array =REARRANGE(Controls,'array',u,gz);
        case 'u'
            array =REARRANGE(Controls,'array',u);
        case 'vxyz'%
            array =REARRANGE(Controls,'array',v,gx,gy,gz);
        case 'vyz'
            array =REARRANGE(Controls,'array',v,gy,gz);
        case 'vxz'
            array =REARRANGE(Controls,'array',v,gx,gz);
        case 'vxy'
            array =REARRANGE(Controls,'array',v,gx,gy);
        case 'vy'
            array =REARRANGE(Controls,'array',v,gy);
        case 'vx'
            array =REARRANGE(Controls,'array',v,gx);
        case 'vz'
            array =REARRANGE(Controls,'array',v,gz);
        case 'v'
            array =REARRANGE(Controls,'array',v);
        case 's'
            array =REARRANGE(Controls,'array',s);
        case 'xyz'
            array =REARRANGE(Controls,'array',gx,gy,gz);
        case 'yz'
            array =REARRANGE(Controls,'array',gy,gz);
        case 'xz'
            array =REARRANGE(Controls,'array',gx,gz);
        case 'xy'
            array =REARRANGE(Controls,'array',gx,gy);
        case 'y'
            array =REARRANGE(Controls,'array',gy);
        case 'x'
            array =REARRANGE(Controls,'array',gx);
        case 'z'
            array =REARRANGE(Controls,'array',gz);
            
        case 'uvxyzs'
            array =REARRANGE(Controls,'array',u,v,gx,gy,gz,s);
        case 'uvyzs'
            array =REARRANGE(Controls,'array',u,v,gy,gz,s);
        case 'uvxzs'
            array =REARRANGE(Controls,'array',u,v,gx,gz,s);
        case 'uvxys'
            array =REARRANGE(Controls,'array',u,v,gx,gy,s);
        case 'uvys'
            array =REARRANGE(Controls,'array',u,v,gy,s);
        case 'uvxs'
            array =REARRANGE(Controls,'array',u,v,gx,s);
        case 'uvzs'
            array =REARRANGE(Controls,'array',u,v,gz,s);
        case 'uvs'
            array =REARRANGE(Controls,'array',u,v,s);
        case 'uxyzs'
            array =REARRANGE(Controls,'array',u,gx,gy,gz,s);
        case 'uyzs'
            array =REARRANGE(Controls,'array',u,gy,gz,s);
        case 'uxzs'
            array =REARRANGE(Controls,'array',u,gx,gz,s);
        case 'uxys'
            array =REARRANGE(Controls,'array',u,gx,gy,s);
        case 'uys'
            array =REARRANGE(Controls,'array',u,gy,s);
        case 'uxs'
            array =REARRANGE(Controls,'array',u,gx,s);
        case 'uzs'
            array =REARRANGE(Controls,'array',u,gz,s);
        case 'us'
            array =REARRANGE(Controls,'array',u,s);
        case 'vxyzs'%
            array =REARRANGE(Controls,'array',v,gx,gy,gz,s);
        case 'vyzs'
            array =REARRANGE(Controls,'array',v,gy,gz,s);
        case 'vxzs'
            array =REARRANGE(Controls,'array',v,gx,gz,s);
        case 'vxys'
            array =REARRANGE(Controls,'array',v,gx,gy,s);
        case 'vys'
            array =REARRANGE(Controls,'array',v,gy,s);
        case 'vxs'
            array =REARRANGE(Controls,'array',v,gx,s);
        case 'vzs'
            array =REARRANGE(Controls,'array',v,gz,s);
        case 'vs'
            array =REARRANGE(Controls,'array',v,s);
            
            
        case 'xyzs'
            array =REARRANGE(Controls,'array',gx,gy,gz,s);
        case 'yzs'
            array =REARRANGE(Controls,'array',gy,gz,s);
        case 'xzs'
            array =REARRANGE(Controls,'array',gx,gz,s);
        case 'xys'
            array =REARRANGE(Controls,'array',gx,gy,s);
        case 'ys'
            array =REARRANGE(Controls,'array',gy,s);
        case 'xs'
            array =REARRANGE(Controls,'array',gx,s);
        case 'zs'
            array =REARRANGE(Controls,'array',gz,s);
    end
else
    
    array = [];
end

end

function [u,v,gx,gy,gz,s] = iARRAY(Controls,array,N,pTx,NS,u_,v_,gx_,gy_,gz_,s_)
if ~isempty(Controls)
    switch Controls
        case 'uvxyz'
            [u,v,gx,gy,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvyz'
            [u,v,gy,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvxz'
            [u,v,gx,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvxy'
            [u,v,gx,gy] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvy'
            [u,v,gy] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvx'
            [u,v,gx] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvz'
            [u,v,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uv'
            [u,v] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uxyz'
            [u,gx,gy,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uyz'
            [u,gy,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uxz'
            [u,gx,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uxy'
            [u,gx,gy] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uy'
            [u,gy] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'ux'
            [u,gx] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uz'
            [u,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'u'
            [u] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'vxyz'%
            [v,gx,gy,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vyz'
            [v,gy,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vxz'
            [v,gx,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vxy'
            [v,gx,gy] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vy'
            [v,gy] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vx'
            [v,gx] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vz'
            [v,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'v'
            [v] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
            
        case 'xyz'
            [gx,gy,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'yz'
            [gy,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'xz'
            [gx,gz] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'xy'
            [gx,gy] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'y'
            [gy] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'x'
            [gx] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'z'
            [gz] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 's'
            [s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvxyzs'
            [u,v,gx,gy,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvyzs'
            [u,v,gy,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvxzs'
            [u,v,gx,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvxys'
            [u,v,gx,gy,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvys'
            [u,v,gy,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvxs'
            [u,v,gx,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvzs'
            [u,v,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uvs'
            [u,v,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'uxyzs'
            [u,gx,gy,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uyzs'
            [u,gy,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uxzs'
            [u,gx,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uxys'
            [u,gx,gy,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uys'
            [u,gy,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uxs'
            [u,gx,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'uzs'
            [u,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'us'
            [u,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            v = v_;
        case 'vxyzs'%
            [v,gx,gy,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vyzs'
            [v,gy,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vxzs'
            [v,gx,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vxys'
            [v,gx,gy,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vys'
            [v,gy,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vxs'
            [v,gx,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vzs'
            [v,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
        case 'vs'
            [v,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            u = u_;
            
        case 'xyzs'
            [gx,gy,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'yzs'
            [gy,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'xzs'
            [gx,gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'xys'
            [gx,gy,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'ys'
            [gy,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'xs'
            [gx,s] = REARRANGE('array',Controls,array,N,pTx,NS);
        case 'zs'
            [gz,s] = REARRANGE('array',Controls,array,N,pTx,NS);
            
    end
end
if ~exist('u','var')
    u = u_;
end
if ~exist('v','var')
    v = v_;
end

if ~exist('gx','var')
    gx = gx_;
end
if ~exist('gy','var')
    gy = gy_;
end
if ~exist('gz','var')
    gz = gz_;
end
if ~exist('s','var')
    s = s_;
end
end

function array_scaled = ARRAY2PHYS(array,arraysc)

array_scaled = array.*arraysc;

end

function array_scaled = ARRAY2NORM(array,arraysc)

array_scaled = array./arraysc;

end

function opt = ALLOCATE(spc,khr,opt)

if opt.OptNum >= 1
    
    if ~isreal(spc.B1map)
        opt.sr = real(spc.B1map);
        opt.si = imag(spc.B1map);
    else
        opt.sr = spc.B1map;
        opt.si = zeros(size(spc.B1map));
    end
    
    
    
    opt.Durations = zeros(opt.MaxIter+1,1);
    opt.Fun = zeros(opt.MaxIter+1,1);
    opt.dFun = nan(opt.MaxIter+1,1);
    
    
    opt.w0 = spc.w0map+spc.V*2*pi;
    opt.yX = spc.X.*khr.gamma;
    opt.yY = spc.Y.*khr.gamma;
    opt.yZ = spc.Z.*khr.gamma;
    if isfield(spc,'Smap')
        opt.yS = spc.Smap.*khr.gamma;
    else
        opt.yS = [];
    end
    
    opt.Gm = khr.GmHW*khr.Gmpct/100;
    opt.Sm = khr.SmHW*khr.Smpct/100;
    
    
    if isfield(opt,'mem') % the user can avoid a memory problem with huge array by specifying 'opt.mem = 'low''. Then only the start and final conditions in M_t and L_t, respectively, will be output.
        if strcmp(opt.mem,'low')
            opt.M_t = spc.M0;
            opt.L_t = spc.Md;
        else
            opt.M_t = zeros(3*spc.P,opt.N+1);
            opt.M_t(:,1:opt.mon(1)) = repmat(spc.M0,1,opt.mon(1));
            opt.L_t = zeros(3*spc.P,opt.N+1);
            opt.L_t(:,end-opt.mon(end)+1:end) = repmat(spc.Md,1,opt.mon(end));
            
            
        end
    else
        opt.M_t = zeros(3*spc.P,opt.N+1);
        opt.M_t(:,1:opt.mon(1)) = repmat(spc.M0,1,opt.mon(1));
        opt.L_t = zeros(3*spc.P,opt.N+1);
        opt.L_t(:,end-opt.mon(end)+1:end) = repmat(spc.Md,1,opt.mon(end));
        
    end
    
    
    
    
    
    
    
    
    
    
end
end

function PRINT(opt,Type,Fun,con)

if nargin == 2
    Type = 'all';
end
ComStrLen = 8;
Iterlen = max(4,length(sprintf('%i',opt.MaxIter)));
Str_Iter = PADSTR('Iter',opt.k-1,Iterlen,[]);
Str_Fun = PADSTR('Fun',Fun,20,[]);

Str_peak_u = PADSTR('',con.rat_peak_u,ComStrLen,opt.Par.Constr.peak_uv.Violu);
Str_peak_v = PADSTR('',con.rat_peak_v,ComStrLen,opt.Par.Constr.peak_uv.Violv);
Str_ave_uv = PADSTR('',con.rat_ave_uv,ComStrLen,opt.Par.Constr.ave_uv.Violuv);
Str_peak_x = PADSTR('',con.rat_peak_x,ComStrLen,opt.Par.Constr.peak_xyz.Violx);
Str_peak_y = PADSTR('',con.rat_peak_y,ComStrLen,opt.Par.Constr.peak_xyz.Violy);
Str_peak_z = PADSTR('',con.rat_peak_z,ComStrLen,opt.Par.Constr.peak_xyz.Violz);
Str_slew_x = PADSTR('',con.rat_slew_x,ComStrLen,opt.Par.Constr.slew_xyz.Violx);
Str_slew_y = PADSTR('',con.rat_slew_y,ComStrLen,opt.Par.Constr.slew_xyz.Violy);
Str_slew_z = PADSTR('',con.rat_slew_z,ComStrLen,opt.Par.Constr.slew_xyz.Violz);

Lab_Iter = PADSTR('Label','Iter',length(sprintf('%i',opt.MaxIter)),[]);
Lab_Fun = PADSTR('Label','J',20,[]);

Lab_peak_u = PADSTR('Label',[char(hex2dec('22BC')),'u'],ComStrLen,[]);
Lab_peak_v = PADSTR('Label',[char(hex2dec('22BC')),'v'],ComStrLen,[]);
Lab_ave_uv = PADSTR('Label',[char(hex2dec('03BC')),'uv'],ComStrLen,[]);
Lab_peak_x = PADSTR('Label',[char(hex2dec('22BC')),'Gx'],ComStrLen,[]);
Lab_peak_y = PADSTR('Label',[char(hex2dec('22BC')),'Gy'],ComStrLen,[]);
Lab_peak_z = PADSTR('Label',[char(hex2dec('22BC')),'Gz'],ComStrLen,[]);
Lab_slew_x = PADSTR('Label',[char(hex2dec('0394')),'Gx/',char(hex2dec('0394')),'t'],ComStrLen,[]);
Lab_slew_y = PADSTR('Label',[char(hex2dec('0394')),'Gy/',char(hex2dec('0394')),'t'],ComStrLen,[]);
Lab_slew_z = PADSTR('Label',[char(hex2dec('0394')),'Gz/',char(hex2dec('0394')),'t'],ComStrLen,[]);

% here it would be nice to have everything printed with constant table
% column widths
Tab = ' ';

switch Type
    case 'all'
        
        if opt.k == 1 || rem(opt.k-1,50) == 0
            
            fprintf('%s\n',[Lab_Iter,Tab,Lab_Fun,Tab,Tab,Tab,Lab_peak_u,Tab,Lab_peak_v,Tab,Lab_ave_uv,Lab_peak_x,Tab,Lab_peak_y,Tab,Lab_peak_z,Tab,Lab_slew_x,Tab,Lab_slew_y,Tab,Lab_slew_z,Tab]);
            
            fprintf('%s\n',[Str_Iter,Tab,Str_Fun,Tab,Tab,Tab,Str_peak_u,Tab,Str_peak_v,Tab,Str_ave_uv,Tab,Str_peak_x,Tab,Str_peak_y,Tab,Str_peak_z,Tab,Str_slew_x,Tab,Str_slew_y,Tab,Str_slew_z,Tab]);
            
            
            
        else
            fprintf('%s\n',[Str_Iter,Tab,Str_Fun,Tab,Tab,Tab,Str_peak_u,Tab,Str_peak_v,Tab,Str_ave_uv,Tab,Str_peak_x,Tab,Str_peak_y,Tab,Str_peak_z,Tab,Str_slew_x,Tab,Str_slew_y,Tab,Str_slew_z,Tab]);
            
        end
        
        
    case 'con'
        fprintf('%s\n',[repmat(' ',1,length(Tab)*3+length(Lab_Iter)+length(Lab_Fun)),Tab,Str_peak_u,Tab,Str_peak_v,Tab,Str_ave_uv,Tab,Str_peak_x,Tab,Str_peak_y,Tab,Str_peak_z,Tab,Str_slew_x,Tab,Str_slew_y,Tab,Str_slew_z,Tab]);
        
end

end

function nP = NEWPAR(uP,dP)

if isempty(uP)
    nP = dP;
else
    dPnames = RECURSIVEFIELDNAMES(dP);
    uPnames = RECURSIVEFIELDNAMES(uP);
    nP = struct;
    
    % first find the user parameters amongst parameters with default
    % parameters and write those into the new parameters struct.
    % those default parameters that are not user specified will be adopted
    % to the new parameters as well
    
    N = length(dPnames);
    for n = 1:N
        curName = dPnames{n};
        test = find(strcmp(curName,uPnames), 1);
        %         test = find(strcmp(curName,uPnames));
        if isempty(test)
            eval(['nP.',curName,'=dP.',curName,';'])
        else
            eval(['nP.',curName,'=uP.',curName,';'])
        end
        
    end
    
    % secondly, keep those user parameters that do not have a default
    % parameter in this function call, and adopt them to the new parameter
    % struct too.
    N = length(uPnames);
    
    for n = 1:N
        curName = uPnames{n};
        test = find(strcmp(curName,dPnames), 1);
        if isempty(test)
            eval(['nP.',curName,'=uP.',curName,';'])
        end
        
    end
end
end

function O = RECURSIVEFIELDNAMES(I)
f=0; q=1; U={};
O=cellfun(@(x)strcat('I.',x),fieldnames(I),'UniformOutput',0);

while q~=0
    H={}; q=length(O);
    for i=1:length(O)
        if isstruct(eval(O{i}))==1
            if f~=1
                A=fieldnames(eval(O{i}));
                A=cellfun(@(x)strcat(sprintf('%s.',O{i}),x),A,'UniformOutput',0);
                H=cat(1,H,A);
            elseif f==1
                H=cat(1,H,O{i});
                q=q-1;
            end
            U=cat(1,U,O{i});
        else
            H=cat(1,H,O{i}); q=q-1;
        end
    end
    O = H;
end
O=cellfun(@(x)sprintf('%s',x(3:end)),O,'UniformOutput',0);
end

function [Controls,PosControls,FixControls] = CONTROLS(In)
p = inputParser;
test = false;
try p.addRequired('Controls',@(x)VALIDATECONTROLS(x));
    test = true;
catch  me;fprintf(2,['CONTROLS: ',me.message,'\n']);
end

if test
    
    test = false; %#ok<*NASGU>
    try p.parse(In);
        test = true;
    catch  me; fprintf(2,['CONTROLS: ',me.message,'\n']);
    end
end

[test,Controls,PosControls,FixControls]= VALIDATECONTROLS(In); %#ok<*ASGLU>

end

function [test,y,z,h]= VALIDATECONTROLS(x)


PosControls = {'uvxyz','uvyz','uvxz','uvxy','uvy','uvx','uvz','uv','xyz','yz','xz','xy','y','x','z','',...
    };


test = false;
if ischar(x)
    switch x
        case PosControls
            test = true;
            
            
            
            if nargout >=2
                y = x;
                if strcmp(y,'')
                    y = PosControls{1};
                end
            end
            if nargout >= 3
                z = PosControls;
            end
            if nargout >= 4
                h = CONTROLSFIX(x);
            end
            
            
        otherwise
            if nargout == 0
                fprintf(2,'The Controls parameter should be any of:\n')
                for n = 1:length(PosControls)-1
                    fprintf(2,'%s,',PosControls{n})
                end
                fprintf(2,'%s\n',PosControls{length(PosControls)})
                
            end
            if nargout >= 2
                y = [];
            end
            if nargout >= 3
                z = [];
            end
            if nargout >= 4
                h = [];
            end
    end
end




end

function FixControls = CONTROLSFIX(Controls)
% Fixed for s by simply replicating the previous without s and then adding
% s. Don't know if this works actually.
switch Controls
    case {'xyz','yz','xz','xy','y','x','z'}
        % a combination where some gradient channel is not needed doesn't
        % put anything to the fixed array because no limits or constraint
        % will need the missing gradient channel.
        FixControls = [];
        
    case {'uvxyz','uvxz','uvxy','uvy','uvx','uvz','uv','uvyz'}
        
        FixControls = [];
        
        
        
end
end

function Str = PADSTR(Type,Val,TotLen,Viol)
Tablespacing = '  ';
switch Type
    case 'Iter'
        Strtemp = sprintf('%i',Val);
        Len = length(Strtemp);
        if Len < TotLen
            Str = [sprintf('%s',repmat(' ',1,TotLen-Len)),Strtemp];
        elseif Len == TotLen
            
            Str = Strtemp;
        else
            Str = Strtemp;
        end
        %         Str = [Str];
    case 'Fun'
        Strtemp = sprintf('%0.16f',Val);
        Len = length(Strtemp);
        if Len < TotLen
            Str = [sprintf('%s',repmat(' ',1,TotLen-Len)),Strtemp];
        elseif Len == TotLen
            
            Str = Strtemp;
        else
            Str = Strtemp;
        end
        Str = [Str];
    case 'Label'
        if ~ischar(Val)
            error('Label must be char')
        end
        
        
        
        Len = length(Val);
        Strtemp = Val;
        if Len < TotLen
            Str = [sprintf('%s',repmat(' ',1,TotLen-Len)),Strtemp];
        elseif Len == TotLen
            
            Str = Strtemp;
        else
            Str = Strtemp;
        end
        
        %         Str = [Str];
        
    otherwise
        Strtemp = sprintf('%0.3f',Val);
        Len = length(Strtemp);
        
        if Len < TotLen
            if Viol == 1
                Strtemp = ['*', Strtemp];
                Str = [repmat(' ',1,TotLen-Len-1),Strtemp];
            else
                Str = [repmat(' ',1,TotLen-Len),Strtemp];
            end
            
            
        elseif Len == TotLen
            Strtemp = sprintf('%0.3f',Val);
            Len = length(Strtemp);
            if Viol == 1
                Strtemp = ['*', Strtemp];
                Str = [repmat(' ',1,TotLen-Len-1),Strtemp];
            else
                Str = [repmat(' ',1,TotLen-Len),Strtemp];
            end
            Str = Strtemp;
        else
            Str = repmat('X',1,TotLen);
        end
        
        
        
end
end

function ksafe = KSAFE(ksafe,Copes,Viols)

N = length(Copes);
ksafe_ = ones(N,1).*ksafe;

for n = 1:N
    if strcmp(Copes{n},'Stop') && Viols(n) > 0
        ksafe_(n) = ksafe_(n)-1;
        if ksafe_(n) < 1
            ksafe_(n) = 1;
        end
    end
end
ksafe = min(ksafe_);
end

function [wx,wy,wxo,wyo] = RF0s(pTx,Non,mon,N,K)

wx_ = zeros(pTx,Non);
wy_ = zeros(pTx,Non);

wx = zeros(pTx,N);
wy = zeros(pTx,N);

wx(:,mon) = wx_;
wy(:,mon) = wy_;


wxo = zeros(pTx,N,K+1);
wyo = zeros(pTx,N,K+1);
wxo(:,:,1) = wx;
wyo(:,:,1) = wy;

end

function [val,wei] = RFROBUSTNESS(Num,Min,Max,Dist,FWHM)


if abs(abs(Max)-abs(Min)) <= eps
    fprintf('RFrobu: RF robustness limits equal, applying eps\n');
    Max = Max + eps;
    Min = Min - eps;
end
if Num == 1
    val = 1;
    wei = 1;
    %     return
else
    val = linspace(Min,Max,Num);
end

switch Dist
    case 'lorentzian'
        % 1./(pi./FWHM/2).*(FWHM/2)^2.
        wei = 1./(1+((val-1)/(FWHM/2)).^2);
        wei = wei./sum(wei);
        
        fprintf('RFrobu:');
        fprintf('    RF inhomogeneity distribution:     %s\n',Dist);
        fprintf('    Number of RF inhomogeneity points: %i\n',Num);
        fprintf('    RF inhomogeneity FWHM [%%]:         %f\n',FWHM*100);
        
        
    case 'evenly'
        
        wei = ones(1,Num)./Num;
        
        fprintf('RFrobu:');
        fprintf('    RF inhomogeneity distribution:     %s\n',Dist);
        fprintf('    Number of RF inhomogeneity points: %i\n',Num);
        fprintf('    RF inhomogeneity minimum [%%]:      %f\n',Min*100);
        fprintf('    RF inhomogeneity maximum [%%]:      %f\n',Max*100);
    otherwise
        wei = ones(1,Num)./Num;
        fprintf('RFrobu: RF robustness weighting distribution unknown, assigning a distribution with even weights\n');
        fprintf('RFrobu:');
        fprintf('    RF inhomogeneity distribution:     %s\n',Dist);
        fprintf('    Number of RF inhomogeneity points: %i\n',Num);
        fprintf('    RF inhomogeneity minimum [%%]:      %f\n',Min*100);
        fprintf('    RF inhomogeneity maximum [%%]:      %f\n',Max*100);
end
fprintf('    Value [%%%%]\t\tWeight\n');
for n=1:Num
    fprintf('    %1.4f\t\t%1.4f\n',val(n)*100,wei(n));
end


end
%% Math functions and array functions
function [res1,res2] = SPDIAGS_old(arg1,arg2,arg3,arg4)
%SPDIAGS Sparse matrix formed from diagonals.
%   SPDIAGS, which generalizes the function "diag", deals with three
%   matrices, in various combinations, as both input and output.
%
%   [B,d] = SPDIAGS(A) extracts all nonzero diagonals from the m-by-n
%   matrix A.  B is a min(m,n)-by-p matrix whose columns are the p
%   nonzero diagonals of A.  d is a vector of length p whose integer
%   components specify the diagonals in A.
%
%   B = SPDIAGS(A,d) extracts the diagonals specified by d.
%   A = SPDIAGS(B,d,A) replaces the diagonals of A specified by d with
%       the columns of B.  The output is sparse.
%   A = SPDIAGS(B,d,m,n) creates an m-by-n sparse matrix from the
%       columns of B and places them along the diagonals specified by d.
%
%   Roughly, A, B and d are related by
%       for k = 1:p
%           B(:,k) = diag(A,d(k))
%       end
%
%   Example: These commands generate a sparse tridiagonal representation
%   of the classic second difference operator on n points.
%       e = ones(n,1);
%       A = spdiags([e -2*e e], -1:1, n, n)
%
%   Some elements of B, corresponding to positions "outside" of A, are
%   not actually used.  They are not referenced when B is an input and
%   are set to zero when B is an output.  If a column of B is longer than
%   the diagonal it's representing, elements of super-diagonals of A
%   correspond to the lower part of the column of B, while elements of
%   sub-diagonals of A correspond to the upper part of the column of B.
%
%   Example: This uses the top of the first column of B for the second
%   sub-diagonal and the bottom of the third column of B for the first
%   super-diagonal.
%       B = repmat((1:n)',1,3);
%       S = spdiags(B,[-2 0 1],n,n);
%
%   See also DIAG, SPEYE.

%   Rob Schreiber
%   Copyright 1984-2016 The MathWorks, Inc.


if nargin <= 2
    % Extract diagonals
    A = arg1;
    if nargin == 1
        % Find all nonzero diagonals
        [i,j] = find(A);
        % Compute d = unique(d) without extra function call
        d = sort(j-i);
        d = d(diff([-inf; d(:)])~=0);
        d = d(:);
    else
        % Diagonals are specified
        d = arg2(:);
    end
    [m,n] = size(A);
    p = length(d);
    B = zeros(min(m,n),p,class(A));
    for k = 1:p
        if m >= n
            i = max(1,1+d(k)):min(n,m+d(k));
        else
            i = max(1,1-d(k)):min(m,n-d(k));
        end
        B(i,k) = DIAGK_old(A,d(k));
    end
    res1 = B;
    res2 = d;
end

if nargin >= 3
    B = arg1;
    d = arg2(:);
    p = length(d);
    if nargin == 3 % Replace specified diagonals
        A = arg3;
    else           % Create new matrix with specified diagonals
        A = sparse(arg3, arg4);
    end
    [m,n] = size(A);
    
    % Check size of matrix B (should be min(m,n)-by-p)
    % For backwards compatibility, only error if the code would
    % previously have errored out in the indexing expression.
    maxIndexRows = max(max(1,1-d), min(m,n-d)) + (m>=n)*d;
    maxIndexRows(max(1,1-d) > min(m,n-d)) = 0;
    if any(maxIndexRows > size(B, 1)) || p > size(B, 2)
        if nargin == 3
            error(message('MATLAB:spdiags:InvalidSizeBThreeInput'));
        else
            error(message('MATLAB:spdiags:InvalidSizeBFourInput'));
        end
    end
    
    % Compute indices and values of sparse matrix with given diagonals
    
    % Compute lengths of diagonals:
    len = max(0, min(m, n-d) - max(1, 1-d) + 1);
    len = [0; cumsum(len)];
    
    a = zeros(len(p+1), 3);
    for k = 1:p
        % Append new d(k)-th diagonal to compact form
        i = (max(1,1-d(k)):min(m,n-d(k)))';
        a((len(k)+1):len(k+1),:) = [i i+d(k) B(i+(m>=n)*d(k),k)];
    end
    
    % Remove diagonal elements in old matrix if necessary
    if nnz(A) ~= 0
        % Process A in compact form
        [i,j,aold] = find(A);
        aold = [i(:) j(:) aold(:)]; % need (:) if A is row vector
        
        % Delete current d(k)-th diagonal, k=1,...,p
        i = any((aold(:, 2) - aold(:, 1)) == d', 2);
        aold(i, :) = [];
        
        % Combine new diagonals and non-diagonal entries of original matrix
        a = [a; aold];
    end
    
    res1 = sparse(a(:,1),a(:,2),a(:,3),m,n);
    if islogical(A) || islogical(B)
        res1 = (res1~=0);
    end
end
end

function D = DIAGK_old(X,k)
% DIAGK  K-th matrix diagonal.
% DIAGK(X,k) is the k-th diagonal of X, even if X is a vector.
if ~isvector(X)
    D = diag(X,k);
    D = D(:);  %Ensure column vector is returned for empty X.
else
    if ~isempty(X) && 0 <= k && 1+k <= size(X,2)
        D = X(1+k);
    elseif ~isempty(X) && k < 0 && 1-k <= size(X,1)
        D = X(1-k);
    else
        D = zeros(0,1,'like',X);
    end
end
end

function rho=EXPV_IK(L,rho,dt)
% IK: Ilya Kuprov, SOTON
% Make sure rho is full
if issparse(rho), rho=full(rho); end

% Determine the number of time steps
norm_mat=norm(L,1)*abs(dt);
nsteps=ceil(norm_mat/2);

% Estimate the scaling coefficient
scaling=max([1 norm(rho,1)]);

% Scale the vector
rho=rho/scaling;

% Run the Krylov procedure
for n=1:nsteps
    next_term=rho; k=1;
    while nnz(abs(next_term)>eps)>0
        next_term=(dt/(k*nsteps))*(L*next_term);
        rho=rho+next_term; k=k+1;
    end
    if k>32, warning(['loss of accuracy in EXPV_IK(), k=' num2str(k) ...
            ', use evolution() instead']); end
end

% Scale the vector back
rho=scaling*rho;

end

%% Constraints

function opt = CONSTRAINHOW(opt)




opt.Par.Constr.sar_l.LimCor = NaN;
opt.Par.Constr.sar_l.UnitCor = 'W/kg';
opt.Par.Constr.sar_l.NQ = 0;
opt.Par.Constr.sar_l.Nvop = 0;

opt.Par.Constr.sar_g.LimCor = NaN;
opt.Par.Constr.sar_g.UnitCor = 'W/kg';



%%
% RF peak
if ~strcmp(opt.Par.Constr.peak_uv.Cope,'Ignore')
    switch opt.Par.Constr.peak_uv.Unit
        
        case 'W'
            
            opt.Par.Constr.peak_uv.Power = 2;
            
            if strcmp(opt.B1_nom_amp_unit,'T/V')
                
                opt.Par.Constr.peak_uv.Unit_cor_factor = 8*50*opt.gamma^2/opt.f_B1_val^2; % (rad/s)^2/W
                
                opt.Par.Constr.peak_uv.LimCor = opt.Par.Constr.peak_uv.Lim*opt.Par.Constr.peak_uv.Unit_cor_factor;
                opt.Par.Constr.peak_uv.UnitCor = '(rad/s)^2';
                
                
                opt.Par.Constr.peak_uv.Unit_cor_factor_VV_2_W = opt.gamma^2*opt.B1_nom_amp_val^2./opt.Par.Constr.peak_uv.Unit_cor_factor;
            else
                fprintf('%s: Par.Constr.peak_uv.Unit''s value (%s) is incompatible with B1 map unit (%s)\n',mfilename,opt.Par.Constr.peak_uv.Unit,opt.B1_nom_amp_unit);
                
            end
            
        case 'V'
            
            opt.Par.Constr.peak_uv.Power = 1;
            if strcmp(opt.B1_nom_amp_unit,'T/V')
                
                opt.Par.Constr.peak_uv.Unit_cor_factor = opt.gamma/opt.f_B1_val;
                
                opt.Par.Constr.peak_uv.LimCor = opt.Par.Constr.peak_uv.Lim*opt.Par.Constr.peak_uv.Unit_cor_factor;
                opt.Par.Constr.peak_uv.UnitCor = 'rad/s';
            else
                fprintf('%s: Par.Constr.peak_uv.Unit''s value (%s) is incompatible with B1 map unit (%s)\n',mfilename,opt.Par.Constr.peak_uv.Unit,opt.B1_nom_amp_unit);
                
            end
        case 'rad/s'
            
            opt.Par.Constr.peak_uv.Power = 1;
            
            opt.Par.Constr.peak_uv.Unit_cor_factor = 1;
            
            opt.Par.Constr.peak_uv.LimCor = opt.Par.Constr.peak_uv.Lim*opt.Par.Constr.peak_uv.Unit_cor_factor;
            opt.Par.Constr.peak_uv.UnitCor = 'rad/s';
        case '(rad/s)^2'
            
            opt.Par.Constr.peak_uv.Power = 2;
            
            opt.Par.Constr.peak_uv.Unit_cor_factor = 1;
            
            opt.Par.Constr.peak_uv.LimCor = opt.Par.Constr.peak_uv.Lim*opt.Par.Constr.peak_uv.Unit_cor_factor;
            opt.Par.Constr.peak_uv.UnitCor = '(rad/s)^2';
        case 'T'
            
            opt.Par.Constr.peak_uv.Power = 1;
            
            opt.Par.Constr.peak_uv.Unit_cor_factor = opt.gamma;
            opt.Par.Constr.peak_uv.LimCor = opt.Par.Constr.peak_uv.Lim*opt.Par.Constr.peak_uv.Unit_cor_factor;
            opt.Par.Constr.peak_uv.UnitCor = 'rad/s';
        case 'Hz'
            
            opt.Par.Constr.peak_uv.Power = 1;
            opt.Par.Constr.peak_uv.Unit_cor_factor = 2*pi;
            opt.Par.Constr.peak_uv.LimCor = opt.Par.Constr.peak_uv.Lim*opt.Par.Constr.peak_uv.Unit_cor_factor;
            opt.Par.Constr.peak_uv.UnitCor = 'rad/s';
    end
    if strcmp(opt.Par.Constr.peak_uv.Cope,'Constr')
        %         if opt.Par.Constr.peak_uv.Power == 1
        %             opt.w1m = opt.Par.Constr.peak_uv.LimCor./sqrt(2);
        %         else
        %             opt.w1m = sqrt(opt.Par.Constr.peak_uv.LimCor)./sqrt(2);
        %         end
        if opt.Par.Constr.peak_uv.Power == 1
            opt.w1m = opt.Par.Constr.peak_uv.LimCor;
        else
            opt.w1m = sqrt(opt.Par.Constr.peak_uv.LimCor);
        end
    end
else
    opt.Par.Constr.peak_uv.Power = NaN;
    opt.Par.Constr.peak_uv.Unit_cor_factor = NaN;
    opt.Par.Constr.peak_uv.LimCor = NaN;
    opt.Par.Constr.peak_uv.UnitCor = 'NaN';
end


if ~isfield(opt,'w1m')
    opt.w1m = 2*pi*10e3;
    
end

if ~strcmp(opt.Par.Constr.ave_uv.Cope,'Ignore')
    if strcmp(opt.Par.Constr.ave_uv.Type,'nlc')
        
        switch opt.Par.Constr.ave_uv.Unit
            
            case 'W'
                
                opt.Par.Constr.ave_uv.Power = 2;
                
                if strcmp(opt.B1_nom_amp_unit,'T/V')
                    
                    
                    opt.Par.Constr.ave_uv.Unit_cor_factor = 8*50*opt.gamma^2/opt.f_B1_val^2;
                    
                    opt.Par.Constr.ave_uv.LimCor = opt.Par.Constr.ave_uv.Lim*opt.Par.Constr.ave_uv.Unit_cor_factor./opt.Par.Constr.Dutycycle;
                    opt.Par.Constr.ave_uv.UnitCor = '(rad/s)^2';
                    
                else
                    fprintf('%s: AvLim''s unit (%s) is incompatible with B1 map unit (%s)\n',mfilename,opt.Par.Constr.ave_uv.Unit,opt.B1_nom_amp_unit);
                    
                end
                
            case '(rad/s)^2'
                
                opt.Par.Constr.ave_uv.Power = 2;
                
                opt.Par.Constr.ave_uv.Unit_cor_factor = 1;
                
                opt.Par.Constr.ave_uv.LimCor = opt.Par.Constr.ave_uv.Lim*opt.Par.Constr.ave_uv.Unit_cor_factor./opt.Par.Constr.Dutycycle;
                opt.Par.Constr.ave_uv.UnitCor = '(rad/s)^2';
            case 'T^2'
                
                opt.Par.Constr.ave_uv.Power = 2;
                opt.Par.Constr.ave_uv.Unit_cor_factor = (opt.gamma)^2;
                opt.Par.Constr.ave_uv.LimCor = opt.Par.Constr.ave_uv.Lim*opt.Par.Constr.ave_uv.Unit_cor_factor./opt.Par.Constr.Dutycycle;
                opt.Par.Constr.ave_uv.UnitCor = '(rad/s)^2';
                
            case 'Hz^2'
                
                opt.Par.Constr.ave_uv.Power = 2;
                opt.Par.Constr.ave_uv.Unit_cor_factor = (2*pi)^2;
                opt.Par.Constr.ave_uv.LimCor = opt.Par.Constr.ave_uv.Lim*opt.Par.Constr.ave_uv.Unit_cor_factor./opt.Par.Constr.Dutycycle;
                opt.Par.Constr.ave_uv.UnitCor = '(rad/s)^2';
            otherwise
                fprintf('%s: The ave_uv unit %s is not supported currently\n',mfilename,opt.Par.Constr.ave_uv.Unit);
                return
        end
    else
        opt.Par.Constr.ave_uv.LimCor = NaN;
        opt.Par.Constr.ave_uv.UnitCor = '(rad/s)^2';
        opt.Par.Constr.ave_uv.Unit_cor_factor = 1;
        opt.Par.Constr.ave_uv.Power = 2;
    end
else
    opt.Par.Constr.ave_uv.LimCor = NaN;
    opt.Par.Constr.ave_uv.UnitCor = '(rad/s)^2';
    opt.Par.Constr.ave_uv.Unit_cor_factor = 1;
    opt.Par.Constr.ave_uv.Power = 2;
end

%%
% G peak
if strcmp(opt.Par.Constr.peak_xyz.Type,'nlc') || strcmp(opt.Par.Constr.peak_xyz.Type,'bnd')
    
    switch opt.Par.Constr.peak_xyz.Unit
        
        
        case 'T/m'
            
            opt.Par.Constr.peak_xyz.Power = 2;
            
            opt.Par.Constr.peak_xyz.Unit_cor_factor = opt.Par.Constr.peak_xyz.Lim;
            opt.Par.Constr.peak_xyz.LimCor = opt.Par.Constr.peak_xyz.Lim*opt.Par.Constr.peak_xyz.Unit_cor_factor;
            opt.Par.Constr.peak_xyz.UnitCor = '(T/m)^2';
        otherwise
            fprintf('%s: The peak_xyz unit %s is not supported currently\n',mfilename,opt.Par.Constr.peak_xyz.Unit);
            return
            
    end
    
    
else
    opt.Par.Constr.peak_xyz.LimCor = NaN;
    opt.Par.Constr.peak_xyz.UnitCor = '(T/m)^2';
end

%%
% G slewrate
if strcmp(opt.Par.Constr.slew_xyz.Type,'nlc')
    
    switch opt.Par.Constr.slew_xyz.Unit
        
        
        case 'T/m/s'
            
            opt.Par.Constr.slew_xyz.Power = 2;
            
            opt.Par.Constr.slew_xyz.Unit_cor_factor = opt.Par.Constr.slew_xyz.Lim;
            opt.Par.Constr.slew_xyz.LimCor = opt.Par.Constr.slew_xyz.Lim*opt.Par.Constr.slew_xyz.Unit_cor_factor;
            opt.Par.Constr.slew_xyz.UnitCor = '(T/m/s)^2';
        otherwise
            fprintf('%s: The slew_xyz unit %s is not supported currently\n',mfilename,opt.Par.Constr.slew_xyz.Unit);
            return
    end
    
    
else
    opt.Par.Constr.slew_xyz.LimCor = NaN;
    opt.Par.Constr.slew_xyz.UnitCor = '(T/m/s)^2';
end


%%




%%

fprintf('\n%s:\n\n',mfilename);

fprintf('\tPeak uv (%s)\n',...
    opt.Par.Constr.peak_uv.Cope);

fprintf('\t\tLimit: %1.2e %s ~ %1.2e %s\n',...
    opt.Par.Constr.peak_uv.Lim,opt.Par.Constr.peak_uv.Unit,opt.Par.Constr.peak_uv.LimCor,opt.Par.Constr.peak_uv.UnitCor);

fprintf('\tAve uv (%s)\n',...
    opt.Par.Constr.ave_uv.Cope);

fprintf('\t\tLimit w. dutycycle (%f): %1.2e/%f %s ~ %1.2e %s ~ %1.2e %s\n',...
    opt.Par.Constr.Dutycycle,...
    opt.Par.Constr.ave_uv.Lim,opt.Par.Constr.Dutycycle,opt.Par.Constr.ave_uv.Unit,...
    opt.Par.Constr.ave_uv.Lim/opt.Par.Constr.Dutycycle,opt.Par.Constr.ave_uv.Unit,...
    opt.Par.Constr.ave_uv.LimCor,opt.Par.Constr.ave_uv.UnitCor);


fprintf('\tPeak g_xyz (%s)\n',...
    opt.Par.Constr.peak_xyz.Cope);

fprintf('\t\tLimit: %1.2e %s ~ %1.2e %s\n',...
    opt.Par.Constr.peak_xyz.Lim,opt.Par.Constr.peak_xyz.Unit,opt.Par.Constr.peak_xyz.LimCor,opt.Par.Constr.peak_xyz.UnitCor);

fprintf('\tSlewrate g_xyz (%s)\n',...
    opt.Par.Constr.slew_xyz.Cope);

fprintf('\t\tLimit: %1.2e %s ~ %1.2e %s\n',...
    opt.Par.Constr.slew_xyz.Lim,opt.Par.Constr.slew_xyz.Unit,opt.Par.Constr.slew_xyz.LimCor,opt.Par.Constr.slew_xyz.UnitCor);


end

function opt = CONSTRAINTSANITY(opt)

opt.Par.Constr.sar_l.PossibleType = {'nlc'};
opt.Par.Constr.sar_l.PossibleCope = {'Ignore'};

opt.Par.Constr.sar_g.PossibleType = {'nlc'};
opt.Par.Constr.sar_g.PossibleCope = {'Ignore'};


opt.Par.Constr.peak_uv.PossibleType = {'bnd'};
opt.Par.Constr.peak_uv.PossibleCope = {'Ignore','Monitor','Stop','Constr'};




opt.Par.Constr.ave_uv.PossibleType = {'nlc'};
opt.Par.Constr.ave_uv.PossibleCope = {'Ignore','Monitor','Stop','Constr'};
opt.Par.Constr.jag_uv.PossibleType = {'nlc'};
opt.Par.Constr.jag_uv.PossibleCope = {'Ignore'};

opt.Par.Constr.peak_xyz.PossibleType = {'bnd'};
opt.Par.Constr.peak_xyz.PossibleCope = {'Ignore','Monitor','Stop','Constr'};
opt.Par.Constr.slew_xyz.PossibleType = {'nlc'};
opt.Par.Constr.slew_xyz.PossibleCope = {'Ignore','Monitor','Stop','Constr'};

opt.Par.Constr.peak_s.PossibleType = {'bnd','nlc'};
opt.Par.Constr.peak_s.PossibleCope = {'Ignore'};
opt.Par.Constr.slew_s.PossibleType = {'nlc'};
opt.Par.Constr.slew_s.PossibleCope = {'Ignore'};

opt.Par.Constr.sumabs_s.PossibleType = {'nlc'};
opt.Par.Constr.sumabs_s.PossibleCope = {'Ignore'};

opt.Par.Constr.edges_uv.PossibleType = {'bnd'};
opt.Par.Constr.edges_uv.PossibleCope = {'Ignore'};
opt.Par.Constr.edges_xyz.PossibleType = {'bnd'};
opt.Par.Constr.edges_xyz.PossibleCope = {'Ignore'};
opt.Par.Constr.edges_s.PossibleType = {'bnd'};
opt.Par.Constr.edges_s.PossibleCope = {'Ignore'};



opt.Par.Constr.Sanity = false;

SARSanity = ones(4,1);
RFSanity = zeros(7,1);
GradSanity = zeros(5,1);
% ShimSanity = zeros(6,1);



if contains( opt.Par.Controls,'u') || contains( opt.Par.Controls,'v')
    switch opt.Par.Constr.peak_uv.Type
        case opt.Par.Constr.peak_uv.PossibleType
            RFSanity(1) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.peak_uv.Type is not right')
    end
    switch opt.Par.Constr.ave_uv.Type
        case opt.Par.Constr.ave_uv.PossibleType
            RFSanity(2) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.ave_uv.Type is not right')
    end
    
    switch opt.Par.Constr.jag_uv.Type
        case opt.Par.Constr.jag_uv.PossibleType
            RFSanity(3) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.jag_uv.Type is not right')
    end
    switch opt.Par.Constr.peak_uv.Cope
        case opt.Par.Constr.peak_uv.PossibleCope
            RFSanity(4) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.peak_uv.Cope is not right')
    end
    
    
    
    switch opt.Par.Constr.ave_uv.Cope
        case opt.Par.Constr.ave_uv.PossibleCope
            RFSanity(5) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.ave_uv.Cope is not right')
    end
    
    switch opt.Par.Constr.jag_uv.Cope
        case opt.Par.Constr.jag_uv.PossibleCope
            RFSanity(6) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.jag_uv.Cope is not right')
    end
    
    
    switch opt.Par.Constr.edges_uv.Cope
        case opt.Par.Constr.edges_uv.PossibleCope
            RFSanity(7) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.edges_uv.Cope is not right')
    end
    
    
    
end

if contains( opt.Par.Controls,'x') || contains( opt.Par.Controls,'y') || contains( opt.Par.Controls,'z')
    
    switch opt.Par.Constr.peak_xyz.Type
        case opt.Par.Constr.peak_xyz.PossibleType
            GradSanity(1) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.peak_xyz.Type is not right')
    end
    
    switch opt.Par.Constr.peak_xyz.Cope
        case opt.Par.Constr.peak_xyz.PossibleCope
            GradSanity(2) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.peak_xyz.Cope is not right')
    end
    
    switch opt.Par.Constr.slew_xyz.Type
        case opt.Par.Constr.slew_xyz.PossibleType
            GradSanity(3) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.slew_xyz.Type is not right')
    end
    
    switch opt.Par.Constr.slew_xyz.Cope
        case opt.Par.Constr.slew_xyz.PossibleCope
            GradSanity(4) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.slew_xyz.Cope is not right')
    end
    
    
    switch opt.Par.Constr.edges_xyz.Cope
        case opt.Par.Constr.edges_xyz.PossibleCope
            GradSanity(5) = 1;
        otherwise
            error('Sanity_Check_Constraints: Par.Constr.edges_xyz.Cope is not right')
    end
    
else
    GradSanity = ones(5,1);
end

ShimSanity = ones(6,1);



if sum(GradSanity) == 5 && sum(RFSanity) == 7 && sum(SARSanity) == 4 && sum(ShimSanity) == 6
    opt.Par.Constr.Sanity = true;
    
end
end

function [con,opt,Grad] = CONSTR_dummyfunc(opt)

%

array = ARRAY(opt.Par.Controls,opt.u,opt.v,opt.gx,opt.gy,opt.gz,opt.s);

array_scaled = ARRAY2NORM(array,opt.arraysc);
[c,ceq,dc,dceq,con,opt]=CONSTR(array_scaled,opt);


if nargout == 3
    Grad = dc;
end
end

function [c,ceq,dc,dceq,con,opt] = CONSTR(array,opt)

ceq = [];
dceq   = [];
array_scaled = ARRAY2PHYS(array,opt.arraysc);


[u,v,gx,gy,gz,s] = iARRAY(opt.Par.Controls,array_scaled,opt.N,opt.pTx,opt.NS,opt.u,opt.v,opt.gx,opt.gy,opt.gz,opt.s);


rf = complex(u,v);
rf_long = rf(:);

% SAR
% local

c_lsar = [];
dc_lsar = [];
con.sar_l = 0;
con.rat_sar_l = 0;
opt.Par.Constr.sar_l.Viol = 0;

c_gsar = [];
dc_gsar = [];
con.sar_g = 0;
con.rat_sar_g = 0;
opt.Par.Constr.sar_g.Viol = 0;

%
% peak estimations
if contains(opt.Par.Controls,'u')
    [c_peak_u,dc_peak_u,con.peak_u,opt.Par.Constr.peak_uv.Violu,con.rat_peak_u] = PEAK(opt.Par.Constr.peak_uv.Cope,opt.Par.Constr.peak_uv.Type,opt.Par.Constr.peak_uv.Power,opt.Par.Constr.peak_uv.Lim,opt.Par.Constr.peak_uv.LimCor,opt.Par.Constr.peak_uv.Unit_cor_factor,u,opt.max_uv);
    c_peak_u = [];
    dc_peak_u = [];
else
    c_peak_u = [];
    dc_peak_u = [];
    con.peak_u = 0;
    con.rat_peak_u = 0;
    opt.Par.Constr.peak_uv.Violu = 0;
end
%
if contains(opt.Par.Controls,'v')
    [c_peak_v,dc_peak_v,con.peak_v,opt.Par.Constr.peak_uv.Violv,con.rat_peak_v] = PEAK(opt.Par.Constr.peak_uv.Cope,opt.Par.Constr.peak_uv.Type,opt.Par.Constr.peak_uv.Power,opt.Par.Constr.peak_uv.Lim,opt.Par.Constr.peak_uv.LimCor,opt.Par.Constr.peak_uv.Unit_cor_factor,v,opt.max_uv);
    c_peak_v = [];
    dc_peak_v = [];
else
    c_peak_v = [];
    dc_peak_v = [];
    con.peak_v = 0;
    con.rat_peak_v = 0;
    opt.Par.Constr.peak_uv.Violv = 0;
end
%
if contains(opt.Par.Controls,'x')
    [c_peak_x,dc_peak_x,con.peak_x,opt.Par.Constr.peak_xyz.Violx,con.rat_peak_x] = PEAK(opt.Par.Constr.peak_xyz.Cope,opt.Par.Constr.peak_xyz.Type,opt.Par.Constr.peak_xyz.Power,opt.Par.Constr.peak_xyz.Lim,opt.Par.Constr.peak_xyz.LimCor,opt.Par.Constr.peak_xyz.Unit_cor_factor,gx,opt.max_g);
    c_peak_x = [];
    dc_peak_x = [];
    [c_slew_x,dc_slew_x,con.slew_x,opt.Par.Constr.slew_xyz.Violx,con.rat_slew_x] = SLEW(opt.Par.Constr.slew_xyz.Cope,opt.Par.Constr.slew_xyz.Power,opt.Par.Constr.slew_xyz.Lim,opt.Par.Constr.slew_xyz.LimCor,opt.Par.Constr.slew_xyz.Unit_cor_factor,gx,opt.dt,opt.max_g,opt.N);
else
    c_peak_x = [];
    dc_peak_x = [];
    con.peak_x = 0;
    con.rat_peak_x = 0;
    opt.Par.Constr.peak_xyz.Violx = 0;
    c_slew_x = [];
    dc_slew_x = [];
    con.slew_x = 0;
    con.rat_slew_x = 0;
    opt.Par.Constr.slew_xyz.Violx = 0;
end
%
if contains(opt.Par.Controls,'y')
    [c_peak_y,dc_peak_y,con.peak_y,opt.Par.Constr.peak_xyz.Violy,con.rat_peak_y] = PEAK(opt.Par.Constr.peak_xyz.Cope,opt.Par.Constr.peak_xyz.Type,opt.Par.Constr.peak_xyz.Power,opt.Par.Constr.peak_xyz.Lim,opt.Par.Constr.peak_xyz.LimCor,opt.Par.Constr.peak_xyz.Unit_cor_factor,gy,opt.max_g);
    c_peak_y = [];
    dc_peak_y = [];
    [c_slew_y,dc_slew_y,con.slew_y,opt.Par.Constr.slew_xyz.Violy,con.rat_slew_y] = SLEW(opt.Par.Constr.slew_xyz.Cope,opt.Par.Constr.slew_xyz.Power,opt.Par.Constr.slew_xyz.Lim,opt.Par.Constr.slew_xyz.LimCor,opt.Par.Constr.slew_xyz.Unit_cor_factor,gy,opt.dt,opt.max_g,opt.N);
else
    c_peak_y = [];
    dc_peak_y = [];
    con.peak_y = 0;
    con.rat_peak_y = 0;
    opt.Par.Constr.peak_xyz.Violy = 0;
    c_slew_y = [];
    dc_slew_y = [];
    con.slew_y = 0;
    con.rat_slew_y = 0;
    opt.Par.Constr.slew_xyz.Violy = 0;
end
%
if contains(opt.Par.Controls,'z')
    [c_peak_z,dc_peak_z,con.peak_z,opt.Par.Constr.peak_xyz.Violz,con.rat_peak_z] = PEAK(opt.Par.Constr.peak_xyz.Cope,opt.Par.Constr.peak_xyz.Type,opt.Par.Constr.peak_xyz.Power,opt.Par.Constr.peak_xyz.Lim,opt.Par.Constr.peak_xyz.LimCor,opt.Par.Constr.peak_xyz.Unit_cor_factor,gz,opt.max_g);
    c_peak_z = [];
    dc_peak_z = [];
    [c_slew_z,dc_slew_z,con.slew_z,opt.Par.Constr.slew_xyz.Violz,con.rat_slew_z] = SLEW(opt.Par.Constr.slew_xyz.Cope,opt.Par.Constr.slew_xyz.Power,opt.Par.Constr.slew_xyz.Lim,opt.Par.Constr.slew_xyz.LimCor,opt.Par.Constr.slew_xyz.Unit_cor_factor,gz,opt.dt,opt.max_g,opt.N);
else
    c_peak_z = [];
    dc_peak_z = [];
    con.peak_z = 0;
    con.rat_peak_z = 0;
    opt.Par.Constr.peak_xyz.Violz = 0;
    c_slew_z = [];
    dc_slew_z = [];
    con.slew_z = 0;
    con.rat_slew_z = 0;
    opt.Par.Constr.slew_xyz.Violz = 0;
end


if contains(opt.Par.Controls,'u') || contains(opt.Par.Controls,'v')
    [c_ave_uv,dc_ave_uv,con.ave_uv,opt.Par.Constr.ave_uv.Violuv,con.rat_ave_uv] = AVE('uv',opt.Par.Constr.ave_uv.Cope,opt.Par.Constr.ave_uv.Power,opt.Par.Constr.ave_uv.Lim,opt.Par.Constr.ave_uv.LimCor,opt.Par.Constr.ave_uv.Unit_cor_factor,u,v,opt.max_uv,opt.N,opt.pTx);
    
else
    c_ave_uv = [];
    dc_ave_uv = [];
    con.ave_uv = 0;
    con.rat_ave_uv = 0;
    opt.Par.Constr.ave_uv.Violuv = 0;
end


c_jag_uv = [];
dc_jag_uv = [];
con.jag_uv = 0;
opt.Par.Constr.jag_uv.Viol = 0;


if contains(opt.Par.Controls,'s')
    
else
    c_peak_s = [];
    dc_peak_s = [];
    con.peak_s = 0;
    con.rat_peak_s = 0;
    opt.Par.Constr.peak_s.Viol = 0;
    c_slew_s = [];
    dc_slew_s = [];
    con.slew_s = 0;
    con.rat_slew_s = 0;
    opt.Par.Constr.slew_s.Viol = 0;
    
    c_sumabs_s = [];
    dc_sumabs_s = [];
    con.sumabs_s = 0;
    con.rat_sumabs_s = 0;
    opt.Par.Constr.sumabs_s.Viol = 0;
    
end

% sizes
[A_lsar,B_lsar] = size(dc_lsar);
[A_gsar,B_gsar] = size(dc_gsar);
[A_ave_uv,B_ave_uv] = size(dc_ave_uv);
[A_jag_uv,B_jag_uv] = size(dc_jag_uv);
[A_slew_x,B_slew_x] = size(dc_slew_x);
[A_slew_y,B_slew_y] = size(dc_slew_y);
[A_slew_z,B_slew_z] = size(dc_slew_z);
[A_slew_s,B_slew_s] = size(dc_slew_s);
[A_sumabs_s,B_sumabs_s] = size(dc_sumabs_s);

c = [c_lsar;c_gsar;c_ave_uv;c_jag_uv;c_slew_x;c_slew_y;c_slew_z;c_slew_s;c_sumabs_s];




Arf = 0;
if contains(opt.Par.Controls,'u')
    Arf = Arf+opt.N*opt.pTx;
end
if contains(opt.Par.Controls,'v')
    Arf = Arf+opt.N*opt.pTx;
end

Ag = 0;
if contains(opt.Par.Controls,'x')
    Ag = Ag+opt.N;
end
if contains(opt.Par.Controls,'y')
    Ag = Ag+opt.N;
end
if contains(opt.Par.Controls,'z')
    Ag = Ag+opt.N;
end
As = 0;
if contains(opt.Par.Controls,'s')
    As = As + opt.N*opt.NS;
end

Atot = Arf+Ag+As;

if ~isempty(dc_lsar)
    dc_lsar = [dc_lsar;zeros(Atot-Arf,B_lsar)];
end

if ~isempty(dc_gsar)
    dc_gsar = [dc_gsar;zeros(Atot-Arf,B_gsar)];
end
if ~isempty(dc_ave_uv)
    dc_ave_uv = [dc_ave_uv;zeros(Atot-Arf,B_ave_uv)];
end


if strcmp(opt.Par.Constr.slew_xyz.Cope,'Constr') && (contains(opt.Par.Controls,'x') || contains(opt.Par.Controls,'y') || contains(opt.Par.Controls,'z') )
    if contains(opt.Par.Controls,'x') && ~contains(opt.Par.Controls,'y') && ~contains(opt.Par.Controls,'z') % X
        
        dc_slew_x_ = [zeros(Arf,B_slew_x);dc_slew_x;zeros(As,B_slew_x)];
        dc_slew_y_ = [];
        dc_slew_z_ = [];
    elseif ~contains(opt.Par.Controls,'x') && contains(opt.Par.Controls,'y') && ~contains(opt.Par.Controls,'z')% Y
        dc_slew_x_ = [];
        dc_slew_y_ = [zeros(Arf,B_slew_y);dc_slew_y;zeros(As,B_slew_y)];
        dc_slew_z_ = [];
    elseif ~contains(opt.Par.Controls,'x') && ~contains(opt.Par.Controls,'y') && contains(opt.Par.Controls,'z') % Z
        dc_slew_x_ = [];
        dc_slew_y_ = [];
        dc_slew_z_ = [zeros(Arf,B_slew_z);dc_slew_z;zeros(As,B_slew_z)];
        
    elseif contains(opt.Par.Controls,'x') && contains(opt.Par.Controls,'y') && ~contains(opt.Par.Controls,'z')% XY
        
        dc_slew_x_ = [zeros(Arf,B_slew_x);dc_slew_x;zeros(size(dc_slew_y));zeros(As,B_slew_x)];
        dc_slew_y_ = [zeros(Arf,B_slew_y);zeros(size(dc_slew_x));dc_slew_y;zeros(As,B_slew_y)];
        dc_slew_z_ = [];
        
    elseif contains(opt.Par.Controls,'x') && ~contains(opt.Par.Controls,'y') && contains(opt.Par.Controls,'z')% XZ
        dc_slew_x_ = [zeros(Arf,B_slew_x);dc_slew_x;zeros(size(dc_slew_z));zeros(As,B_slew_x)];
        dc_slew_y_ = [];
        dc_slew_z_ = [zeros(Arf,B_slew_z);zeros(size(dc_slew_x));dc_slew_z;zeros(As,B_slew_z)];
    elseif ~contains(opt.Par.Controls,'x') && contains(opt.Par.Controls,'y') && contains(opt.Par.Controls,'z')% YZ
        dc_slew_x_ = [];
        dc_slew_y_ = [zeros(Arf,B_slew_y);dc_slew_y;zeros(size(dc_slew_z));zeros(As,B_slew_y)];
        dc_slew_z_ = [zeros(Arf,B_slew_z);zeros(size(dc_slew_y));dc_slew_z;zeros(As,B_slew_z)];
    elseif contains(opt.Par.Controls,'x') && contains(opt.Par.Controls,'y') && contains(opt.Par.Controls,'z') % XYZ
        dc_slew_x_ = [zeros(Arf,B_slew_x);dc_slew_x;zeros(size(dc_slew_y));zeros(size(dc_slew_z));zeros(As,B_slew_x)];
        dc_slew_y_ = [zeros(Arf,B_slew_y);zeros(size(dc_slew_x));dc_slew_y;zeros(size(dc_slew_z));zeros(As,B_slew_y)];
        dc_slew_z_ = [zeros(Arf,B_slew_z);zeros(size(dc_slew_x));zeros(size(dc_slew_y));dc_slew_z;zeros(As,B_slew_z)];
    else
        dc_slew_x_ = [];
        dc_slew_y_ = [];
        dc_slew_z_ = [];
    end
else
    dc_slew_x_ = [];
    dc_slew_y_ = [];
    dc_slew_z_ = [];
end


dc = sparse([dc_lsar,dc_gsar,dc_ave_uv,dc_slew_x_,dc_slew_y_,dc_slew_z_,dc_slew_s,dc_sumabs_s]);

if ...
        (strcmp(opt.Par.Constr.peak_uv.Cope,'Constr') && strcmp(opt.Par.Constr.peak_uv.Type,'nlc') && (contains(opt.Par.Controls,'u') || contains(opt.Par.Controls,'v') )) || ...
        (strcmp(opt.Par.Constr.peak_xyz.Cope,'Constr') && strcmp(opt.Par.Constr.peak_xyz.Type,'nlc') && (contains(opt.Par.Controls,'u') || contains(opt.Par.Controls,'v'))) ||...
        (strcmp(opt.Par.Constr.ave_uv.Cope,'Constr') && (contains(opt.Par.Controls,'u') || contains(opt.Par.Controls,'v') ))||...
        (strcmp(opt.Par.Constr.jag_uv.Cope,'Constr') && (contains(opt.Par.Controls,'u') || contains(opt.Par.Controls,'v') ))||...
        (strcmp(opt.Par.Constr.slew_xyz.Cope,'Constr') && (contains(opt.Par.Controls,'x') || contains(opt.Par.Controls,'y') || contains(opt.Par.Controls,'z')))
    opt.con = con;
end


end

function [c,dc,Val,Viol,Rat] = PEAK(Cope,Type,Power,Lim,LimCor,Unit_cor_factor,Candidate,Max)

if strcmp(Cope,'Ignore')
    
    c = [];
    dc = [];
    Viol = 0;
    Val = 0;
else
    Candidate_ = permute(Candidate,[2,1]); Candidate_ = Candidate_(:).';
    if Power == 2
        c = abs(Candidate_).^2 - LimCor;
        if strcmp(Cope,'Constr') && strcmp(Type,'nlc')
            dc = [2*diag(Candidate_)].*Max;
        else
            dc = [];
        end
    else
        
        c = abs(Candidate_).^2 - LimCor^2;
        
        if strcmp(Cope,'Constr') && strcmp(Type,'nlc')
            dc = [2*diag(Candidate_)].*Max;
        else
            dc = [];
        end
    end
    c = c(:);
    
    
    if Power == 2
        Val = max(c)+LimCor;
        Val = Val./Unit_cor_factor;
    else
        %         peak_uv = max(c)+LimCor^2;
        Val = max(max(Candidate_(:)));
        Val = Val./Unit_cor_factor;
    end
    
    if Lim-Val < 0
        Viol = 1;
    else
        Viol = 0;
    end
    if strcmp(Cope,'Stop') || strcmp(Cope,'Monitor')
        c = [];
    end
    if nargout == 5
        Rat = Val/Lim;
    end
end
end

function [c,dc,Val,Viol,Rat] = AVE(uorv,Cope,Power,Lim,LimCor,Unit_cor_factor,u,v,maxRF,N,pTx)

if strcmp(Cope,'Ignore')
    
    Viol = 0;
    Val = 0;
    c = [];
    dc = [];
else
    switch uorv
        case 'uv'
            rf = complex(u,v);
            if Power == 2
                c = sum(abs(rf).^2,2)./N - LimCor;
                
                if strcmp(Cope,'Constr')
                    tmp = zeros(N*pTx,pTx);
                    k = 1;
                    for n = 1:N
                        tmp(k:pTx+k-1,:) = diag(rf(:,n));
                        k = k+pTx;
                    end
                    
                    
                    dc_ = 2.*[real(tmp);imag(tmp)];
                    dc = zeros(size(dc_));
                    for n = 1:pTx
                        dc(:,n) = b6_Rearrange_controls(dc_(:,n).',N,pTx);
                    end
                    dc = dc.*maxRF./N;
                    
                else
                    dc =[];
                end
                Val = max(c)+LimCor;
                
                if LimCor-Val < 0
                    
                    Viol = 1;
                else
                    Viol = 0;
                end
                
                
                Val = Val./Unit_cor_factor; % This conversion must happen after Viol is made to take dutycycle into account
                
                if strcmp(Cope,'Stop') || strcmp(Cope,'Monitor')
                    c = [];
                end
            end
        otherwise
            fprintf(2,'%s: only implemented for uv',mfilename)
            Viol = 0;
            Val = 0;
            c = [];
            dc = [];
    end
    
end
if nargout == 5
    Rat = Val/Lim;
end
end

function [c,dc,Val,Viol,Rat] = SLEW(Cope,Power,Lim,LimCor,Unit_cor_factor,Candidate,dt,Max,N)

[A,B] = size(Candidate);

if strcmp(Cope,'Ignore')
    Viol = 0;
    Val = 0;
    c = [];
    dc = [];
else
    
    if Power == 2
        
        
        Candidate = permute(Candidate,[2,1]);
        %         Candidate = Candidate(:);
        g0 = [zeros(1,A);Candidate;zeros(1,A)];
        
        c_ = (diff(g0,1,1)./dt).^2;
        
        c = c_(:)-LimCor;
        if strcmp(Cope,'Constr')
            
            d0 = 2*Candidate-2*[zeros(1,A);Candidate(1:end-1,:)];
            d1 = 2*Candidate-2*[Candidate(2:end,:);zeros(1,A)];
            
            %             d0 = [d0,zeros(1,N*2)];
            %             d1 = [d1,d0y,zeros(1,N)];
            %             d2 = [zeros(1,N),d1y,d0z];
            %             d3 = [zeros(1,N*2),d1z];
            
            factor = 1./dt^2.*Max;
            
            dc = spdiags([d0(:),d1(:)],[0,1],N*A,(N+1)*A);
            dc = dc.*factor;
            
            dc(isnan(dc))  = 0;
            dc(isinf(-dc)) = 0;
            dc(isinf(dc))  = 0;
        else
            dc = [];
        end
        
    end
    c = c(:);
    if Power == 2
        Val = max(c)+LimCor;
        Val = Val./Unit_cor_factor;
        
    end
    
    if Lim-Val < 0 %#ok<*BDSCI>
        Viol = 1;
    else
        Viol = 0;
    end
    if strcmp(Cope,'Stop') || strcmp(Cope,'Monitor')
        c = [];
    end
    if nargout == 5
        Rat = Val/Lim;
    end
end
end

