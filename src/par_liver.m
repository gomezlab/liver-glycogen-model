%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% Uses UNC's Copyright and Permission Notice:	
% COPYRIGHT AND PERMISSION NOTICE	
% UNC Software:  << Liver_glycogen_model >>
% Copyright (C) 2009 The University of North Carolina at Chapel
% All rights reserved.
% The University of North Carolina at Chapel Hill (?UNC?) and the
% developers (?Developers?) of << Liver_glycogen_model >> (?Software?) give
% recipient (?Recipient?) and Recipient?s Institution (?Institution?)
% permission to use and copy the software in source and binary forms, with
% or without modification for non-commercial purposes only provided that
% the following conditions are met:
% 1)	All copies of Software in binary form and/or source code, related
% documentation and/or other materials provided with the Software must
% reproduce and retain the above copyright notice, this list of conditions
% and the following disclaimer.    
% 2)	Recipient and Institution shall not distribute Software to any
% third parties. 
% 3)	The Software is provided ?As Is.? The Developers can not guarantee
% the provision of technical support or consultation for the Software. The
% Developers may provide a location on a UNC Web Site for Recipients to
% post comments, questions, and suggestions at some time in the future.
% Recipient may provide the Developers with feedback on the use of the
% Software in their research at that time.  The Developers and UNC are
% permitted to use any information Recipient provides in making changes to
% the Software.        
% 4)	Recipient acknowledges that the Developers, UNC and its licensees
% may develop modifications to Software that may be substantially similar
% to Recipient?s modifications of Software, and that the Developers, UNC
% and its licensees shall not be constrained in any way by Recipient in
% UNC?s or its licensees? use or management of such modifications.
% Recipient acknowledges the right of the Developers and UNC to prepare and
% publish modifications to Software that may be substantially similar or
% functionally equivalent to your modifications and improvements, and if
% Recipient or Institution obtains patent protection for any modification
% or improvement to Software, Recipient and Institution agree not to allege
% or enjoin infringement of their patent by the Developers, UNC or any of
% UNC?s licensees obtaining modifications or improvements to Software from
% the UNC or the Developers.            
% 5)	Recipient and Developer will acknowledge in their respective
% publications the contributions made to each other?s research involving or
% based on the Software. The current citations for Software are:
% <<<PENDING LIST>>>	
% 6)	Any party desiring a license to use the Software for commercial
% purposes shall contact The Office of Technology Development at UNC at
% 919-966-3929.  
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS, CONTRIBUTORS, AND THE
% UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL "AS IS" AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT OWNER, CONTRIBUTORS OR THE UNIVERSITY OF
% NORTH CAROLINA AT CHAPEL HILL BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
%
%Description:
% par_liver.m 
% called by ode_liver.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function par = par_liver()

% LIVER:
k1 = 3;.2;        % GK, gluc -> g6p
km1 = 7.7;       %   k1 * gluc /(km1 + gluc) * (1 + (B_ins/kDins)^ep1)
ep1 = 10;

k_1 = 4;.2;.1;.2;       % G6Pase, g6p -> gluc
km_1 = 1.3;      %   k_1 * g6p /(km_1 + g6p) * (1 + (B_glucgn/kDglucgn)^ep9)
ep9 = 10;

% Mutalik
k2 = 200;        % GS, g6p -> glyc
km2 = .57;         %   k2 * GSa * g6p / (km2 + g6p)
k_2 = 20;      % GP, glyc -> g6p
km_2 = 1.4;        %   k_2 * GPa * glycgn / (km_2 + glycgn)

k3 = .1;.05;.1;    % g6p -> pep
km3 = .01;   %   k3 * g6p/(km3+g6p) * (1 + (B_ins/kDins)^ep2) * kDglucgn^en1/(kDglucgn^en1 + B_glucgn^en1)
ep2 = 10;
en1 = 10;

k_3 = .3;   % pep -> g6p
km_3 = .0034;    %   k_3 * pep/(km_3+pep) * (1 + (B_glucgn/kDglucgn)^ep8) * kDins^en6/(kDins^en6 + B_ins^en6)
ep8 = 10;
en6 = 10;

k4 = 2;    % pep -> pyr
km4 = .18;   %   k4 * pep  * (1 + (B_ins/kDins)^ep3) * kDglucgn^en2/(kDglucgn^en2 + B_glucgn^en2) * ki13/(ki13 + alan)
ki13 = 2;
ep3 = 10;
en2 = 10;

k5 = .0;    % LDH, pyr -> lac
km5 = .03;       % k5 * pyr / (km5 + pyr)
k_5 = .1;   % LDH, lac -> pyr
km_5 = .8;       % k_5 * lac / (km_5 + lac)

k6 = 1;    % PC, pyr -> oa_m
km6 = 0.22;
p2 = 1;     %   k6 * pyr /(km6+pyr)* (1 + p2*acet_m/(kp2 + acet_m)) * (1 + (B_glucgn/kDglucgn)^ep4)
kp2 = 2;     
p7 = 10;   
ep4 = 10;
k7 = 1;    % PDH, pyr -> acet_m
ki8 = 2;    %   k7 * pyr / (km7 + pyr) * ki8/(ki8 + acet_m)
km7 = .0204;
k8 = .1;    % oa_m, acet_m -> cit
ki4 = 3;    %   k8 * oa_m * acet_m * ki4/(ki4 + palmCoA)
k9 = .1;    % cit -> aK
k10 = .1;   % aK -> malate
k11 = .6;   % oa_m -> malate
ep5 = 10;        % k11 * oa_m * (1 + (B_glucgn/kDglucgn)^ep5)

k_11 = .01;   % malate -> oa_m
k12 = .6;   % malate -> oa_c
ep6 = 10;        % k12 * malate * kDins^en4/(kDins^en4 + B_ins^en4)* (1 + (B_glucgn/kDglucgn)^ep6)
en4 = 10;
k13 = .5;   % PEPCK, oa_c -> pep
ep7 = 10;    %   k13 * oa_c * (1 + (B_glucgn/kDglucgn)^ep7) * kDins^en5/(kDins^en5 + B_ins^en5)
en5 = 10;
k14 = .01;   % cit -> acet_c, oa_c
ep10 = 10;     %  k14 * citrate * (1 + (B_ins/kDins)^ep10) * kDglucgn^en7/(kDglucgn^en7 + B_glucgn^en7)
en7 = 10;
k15 = .01;   % acet_c -> malonyl
p1 = 1;     %   k15 * acet_c * (1 + p1*citrate/(kp1 + citrate)) * ki5/(ki5 + palmCoA)
kp1 = .5;     
ki5 = 1;    
k16 = .01;   % acet_c, malonyl -> palm
                 % k16 * acet_c * malonyl^7
k17 = .01;   % palm -> palmCoA
ki1 = .1;   %   k17 * palm * ki1/(ki1 + malonyl)
k18 = .01;   % palmCoA -> acet_m
ki2 = 1;    %   k18 * palmCoA * ki2/(ki2 + malonyl)
k19 = .01;   % acet_m -> ket
ep11 = 10;     %   k19 * acet_m^2 * (1 + (cAMP/kdcAMP)^ep11) * kDins^en8/(kDins^en8 + B_ins^en8)
en8 = 10;
k20 = .2;   % alan -> pyr, non-ALT pathway
k21 = .001;    % pyr,glutamate -> alan,aK
en3 = 10;                % k21 * pyr * glutamate * kDins^en3/(kDins^en3 + B_ins^en3)
km21p = .4;
km21g = 4.3;
k_21 = .2;   % alan,aK -> pyr,glutamate
                 %  k_21 * alan * aK * acet_m * kDins^en3/(kDins^en3 + B_ins^en3) 
km_21a = 21;
km_21k = .22;
k22 = .01;  % glutamate sink


kc1 = 1;   % cAMP growth due to B_glucgn 
kcm1 = 40*10^(-9); % kc1 * B_glucgn/(kcm1 + B_glucgn)
kc2 = kc1 * 10^(5.5);   % cAMP decay due to B_ins
kcm2 = 1*10^(-6);
% dcAMP = kc1 * B_glucgn/(kcm1 + B_glucgn) ...
%       - kc2 * B_ins/(kcm2 + B_ins) * cAMP;

kDins = 1 * 10^(-6);
kDins2 = .75*10^(-6);
kDglucgn = 40 * 10^(-9);
kdcAMP = 1*10^(-5.5);

par=[k1; km1; k_1; km_1; k2; km2; k_2; km_2; k3; km3; k_3; km_3; 
    k4; km4; k5; km5; k_5; km_5; k6; k7; km7; k8; k9; k10; 
    k11; k_11; k12; k13; k14; k15; k16; k17; k18; k19; k21; k_21;
    ki1; ki2; ki4; ki5; ki8; ki13; p1; p2; p7; 
    kc1; kc2; kcm1; kcm2; kDins; kDglucgn; kdcAMP; kp1; kp2; 
    ep1; ep2; ep3; ep4; ep5; ep6; ep7; ep8; ep9; ep10; ep11; 
    en1; en2; en3; en4; en5; en6; en7; en8; kDins2; k20; k22; km6; km21p;
    km21g; km_21a; km_21k];

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



