%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% Uses UNC's Copyright and Permission Notice:	
% COPYRIGHT AND PERMISSION NOTICE	
% UNC Software:  << Liver_glycogen_model >>
% Copyright (C) 2009 The University of North Carolina at Chapel
% All rights reserved.
% The University of North Carolina at Chapel Hill (?UNC?) and the
% developers (?Developers?) of <<Liver_glycogen_model>> (?Software?) give
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
%Description 
% ode_liver.m
% Called by ode_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = ode_liver(t,y,B_gluc,B_ffa,B_ins,B_alan,B_lac,B_glucgn,GSa,GPa,run)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters/variables:

par=par_liver();
k1=par(1); km1=par(2); k_1=par(3); km_1=par(4); k2=par(5); 
km2=par(6); k_2=par(7); km_2=par(8); k3=par(9); km3=par(10); 
k_3=par(11); km_3=par(12); k4=par(13); km4=par(14); k5=par(15); 
km5=par(16); k_5=par(17); km_5=par(18); k6=par(19); k7=par(20); 
km7=par(21); k8=par(22); k9=par(23); k10=par(24); k11=par(25); 
k_11=par(26); k12=par(27); k13=par(28); k14=par(29); k15=par(30); 
k16=par(31); k17=par(32); k18=par(33); k19=par(34); k21=par(35); 
k_21=par(36); ki1=par(37); ki2=par(38); ki4=par(39); ki5=par(40); 
ki8=par(41); ki13=par(42); p1=par(43); p2=par(44); p7=par(45); 
kc1=par(46); kc2=par(47); kcm1=par(48); kcm2=par(49); kDins=par(50); 
kDglucgn=par(51); kdcAMP=par(52); kp1=par(53); kp2=par(54); ep1=par(55); 
ep2=par(56); ep3=par(57); ep4=par(58); ep5=par(59); ep6=par(60); 
ep7=par(61); ep8=par(62); ep9=par(63); ep10=par(64); ep11=par(65); 
en1=par(66); en2=par(67); en3=par(68); en4=par(69); en5=par(70); 
en6=par(71); en7=par(72); en8=par(73); kDins2=par(74);
k20=par(75); k22=par(76); km6=par(77); km21p=par(78);
km21g=par(79); km_21a=par(80); km_21k=par(81);

par=par_blood();
kf=par(1); ktL1=par(2); ktL2=par(3); ktL3=par(4); ktL4=par(5); 
ktL5=par(6); ktL6=par(7); ktF1=par(8); ktF3=par(9);
ep13=par(34); en9=par(35); en10=par(36); en11=par(37);

% transcriptional regulation
%     tr_par=par_transcripts(t,k1,k2,k_2,k3,k4,k5,k_5,k6,k7,k11,k12,k13,k14,k16,k17,k18,k19,trans);
%     k1 = tr_par(1); k2 = tr_par(2); k_2 = tr_par(3); k3 = tr_par(4);
%     k4 = tr_par(5); k5 = tr_par(6); k_5 = tr_par(7); k6 = tr_par(8);
%     k7 = tr_par(9); k11 = tr_par(10); k12 = tr_par(11);
%     k13 = tr_par(12); k14 = tr_par(13); k16 = tr_par(14); k17 = tr_par(15);
%     k18 = tr_par(16); k19 = tr_par(17);

gluc=y(1); g6p=y(2); glycgn=y(3); pep=y(4); pyr=y(5); 
lac=y(6); oa_m=y(7); acet_m=y(8); citrate=y(9); aK=y(10); 
malate=y(11); oa_c=y(12); acet_c=y(13); malonyl=y(14); palm=y(15); 
palmCoA=y(16); ket=y(17); alan=y(18); cAMP=y(19); glutamate=y(20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIVER ode eqns:

dgluc = ktL1 * (B_gluc - gluc) ...      % gluc transport, B_gluc <-> gluc
      - k1 * gluc /(km1 + gluc) * (1 + B_ins^ep1/(kDins^ep1+B_ins^ep1)) ...   % gluc -> g6p
      + k_1 * g6p /(km_1 + g6p) * (1 + B_glucgn^ep9/(kDglucgn^ep9+B_glucgn^ep9));      % g6p -> gluc

dg6p = k1 * gluc /(km1 + gluc) * (1 + B_ins^ep1/(kDins^ep1+B_ins^ep1)) ...        % gluc -> g6p
     - k_1 * g6p /(km_1 + g6p) * (1 + B_glucgn^ep9/(kDglucgn^ep9+B_glucgn^ep9)) ...        % g6p -> gluc
     - k2 * GSa * g6p / (km2 + g6p) ...         % g6p -> glycgn
     + k_2 * GPa * glycgn / (km_2 + glycgn) ...     % glycgn -> g6p
     - k3 * g6p/(km3+g6p) * (1 + B_ins^ep2/(kDins^ep2+B_ins^ep2)) * kDglucgn^en1/(kDglucgn^en1 + B_glucgn^en1)  ...         % g6p -> pep
     + k_3 * pep/(km_3+pep) * (1 + B_glucgn^ep8/(kDglucgn^ep8+B_glucgn^ep8)) * kDins2^en6/(kDins2^en6 + B_ins^en6);          % pep -> g6p
 
dglycgn = k2 * GSa * g6p / (km2 + g6p) ...         % g6p -> glycgn
        - k_2 * GPa * glycgn / (km_2 + glycgn);     % glycgn -> g6p
    
dpep = 2 * k3 * g6p/(km3+g6p) * (1 + B_ins^ep2/(kDins^ep2+B_ins^ep2)) * kDglucgn^en1/(kDglucgn^en1 + B_glucgn^en1) ...        % g6p -> pep
     - 2 * k_3 * pep/(km_3+pep) * (1 + B_glucgn^ep8/(kDglucgn^ep8+B_glucgn^ep8)) * kDins2^en6/(kDins2^en6 + B_ins^en6) ...       % pep -> g6p
     - k4 * pep / (km4 + pep) * (1 + B_ins^ep3/(kDins^ep3+B_ins^ep3)) * kDglucgn^en2/(kDglucgn^en2 + B_glucgn^en2) * ki13/(ki13 + alan) ...
     + k13 * oa_c * (1 + B_glucgn^ep7/(kDglucgn^ep7+B_glucgn^ep7)) * kDins^en5/(kDins^en5 + B_ins^en5);         % pep carboxylase, oa_c --> pep

dpyr = k4 * pep / (km4 + pep) * (1 + B_ins^ep3/(kDins^ep3+B_ins^ep3)) * kDglucgn^en2/(kDglucgn^en2 + B_glucgn^en2) * ki13/(ki13 + alan) ...        % pyruvate kinase, pep -> pyr
     - k5 * pyr / (km5 + pyr)...        % pyr -> lac
     + k_5 * lac / (km_5 + lac)...       % lac -> pyr
     - k6 * pyr / (km6 + pyr) * (1 + p2*acet_m/(kp2 + acet_m)) * (1 + B_glucgn^ep4/(kDglucgn^ep4+B_glucgn^ep4)) ...         % pyruvate carboxylase, pyr -> oa_m
     - k7 * pyr / (km7 + pyr)  * ki8/(ki8 + acet_m) ...            % pyruvate dehydrogenase, pyr -> acet_m
     - k21 ./(km21p + pyr) ./ (km21g + glutamate) * pyr * glutamate * kDins^en3/(kDins^en3 + B_ins^en3)  ...           % pyr,glutamate -> alan,aK
     + k_21 ./(km_21a + alan) ./ (km_21k + aK) * alan * aK * acet_m * kDins^en3/(kDins^en3 + B_ins^en3) ...   % alan,aK -> pyr,glutamate
     + k20 * alan * kDins2^en3/(kDins2^en3 + B_ins^en3);           % alternate pathway, alan -> pyr
 
dlac = k5 * pyr / (km5 + pyr)...        % pyr -> lac
     - k_5 * lac / (km_5 + lac)...       % lac -> pyr
     + ktL2 * B_lac;     % lac transport (bl to liv), B_lac --> lac

doa_m = k6 * pyr / (km6 + pyr) * (1 + p2*acet_m/(kp2 + acet_m)) * (1 + B_glucgn^ep4/(kDglucgn^ep4+B_glucgn^ep4)) ...            % pyruvate carboxylase, pyr -> oa_m
      - k8 * oa_m * acet_m * ki4/(ki4 + palmCoA) ...  % cit synthase, oa_m + acet_m -> citrate
      - k11 * oa_m * (1 + B_glucgn^ep5/(kDglucgn^ep5+B_glucgn^ep5))...              % malate dehydrogenase, oa_m -> malate
      + k_11 * malate;      % malate dehydrogenase, malate -> oa_m
  
dacet_m = k7 * pyr / (km7 + pyr) * ki8/(ki8 + acet_m) ...          % pyruvate dehydrogenase, pyr -> acet_m
        - k8 * oa_m * acet_m * ki4/(ki4 + palmCoA) ...  % cit synthase, oa_m + acet_m -> citrate
        + 8 * k18 * palmCoA * ki2/(ki2 + malonyl) ...     % beta-oxidation, palmCoA -> acet_m
        - 2 * k19 * acet_m^2 * (1 + cAMP^ep11/(kdcAMP^ep11+cAMP^ep11)) * kDins^en8/(kDins^en8 + B_ins^en8);         % HMG synthase/lyase, acet_m -> ket

dcitrate = k8 * oa_m * acet_m * ki4/(ki4 + palmCoA) ... % cit synthase, oa_m + acet_m -> citrate
         - k9 * citrate ...     % parital TCA cycle, citrate -> alpha-K
         - k14 * citrate * (1 + B_ins^ep10/(kDins^ep10+B_ins^ep10)) * kDglucgn^en7/(kDglucgn^en7 + B_glucgn^en7);     % citrate lyase, citrate -> acet_c, oa_c

daK = k9 * citrate ...          % parital TCA cycle, citrate -> alpha-K
    - k10 * aK ...                  % parital TCA cycle, alpha-K -> malate
    + k21 ./(km21p + pyr) ./ (km21g + glutamate) * pyr * glutamate * kDins^en3/(kDins^en3 + B_ins^en3) ...           % pyr,glutamate -> alan,aK
    - k_21 ./(km_21a + alan) ./ (km_21k + aK) * alan * aK * acet_m * kDins^en3/(kDins^en3 + B_ins^en3);   % alan,aK -> pyr,glutamate

dmalate = k10 * aK ...           % parital TCA cycle, alpha-K -> malate
        + k11 * oa_m * (1 + B_glucgn^ep5/(kDglucgn^ep5+B_glucgn^ep5))...         % malate dehydrogenase, oa_m -> malate
        - k_11 * malate ...      % malate dehydrogenase, malate -> oa_m
        - k12 * malate * kDins^en4/(kDins^en4 + B_ins^en4)* (1 + B_glucgn^ep6/(kDglucgn^ep6+B_glucgn^ep6));          % malate shuttle/DH, malate -> oa_c
    
doa_c = k12 * malate * kDins^en4/(kDins^en4 + B_ins^en4) * (1 + B_glucgn^ep6/(kDglucgn^ep6+B_glucgn^ep6))...         % malate shuttle/DH, malate -> oa_c
      - k13 * oa_c * (1 + B_glucgn^ep7/(kDglucgn^ep7+B_glucgn^ep7)) * kDins^en5/(kDins^en5 + B_ins^en5) ...     % pep carboxylase, oa_c --> pep
      + k14 * citrate * (1 + B_ins^ep10/(kDins^ep10+B_ins^ep10)) * kDglucgn^en7/(kDglucgn^en7 + B_glucgn^en7);     % citrate lyase, citrate -> acet_c, oa_c

dacet_c = k14 * citrate * (1 + B_ins^ep10/(kDins^ep10+B_ins^ep10)) * kDglucgn^en7/(kDglucgn^en7 + B_glucgn^en7) ...     % citrate lyase, citrate -> acet_c, oa_c
        - k15 * acet_c * (1 + p1*citrate/(kp1 + citrate)) * ki5/(ki5 + palmCoA) ...      % acet carboxylase, acet_c -> malonyl
        - k16 * acet_c * malonyl^7;   % fas, acet_c + malonyl -> palm
    
dmalonyl = k15 * acet_c * (1 + p1*citrate/(kp1 + citrate)) * ki5/(ki5 + palmCoA) ...     % acet carboxylase, acet_c -> malonyl
         - 7 * k16 * acet_c * malonyl^7;  % fas, acet_c + malonyl -> palm
     
dpalm = k16 * acet_c * malonyl^7 ...  % fas, acet_c + malonyl -> palm
      - k17 * palm * ki1/(ki1 + malonyl) ...          % acyl synth/CPT, palm --> palmCoA
      + ktL5 * B_ffa;       % palm transport (bl to liv), B_palm -> palm
  
dpalmCoA = k17 * palm * ki1/(ki1 + malonyl) ...       % acyl synth/CPT, palm --> palmCoA
         - k18 * palmCoA * ki2/(ki2 + malonyl);       % beta-oxidation, palmCoA -> acet_m
     
dket = k19 * acet_m^2 * (1 + cAMP^ep11/(kdcAMP^ep11+cAMP^ep11)) * kDins^en8/(kDins^en8 + B_ins^en8) ...         % HMG synthase/lyase, acet_m -> ket
     - ktL3 * ket;               % ket transport (liv to bl) ket -> B_ket

dalan = ktL6 * B_alan * (kDins2)^en10/((kDins2)^en10 + B_ins^en10) ...   % alan transport, B_alan -> alan
      + k21 ./(km21p + pyr) ./ (km21g + glutamate) * pyr * glutamate * kDins^en3/(kDins^en3 + B_ins^en3) ...    % pyr,glutamate -> alan,aK
      - k_21 ./(km_21a + alan) ./ (km_21k + aK) * alan * aK * acet_m * kDins^en3/(kDins^en3 + B_ins^en3) ...    % alan,aK -> pyr,glutamate
      - k20 * alan * kDins2^en3/(kDins2^en3 + B_ins^en3);           % alternate pathway, alan -> pyr
 
dcAMP = kc1 * B_glucgn^10/(kcm1^10 + B_glucgn^10) ...
      - kc2 * B_ins^10/(kcm2^10+B_ins^10) * cAMP;
 
dglutamate = - k21 ./(km21p + pyr) ./ (km21g + glutamate) * pyr * glutamate * kDins^en3/(kDins^en3 + B_ins^en3) ...           % pyr,glutamate -> alan,aK
             + k_21 ./(km_21a + alan) ./ (km_21k + aK) * alan * aK * acet_m * kDins^en3/(kDins^en3 + B_ins^en3) ...   % alan,aK -> pyr,glutamate
             - k22 * glutamate;
         
%dco2 = 0;

f = [dgluc; dg6p; dglycgn; dpep; dpyr; dlac; doa_m; dacet_m; dcitrate; 
    daK; dmalate; doa_c; dacet_c; dmalonyl; dpalm; 
    dpalmCoA; dket; dalan; dcAMP; dglutamate];
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
