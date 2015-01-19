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
% Description:
% ode_blood.m
% Called by ode_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = ode_blood(t,y,gluc,lac,ket,F_ffa,M_alan,M_lac,M_ket,run,glucinput)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters/variables:

par=par_blood();
kf=par(1); ktL1=par(2); ktL2=par(3); ktL3=par(4); ktL4=par(5); 
ktL5=par(6); ktL6=par(7); ktF1=par(8); ktF3=par(9);
ktS1=par(10); ktS2=par(11); ktS3=par(12); ktS4=par(13);
kd_Bgluc=par(14); kd_Bins=par(15); kd_Bglucgn=par(16); kd_Bffa=par(17);
kd_Blac=par(18); kd_Bket=par(19); kd_Balan=par(20);
k1ins=par(21); k_mIns=par(22); ni=par(23); k1glucgn=par(24); 
k_mGlgn=par(25); ng=par(26); k_ins=par(27); k_glucgn=par(28);
Imax=par(29); Gmax=par(30); Imin=par(31); Gmin=par(32); ep12=par(33);
ep13=par(34); en9=par(35); en10=par(36); en11=par(37);

par=par_liver();
kDins=par(50); kDglucgn=par(51); kDins2=par(74);

B_gluc=y(1); B_ins=y(2); B_glucgn=y(3);  
B_ffa=y(4); B_lac=y(5); B_ket=y(6); B_alan=y(7);
%sink=y(8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEDING pattern:
switch  run
    case 1
        kfeed=kf; % constant feeding
    case 2 
        if t<0
           kfeed = kf;     % constant feeding
        else
           u = 0;
           sd = 1/(kf*sqrt(2*pi));  % gives range of [0,kf]
           kfeed = exp(-(.014*t-u).^2./(2.*sd^2))./(sd.*sqrt(2.*pi));  % decrease feeding
           %kfeed = .25*cos(.01*t) + .25;       % oscillatory feeding
           %kfeed = .25*cos(.005*t) + .25;            % oscillatory feeding

        end
    case 3 
        kfeed = glucinput;  % resume constant feeding, input rate from call function
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLOOD ode eqns:

dB_gluc = kfeed ...        % feeding
        - ktL1 * (B_gluc - gluc) ...   % gluc transport, B_gluc <-> gluc
        - ktF1 * B_gluc * (1 + B_ins^ep12/(kDins^ep12+B_ins^ep12)) ...    % transport, B_gluc -> F_g6p
        - ktS1 * B_gluc * (1 + B_ins^ep13/(kDins^ep13+B_ins^ep13)) ...    % gluc transport, B_gluc -> M_gluc
        - kd_Bgluc * B_gluc;           % output of B_gluc

dB_ins = k_ins ...              % basal insulin secretion
       + k1ins*B_gluc^ni / (k_mIns^ni + B_gluc^ni)...    % pos reg by B_gluc
       - kd_Bins * B_ins;                                    % insulin decay

dB_glucgn = k_glucgn - k1glucgn*B_gluc^ng / (k_mGlgn^ng + B_gluc^ng)...     % neg reg by gluc
        - kd_Bglucgn * B_glucgn;                                    % glucgn decay

dB_ffa = - ktL5 * B_ffa ...         % palm transport (bl to liv), B_ffa -> palm
         + ktF3 * F_ffa^8 * (kDins2)^en9/(kDins2^en9 + B_ins^en9) ...    % transport, F_ffa -> B_ffa
        - kd_Bffa * B_ffa;          % output of B_ffa
        
dB_lac = ktS3 * M_lac ...     % lac transport, mus to blood
       - ktL2 * B_lac ...     % lac transport (bl to liv), B_lac --> lac
       - kd_Blac * B_lac;     % output of B_lac

dB_ket = ktL3 * ket ...         % ket transport (liv to bl) ket -> B_ket
       - ktS2 * B_ket ...       % ket transport, blood to mus
       - kd_Bket * B_ket;       % output of B_ket

dB_alan = ktS4 * M_alan* (.75.*10^(-6))^en11/((.75.*10^(-6))^en11 + B_ins^en11)...     % alan transport, M_alan -> B_alan
        - ktL6 * B_alan * (.75.*10^(-6))^en10/((.75.*10^(-6))^en10 + B_ins^en10) ...             % alan transport, B_alan -> alan
        - kd_Balan * B_alan;            % output of B_alan

%dsink = 0;
%       6* kd_Bgluc * B_gluc ...            % output of B_gluc
%        + 16* kd_Bffa * B_ffa ...          % output of B_ffa
%        + 3* kd_Blac * B_lac ...     % output of B_lac
%        + 4* kd_Bket * B_ket ...       % output of B_ket
%        + 3* kd_Balan * B_alan;            % output of B_alan
   
f = [dB_gluc; dB_ins; dB_glucgn; dB_ffa; dB_lac; dB_ket; dB_alan];

%f = [dB_gluc; dB_ins; dB_glucgn; dB_ffa; dB_lac; dB_ket; dB_alan; dsink];

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
