%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% Uses UNC's Copyright and Permission Notice:	
% COPYRIGHT AND PERMISSION NOTICE	
% UNC Software:  << Liver_glycogen_model >>
% Copyright (C) 2009 The University of North Carolina at Chapel
% All rights reserved.
% The University of North Carolina at Chapel Hill (?UNC?) and the
% developers (?Developers?) of << Liver_glycogen_model>> (?Software?) give
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
% par_blood.m 
% called by ode_blood.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function par = par_blood()

% BLOOD
kf = .5;     % feeding to bloodstream, kf
Imax = 1.3 * 10^(-6); 
Imin = .7 * 10^(-6);
k_mIns = 8; 
ni = 10; 
Gmax = 50 * 10^(-9); 
Gmin = 30 * 10^(-9);
k_mGlgn = 8; 
ng = 10;

% LIVER/BLOOD transport:
ktL1 = 100;    % gluc <-> B_gluc
ktL2 = .1;    % B_lac -> lac
ktL3 = .1;    % ket -> B_ket
ktL4 = .1;    % %%%%%%%%%%%%%%%%%%%%%  DELETE LATER
ktL5 = .1;    % B_ffa -> palm
ktL6 = 1;    % B_alan -> alan
en10 = 10;          % ktL6 * B_alan * kDins^en10/(kDins^en10 + B_ins^en10)

% FAT/BLOOD transport
ktF1 = .01;       % B_gluc -> F_g6p
ep12 = 10;       %   ktF1 * B_gluc * (1 + B_ins^ep12/(kDins^ep12+B_ins^ep12))
%ktF3 = .01;      % F_ffa -> B_ffa
ktF3 = .008;
en9 = 20;      %   ktF3 * F_ffa^8 * kDins^en9/(kDins^en9 + B_ins^en9)

% MUSCLE/BLOOD transport
ktS1 = .01;     % B_gluc -> M_g6p
ep13 = 10;      %   ktS1 * B_gluc * (1 + B_ins^ep13/(kDins^ep13+B_ins^ep13))
ktS2 = 1;     % B_ket -> M_ket
ktS3 = .01;     % M_lac -> B_lac
ktS4 = 3;     % M_alan -> B_alan
en11 = 10;     %  ktS4 * M_alan* kDins^en11/(kDins^en11 + B_ins^en11)

% OUTFLOW of blood concentrations
kd_all = 0.015;
kd_Bgluc = kd_all;       % outflow of B_gluc
kd_Bins = kd_all;       % outflow of B_ins
kd_Bglucgn = kd_all;       % outflow of B_glucgn
kd_Bffa = kd_all;       % outflow of B_ffa
kd_Blac = kd_all;       % outflow of B_lac
kd_Bket = kd_all;       % outflow of B_ket
kd_Balan = kd_all;       % outflow of B_alan

k_ins = kd_Bins * Imin;
k1ins = kd_Bins * (Imax-Imin);
k_glucgn = kd_Bglucgn * Gmax;
k1glucgn = kd_Bglucgn * (Gmax-Gmin);

par=[kf; ktL1; ktL2; ktL3; ktL4; ktL5; ktL6;
    ktF1; ktF3;
    ktS1; ktS2; ktS3; ktS4;
    kd_Bgluc; kd_Bins; kd_Bglucgn; kd_Bffa;
    kd_Blac; kd_Bket; kd_Balan;
    k1ins; k_mIns; ni; k1glucgn; k_mGlgn; ng; k_ins; k_glucgn;
    Imax; Gmax; Imin; Gmin; ep12; ep13; en9; en10; en11];

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% kf=par(1); ktL1=par(2); ktL2=par(3); ktL3=par(4); ktL4=par(5); 
% ktL5=par(6); ktL6=par(7); ktF1=par(8); ktF3=par(9);
% ktS1=par(10); ktS2=par(11); ktS3=par(12); ktS4=par(13);
% kd_Bgluc=par(14); kd_Bins=par(15); kd_Bglucgn=par(16); kd_Bffa=par(17);
% kd_Blac=par(18); kd_Bket=par(19); kd_Balan=par(20);
% k1ins=par(21); k_mIns=par(22); ni=par(23); k1glucgn=par(24); 
% k_mGlgn=par(25); ng=par(26); k_ins=par(27); k_glucgn=par(28);
% Imax=par(29); Gmax=par(30); Imin=par(31); Gmin=par(32); ep12=par(33);
% ep13=par(34); en9=par(35); en10=par(36); en11=par(37);