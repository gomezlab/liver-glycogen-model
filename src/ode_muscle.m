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
% ode_muscle.m
% Called by ode_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = ode_muscle(t,y,B_gluc,B_ket,B_ins)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters/variables:

par=par_muscle();
ks1=par(1); ks_1=par(2); ks2=par(3); ks3=par(4); 
ks_3=par(5); ks4=par(6); ks_4=par(7); ks_dket=par(8);
ks5=par(9);

par=par_liver();
kDins=par(50);

par=par_blood();
ktS1=par(10); ktS2=par(11); ktS3=par(12); ktS4=par(13);
ep13=par(34); en9=par(35); en10=par(36); en11=par(37);

M_g6p=y(1); M_glycgn=y(2); M_pyr=y(3); M_lac=y(4);
M_ket=y(5); M_alan=y(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUSCLE ode eqns:

dM_g6p = ktS1 * B_gluc * (1 + B_ins^ep13/(kDins^ep13+B_ins^ep13)) ...     % gluc transport, B_gluc -> M_gluc
       - ks1 * M_g6p ...       % M_g6p -> M_glycgn
       + ks_1 * M_glycgn ...   % M_glycgn -> M_g6p
       - ks2 * M_g6p;          % M_g6p -> M_pyr

dM_glycgn = ks1 * M_g6p ...     % M_g6p -> M_glycgn
          - ks_1 * M_glycgn;    % M_glycgn -> M_g6p
      
dM_pyr = 2 * ks2 * M_g6p ...    % M_g6p -> M_pyr
       - ks3 * M_pyr ...        % M_pyr -> M_lac
       + ks_3 * M_lac ...       % M_lac -> M_pyr
       - ks4 * M_pyr ...        % M_pyr -> M_alan
       + ks_4 * M_alan * M_ket;    % M_alan -> M_pyr

dM_lac = ks3 * M_pyr ...        % M_pyr -> M_lac
       - ks_3 * M_lac ...       % M_lac -> M_pyr
       - ktS3 * M_lac;          % lac transport, mus to blood

dM_ket = ktS2 * B_ket ...       % ket transport, blood to mus
       ... - ktS4 * M_alan* M_ket * kis1/(kis1 + (B_ins-Imin)/Imax) ...     % alan transport, M_alan -> B_alan
       - ks_dket * M_ket;       % decay term for M_ket (no other output)

dM_alan = 0;    % hold alanine levels constant
        ks5 ...               % --> M_alan
        + ks4 * M_pyr ...        % M_pyr -> M_alan
        - ks_4 * M_alan * M_ket...      % M_alan -> M_pyr
        - ktS4 * M_alan* kDins^en11/(kDins^en11 + B_ins^en11);     % alan transport, M_alan -> B_alan

f = [dM_g6p; dM_glycgn; dM_pyr; dM_lac; dM_ket; dM_alan];

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
