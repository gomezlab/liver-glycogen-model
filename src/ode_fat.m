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
% ode_fat.m
% Called by ode_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = ode_fat(t,y,B_gluc,B_ins)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters/variables:

par = par_fat();
kf1=par(1); kf3=par(2); kf4=par(3); kf5=par(4);
ep14=par(5); ep15=par(6); en12=par(7);

par=par_liver();
kDins=par(50); kDglucgn=par(51);

par=par_blood();
ktL5=par(6); ktL6=par(7); ktF1=par(8); ktF3=par(9);
Imax=par(29); Gmax=par(30); Imin=par(31); Gmin=par(32); ep12=par(33);
ep13=par(34); en9=par(35); en10=par(36); en11=par(37);

F_g6p=y(1); F_acyl=y(2); F_TG=y(3); F_ffa=y(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAT ode eqns:

dF_g6p = ktF1 * B_gluc * (1 + B_ins^ep12/(kDins^ep12+B_ins^ep12))...       % transport, B_gluc -> F_g6p
       - kf1 * F_g6p * (1 + B_ins^ep14/(kDins^ep14+B_ins^ep14));         % F_g6p -> F_acyl
   
dF_acyl = kf1 * F_g6p * (1 + B_ins^ep14/(kDins^ep14+B_ins^ep14))...        % F_g6p -> F_acyl
        - kf3 * F_acyl * (1 + B_ins^ep15/(kDins^ep15+B_ins^ep15))...  % F_glyc3p, F_acyl -> F_TG
        + kf5 * F_ffa^3;         % F_ffa2 -> F_acyl
    
dF_TG = kf3 * F_acyl  * (1 + B_ins^ep15/(kDins^ep15+B_ins^ep15))...  % F_glyc3p, F_acyl -> F_TG
      - kf4 * F_TG * kDins^en12/(kDins^en12 + B_ins^en12);             % F_TG -> F_glyc, F_ffa
  
dF_ffa = 0;   % hold fat levels constant
       3* kf4 * F_TG * kDins^en12/(kDins^en12 + B_ins^en12)...          % F_TG -> F_glyc, F_ffa
       - 3 * kf5 * F_ffa^3 ...        % F_ffa2 -> F_acyl
       - 8 * ktF3 * F_ffa^8 * (.75*10^(-6))^en9/((.75*10^(-6))^en9 + B_ins^en9);         % transport, F_ffa2 -> B_ffa

f = [dF_g6p; dF_acyl; dF_TG; dF_ffa];

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
