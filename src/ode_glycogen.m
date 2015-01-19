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
% ode_glycogen.m
% Called by ode_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f=ode_glycogen(t,y,gluc,g6p,glycgn,run)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters/variables:

par=par_glycogen();
kg1=par(1); kg2=par(2); kg3=par(3); kg4=par(4); kg5=par(5);
kg6=par(6); kg7=par(7); kg8=par(8); kmg1=par(9); kmg2=par(10);
kmg3=par(11); kmg4=par(12); kmg5=par(13); kmg6=par(14); kmg7=par(15);
kmg8=par(16); k_a=par(17); ka=par(18); kgc_1=par(19); kgc1=par(20); 
kgc_2=par(21); kgc2=par(22); s1=par(23); kg2=par(24); s2=par(25);
kgi=par(26); capkt = par(27); kt = par(28); pt = par(29); 
st = par(30); PP1t = par(31);

%%%%%%%%%%%%%%%%%%%%%
%  assign values

PP1=y(1); PP1_GPa=y(2); PKa=y(3); GPa=y(4); GSa=y(5); 
R2C2=y(6); C=y(7); R2_C_cAMP2=y(8); R2_cAMP4=y(9); cAMP=y(10);

%% %%%%%%%%%%%%%%%%%%%
 %kd is regulated by glycogen

    c0 = 5.0;
    index = 5.;
    cmax = 3.2/1000;
    cmin =0.002/1000;
    kd = (cmax - cmin)*c0^index/(c0^index+glycgn^index) + cmin;   
    k_a = .001 *60*1000;    % min-1 mM-1 dissociation rate constant for PP1_GPa
    ka = k_a/kd;    % min-1 association rate constant for PP1_GPa
 
    kmg5s = kmg5.*(1+s1*g6p./kg2); 
    kmg6s = kmg6./(1+s2*gluc./kgi);
    kmg7s = kmg7.*(1+s1*g6p./kg2);
    kmg8s = kmg8./(1+s1*g6p./kg2);
    
    Jg1 = 	0;
    Jg2 = 	0;
    Jg3 = 	kg3 .* C .* (kt - PKa) ./ (kmg3 + (kt - PKa));
    Jg4 = 	kg4 .* (PP1 + PP1_GPa) .* PKa ./ (kmg4 + PKa);
    Jg5 = 	kg5 .* PKa .* (pt - GPa) ./ (kmg5s + (pt - GPa));
    Jg6 = 	kg6 .* (PP1 + PP1_GPa) .* GPa ./ (kmg6s + GPa);
    Jg7 = 	kg7 .* (PKa + C) .* GSa ./ (kmg7s + GSa);
    Jg8 = 	kg8 .* (PP1) .* (st - GSa) ./ (kmg8s + (st - GSa));
    Jg9 = 	ka .* PP1 .* GPa  -  k_a .* PP1_GPa;
    Jg10 = 	kgc1 * R2C2 .* cAMP.^2 - kgc_1 * R2_C_cAMP2 .* C;
    Jg11 = 	kgc2 * R2_C_cAMP2 .* cAMP.^2 - kgc_2 * R2_cAMP4 .* C;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLYCOGEN ode eqns:

dPP1 = - Jg9; 
dPP1_GPa = Jg9; 
dPKa = Jg3 - Jg4; 
dGPa = Jg5 - Jg6 - Jg9; 
dGSa = Jg8 - Jg7; 
dR2C2 = - Jg10; 
dC = Jg10 + Jg11; 
dR2_C_cAMP2 = Jg10 - Jg11; 
dR2_cAMP4 = Jg11; 
dcAMP = - 2*Jg10 - 2*Jg11;

f = [dPP1; dPP1_GPa; dPKa; dGPa; dGSa; 
    dR2C2; dC; dR2_C_cAMP2; dR2_cAMP4; dcAMP];

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%