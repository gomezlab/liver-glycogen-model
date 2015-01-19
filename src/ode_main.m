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
% ode_main.m
% Called by main.m
% All ode's are grouped according to tissues.  All are called here in the
% ode_main.m file and output is returned in one main vector, f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = ode_main(t,y,run,gluc_input)

% liver:
gluc=y(1); g6p=y(2); glycgn=y(3); pep=y(4); pyr=y(5); 
lac=y(6); oa_m=y(7); acet_m=y(8); citrate=y(9); aK=y(10); 
malate=y(11); oa_c=y(12); acet_c=y(13); malonyl=y(14); palm=y(15); 
palmCoA=y(16); ket=y(17); alan=y(18); cAMP=y(19); glutamate=y(20);


% blood:
B_gluc=y(21); B_ins=y(22); B_glucgn=y(23); 
B_ffa=y(24); B_lac=y(25); B_ket=y(26); B_alan=y(27);


% fat:
F_g6p=y(28); F_acyl=y(29); F_TG=y(30); F_ffa=y(31);

% muscle:
M_g6p=y(32); M_glyc=y(33); M_pyr=y(34); M_lac=y(35);
M_ket=y(36); M_alan=y(37);

% glycogen:
PP1=y(38); PP1_GPa=y(39); PKa=y(40); GPa=y(41); GSa=y(42); 
R2C2=y(43); C=y(44); R2_C_cAMP2=y(45); R2_cAMP4=y(46);

% call ode files for each tissue:
f_liver = ode_liver(t,y(1:20),B_gluc,B_ffa,B_ins,B_alan,B_lac,B_glucgn,GSa,GPa,run);
f_blood = ode_blood(t,y(21:27),gluc,lac,ket,F_ffa,M_alan,M_lac,M_ket,run,gluc_input);
f_fat = ode_fat(t,y(28:31),B_gluc,B_ins);
f_muscle = ode_muscle(t,y(32:37),B_gluc,B_ket,B_ins);

% glycogen regulatory circuit ode:
dcAMP1 = f_liver(19);   % involved in ode_liver.m and ode_glycogen.m, combined on line 38
f_glycogen = ode_glycogen(t,[y(38:46);cAMP],gluc,g6p,glycgn,run);
dcAMP2 = f_glycogen(end);

f_glycogen=f_glycogen(1:end-1);
f_liver(19) = dcAMP1 + dcAMP2;

f = [f_liver; f_blood; f_fat; f_muscle; f_glycogen];    % save all output in one vector

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%