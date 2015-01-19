
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
%
% Description:  
% par_glycogen.m 
% called by ode_glycogen.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function par=par_glycogen()

k1 = 1.4 *60; % min-1 rate constant for phosphorylation of inhibitor [48]
k2 = 0.01 *60; % min-1 rate constant for dephosphorylation of inhibitor [assumed]
k3 = 20 *60; % min-1 rate constant for phosphorylation of phosphorylase kinase [assumed]
k4 = 5 *60; % min-1 rate constant for dephosphorylation of phosphorylase kinase [assumed]
k5 = 20 *60; % min-1 rate constant for phosphorylation of Phosphorylase [42]
k6 = 5 *60; % min-1 rate constant for dephosphorylation of Phosphorylase [49]
k7 = 20 *60; % min-1 rate constant for phosphorylation of glycogen synthase [assumed]
k8 = 0.05 *60; % min-1 rate constant for dephosphorylation of glycogen synthase [assumed]
k8 = 5*60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Michaelis-Menten constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
km1 = 5 /1000; % mM for inhibitor phosphorylation [48]
km2 = 0.7 /1000; % mM for dephosphorylation of Inhibitor [52]
km3 = 0.4 /1000; % mM for Phosphorylation of phosphorylase kinase [assumed]
km4 = 1.1 /1000; % mM for dephosphorylation of phosphorylase kinase [52]
km5 = 10 /1000; % mM for phosphorylation of phosphorylase [25]
km6 = 5 /1000; % mM for dephosphorylation of phosphorylase [47]
km7 = 15 /1000; % mM for phosphorylation of glycogen synthase [assumed]
km8 = 0.12 /1000; % mM for dephosphorylation of glycogen synthase [50]
%km8 =1.2/1000;
kd = 0.002 /1000; % dissociation of PP1 and phosphorylated PP1 Inhibitor, and also phosphorylase a with synthase PP1 [47]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Concentrations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
capkt = 0.25 /1000; % mM total R2C2 ie. cAMP dependent protein kinase, CAPK [3]
It = 1.8 /1000; % mM total Inhibitor concentration [3]
kt = 2.5 /1000; % mM total Phosphorylase kinase [35]
pt = 70 /1000; % mM total Glycogen Phosphorylase [3]
%pt = 80 /1000; % mM total Glycogen Phosphorylase [3]
st = 3 /1000; % mM total Glycogen synthase [3]
PP1 = 0.25 /1000; % mM PTPase 1 [33]
PP2A = 0.025 /1000; % mM PTPase 2 [3]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other parameters: (chosen as per various qualitative observations are in physiological ranges as given in
% [3,9,12,17,25,27,29,33,35,42,46-54])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k11 = 0.043 /1000; % mM Dissociation constant of cAMP [35]
k22 = 0.7 /1000; % mM Dissociation constant of cAMP [35]
ki = 100 /1000; % mM cAMP inhibition constant 
campt = 10 /1000; % mM maximum cAMP [3]
kg = 349500 /1000; % mM activation constant of glucose-6-phosphate for synthase PP1
g6pt = 700 /1000; % mM maximum glucose-6-Phosphate [33]
kgi = 10000 /1000; % mM activation constant of glucose for phosphorylase phosphatase
s1 = 100; % a multiplicative factor for glucose-6-phosphate effect on glycogen synthase dephosphorylation
kg2 = 500 /1000; % mM inhibition due to glucose-6-phosphate = 0.05 mM
s2 = 0.0010; % a multiplicative factor for glucose effect on phosphorylase phosphatase


k_a = .001 *60*1000;    % min-1 mM-1 dissociation rate constant for PP1_GPa
ka = k_a/kd;    % min-1 association rate constant for PP1_GPa

%k_gc1 = 1000 *60*1000^2;           % dissociation rate constant for [R2C(cAMP)2]+[C]
k_gc1 = 60*1000;
kgc1 = k_gc1/k11;    % association rate constant for [R2C(cAMP)2]+[C]
%k_gc2 = 1000 *60*1000^2;          % dissociation rate constant for [R2(cAMP)4]+[C]
k_gc2 = 1000*60;
kgc2 = k_gc2/k22;    % association rate constant for [R2(cAMP)4]+[C]

par = [k1; k2; k3; k4; k5; k6; k7; k8; km1; km2; km3; km4; km5; km6; km7; km8;
    k_a; ka; k_gc1; kgc1; k_gc2; kgc2; s1; kg2; s2; kgi; capkt; kt; pt; st; PP1];

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%