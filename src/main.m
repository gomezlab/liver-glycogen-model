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
% Called by runmain.m. 1st run is for the system to reach fed
% steady state. The values at steady state are then used in the 2nd run as
% initial conditions. 2nd run is the post-absorption state with a decreasing
% glucose input rate (drops below 5% in 140 mins). 3rd run is the
% refeeding period. Glucose feeding rate is set in 'gluc_input' variable
% called by 'runmain.m'.'state' variable determines the fasting time in
% the 2nd run. If state = 'Fed', fasting time is set to be 250 mins.
% Otherwise, it is 1200 mins. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main(state,gluc_input)

    
    y0 = init_main;     % load initial conditions from init_main.m

% 1st run -- to steady state under constant feeding
    tmax = 20000;      
    run=1;      % determines what feeding level will be (see ode_blood.m)
    options = odeset('AbsTol',10^-12,'RelTol',10^-12);
    [t,y] = ode15s(@ode_main,[0:tmax/200:tmax],y0,options,run,gluc_input);

    yy = y;
    tt = t;

% 2nd run -- see ode_blood.m for feeding pattern (decreasing feeding is chosen)
    y0 = y(end,:);      % save output from last run as updated init conds
    y0(3)=120;    %glycgn, start at same glycgn concentration
    if state(1:3) == 'fed'
        tmax = 250;       % fasting time in fed livers
    else
        tmax = 1200;       % fasting time in fasted livers
    end
    run=2;              % determines what feeding level will be (see ode_blood.m)
    [t,y] = ode15s(@ode_main,[0:tmax/200:tmax],y0,options,run,gluc_input);
    yy = y;
    tt = t;

% 3rd run -- REFEED with constant feeding
    y0 = y(end,:);      % save output from last run as updated init conds
    tmax = 20;
    run=3;              % determines what feeding level will be (see ode_blood.m)
    [t,y] = ode15s(@ode_main,[0:tmax/400:tmax],y0,options,run,gluc_input);
    yy = [yy(1:end-1,:);y];
    y = yy;
    tt = [tt(1:end-1);t+tt(end)];
    t = tt;  % run 2 and 3 are saved in t and y 
    
    


%% %%%%%%%%%%%%%%  save  data files  %%%%%%%%%%%%%%%%%%%%%%%%%%%

chartemp = num2str(gluc_input);
chartemp(find(chartemp=='.')) = 'p';
filename = ['./results/y_',state,'_',chartemp];
save (filename,'t','y');
display(['gluc input rate = ',num2str(gluc_input),' in ', state, ' liver...finished']);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

