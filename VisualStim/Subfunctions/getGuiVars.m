function [VarNames,VarVals] = getGuiVars(handles,nrTrials)
% short code to isolate variables from vStim GUI. This finds both statice
% and flexible variables are returns their names in a single cell array. 
% Also returns a matrix with variables values. Matrix size is of size
% variables x trials with equal amounts of all possible flexible case combinations.

if ~exist('nrTrials','var') || isempty(nrTrials)
    nrTrials = str2double(handles.NrTrials.String);
end

%% Isolate all relevant variables from MainGUI
StatNames = {}; FlexNames = {};
StatVals = []; FlexVals = {};

for iVars = 1:size(handles.StaticVariableNames.String,1)
    temp = textscan(handles.StaticVariableNames.String{iVars},'%s%f');
    StatNames{iVars} = temp{1}{1};
    StatVals(iVars) = temp{2};
    clear temp
end
for iVars = 1:size(handles.FlexibleVariableNames.String,1)
    temp = textscan(handles.FlexibleVariableNames.String{iVars},'%s%f');
    FlexNames{iVars} = temp{1}{1};
    FlexVals{iVars} = str2num(handles.FlexibleVariableNames.String{iVars}(length(temp{1}{1})+1:end));
    clear temp
end

VarNames = [StatNames FlexNames];  %all names  variables that are required for correct function
FlexCases = combvec(FlexVals{:}); %get combinations for flexible variables
if isempty(FlexCases)
    FlexCases = 1; %at least one case, even if there are no flexible variables
end
VarVals = zeros(length(VarNames),size(FlexCases,2)); %values of each  variable that is required. these can change based on the amount of cases that are produced.
    
for x = 1:length(FlexVals)
    VarVals(ismember(VarNames,FlexNames{x}),:) = FlexCases(x,:); % fill flexibe variable values from flexcases
end
for x = 1:length(StatNames)
    VarVals(ismember(VarNames,StatNames{x}),:) = repmat(StatVals(x),1,size(FlexCases,2)); % fill static variable values from StatVals
end

VarVals = repmat(VarVals,1,ceil(nrTrials/size(VarVals,2))); %produce enough cases to cover all trials
if handles.RandomTrials.Value %randomize order of trials in each block of cases
    ind = [];
    for x = 1:ceil(nrTrials/size(FlexCases,2))
        ind = [ind randperm(size(FlexCases,2))+size(FlexCases,2)*(x-1)];
    end
    VarVals = VarVals(:,ind);
end
