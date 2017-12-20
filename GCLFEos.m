
%% GCLF-EoS main program
polymerchoice =0;
phasechoice ='';
predictionchoice = input('Would you like to use a binary (B) prediction or a ternary(T) prediction?','s');
if predictionchoice == 'B' || predictionchoice == 'b'
predictionchoice = 'B';

elseif predictionchoice =='T' || predictionchoice == 't'
predictionchoice = 'T';
end


polymerchoice = input('Please enter the id number of the polymer. Press 0 for a list of polymers.');
  
if polymerchoice == 0
      %polymerchoice = input('Please enter the id number of the polymer. Press 0 for a list of polymers.');
  [num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','PolymerNames','A1:A47');%read in list of polymer names from excel sheet
  disp(txt)%display only names
end
polymerchoice = input('Please enter the id number of the polymer. Press 0 for a list of polymers.');
if polymerchoice ~=0
  %polymerchoice = input('Please enter the id number of the polymer. Press 0 for a list of polymers.');
  phasechoice = input('Please specify whether you would like to use Vapor(V) data or Liquid(L) data?','s');
    if phasechoice == 'V'
      disp('You have chosen vapor data.');
      phasechoice ='V';
    elseif phasechoice =='L'
       disp('You have chosen liquid data.');
       phasechoice ='L';
    end
end

%User decides solvent
solventchoice = input('Please enter the id number of the solvent. Press 0 for a list of solvents.');
if solventchoice == 0
  [num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','SolventNames','A1:A325');
  disp(txt) %display solvent names
  solventchoice = input('Please enter the id number of the solvent. Press 0 for a list of solvents.');
elseif solventchoice ~=0
  solventchoice = input('Please enter the id number of the solvent. Press 0 for a list of solvents.');
end
  
  
  %Send data for Binary Prediction and Ternary Prediction functions
  if predictionchoice == 'B'
    %disp('In if.');
    %BinaryPrediction; %calls script file
    BPrediction(polymerchoice,phasechoice,solventchoice);
  elseif predictionchoice =='T'
    TPrediction(polymerchoice,phasechoice,solventchoice)
  end
    
  %% Binary Prediction
function BPrediction(polymerchoice,phasechoice,solventchoice)
  %disp('In BPrediction')

%disp(phasechoice)
%display user choice of polymer based on id chosen using xlsread
[num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','PolymerNames','A1:A47');
for i = 1:47
  if i == polymerchoice
    X = sprintf('The polymer you have chosen is:%s', txt{i});
    disp(X)
  end
end
%display user choice of solvent based on id chosen using xlsread
[num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','SolventNames','A1:A325');
for j = 1:325
  if j == solventchoice
    Y = sprintf('The solvent you have chosen is:   %s', txt{j});
    disp(Y)
    end
end

%Ask user to define temperature
temperature = input('Please enter the temperature of your system in Kelvin. The reference temperature used will be 273.15 K: ');

%send subgroups and id's to vi function
disp('*****SUBGROUP DETERMINATION*****')
totalnumberofsubgroups = input('Please enter the number of subgroups that are present in your molecule:');
idsubgroups = zeros(1,totalnumberofsubgroups); %creates 1D rowvector that contains id# of each subgroup
specificnumberofsubgroups = ones(1,totalnumberofsubgroups);%creates another 1D row vector that contains the number of each subgroup present within the molecule

%ask user what subgroups are present, and how many there are of each subgroup in the molecule
idinput = input('If you want to see the list of subgroups, enter 0: ');
subgroupcounter =0;
if idinput == 0
  [num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','Subgroup Parameters','B2:B10');
  disp(txt);

else
    for i = 1:totalnumberofsubgroups
    idsubgroups(i) = input(sprintf('Please enter the id number of subgroup #%d: ',i));
    specificnumberofsubgroups(i) = input(sprintf('Please enter the number of times that subgroup appears in your polymer molecule:'));
    %disp(specificnumberofsubgroups(i))
    subgroupcounter = subgroupcounter+1; %determines how many subgroups user entered to determine size of 2D array in referencevolumeparameter function
    end
end
 %disp(subgroupcounter);
subgrouparray = cat(1,idsubgroups,specificnumberofsubgroups); %concatenates these two arrays into one array
%disp(subgrouparray)


%sends array to function that calculates vi*
%ReferenceVolumeParameter;
refvolumeparameter(subgrouparray,temperature,subgroupcounter)
epsilonii(subgrouparray,temperature,subgroupcounter)
%endfunction

end
%% calculating vi*
function refvolumeparameter(subgrouparray,temperature,subgroupcounter)
  reftemp = 273.15;
  %disp(length(subgrouparray));
  Rkvaluescell = {length(subgrouparray),1};
  %Rkvaluesarray = zeros(subgroupcounter,3); %preallocate space for array
  %selection of correct subgroup numbers
  [num,txt,excelvalues] = xlsread('Interaction Parameters Sheet.xlsx','Subgroup Parameters','F1:H4');
  %disp(size(subgrouparray))
 %disp(size(Rkvalues))
 %select specific Rk values for the user specified subgroup and store in array
 
 for i = 1:10 %loop through total subgroup list
   for j = 1:length(subgrouparray) %loop through user defined subgroup list
     if i == subgrouparray(1,j)% if a number from the total list matches one of the numbers from the user defined subgroup list
       Rkvaluescell{j}= cell2mat(excelvalues(subgrouparray(1,j),:)); %accumulate all excel values corresponding to the selected subgroup in this cell array 
     end
   end
 end
Rk = zeros(1,length(Rkvaluescell));
%Rkvaluescell is now a matrix

%Perform Rk calculation
 for i = 1:length(Rkvaluescell)
     tempRk = Rkvaluescell{i}; %temp assign each cell of the cell array to a variable
     Rk(i) = tempRk(1) + tempRk(2)*(temperature/reftemp) + tempRk(3)*(temperature/reftemp)^2; %loop through the specific cell and calculate its rk value  
 end
 Rk = Rk ./1000;
 disp('The individual subgroup parameters will appear below')
 disp('')
 fprintf('Rk = %f\n',Rk)
 disp('Note: Each line of the above output corresponds to the Rk values for each subgroup')
 
 %calculation of vi*
 vi = sum(subgrouparray(2,:) .* Rk);
 fprintf('vi* = %f\n',vi)
 
end
%% calculating epsilonii
function epsilonii(subgrouparray,temperature,specificnumberofsubgroups)
reftemp = 273.15;
%calculate ekk values first
ekvaluescell = {length(subgrouparray),1};
[num,txt,excelvalues] = xlsread('Interaction Parameters Sheet.xlsx','Subgroup Parameters','C1:E4');
%disp(excelvalues)
for i = 1:10 %loop through subgroupvalues
    for j = 1:length(subgrouparray)
        if i == subgrouparray(1,j)
            ekvaluescell{j} = cell2mat(excelvalues(subgrouparray(1,j),:));
        end
    end
end

ek = zeros(1,length(ekvaluescell));
%perform ek calculation
for i=1:length(ekvaluescell)
    tempek = ekvaluescell{i};
    %disp(tempek)
    ek(i) = tempek(1)+ tempek(2) * (temperature/reftemp) + tempek(3) * (temperature/reftemp)^2;
end
disp('')
fprintf('ek = %f\n',ek)
disp('Note: Each line of the above output corresponds to the ek values for each subgroup')
%calculating Theta(k,i)
Qkvalues = {length(subgrouparray),1};
[num,txt,excelvalues] = xlsread('Interaction Parameters Sheet.xlsx','Subgroup Parameters','I1:I23');
%disp(excelvalues)
for i = 1:10
    for j = 1:length(subgrouparray)
        if i == subgrouparray(1,j)
            Qkvalues(j) = excelvalues(subgrouparray(1,j),1);
        end
    end
end
Qkvalues = cell2mat(Qkvalues);
Thetak = (subgrouparray(2,:) .* Qkvalues)/(sum(subgrouparray(2,:) .* Qkvalues));
disp('')
fprintf('Theta(i,k) = %f\n',Thetak)
disp('Note: Each line of the above output corresponds to the Theta(i,k) values for each subgroup')
epsilon = sum(Thetak .* ek);
fprintf('epsilonii = %f \n',epsilon)

end

  