
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
  [num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','PolymerNames','A2:A47');%read in list of polymer names from excel sheet
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
  [num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','SolventNames','A2:A325');
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
[num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','PolymerNames','A2:A48');
for i = 1:47
  if i == polymerchoice
    X = sprintf('The polymer you have chosen is:%s', txt{i});
    disp(X)
  end
end
%display subgroups of chosen polymer using xlsread
polymergroups = zeros(1,10);
[num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','PolymerNames','B2:K20');
for j = 1:10
    polymergroups(1,j) = num(polymerchoice,j);
end
disp(polymergroups)
%display user choice of solvent based on id chosen using xlsread

[num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','SolventNames','A2:A325');
for j = 1:325
  if j == solventchoice
    Y = sprintf('The solvent you have chosen is:   %s', txt{j});
    disp(Y)
    end
end
%displays subgroup numbers of chosen solvent
solventgroups = zeros(1,10);
[num,txt,raw] = xlsread('Interaction Parameters Sheet.xlsx','SolventNames','B2:E13');
for j = 1:4
    solventgroups(j) = num(solventchoice,j);
end
disp(solventgroups)

    
%Ask user to define temperature
temperature = input('Please enter the temperature of your system in Kelvin. The reference temperature used will be 273.15 K: ');

%send subgroups and id's to vi function
disp('*****SUBGROUP DETERMINATION*****')
%pull up specific interaction values for each subgroup based on id of subgroup
k=1;
for i = 1:2:length(polymergroups) %parse through all subgroups for polymer
    polymertemparray{k} = polymergroups(1,i:i+1);%break polymergroup array into pairs of 2 arrays
    k = k+1;
end

l=1;
for i = 1:2:length(solventgroups)
    solventtemparray{l} = solventgroups(1,i:i+1);
    l=l+1;
end


%separates num of each subgroup from the temp array subgroup
for j = 1:5
     polymernumofsubgrouparray(j) = polymertemparray{1,j}(1,1); %move the numbers of each subgroup into a separate array
end
 
for k = 1:5
    solventnumofsubgrouparray(k) = solventtemparray{1,k}(1,1);
end
%sends array to function that calculates vi*
%ReferenceVolumeParameter;
refvolumeparameter(polymertemparray,temperature,polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray)
epsilonii(polymertemparray,temperature,polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray)
%endfunction


end
%% calculating vi*
function refvolumeparameter(polymertemparray,temperature,polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray)
  reftemp = 273.15;
  %disp(length(subgrouparray));
  polymerRkvaluescell = cell([1 5]);
  solventRkvaluescell = cell([1,5]);
  %selection of correct subgroup numbers
  [num,txt,excelvalues] = xlsread('Interaction Parameters Sheet.xlsx','Subgroup Parameters','F1:H23');
 %select specific Rk values for the user specified subgroup and store in array
  subgrouplist = 1:47;
  %loop through total subgroup list
   for j = 1:5 %loop through user defined subgroup list
     for i = 1:47
        if subgrouplist(i) == polymertemparray{1,j}(1,2)% if a number from the total list matches one of the numbers from the user defined subgroup list
            polymerRkvaluescell{1,j}= cell2mat(excelvalues(polymertemparray{1,j}(1,2),:));%accumulate all excel values corresponding to the selected subgroup in this cell array
        end
     end 
   end
 
[num,txt,excelvalues] = xlsread('Interaction Parameters Sheet.xlsx','Subgroup Parameters','F1:H23');
for j = 1:5
    for i = 1:7
        if i == solventtemparray{1,j}(1,2)
            solventRkvaluescell{1,j} = cell2mat(excelvalues(solventtemparray{1,j}(1,2),:));
        end
    end
end

Rkpolymers = zeros(1,length(polymerRkvaluescell));
Rksolvent = zeros(1,length(solventRkvaluescell));
%Rkvaluescell is now a matrix


%Perform Rk calculation
 for i = 1:3
     tempRk = polymerRkvaluescell{i}; %temp assign each cell of the cell array to a variable
     Rkpolymers(i) = tempRk(1) + tempRk(2)*(temperature/reftemp) + tempRk(3)*(temperature/reftemp)^2; %loop through the specific cell and calculate its rk value  
 end
 Rkpolymers = Rkpolymers ./1000;
 disp('The individual subgroup parameters for the polymers will appear below')
 fprintf('Rk = %f\n',Rkpolymers)
%  disp('Note: Each line of the above output corresponds to the Rk values for each subgroup')
%  
 for i = 1:2
     tempRk=solventRkvaluescell{i};
     Rksolvent(i) = tempRk(1) + tempRk(2)*(temperature/reftemp) + tempRk(3)*(temperature/reftemp)^2;
 end
 
 Rksolvent = Rksolvent ./1000;
 disp('The individual subgroup parameters for the solvent will appear below')
 fprintf('Rk = %f\n',Rksolvent)
%  disp('Note: Each line of the above output corresponds to the Rk values for each subgroup')
 
 %calculation of vi*
 
vipolymer = sum(polymernumofsubgrouparray .* Rkpolymers); %multiply each subgroup number by each element of the Rk array
fprintf('vi*polymer = %.4f\n',vipolymer) %display the final value

visolvent = sum(solventnumofsubgrouparray .* Rksolvent);
fprintf('vi*solvent = %.4f\n',visolvent)
 
end
%% calculating epsilonii
function epsilonii(polymertemparray,temperature,polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray)
reftemp = 273.15;
%calculate ekk values first
polymerekvaluescell = {length(polymertemparray),1};
solventekvaluescell = {length(solventtemparray),1};
[num,txt,excelvalues] = xlsread('Interaction Parameters Sheet.xlsx','Subgroup Parameters','C1:E23');
%disp(excelvalues)
for i = 1:47 %loop through subgroupvalues
    for j = 1:length(polymertemparray)
        if i == polymertemparray{1,j}(1,2)
            polymerekvaluescell{j} = cell2mat(excelvalues(polymertemparray{1,j}(1,2),:));
        end
    end
end

for i = 1:47 %loop through subgroupvalues
    for j = 1:length(solventtemparray)
        if i == solventtemparray{1,j}(1,2)
            solventekvaluescell{j} = cell2mat(excelvalues(solventtemparray{1,j}(1,2),:));
        end
    end
end

ekpolymer = zeros(1,length(polymerekvaluescell));
eksolvent = zeros(1,length(solventekvaluescell));
%perform ek calculation
for i=1:length(polymerekvaluescell)
    tempek = polymerekvaluescell{i};
    %disp(tempek)
    ekpolymer(i) = tempek(1)+ tempek(2) * (temperature/reftemp) + tempek(3) * (temperature/reftemp)^2;
end

for i=1:length(solventekvaluescell)
    tempek = solventekvaluescell{i};
    %disp(tempek)
    eksolvent(i) = tempek(1)+ tempek(2) * (temperature/reftemp) + tempek(3) * (temperature/reftemp)^2;
end


polyzeroarray = zeros(1,2);
ekpoly = [ekpolymer, polyzeroarray];
disp('')
 fprintf('ek polymer = %f\n',ekpoly)
solvzeroarray = zeros(1,3);
disp('Note: Each line of the above output corresponds to the ek values for each subgroup in the polymer')
eksolv = [eksolvent,solvzeroarray];
fprintf('ek solvent = %f\n',eksolv)
disp('Note: Each line of the above output corresponds to the ek values for each subgroup in the solvent')
%calculating Theta(k,i)

polymerQkvalues = {1,5};
solventQkvalues = {1,5};
[num,txt,excelvalues] = xlsread('Interaction Parameters Sheet.xlsx','Subgroup Parameters','I1:I23');
%disp(excelvalues)
for i = 1:23
    for j = 1:length(polymertemparray)
        if i == polymertemparray{1,j}(1,2)
            if polymertemparray{1,j}(1,2) == 0.00
                polymerQkvalues{j} = 0;
            else
            polymerQkvalues(j) = excelvalues(polymertemparray{1,j}(1,2),1);
            end
        end
    end
end

for i = 1:23
    for j = 1:length(solventtemparray)
        if i == solventtemparray{1,j}(1,2)
            if solventtemparray{1,j}(1,2) == 0.00
                solventQkvalues(j) = 0;
            else
            solventQkvalues(j) = excelvalues(solventtemparray{1,j}(1,2),1);
            end
        end
    end
end

%zeroarraypolymer = zeros(1,2); %since size of Qkvalues and numofsubgrouparray are different,
%concatenate 2 empty spaces onto end of Qkvalues to allow for matrix operations
polymerQkvalues = [polymerQkvalues zeros(1,2)];
%disp(polymerQkvalues)
polymerQkvalues = cell2mat(polymerQkvalues);
%polymerQkvalues = [polymertemparray polymerQkvalues]; %concatenate these two vectors
zeroarraysolvent = zeros(1,3);
solventQkvalues = [solventQkvalues zeroarraysolvent]; %concatenate these two vectors
solventQkvalues = cell2mat(solventQkvalues);

Thetakpolymer = (polymernumofsubgrouparray .* polymerQkvalues)/(sum(polymernumofsubgrouparray .* polymerQkvalues));
fprintf('Theta(i,k) polymer = %f\n',Thetakpolymer)
% disp('Note: Each line of the above output corresponds to the Theta(i,k) values for each subgroup in the polymer')
epsilonpolymer = sum(Thetakpolymer .* ekpoly);
fprintf('epsiloniipolymer = %f \n',epsilonpolymer)

Thetaksolvent = (solventnumofsubgrouparray .* solventQkvalues)/(sum(solventnumofsubgrouparray .* solventQkvalues));
fprintf('Theta(i,k) solvent = %f\n',Thetaksolvent)
epsilonsolvent = sum(Thetaksolvent .* eksolv);
fprintf('epsiloniisolvent = %f \n',epsilonsolvent)



Thetaksolvent = (solventnumofsubgrouparray .* solventQkvalues)/(sum(solventnumofsubgrouparray .* solventQkvalues));
% fprintf('Theta(i,k) = %f\n',Thetaksolvent)
% disp('Note: Each line of the above output corresponds to the Theta(i,k) values for each subgroup in the solvent');
epsilonsolvent = sum(Thetaksolvent .* eksolv);
% fprintf('epsiloniisolvent = %f \n',epsilonsolvent)

%Calculate kij
kij(polymertemparray,temperature,polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray,polymerQkvalues,solventQkvalues);


end
%% determining which method to calculate kij
function kij(polymertemparray,temperature,polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray,polymerQkvalues,solventQkvalues)
disp('There are two methods to calculate kij')
disp('1. kij is based on surface area fractions of group m in the mixture (Lee-Danner method).')
disp('2. kij is based solely on the interactions of groups of unlike species (Hamedi et al. method)')

kijmethodselection = input('Which method would you like to use? Please enter 1 or 2:');

if(kijmethodselection == 1)
    kijmethod1(polymertemparray, polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray,temperature,polymerQkvalues,solventQkvalues);
elseif(kijmethodselection==2)
    kijmethod2(polymertemparray, polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray,temperature,polymerQkvalues,solventQkvalues);
end

end

%% calculating kij using Lee-Danner Method
function kijmethod1(polymertemparray, polymernumofsubgrouparray,solventtemparray,solventnumofsubgrouparray,temperature,polymerQkvalues,solventQkvalues)
disp('You have chosen the Lee-Danner method for calculating kij')

combinedsolpolyarray = [polymertemparray; solventtemparray];
disp('combined sol polymer array')
% for i = 1:10
%     disp(combinedsolpolyarray{i})
% end

%disp(polymerQkvalues)


for j = 1:5
     polymeridarray(j) = polymertemparray{1,j}(1,2); %move the subgroup id of each subgroup into a separate array
end

polyQkarray = vertcat(polymeridarray,polymernumofsubgrouparray,polymerQkvalues);
%disp(polyQkarray)

for j = 1:5
     solventidarray(j) = solventtemparray{1,j}(1,2); %move the numbers of each subgroup into a separate array
end

solvQkarray = vertcat(solventidarray,solventnumofsubgrouparray,solventQkvalues);
%disp(solvQkarray)


combinedarray = horzcat(polyQkarray,solvQkarray);
disp(combinedarray)
%temp = zeros(3,1);
index=0;
x = combinedarray(1,:);
disp(length(combinedarray))
disp(x)
 for i = 1:23   
    for k = 1:length(combinedarray)
        %disp(combinedarray(1,k))
        if i == x(k)
            %disp('in')
            index = find(x == i);
            %disp(index)
            temp{i} = combinedarray(1:3,index);
            
                
        end
            %disp('not in if')
    end
 end

    
end
%% calculating kij using Hamedi et al. method
function kijmethod2(polymertemparray, polymernumofsubgrouparray,solventtemparray,solventnumofsubggrouparray,temperature,polymerQkvalues,solventQkvalues)
disp('You have chosen the Hamedi et al. method for calculating kij')

end



  