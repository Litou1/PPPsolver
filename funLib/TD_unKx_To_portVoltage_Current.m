function [Vin_t,Is_t]=TD_unKx_To_portVoltage_Current(sourceType,unkX,numVia)
% do addtional work to extract Vin and get the correct sign of current
if strcmp(sourceType,'voltage')
    % voltage at the port Vin
    Vin_t= squeeze(unkX(numVia , 1 , :)-unkX(end-numVia , 1 , :));
    % current at the port (-Is), take the negative to get Is
    Is_t= -squeeze(unkX(end , 1 , :));
else
    Vin_t=-squeeze(unkX(end , 1 , :));
    Is_t= squeeze(unkX(end , 1 , :));
end