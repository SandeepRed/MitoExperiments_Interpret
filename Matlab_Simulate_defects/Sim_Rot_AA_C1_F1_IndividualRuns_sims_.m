%%%Also change the number of simulations to 50 or 1000 in Beard_NC_simulatePopulation
%%%Change the outputs in Beard_NC_simulatePopulation to median or all
%%%outputs

load("All_Par_with_ROS_2_Percent.mat")

%Define Paramters to vary, 6 different combinations, 8 paramters

Par_vary = ones(6,8);
Par_vary(2,2) = All_Par_with_ROS_2_Percent(15,2) %CI 70 Percent
Par_vary(3,5) = All_Par_with_ROS_2_Percent(5,5)	 %FI 90 Percent
Par_vary(4,6) = All_Par_with_ROS_2_Percent(30,6) %Hle 40 Percent



Par_vary(5,2) = All_Par_with_ROS_2_Percent(15,2)
Par_vary(5,5) = All_Par_with_ROS_2_Percent(5,5)

Par_vary(6,2) = All_Par_with_ROS_2_Percent(15,2)
Par_vary(6,5) = All_Par_with_ROS_2_Percent(5,5)
Par_vary(6,6) = All_Par_with_ROS_2_Percent(30,6)





percentage_cells = 1; %All cells with defect

k = 1;

for  i = [1,2] %Appropriate experiment, For simulating OCR, AA, Rot
    
    
    for j = 1:6

%Appropriate experiment, For simulating OCR, AA, Rot
%Save median OCR for each defect

  %%%Change the outputs in Beard_NC_simulatePopulation to median      

%     [OCR(:,k), delPsi_basal, ATP_m_basal, NADH_m_basal, ATP_c_basal,H2O2_basal,...
%     delPsi_drug1, ATP_m_drug1, NADH_m_drug1, ATP_c_drug1,...
%     delPsi_drug2, ATP_m_drug2, NADH_m_drug2, ATP_c_drug2,...
%     delPsi_FC_drug1, delPsi_FC_drug2, ATP_m_FC_drug1, ATP_m_FC_drug2,...
%     NADH_m_FC_drug1, NADH_m_FC_drug2] = Beard_NC_simulatePopulation(i, Par_vary(j,:), percentage_cells)


%Save all the simulations basal and drug delpSi for each defect
%%%Change the outputs in Beard_NC_simulatePopulation to all outputs

    [OCR, delPsi_basal(:,k), ATP_m_basal, NADH_m_basal, ATP_c_basal,H2O2_basal,...
    delPsi_drug1(:,k), ATP_m_drug1, NADH_m_drug1, ATP_c_drug1,...
    delPsi_drug2, ATP_m_drug2, NADH_m_drug2, ATP_c_drug2,...
    delPsi_FC_drug1(:,k), delPsi_FC_drug2, ATP_m_FC_drug1, ATP_m_FC_drug2,...
    NADH_m_FC_drug1, NADH_m_FC_drug2] = Beard_NC_simulatePopulation(i, Par_vary(j,:), percentage_cells)

    
    k = k+1
    
    
    end
end

