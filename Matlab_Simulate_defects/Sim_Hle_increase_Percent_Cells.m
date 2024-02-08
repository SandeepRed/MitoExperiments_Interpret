%For Baseline Simulations, Set all paramters to 1 and Percentage cells to 1

%%%Also change the number of simulations to 50 in Beard_NC_simulatePopulation
%%%Change the outputs in Beard_NC_simulatePopulation to median Serch median
%%%outputs



Par_vary = ones(1,8);
percentage_cells = 1;

load("All_Par_with_ROS_2_Percent.mat")

k = 1;


for  i = 2:10 
    
    Par_vary(6) = i


            


        [OCR_Sea(:,k),  delPsi_basal(k), ATP_m_basal(k), NADH_m_basal(k), ATP_c_basal(k), H2O2_c_basal(k),...
        delPsi_Oligo(k), ATP_m_Oligo(k), NADH_m_Oligo(k), ATP_c_Oligo(k),...
        delPsi_Oligo_FCCP(k), ATP_m_Oligo_FCCP(k), NADH_m_Oligo_FCCP(k), ATP_c_Oligo_FCCP(k),...
        delPsi_FC_Oligo(k), delPsi_FC_Oligo_FCCP(k), ATP_m_FC_Oligo(k), ATP_m_FC_Oligo_FCCP(k),...
        NADH_m_FC_Oligo(k), NADH_m_FC_Oligo_FCCP(k)...
        ] = Beard_NC_simulatePopulation(3, Par_vary, percentage_cells)
    
  
        
                            [Oxygen_1(:,k), delPsi_basal_1(k), ATP_m_basal_1(k), NADH_m_basal_1(k), ATP_c_basal_1(k),H2O2_c_basal_1(k),...
        delPsi_Rot(k), ATP_m_Rot(k), NADH_m_Rot(k), ATP_c_Rot(k),...
        delPsi_Rot_Oligo(k), ATP_m_Rot_Oligo(k), NADH_m_Rot_Oligo(k), ATP_c_Rot_Oligo(k),...
        delPsi_FC_Rot(k), delPsi_FC_Rot_Oligo(k), ATP_m_FC_Rot(k), ATP_m_FC_Rot_Oligo(k),...
        NADH_m_FC_Rot(k), NADH_m_FC_Rot_Oligo(k)...
        ] = Beard_NC_simulatePopulation(1, Par_vary, percentage_cells)
    
                                [Oxygen_2(:,k), delPsi_basal_2(k), ATP_m_basal_2(k), NADH_m_basal_2(k), ATP_c_basal_2(k),H2O2_c_basal_2(k),...
        delPsi_AA(k), ATP_m_AA(k), NADH_m_AA(k), ATP_c_AA(k),...
        delPsi_AA_Oligo(k), ATP_m_AA_Oligo(k), NADH_m_AA_Oligo(k), ATP_c_AA_Oligo(k),...
        delPsi_FC_AA(k), delPsi_FC_AA_Oligo(k), ATP_m_FC_AA(k), ATP_m_FC_AA_Oligo(k),...
        NADH_m_FC_AA(k), NADH_m_FC_AA_Oligo(k)...
        ] = Beard_NC_simulatePopulation(2, Par_vary, percentage_cells)


    
                                [Oxygen_5(:,k),  delPsi_basal_5(k), ATP_m_basal_5(k), NADH_m_basal_5(k), ATP_c_basal_5(k), H2O2_c_basal_5(k),...
        delPsi_FCCP(k), ATP_m_FCCP(k), NADH_m_FCCP(k), ATP_c_FCCP(k),...
        delPsi_FCCP_Rot(k), ATP_m_FCCP_Rot(k), NADH_m_FCCP_Rot(k), ATP_c_FCCP_Rot(k),...
        delPsi_FC_FCCP(k), delPsi_FC_FCCP_Rot(k), ATP_m_FC_FCCP(k), ATP_m_FC_FCCP_Rot(k),...
        NADH_m_FC_FCCP(k), NADH_m_FC_FCCP_Rot(k)...
        ] = Beard_NC_simulatePopulation(5, Par_vary, percentage_cells)
    
    
    
    
    
    


    
          k = k + 1;
          
              
       
    end


