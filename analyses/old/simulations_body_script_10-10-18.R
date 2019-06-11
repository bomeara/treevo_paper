#############################################################################
# analyses


#############################################
   # runLabel = An_Emp-DispBound 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 


intrinsicModel<-
extrinsicModel<-
   
   



#############################################
   # runLabel = An_Emp-Bound 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 






#############################################
   # runLabel = An_SimDispBound-DispBound 
   # treeSet = empirical-anolis_tree 
   # simTrait.Intrinsic = An_Emp-DispBound 
   # simTrait.Extrinsic = An_Emp-DispBound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # nSimTrait = 10 






#############################################
   # runLabel = An_SimBound-DispBound 
   # treeSet = empirical-anolis_tree 
   # simTrait.Intrinsic = An_Emp-Bound 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # nSimTrait = 10 






#############################################
   # runLabel = An_Emp-BrownMotion 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 






#############################################
   # runLabel = An_Emp-Disp 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 






#############################################
   # runLabel = An_Emp-TimeReg 
   # treeSet = empirical-anolis_tree 
   # empiricalTraitData = Anolis_Size_Data 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 






#############################################
   # runLabel = Aq_Emp-3Opt 
   # treeSet = empirical-Aquilegia_tree 
   # empiricalTraitData = Aquilegia_Nectar_Spur_Data 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 





#############################################
   # runLabel = Aq_Emp-BrownMotion 
   # treeSet = empirical-Aquilegia_tree 
   # empiricalTraitData = Aquilegia_Nectar_Spur_Data 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 






#############################################
   # runLabel = Aq_Sim3Opt-3Opt 
   # treeSet = empirical-Aquilegia_tree 
   # simTrait.Intrinsic = Aq_Emp-3Opt 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # nSimTrait = 10 






#############################################
   # runLabel = Aq_SimBM-3Opt 
   # treeSet = empirical-Aquilegia_tree 
   # simTrait.Intrinsic = Aq_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBM-BM 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBMpriorBiased 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Null 
   # prior = rexp_with_mean_*not*_at_true_sigmasq 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBM-Disp 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDisp-Disp 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Disp 
   # simTrait.Extrinsic = An_Emp-Disp 
   # doRun.Intrinsic = Pure_BM 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBound-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Bound 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDispBoundNear-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Bound_BoundByStartingState 
   # simTrait.Extrinsic = An_Emp-Bound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDispBoundMod-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Bound_BoundByMinValue 
   # simTrait.Extrinsic = An_Emp-Bound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDispBoundFar-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-Bound_BoundOneRangeAway 
   # simTrait.Extrinsic = An_Emp-Bound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimDispBound-DispBound 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-DispBound 
   # simTrait.Extrinsic = An_Emp-DispBound 
   # doRun.Intrinsic = BM_w/_LowerBound 
   # doRun.Extrinsic = Displacement 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBM-TimeReg 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Time-AutoRegressive_Model 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimTimeReg-TimeReg 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = An_Emp-TimeReg 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = Time-AutoRegressive_Model 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_SimBM-3Opt 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = Aq_Emp-BrownMotion 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 






#############################################
   # runLabel = Ideal_Sim3Opt-3Opt 
   # treeSet = Ideal-Simulated 
   # simTrait.Intrinsic = Aq_Emp-3Opt 
   # simTrait.Extrinsic = Null 
   # doRun.Intrinsic = 3-Optima 
   # doRun.Extrinsic = Null 
   # prior = standard_(uniform) 
   # idealTreeSets = c("Ideal-Balanced", "Ideal-Pectinate", "Ideal-Star") 
   # nTipSets = c(8, 16, 64) 
   # nSimTrait = 10 
