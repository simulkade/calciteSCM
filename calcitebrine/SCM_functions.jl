function hirasaki_model(;log_k=[0.30, 1.74, 1.62, 0.14, 0.5, 2.23, 1.0, 0.09, -0.64], 
                        dh=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return """
    SURFACE_MASTER_SPECIES
        Chalk_a    Chalk_aOH-0.75
        Chalk_c    Chalk_cOH+0.75
    SURFACE_SPECIES
        Chalk_aOH-0.75 = Chalk_aOH-0.75
        log_k   0.0
        Chalk_cOH+0.75 = Chalk_cOH+0.75
        log_k   0.0
        # calcium site
        Chalk_aOH-0.75 + H+ = Chalk_aOH2+0.25
        log_k   $(log_k[1]) # 0.30
        delta_h      $(dh[1])
        Chalk_aOH-0.75 + Ca+2 = Chalk_aOHCa+1.25
        log_k   $(log_k[2]) # 1.74
        delta_h      $(dh[2])
        Chalk_aOH-0.75 + Mg+2 = Chalk_aOHMg+1.25
        log_k   $(log_k[3]) # 1.62
        delta_h      $(dh[3])
        Chalk_aOH-0.75 + Na+ = Chalk_aOHNa+0.25
        log_k   $(log_k[4]) # 0.14
        delta_h      $(dh[4])
        # carbonate site
        Chalk_cOH+0.75 + OH- = Chalk_cO-0.25 + H2O
        log_k   $(log_k[5]) # 0.5
        delta_h      $(dh[5])
        Chalk_cOH+0.75 + CO3-2 = Chalk_cOHCO3-1.25
        log_k   $(log_k[6]) # 2.23
        delta_h      $(dh[6])
        Chalk_cOH+0.75 + SO4-2 = Chalk_cOHSO4-1.25
        log_k   $(log_k[7]) # 1.0
        delta_h      $(dh[7])
        Chalk_cOH+0.75 + HCO3- = Chalk_cOH2CO3-0.25
        log_k   $(log_k[8]) # 0.09
        delta_h      $(dh[8])
        Chalk_cOH+0.75 + Cl- = Chalk_cOHCl-0.25
        log_k   $(log_k[9]) # -0.64
        delta_h      $(dh[9])
    END
    """
end

function modified_heberling(;log_k=[0.501, 0.56, 1.68, 1.48, -0.05, 0.04, -0.07, 0.4, 0.56, 1.68, 1.48], 
                        dh=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return """
    SURFACE_MASTER_SPECIES
    Chalk_a Chalk_aOH1.5
    Chalk_c Chalk_cH0.5
    
    SURFACE_SPECIES
    Chalk_aOH1.5 = Chalk_aOH1.5
    -cd_music 0 0 0
    log_k 0.0
    Chalk_cH0.5 = Chalk_cH0.5
    -cd_music 0 0 0
    log_k 0.0
    
    #dummy reactions producing surface species with fractional charges
    # are these necessary? Is there a way around it?
    Chalk_aOH1.5 = Chalk_aOH-0.5 + 0.5 H+
    -cd_music -0.5 0 0
    log_k 20
    Chalk_cH0.5 = Chalk_c-0.5 + 0.5 H+
    -cd_music -0.5 0 0
    log_k 20
    
    #inner sphere reactions
    Chalk_aOH-0.5 + H+ = Chalk_aOH2+0.5
    -cd_music 1 0 0
    log_K $(log_k[1]) # 0.501
    delta_h      $(dh[1])
    
    #outer sphere reactions
    Chalk_aOH-0.5 + Na+ = Chalk_aOHNa+0.5
    -cd_music 0 1 0
    log_K $(log_k[2]) # 0.56
    delta_h      $(dh[2])
    Chalk_aOH-0.5 + Ca+2 = Chalk_aOHCa+1.5
    -cd_music 0 2 0
    log_K $(log_k[3]) # 1.68
    delta_h      $(dh[3])
    Chalk_aOH-0.5 + Mg+2 = Chalk_aOHMg+1.5
    -cd_music 0 2 0
    log_K $(log_k[4]) # 1.48
    delta_h      $(dh[4])
    Chalk_aOH2+0.5 + Cl- = Chalk_aOH2Cl-0.5
    -cd_music 0 -1 0
    log_K $(log_k[5]) # -0.05
    delta_h      $(dh[5])
    Chalk_aOH2+0.5 + HCO3- = Chalk_aOH3CO3-0.5
    -cd_music 0 -1 0
    log_K $(log_k[6]) # 0.04
    delta_h      $(dh[6])
    Chalk_aOH2+0.5 + HCO3- = Chalk_aOH2CO3-1.5+ H+ 
    -cd_music 0 -2 0
    log_K $(log_k[7]) # -7.07
    delta_h      $(dh[7])
    Chalk_aOH2+0.5 + SO4-2 = Chalk_aOH2SO4-1.5
    -cd_music 0 -2 0
    log_K $(log_k[8]) # 0.4
    delta_h      $(dh[8])
    Chalk_c-0.5 + Na+ = Chalk_cNa+0.5
    -cd_music 0 1 0
    log_K $(log_k[9]) #  0.56
    delta_h      $(dh[9])
    Chalk_c-0.5 + Ca+2 = Chalk_cCa+1.5
    -cd_music 0 2 0
    log_K $(log_k[10]) # 1.68
    delta_h      $(dh[10])
    Chalk_c-0.5 + Mg+2 = Chalk_cMg+1.5
    -cd_music 0 2 0
    log_K $(log_k[11]) # 1.48
    delta_h      $(dh[11])
    """
end

function modified_so(;log_k=[-3.52108, -2.61239, -2.72276, 14.5456, -24.73, 9.26094, 1.37606, 0.205842, -9.8144, -8.8383], 
    dh=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return """
    SURFACE_MASTER_SPECIES
        Chalk_a Chalk_aOH-0.667
        Chalk_c Chalk_cH+0.667
    SURFACE_SPECIES
        Chalk_cH+0.667 = Chalk_cH+0.667
        log_k 0.0 
        -cd_music 0 0 0 0 0
        Chalk_aOH-0.667 = Chalk_aOH-0.667 
        log_k 0.0 
        -cd_music 0 0 0 0 0
        Chalk_cH+0.667 = Chalk_c-0.333 + H+ 
        log_k $(log_k[1]) # -3.58
        delta_h      $(dh[1])
        -cd_music -1 0 0 0 0    
        Chalk_cH+0.667 + Ca+2 = Chalk_cCa+1.667 + H+ 
        log_k $(log_k[2]) # -2.8
        delta_h      $(dh[2])
        -cd_music -1 2 0 0 0
        Chalk_cH+0.667 + Mg+2 = Chalk_cMg+1.667 + H+ 
        log_k $(log_k[3]) # -2.6
        delta_h      $(dh[3])
        -cd_music -1 2 0 0 0
        Chalk_aOH-0.667 + H+ = Chalk_aOH2+0.333 
        log_k $(log_k[4]) # 12.85
        delta_h      $(dh[4])
        -cd_music 1 0 0 0 0
        Chalk_aOH-0.667 = Chalk_aO-1.667 + H+
        log_k $(log_k[5]) # -24.73
        delta_h      $(dh[5])
        -cd_music -1 0 0 0 0
        Chalk_aOH-0.667 + CO3-2 + H+ = Chalk_aHCO3-0.667 + OH-
        log_k $(log_k[6]) # 10.15
        delta_h      $(dh[6])
        -cd_music 0.6 -0.6 0 0 0  
        Chalk_aOH-0.667 + CO3-2 = Chalk_aCO3-1.667 + OH-
        log_k $(log_k[7]) # 1.55
        delta_h      $(dh[7])
        -cd_music 0.6 -1.6 0 0 0
        Chalk_aOH-0.667 + SO4-2 = Chalk_aSO4-1.667 + OH- 
        log_k $(log_k[8]) # 1.55
        delta_h      $(dh[8])
        -cd_music 0.6 -1.6 0 0 0
        Chalk_cH+0.667 + Na+ = Chalk_cNa+0.667+H+        log_k $(log_k[9])  #-9.11        delta_h      $(dh[9])        -cd_music -1 1 0 0 0        Chalk_aOH-0.667 + Cl- = Chalk_aCl-0.667 + OH-         log_k $(log_k[10])  #-5.0        delta_h      $(dh[10])        -cd_music 1 -1 0 0 0 

    """
end

function modified_wolthers(;log_k=[-3.58, -2.8, -2.2, 12.85, -24.73, 10.15, 1.55, 0.35], 
    dh=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    return """
    SURFACE_MASTER_SPECIES
        Chalk_a Chalk_aOH-0.667
        Chalk_c Chalk_cH+0.667
    SURFACE_SPECIES
        Chalk_cH+0.667 = Chalk_cH+0.667
        log_k 0.0 
        -cd_music 0 0 0 0 0
        Chalk_aOH-0.667 = Chalk_aOH-0.667 
        log_k 0.0 
        -cd_music 0 0 0 0 0
        Chalk_cH+0.667 = Chalk_c-0.333 + H+ 
        log_k $(log_k[1]) # -3.58
        delta_h      $(dh[1])
        -cd_music -1 0 0 0 0    
        Chalk_cH+0.667 + Ca+2 = Chalk_cCa+1.667 + H+ 
        log_k $(log_k[2]) # -2.8
        delta_h      $(dh[2])
        -cd_music -1 2 0 0 0
        Chalk_cH+0.667 + Mg+2 = Chalk_cMg+1.667 + H+ 
        log_k $(log_k[3]) # -2.6
        delta_h      $(dh[3])
        -cd_music -1 2 0 0 0
        Chalk_aOH-0.667 + H+ = Chalk_aOH2+0.333 
        log_k $(log_k[4]) # 12.85
        delta_h      $(dh[4])
        -cd_music 1 0 0 0 0
        Chalk_aOH-0.667 = Chalk_aO-1.667 + H+
        log_k $(log_k[5]) # -24.73
        delta_h      $(dh[5])
        -cd_music -1 0 0 0 0
        Chalk_aOH-0.667 + CO3-2 + H+ = Chalk_aHCO3-0.667 + OH-
        log_k $(log_k[6]) # 10.15
        delta_h      $(dh[6])
        -cd_music 0.6 -0.6 0 0 0  
        Chalk_aOH-0.667 + CO3-2 = Chalk_aCO3-1.667 + OH-
        log_k $(log_k[7]) # 1.55
        delta_h      $(dh[7])
        -cd_music 0.6 -1.6 0 0 0
        Chalk_aOH-0.667 + SO4-2 = Chalk_aSO4-1.667 + OH- 
        log_k $(log_k[8]) # 1.55
        delta_h      $(dh[8])
        -cd_music 0.6 -1.6 0 0 0

    """
end