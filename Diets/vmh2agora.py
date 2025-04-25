#!/usr/bin/env python
# coding: utf-8

import json
import sys

####This script requires a VMH diet file downloaded from https://www.vmh.life/#nutrition in tab separated format as argument


raw_vmh_dict = {}

## conversion to BiGG identifier format of the vmh diet
with open(sys.argv[1]) as filer:
    for i in filer:
        content = i.strip().split('\t')
        content[0] = content[0].replace('[e]','_e')
        reac_name = content[0].replace('(e)','_e')
        
        if reac_name.endswith('_L_e'):
            reac_name = reac_name.replace('_L_e','__L_e')
        if reac_name.endswith('_D_e'):
            reac_name = reac_name.replace('_D_e','__D_e')
        if reac_name.endswith('_R_e'):
            reac_name = reac_name.replace('_R_e','__R_e')
                    
        if reac_name != 'Reaction':
            raw_vmh_dict[reac_name] = (float(content[1])*-1)/12.63   ###  FLUX scaling applied


##########################################################################
### ADDING Essential Metabolites to the vmh diet dictionary ##############

essentialMetabolites = ['EX_12dgr180_e', 'EX_26dap__M_e', 'EX_2dmmq8_e', 'EX_2obut_e', 'EX_3mop_e', 'EX_4abz_e', 'EX_4hbz_e', 'EX_ac_e', 'EX_acgam_e', 'EX_acmana_e', 'EX_acnam_e', 'EX_ade_e', 'EX_adn_e', 'EX_adocbl_e', 'EX_adpcbl_e', 'EX_ala__D_e', 'EX_ala__L_e', 'EX_amet_e', 'EX_amp_e', 'EX_arab__D_e', 'EX_arab__L_e', 'EX_arg__L_e', 'EX_asn__L_e', 'EX_btn_e', 'EX_ca2_e', 'EX_cbl1_e', 'EX_cgly_e', 'EX_chor_e', 'EX_chsterol_e', 'EX_cit_e', 'EX_cl_e', 'EX_cobalt2_e', 'EX_csn_e', 'EX_cu2_e', 'EX_cys__L_e', 'EX_cytd_e', 'EX_dad_2_e', 'EX_dcyt_e', 'EX_ddca_e', 'EX_dgsn_e', 'EX_fald_e', 'EX_fe2_e', 'EX_fe3_e', 'EX_fol_e', 'EX_for_e', 'EX_gal_e', 'EX_glc__D_e', 'EX_gln__L_e', 'EX_glu__L_e', 'EX_gly_e', 'EX_glyc_e', 'EX_glyc3p_e', 'EX_gsn_e', 'EX_gthox_e', 'EX_gthrd_e', 'EX_gua_e', 'EX_h_e', 'EX_h2o_e', 'EX_h2s_e', 'EX_his__L_e', 'EX_hxan_e', 'EX_ile__L_e', 'EX_k_e', 'EX_lanost_e', 'EX_leu__L_e', 'EX_lys__L_e', 'EX_malt_e', 'EX_met_L_e', 'EX_mg2_e', 'EX_mn2_e', 'EX_mqn7_e', 'EX_mqn8_e', 'EX_nac_e', 'EX_ncam_e', 'EX_nmn_e', 'EX_no2_e', 'EX_ocdca_e', 'EX_ocdcea_e', 'EX_orn_e', 'EX_phe__L_e', 'EX_pheme_e', 'EX_pi_e', 'EX_pnto__R_e', 'EX_pro__L_e', 'EX_ptrc_e', 'EX_pydx_e', 'EX_pydxn_e', 'EX_q8_e', 'EX_rib__D_e', 'EX_ribflv_e', 'EX_ser__L_e', 'EX_sheme_e', 'EX_so4_e', 'EX_spmd_e', 'EX_thm_e', 'EX_thr__L_e', 'EX_thymd_e', 'EX_trp__L_e', 'EX_ttdca_e', 'EX_tyr__L_e', 'EX_ura_e', 'EX_val__L_e', 'EX_xan_e', 'EX_xyl__D_e', 'EX_zn2_e']

for i in essentialMetabolites:
    if i not in raw_vmh_dict.keys():
        raw_vmh_dict[i] = -0.1


##########################################################################
### ADDING unmapped Compounds to the vmh diet dictionary #################

unmappedCompounds = ['EX_asn__L_e','EX_gln__L_e','EX_crn_e','EX_elaid_e','EX_hdcea_e','EX_dlnlcg_e','EX_adrn_e','EX_hco3_e','EX_sprm_e', 'EX_carn_e','EX_7thf_e','EX_Lcystin_e','EX_hista_e','EX_orn_e','EX_ptrc_e','EX_creat_e','EX_cytd_e','EX_so4_e'];

for i in unmappedCompounds:
    if i not in raw_vmh_dict.keys():
        raw_vmh_dict[i] = -50
        
raw_vmh_dict['EX_chol_e'] = -41.251

##########################################################################
### ADDING Micronutrients to the vmh diet dictionary #####################

micronutrients =['EX_adocbl_e','EX_vitd2_e','EX_vitd3_e','EX_psyl_e','EX_gum_e','EX_bglc_e','EX_phyQ_e','EX_fol_e','EX_5mthf_e','EX_q10_e','EX_retinol_9_cis_e','EX_pydxn_e','EX_pydam_e','EX_pydx_e','EX_pheme_e','EX_ribflv_e','EX_thm_e','EX_avite1_e','EX_pnto_R_e','EX_na1_e','EX_cl_e','EX_k_e','EX_pi_e','EX_zn2_e','EX_cu2_e']


for i in micronutrients:
    if i in raw_vmh_dict.keys():
        if raw_vmh_dict[i] > -0.1:
            raw_vmh_dict[i] = -0.1
    if 'EX_fol_e' in raw_vmh_dict.keys():
        if raw_vmh_dict['EX_fol_e'] > -1:
            raw_vmh_dict['EX_fol_e'] = -1



### Saving dictionary as json ############################################
with open('vmhToAGORA_diet.json','w') as ofile:
    json.dump(raw_vmh_dict, ofile)