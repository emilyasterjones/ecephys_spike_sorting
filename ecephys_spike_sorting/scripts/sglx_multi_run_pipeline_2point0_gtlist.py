import os
import sys
import subprocess
import glob
import pandas as pd
import numpy as np
from datetime import datetime

from helpers import SpikeGLX_utils
from helpers import log_from_json
from helpers import run_one_probe
from create_input_json import createInputJson

# To use this script, you must modify catGT_helper/main.py by
# commenting out lines 60 & 61


# script to run CatGT, KS2, postprocessing and TPrime on data collected using
# SpikeGLX. The construction of the paths assumes data was saved with
# "Folder per probe" selected (probes stored in separate folders) AND
# that CatGT is run with the -out_prb_fld option

# -------------------------------
# -------------------------------
# User input -- Edit this section
# -------------------------------
# -------------------------------

# set for each run
base_dir = '//oak-smb-giocomo.stanford.edu/groups/giocomo/emijones/WT_Sequences/2024_winter/'
###base_dir = 'D:/WT_Sequences/2024_winter/'
prefix = 'all' #eg datetime.today().strftime('%Y%m%d')

# should be the same across all experiments
rec_file_list = '//oak-smb-giocomo.stanford.edu/groups/giocomo/emijones/WT_Sequences/all_sessions_Xmaze_OF.csv'

log_file = os.path.join(base_dir, 'Preprocessed_Data/Provenance',prefix+'_ecephys_log.csv')
log_stream = open(log_file, 'a+')
log_stream.write('Index,File,Sort_Error,Sort_Error_Description\n')
sessions = pd.read_csv(rec_file_list)

# set which parts of script are running
run_CatGT = True
# List of modules to run per probe; CatGT and TPrime are called once for each run.
modules = [
            'kilosort_helper',
            'kilosort_postprocessing',
            'noise_templates',    
            'mean_waveforms',
            'quality_metrics',
            ###'depth_estimation',
            'prePhy_filters'
			]
run_TPrime = True

# brain region specific params
# can add a new brain region by adding the key and value for each param
# can add new parameters -- any that are taken by create_input_json --
# by adding a new dictionary with entries for each region and setting the 
# according to the new dictionary in the loop to that created json files.
refPerMS_dict = {'default': 2.0, 'cortex': 2.0, 'medulla': 1.5, 'thalamus': 1.0}

# threhold values appropriate for KS2, KS2.5
#ksTh_dict = {'default':'[10,4]', 'cortex':'[10,4]', 'medulla':'[10,4]', 'thalamus':'[10,4]'}
# threshold values appropriate for KS3.0
ksTh_dict = {'default':'[9,9]', 'cortex':'[9,9]', 'medulla':'[9,9]', 'thalamus':'[9,9]'}

# ------------
# CatGT params
# ------------
# CAR mode for CatGT. Must be equal to 'None', 'gblcar' or 'loccar'
car_mode = 'gblcar'
# inner and outer radii, in um for local comman average reference, if used
loccar_min = 40
loccar_max = 160

# CatGT commands for bandpass filtering, artifact correction, and zero filling
# Note 1: directory naming in this script requires -prb_fld and -out_prb_fld
# Note 2: this command line includes specification of edge extraction
# see CatGT readme for details
# these parameters will be used for all runs

# notes from Mari's repo:
# gfix=0,0.10,0.02 -- artifact removal; params: |thresh_amp(mV)|,|slope(mV/sample)|,noise
# -t_miss_ok option required to concatenate over missing g or t indices
# -zerofillmax=500 option required to fill gaps only up to 500ms of zeros,
# so kilsort doesn't crash
catGT_cmd_string = '-t_miss_ok -zerofillmax=500 -prb_fld -out_prb_fld -gfix=0.4,0.10,0.02 -apfilter=butter,12,300,10000 -lffilter=butter,12,1,500 -gtlist={1,0,0}{3,0,0}'

ni_present = True

# ----- NIDAQ INPUTS -----
# -- Each XA gets its own word, starting with 0. XA inputs must come first.
# -- XD inputs come next, on separate words from the XA inputs. 
# -- Each XD word contains up to 16 bits (0:15)
# Dec 2023 update: new catGT requires 2 leading zeros + xa must be lowercase
# -- Yggdrasil inputs: --
# XA=0,1,3,500 -- sync channel on nidaq: word 0, thresh 1 V, must stay above 3V, dur 500 ms
# XA=1,1,3,0 -- Arduino: word 1, thresh 1V, must stay above 3V, dur ??? ms
# XA=2,1,3,0 -- camera: word 2, thresh 1V, must stay above 3V, dur ??? ms
# - duration must be within +/-20% of the actual pulse width; 0 ignores pulse width requirement
ni_extract_string = '-xa=0,0,0,1,3,500 -xa=0,0,1,1,3,0 -xa=0,0,2,1,3,0'


# ----------------------
# KS2 or KS25 parameters
# ----------------------
# parameters that will be constant for all recordings
# Template ekmplate radius and whitening, which are specified in um, will be 
# translated into sites using the probe geometry.
ks_remDup = 0
ks_saveRez = 1
ks_copy_fproc = 0
ks_templateRadius_um = 163
ks_whiteningRadius_um = 163
ks_minfr_goodchannels = 0.05 #0.1
ks_CAR = 0          # CAR already done in catGT
ks_nblocks = 1      # for KS2.5 and KS3; 1 for rigid registration in drift correction, 
                    # higher numbers to allow different drift for different 'blocks' of the probe

# If running KS20_for_preprocessed_data:
# (https://github.com/jenniferColonell/KS20_for_preprocessed_data)
# can skip filtering with the doFilter parameter.
# Useful for speed when data has been filtered with CatGT.
# This parameter is not implemented in standard versions of kilosort.
ks_doFilter = 0

ks_output_tag = 'ks'

# ----------------------
# C_Waves snr radius, um
# ----------------------
c_Waves_snr_um = 160

# ----------------------
# psth_events parameters
# ----------------------
# extract param string for psth events -- copy the CatGT params used to extract
# events that should be exported with the phy output for PSTH plots
# If not using, set to empty string
event_ex_param_str = ''

# -----------------
# TPrime parameters
# -----------------
sync_period = 1.0   # true for SYNC wave generated by imec basestation
toStream_sync_params = 'ni'


for a, row in sessions.iloc[np.concatenate((np.arange(82,90), np.arange(100,116)))].iterrows():
    if (not row['Position_Error']) and (not row['Recording_Error']):
        print(row['File'])
        rec_file_stem = row['File'].split('/')[1]
        rec_file_stem = rec_file_stem[:-3]
        animal = row['Animal']
        log_string = f'{a},{rec_file_stem},'
        
        npx_directory = os.path.join(row['Base_Directory'],'Raw_Data/Neural_Traces',animal)
        # all output will be written here.
        # Output will be in the standard SpikeGLX directory structure:
        # run_folder/probe_folder/*.bin
        catGT_dest = os.path.join(row['Base_Directory'], 'Preprocessed_Data/Spikes/g1_3', animal, 'Ecephys', rec_file_stem)
        json_directory = catGT_dest
        
        if not os.path.exists(catGT_dest):
            os.mkdir(catGT_dest)

        # -----------
        # Input data
        # -----------
        # Name for log file for this pipeline run. Log file will be saved in the
        # output destination directory catGT_dest
        # If this file exists, new run data is appended to it
        logName = f'{rec_file_stem}_log.csv'
        
        # run_specs = name, gate, trigger and probes to process
        # Each run_spec is a list of 4 strings:
        #   undecorated run name (no g/t specifier, the run field in CatGT)
        #   gate index range, as a string (e.g. '0,9'); can be greater than actual range
        #       later code will look for all g indices in this range with the same filename stem
        #       EAJ edit: can't currently handle gate indices with >1 digit or multiple triggers
        #   triggers to process/concatenate, as a string e.g. '0,400', '0,0 for a single file
        #           can replace first limit with 'start', last with 'end'; 'start,end'
        #           will concatenate all trials in the probe folder
        #   probes to process, as a string, e.g. '0', '0,3', '0:3'
        #   brain regions, list of strings, one per probe, to set region specific params
        #           these strings must match a key in the param dictionaries above.
        run_specs = [				
                                # [rec_file_stem, 'gate', 'triggers', 'probes', 'regions']
                                [rec_file_stem, '1,3', '0,0', '0', ['cortex']]
                                # [rec_file_stem, '0,0', '0,0', '0', ['cortex']]
        ]
        # nshanks = 4
        # ks_trange_list = ['[0 600]', '[603 1203]', '[1206 1806]', '[1809 Inf]']
        # chanMap_list = []
        nshanks = 1
        ks_trange_list = ['[0 Inf]']
        # chanMap_list = ['_HPC', '_MEC']
        chanMap_list = ['_MEC']
        
        # check for existence of log file, create if not there
        logFullPath = os.path.join(catGT_dest, logName)
        if not os.path.isfile(logFullPath):
            # create the log file, write header
            log_from_json.writeHeader(logFullPath)
            

        for spec in run_specs:
        
            # Make list of probes from the probe string
            prb_list = SpikeGLX_utils.ParseProbeStr(spec[3])
            
            # build path to the first probe folder; look into that folder
            # to determine the range of trials if the user specified t limits as
            # start and end
            [first_gate, last_gate] = SpikeGLX_utils.ParseGateStr(spec[1])
            run_folder_name = spec[0] + '_g' + repr(first_gate)
            prb0_fld_name = run_folder_name + '_imec' + prb_list[0]
            prb0_fld = os.path.join(npx_directory, run_folder_name, prb0_fld_name)
            first_trig, last_trig = SpikeGLX_utils.ParseTrigStr(spec[2], prb_list[0], first_gate, prb0_fld)
            trigger_str = repr(first_trig) + ',' + repr(last_trig)
            
            # from MS fork
            # get list of g-indices to concatenate from data directory
            g_range = '[' + spec[1][0] + '-' + spec[1][-1] + ']'
            g_tocat = sorted(glob.glob(os.path.join(npx_directory,(rec_file_stem + '_g' + g_range))))
            glist = ''.join((x[-1]+'-') for x in g_tocat)[:-1] # g inds separated by dashes, minus the last dash        
        
            print('Concatenating g indices ' + glist)
            
            # loop over all probes to build json files of input parameters
            # initalize lists for input and output json files
            catGT_input_json = []
            catGT_output_json = []
            module_input_json = []
            module_output_json = []
            session_id = []
            catgt_output_dir = []
            data_directory = []
            
            # first loop over probes creates json files containing parameters for
            # both preprocessing (CatGt) and sorting + postprocessing
            for i, prb in enumerate(prb_list):
                    
                #create CatGT command for this probe
                print('Creating json file for CatGT on probe: ' + prb)
                
                # Run CatGT
                catGT_input_json.append(os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_CatGT' + '-input.json'))
                catGT_output_json.append(os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_CatGT' + '-output.json'))
                
                # if this is the first probe proceessed, process the ni stream with it
                if i == 0 and ni_present:
                    catGT_stream_string = '-ni -ap -lf'
                    extract_string = ni_extract_string
                else:
                    catGT_stream_string = '-ap -lf'
                    extract_string = ''

                # build name of first trial to be concatenated/processed;
                # allows reading of the metadata
                run_str = spec[0] + '_g' + str(first_gate)
                run_folder = run_str
                prb_folder = run_str + '_imec' + prb
                input_data_directory = os.path.join(npx_directory, run_folder, prb_folder)
                fileName = run_str + '_t' + repr(first_trig) + '.imec' + prb + '.ap.bin'
                continuous_file = os.path.join(input_data_directory, fileName)
                metaName = run_str + '_t' + repr(first_trig) + '.imec' + prb + '.ap.meta'
                input_meta_fullpath = os.path.join(input_data_directory, metaName)
                
                print(input_meta_fullpath)
                 
                info = createInputJson(catGT_input_json[i], npx_directory=npx_directory, 
                                                continuous_file = continuous_file,
                                                kilosort_output_directory=catGT_dest,
                                                input_meta_path = input_meta_fullpath,
                                                catGT_run_name = spec[0],
                                                gate_string = spec[1],
                                                trigger_string = trigger_str,
                                                probe_string = prb,
                                                catGT_stream_string = catGT_stream_string,
                                                catGT_car_mode = car_mode,
                                                catGT_loccar_min_um = loccar_min,
                                                catGT_loccar_max_um = loccar_max,
                                                catGT_cmd_string = catGT_cmd_string + ' ' + extract_string,
                                                extracted_data_directory = catGT_dest
                                                )      

                #create json files for the other modules
                session_id.append(spec[0] + '_imec' + prb)
                module_input_json.append(os.path.join(json_directory, session_id[i] + '-input.json'))
                
                # location of the binary created by CatGT, using -out_prb_fld
                run_folder = 'catgt_' + run_str
                prb_folder = run_str + '_imec' + prb
                catgt_output_dir = os.path.join(catGT_dest, run_folder)
                data_directory.append(os.path.join(catGT_dest, run_folder, prb_folder))
                fileName = run_str + '_tcat.imec' + prb + '.ap.bin'
                continuous_file = os.path.join(data_directory[i], fileName)
        
                # all post-processing modules alter KS output for Phy, so store the original
                ks_make_copy = True
                
                imro_dir = os.path.join(row['Base_Directory'],'Raw_Data/Neural_Traces/IMRO')
                
                for shank in range(nshanks):
                    if nshanks==1:
                        outputName = 'imec' + prb + '_ks'
                        chanMap = os.path.join(imro_dir, animal+chanMap_list[i]+'.mat')
                    else:
                        outputName = 'imec' + prb + 'shank' + str(shank) + '_ks'
                        chanMap = os.path.join(imro_dir, 'shank'+str(shank)+'bank0.mat')
                    ks_trange = ks_trange_list[shank]
            
                    kilosort_output_dir = os.path.join(data_directory[i], outputName)
                    
                    # get region specific parameters
                    ks_Th = ksTh_dict.get(spec[4][i])
                    refPerMS = refPerMS_dict.get(spec[4][i])
                    print( 'ks_Th: ' + repr(ks_Th) + ' ,refPerMS: ' + repr(refPerMS))
            
                    info = createInputJson(module_input_json[i], npx_directory=npx_directory, 
                                                   continuous_file = continuous_file,
                                                   input_meta_path = input_meta_fullpath,
                                                   kilosort_output_directory=kilosort_output_dir,
                                                   ks_make_copy = ks_make_copy,
                                                   noise_template_use_rf = False,
                                                   catGT_run_name = session_id[i],
                                                   gate_string = spec[1],
                                                   probe_string = spec[3], 
                                                   ks_remDup = ks_remDup,                   
                                                   ks_finalSplits = 1,
                                                   ks_labelGood = 1,
                                                   ks_saveRez = ks_saveRez,
                                                   ks_copy_fproc = ks_copy_fproc,
                                                   ks_minfr_goodchannels = ks_minfr_goodchannels,                  
                                                   ks_whiteningRadius_um = ks_whiteningRadius_um,
                                                   ks_doFilter = ks_doFilter,
                                                   ks_Th = ks_Th,
                                                   ks_CSBseed = 1,
                                                   ks_LTseed = 1,
                                                   ks_templateRadius_um = ks_templateRadius_um,
                                                   ks_nblocks = ks_nblocks,
                                                   ks_CAR = ks_CAR,
                                                   extracted_data_directory = data_directory[i],
                                                   event_ex_param_str = event_ex_param_str,
                                                   c_Waves_snr_um = c_Waves_snr_um,                
                                                   qm_isi_thresh = refPerMS/1000,
                                                   ks_trange = ks_trange,
                                                   ks_chanMap = chanMap,
                                                   halfwidth_max = 0.35,
                                                   )   
            
                    # Run each module --- CatGT and KS are run here ---
                    run_one_probe.runOne( session_id[i],
                         json_directory,
                         data_directory[i],
                         run_CatGT,
                         catGT_input_json[i],
                         catGT_output_json[i],
                         modules,
                         module_input_json[i],
                         logFullPath )
                             
                    
            if run_TPrime:
                # after loop over probes, run TPrime to create files of 
                # event times -- edges detected in auxialliary files and spike times 
                # from each probe -- all aligned to a reference stream.
            
                # create json files for calling TPrime
                input_json = os.path.join(json_directory, spec[0] + '_g' + glist + '_TPrime' + '-input.json') 
                output_json = os.path.join(json_directory, spec[0] + '_g' + glist + '_TPrime' + '-input.json')
                               
                info = createInputJson(input_json, npx_directory=npx_directory, 
                                       continuous_file = continuous_file,
                                       input_meta_path = input_meta_fullpath,
                                       catGT_run_name = spec[0],
                                       kilosort_output_directory=kilosort_output_dir,
                                       extracted_data_directory = catGT_dest,
                                       tPrime_ni_ex_list = ni_extract_string,
                                       event_ex_param_str = event_ex_param_str,
                                       sync_period = 1.0,
                                       toStream_sync_params = toStream_sync_params,
                                       gate_string = spec[1],
                                       ks_output_tag = ks_output_tag
                                                   ) 
                
                command = sys.executable + " -W ignore -m ecephys_spike_sorting.modules." + 'tPrime_helper' + " --input_json " + input_json \
                              + " --output_json " + output_json
                subprocess.check_call(command.split(' '))  

        #    log_string+='False,'
        #except Exception as e:
        #    log_string+=f'True,{e}'
            
        log_stream.write(log_string+'\n')
        #except:
        #    print("Error running current file")
        
log_stream.close()
