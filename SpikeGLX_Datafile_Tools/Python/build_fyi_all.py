# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 17:38:31 2024


Requires python 3

Build fyi_all files for pipeline output.
At the dialog, select the parent folder of the 'catgt_<run_name>' parent folders of the output.

Requires access to functions in DemoReadSGLX module to interpret metadata.

@author: colonellj
"""

import os
import shutil
import fnmatch
from pathlib import Path
from tkinter import Tk
from tkinter import filedialog
from DemoReadSGLXData import readSGLX


def main():
    
    # Get directory of pipeline output from user
    root = Tk()         # create the Tkinter widget
    root.withdraw()     # hide the Tkinter root window

    # Windows specific; forces the window to appear in front
    root.attributes("-topmost", True)

    pipeline_outDir = Path(filedialog.askdirectory(title="Select parent directory of output"))
    root.destroy()      # destroy the Tkinter widget
    
    fileList = os.listdir(pipeline_outDir)
    
    # loop over files, if a directory with name starting with catgt, process
    for nf, fName in enumerate(fileList):
        cp = os.path.join(pipeline_outDir,fName)
        if fName.startswith('catgt') & Path(cp).is_dir():
            # this is an output directory
            # parse for run name
            # check for all_fyi, make a copy
            run_name = fName[6:]
            all_fyi_path =  os.path.join(pipeline_outDir, fName, (run_name + '_all_fyi.txt'))
            if Path(all_fyi_path).is_file():
               # make a backup copy
               all_fyi_backup = os.path.join(pipeline_outDir, fName, (run_name + '_all_fyi_backup.txt'))
               shutil.copyfile(all_fyi_path,all_fyi_backup)
            
            fo = open(all_fyi_path,'w')
            cp_fileList = os.listdir(cp)
            # if nidq present, add lines for those
            nidq_metaName = run_name + '_tcat.nidq.meta'
            
            if cp_fileList.count(nidq_metaName):
                # read ni metadata, get identity of sync signal
                # loop over files in cp, find ni edge files
                nidq_binName = run_name + '_tcat.nidq.bin'
                ni_meta = readSGLX.readMeta(Path(os.path.join(cp,nidq_binName)))
                syncNiChanType = ni_meta['syncNiChanType']
                syncNiChan = ni_meta['syncNiChan']
                MN, MA, XA, DW = readSGLX.ChannelCountsNI(ni_meta)
                if syncNiChanType == '0':
                    # digital channel                
                    DW_index = MN+MA+XA
                    syncNIedgeName = run_name + '_tcat.nidq.xd_' + repr(DW_index) + '_' + syncNiChan + '_500.txt'                
                else:
                    # analog channel
                    xa_index = MN+MA+ int(syncNiChan)
                    syncNIedgeName = run_name + '_tcat.nidq.xa_' + syncNiChan + '_500.txt'
                
                # loop over filelist to find other ni times files
                ni_edge_pat = run_name + '_tcat.nidq.*.txt'
                ni_edge_list = fnmatch.filter(cp_fileList, ni_edge_pat);
                # write ni sync line to fyi
                ni_sync_path = os.path.join(cp,syncNIedgeName)
                ni_sync_line = 'sync_ni=' + ni_sync_path + '\n'
                fo.write(ni_sync_line)
                
                for edgeName in ni_edge_list:
                    if edgeName != syncNIedgeName:
                        #write a 'times' line
                        times_line = 'times_ni=' + os.path.join(cp,edgeName) + '\n'
                        fo.write(times_line)
            
            # find probe folders and create fyi entries for each
            probe_folder_pat = run_name + '_imec*'
            probe_folderList = fnmatch.filter(cp_fileList, probe_folder_pat)
            # these are common for all probes
            outpath_top_line = 'outpath_top=' + cp + '\n'
            fo.write(outpath_top_line)
            cp_parts = os.path.split(cp)
            sc_elem_line = 'supercat_element={' + cp_parts[0] + ',' + cp_parts[1] + '}' + '\n'
            fo.write(sc_elem_line)
            for prb in probe_folderList:
                # get probe string
                prb_str = prb[prb.index('_imec')+5 :]
                # find edge file
                prb_fileList = os.listdir(os.path.join(cp,prb))
                prb_edge_pat = run_name + '_tcat.imec' + prb_str + '.ap.xd_*_6_500.txt'
                prb_edgeName = fnmatch.filter(prb_fileList, prb_edge_pat)
                prb_path_line = 'outpath_probe' + prb_str + '=' + prb + '\n'
                fo.write(prb_path_line)
                prb_sync_line = 'sync_imec' + prb_str + '=' + os.path.join(prb,prb_edgeName[0]) + '\n'
                fo.write(prb_sync_line)
            fo.close()     
                   
    
    
if __name__ == "__main__":
    main()