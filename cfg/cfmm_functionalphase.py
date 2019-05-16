import os
import numpy
from cfmm_base import infotodict as cfmminfodict
from cfmm_base import create_key

def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """

    # call cfmm for general labelling and get dictionary
    info = cfmminfodict(seqinfo)

    task_ge = create_key('{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-{task}_part-mag_run-{item:02d}_bold')
    task_ge_phase = create_key('{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-{task}_part-phase_run-{item:02d}_bold')
    task_sbref = create_key('{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-{task}_run-{item:02d}_sbref')

    gre_diff = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_run-{item:02d}_phasediff')
    gre_magnitude = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_run-{item:02d}_magnitude')

    fmap_diff = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_phasediff')
    fmap_magnitude = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_magnitude')
    del info[fmap_diff]
    del info[fmap_magnitude]

    info[task_ge]=[]
    info[task_ge_phase]=[]
    info[task_sbref]=[]
    info[gre_diff]=[]
    info[gre_magnitude]=[]

    magacq=[]
    phaseacq=[]

    for idx, s in enumerate(seqinfo):

        if ('bold' in (s.series_description).strip()):
            if s.dim4==1 and 'SBRef' in (s.series_description).strip() and 'M' == s.image_type[2].strip():
                info[task_sbref].append({'item': s.series_id, 'task':'rest'})
            elif (s.dim4>1):
                if 'M' == s.image_type[2].strip():
                    acquisition_time=s.series_uid.split('.')[-4][8:14]
                    info[task_ge].append({'item': s.series_id, 'acqtime': acquisition_time, 'task': 'rest'})
                    magacq.append(float(acquisition_time))
                elif 'P' == s.image_type[2].strip():
                    acquisition_time=s.series_uid.split('.')[-4][8:14]
                    info[task_ge_phase].append({'item': s.series_id, 'acqtime': acquisition_time,  'task': 'rest'})
                    phaseacq.append(float(acquisition_time))

        if ('field_mapping' in s.protocol_name):
            if (s.dim4==1 and 'gre_field_mapping' in (s.series_description).strip()):
                if('P' in (s.image_type[2].strip()) ):
                    info[gre_diff].append({'item': s.series_id})
                if('M' in (s.image_type[2].strip()) ):
                    info[gre_magnitude].append({'item': s.series_id})

    magacq=[]
    phaseacq=[]

    # Now we have all the sequences in the right bins we need to link up the magnitudes and phases from the GE to ensure the runs match
    # find the unique magnitude runs and sort by start times
    magtimes = list(sorted(set(magacq)))

    # find the pairs of indicies corresponding to the same run
    magind =[None]*len(magtimes)
    phaseind =[None]*len(magtimes)
    for i in range(len(magtimes)):
        #find last phase file with that acq time (vunerable: assumes last phase recon is correct recon)
        magind[i]=magacq.index(magtimes[i])
        try:
            phaseind[i]=len(phaseacq)-1-(phaseacq[::-1]).index(magtimes[i])
        except:
            phaseind[i]=-1 # no phase run found do not append it
    # put the results sorted back into info
    maginfo = info[task_ge]
    phaseinfo = info[task_ge_phase]
    info[task_ge]=[]
    info[task_ge_phase]=[]

    for i in range(len(magind)):
        info[task_ge].append(maginfo[magind[i]])
        info[task_ge][i]['run']=str(i+1).zfill(2)
        if phaseind[i]>-1: # if no run found phase is not appended
            phaseinfo[phaseind[i]]['run']=str(i+1).zfill(2)
            info[task_ge_phase].append(phaseinfo[phaseind[i]])

    print(info)

    return info
