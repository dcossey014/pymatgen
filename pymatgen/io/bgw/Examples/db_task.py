#!/usr/bin/env python

from fireworks import Firework
from pymatgen.io.bgw.custodian_jobs import BgwDB

db_fw = Firework(BgwDB(config_file='bgw_db.yaml', insert_to_db=True), name="BGW DB Task")
db_fw.spec['PREV_DIRS'] = {'ESPRESSO' : 
        {'scf' : '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-30-21-08-51-280825/ESPRESSO/scf',
        'wfn' : '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-30-21-08-51-280825/ESPRESSO/wfn',
        'wfnq': '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-30-21-08-51-280825/ESPRESSO/wfnq',
        'wfn_co': '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-30-21-08-51-280825/ESPRESSO/wfn_co',
        'wfn_fi': '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-30-21-08-51-280825/ESPRESSO/wfn_fi',
        'wfnq_fi': '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-30-21-08-51-280825/ESPRESSO/wfnq_fi'
        },
    'BGW' : 
        {'epsilon' : '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-30-22-09-23-943604',
        'sigma' : '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-30-23-40-51-634090',
        'kernel': '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-31-16-22-31-084675',
        'absorption' : '/p/work1/workspace/dec014/BGW/GARY/launcher_2018-01-31-17-43-26-019731'
        } 
    } 

db_fw.to_file('db_task.yaml')
