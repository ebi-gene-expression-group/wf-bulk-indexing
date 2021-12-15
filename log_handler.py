import sys
import logging
import os

l = logging.getLogger(__name__)
l.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.environ.get('LOG_PATH'), mode='w')
l.addHandler(fh)

jobs = {}

def log_handler(msg):

    if "jobid" in msg:
        if msg['jobid'] in jobs:
            if 'level' in msg:
                jobs[msg['jobid']]['status'] = msg['level']
        else:
            job = { 'name': msg['name'], 'status': 'created' }
            if 'accession' in msg['wildcards']:
                job['accession'] = msg['wildcards']['accession']
            if 'log' in msg and len(msg['log']) > 0:
                job['log'] = msg['log'][0]
            jobs[msg['jobid']] = job

        #for k, v in msg.items():
        #    print(f"Key: {k}")
        #    print(f"Value: {v}")

        jobs_to_remove = []
        for j_id, j_cont in jobs.items():
            if 'accession' in j_cont and j_cont['status'] != 'created':
                log=""
                if 'log' in j_cont:
                    log=j_cont['log']
                l.info(f"{j_cont['accession']}\t{j_cont['name']}\t{j_cont['status']}\t{log}")
                jobs_to_remove.append(j_id)
                #sys.stdout.flush()

        for j_id in jobs_to_remove:
            del jobs[j_id]

    #if msg.get('level') == 'progress':
    #    cur = msg['done']
    #    total = msg['total']
    #    print(get_bar(cur, total) + f" {cur}/{total}")
    #if 'msg' in msg and msg['msg'] == 'removed all locks':
    #    sys.stdout.close()
    