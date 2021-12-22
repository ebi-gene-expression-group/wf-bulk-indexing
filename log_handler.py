import sys
import logging
import os
import re

l = logging.getLogger(__name__)
l.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.environ.get('LOG_PATH'), mode='w')
l.addHandler(fh)

#jobs = {}
error_list=[]
job_error_list=[]
msg_error_list=[]
chunks=[]
# logging of errors in the workflow

def log_handler(msg):
    #for k, v in msg.items():
    #    l.info(f"Key: {k}")
    #    l.info(f"Value: {v}")

    #l.info(f"new msg")
    if "level" in msg:
        if msg['level'] == 'error':
            l.info(f"---error found")
            error_list.append(msg['level'])
        if msg['level'] == 'job_error':
            l.info(f"---job_error found")
            job_error_list.append(msg['level'])

    if 'wildcards' in msg:
        if msg['wildcards']:
            chunks.append(msg['wildcards']['chunk'])
            l.info(f"---Wildcards found in msg:\t{chunks}")


    if 'msg' in msg:
        A=msg['msg']
        if A != None and isinstance(A, str):
            if re.search("Error executing rule", A, re.IGNORECASE) :
                l.info(f"---\t{A}")
                msg_error_list.append(A)
            # end of the workflow:
            elif re.search("removed all locks", A, re.IGNORECASE) :

                l.info(f"---------------------")
                l.info(f"Log handler- Summary:")
                l.info(f"{error_list}" )
                l.info(f"{job_error_list}" )
                l.info(f"{msg_error_list}" )

                l.info(f"length of error_list = {len(error_list)}" )
                l.info(f"length of job_error_list = {len(job_error_list)}" )
                l.info(f"length of msg_error_list = {len(msg_error_list)}" )

                l.info(f"Some errors shown above may have been solved through retrials")
                l.info(f"---------------------")
                l.info(f"Detect and print content of files with failed or missing accessions:")
                wdir = os.getcwd()
                folders=["analytics_jsonl_files", "load_bulk_analytics_index", "update_coexpressions", "update_experiment_designs"]
                filetype=["missing", "failed"]
                if set(chunks):
                    for folder in folders:
                        for ch in set(chunks):
                            for fltp in filetype:
                                acc_info_file=wdir+"/"+folder+"/"+ch+"/"+fltp+"_accessions.txt"
                                if os.path.isfile(acc_info_file):
                                    l.info(f"File found:  \t{acc_info_file}")
                                    f = open(acc_info_file, 'r')
                                    file_contents = f.read()
                                    l.info(f"\t{file_contents}")
                                    f.close()

