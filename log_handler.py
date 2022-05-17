import sys
import logging
import os
import re

import pandas as pd

# print full data frame
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", -1)


def format_percentage(done, total):
    """Format percentage from given fraction while avoiding superflous precision."""
    if done == total:
        return "100%"
    if done == 0:
        return "0%"
    precision = 0
    fraction = done / total
    fmt_precision = "{{:.{}%}}".format
    fmt = lambda fraction: fmt_precision(precision).format(fraction)
    while fmt(fraction) == "100%" or fmt(fraction) == "0%":
        precision += 1
    return fmt(fraction)


l = logging.getLogger(__name__)
l.setLevel(logging.DEBUG)
if os.environ.get('LOG_PATH'):
    fh = logging.FileHandler(os.environ.get('LOG_PATH'), mode='w')
    l.addHandler(fh)

reit_number = int(os.environ.get("RESTART_TIMES"))

# jobs = {}
error_list = []
job_error_list = []
msg_error_list = []
chunks = []
n = 0
# logging of errors in the workflow


def log_handler(msg):
    # for k, v in msg.items():
    #    l.info(f"Key: {k}")
    #    l.info(f"Value: {v}")

    global n  # counter for errors in rules
    global df  # data frame with error summary
    # l.info(f"new msg")

    # ======================== RUNTIME LOG =========================
    if "level" in msg:
        if msg["level"] == "error":
            # l.info(f"---error found")
            error_list.append(msg["level"])
        if msg["level"] == "job_error":
            # l.info(f"---job_error found")
            job_error_list.append(msg["level"])
        if msg["level"] == "progress":
            done = msg["done"]
            total = msg["total"]
            l.info(
                "{} of {} steps ({}) done".format(
                    done, total, format_percentage(done, total)
                )
            )

    if "wildcards" in msg:
        if msg["wildcards"]:
            chunks.append(msg["wildcards"]["chunk"])
            # l.info(f"---Wildcards found in log msg:\t{chunks}")

    if "msg" in msg:
        A = msg["msg"]
        if A != None and isinstance(A, str):
            if re.search("Error executing rule", A, re.IGNORECASE):
                # lets not print all these on runtime and just print a final summary
                # l.info(f"---\t{A}")
                msg_error_list.append(A)
                # from A, create and fill array with Rule, Chunk, Error_occurrence, Error_out
                n += 1
                rule = re.split(r"on cluster", re.split(r"Error executing", A)[1])[0]
                if re.search("chunk", A, re.IGNORECASE):
                    # for rules than split in chunks
                    chunk = re.split(r"/", re.split(r"chunk=", A)[1])[0]
                else:
                    # for unique rules
                    chunk = "unique"
                error_occurrence = 1
                error_out = re.split(r",", re.split(r" ", A)[10])[0]

                if n == 1:
                    df = pd.DataFrame(
                        [[rule, chunk, error_occurrence, error_out]],
                        columns=["Rule", "Chunk", "Error_occurrence", "Error_out"],
                    )
                    # l.info(f"{df} ")
                else:
                    i = df[(df["Rule"] == rule) & (df["Chunk"] == chunk)].index
                    # l.info(f"{i} ")
                    if len(i) == 1:
                        # update row
                        df["Error_occurrence"].iloc[i] += 1
                        df["Error_out"].iloc[i] = error_out
                    else:
                        # add new row
                        temp_df = pd.DataFrame(
                            [[rule, chunk, error_occurrence, error_out]],
                            columns=["Rule", "Chunk", "Error_occurrence", "Error_out"],
                        )
                        df = pd.concat([df, temp_df], ignore_index=True)
                        del temp_df
                    # l.info(f"{df} ")
                del rule
                del chunk
                del error_occurrence
                del error_out
            # end of the workflow.
            # ======================== LOG SUMMARY =========================
            elif re.search("removed all locks", A, re.IGNORECASE):

                l.info(f"---------------------")
                l.info(f"Log handler - Summary:")
                # l.info(f"{error_list}" )
                # l.info(f"{job_error_list}" )
                # l.info(f"{msg_error_list}" )

                l.info(f"length of error_list = {len(error_list)}")
                l.info(f"length of job_error_list = {len(job_error_list)}")
                l.info(f"length of msg_error_list = {len(msg_error_list)}")

                l.info(f"Some errors shown above may have been solved through retrials")
                l.info(f"---------------------")
                l.info(
                    f"Detect and print content of files with failed or missing accessions:"
                )
                wdir = os.getcwd()
                folders = [
                    "analytics_jsonl_files",
                    "load_bulk_analytics_index",
                    "update_coexpressions",
                    "update_experiment_designs",
                ]
                filetype = ["missing", "failed"]
                if set(chunks):
                    for folder in folders:
                        for ch in set(chunks):
                            for fltp in filetype:
                                acc_info_file = (
                                    wdir
                                    + "/"
                                    + folder
                                    + "/"
                                    + ch
                                    + "/"
                                    + fltp
                                    + "_accessions.txt"
                                )
                                if os.path.isfile(acc_info_file):
                                    l.info(f"File found:  \t{acc_info_file}")
                                    f = open(acc_info_file, "r")
                                    file_contents = f.read()
                                    l.info(f"\t{file_contents}")
                                    f.close()
                if "df" in globals():
                    l.info(
                        f"There have been error/s. If Error_occurrence is == RESTART_TIMES+1 then error should be investigated"
                    )
                    l.info(f"RESTART_TIMES={reit_number}")
                    if len(df.loc[df["Error_occurrence"] >= reit_number + 1].index) > 0:
                        l.info(f"The filtered data frame is:")
                        l.info(f"{df.loc[df['Error_occurrence'] >= reit_number + 1 ]}")
                    else:
                        l.info(f"Completed. No errors after filtering.")
                else:
                    l.info(f"Completed. No errors.")
