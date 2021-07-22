#!/usr/bin/env python3
#upload_execution_run.py
#*
#* --------------------------------------------------------------------------
#* Licensed under MIT (https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/blob/14a1c222e53f59391d96a2a2e1fd4995474c0d15/LICENSE)
#* --------------------------------------------------------------------------
#*

import argparse
from deriva.core import ErmrestCatalog, get_credential, BaseCLI
import sys
import csv

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repRID', help="replicate RID", required=True)
    parser.add_argument('-w', '--workflowRID', help="workflow RID", required=True)
    parser.add_argument('-g', '--referenceRID', help="reference genome RID", required=True)
    parser.add_argument('-i', '--inputBagRID', help="inputBag RID", required=True)
    parser.add_argument('-n', '--notes', help="notes", default="", required=False)
    parser.add_argument('-s', '--status', help="run status", default="", required=False)
    parser.add_argument('-d', '--statusDetail', help="status detail", default="", required=False)
    parser.add_argument('-o', '--host', help="datahub host", required=True)
    parser.add_argument('-c', '--cookie', help="cookie token", required=True)
    parser.add_argument('-u', '--update', help="update?", default="F", required=True)
    args = parser.parse_args()
    return args

def main(hostname, catalog_number, credential):
    catalog = ErmrestCatalog('https', hostname, catalog_number, credential)
    pb = catalog.getPathBuilder()
    run_table = pb.RNASeq.Execution_Run

    if args.update == "F":
        run_data = {
            "Replicate": args.repRID,
            "Workflow": args.workflowRID,
            "Reference_Genome": args.referenceRID,
            "Input_Bag": args.inputBagRID,
            "Notes": args.notes,
            "Execution_Status": args.status,
            "Execution_Status_Detail": args.statusDetail.replace('\\n','\n')
        }
        entities = run_table.insert([run_data], add_system_defaults=False, defaults={"RID", "RCT", "RMT"})
        rid = entities[0]["RID"]
    else:
        run_data = {
            "RID": args.update,
            "Replicate": args.repRID,
            "Workflow": args.workflowRID,
            "Reference_Genome": args.referenceRID,
            "Input_Bag": args.inputBagRID,
            "Notes": args.notes,
            "Execution_Status": args.status,
            "Execution_Status_Detail": args.statusDetail.replace('\\n','\n')
        }
        entities = run_table.update([run_data])
        rid = args.update

    print(rid)


if __name__ == '__main__':
    args = get_args()
    cli = BaseCLI("Custom RNASeq query", None, 1)
    cli.remove_options(["--config-file"])
    host = args.host
    credentials = {"cookie": args.cookie}
    main(host, 2, credentials)
