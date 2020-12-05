import argparse
from deriva.core import ErmrestCatalog, get_credential, BaseCLI
import sys
import csv

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--repRID', help="replicate RID", required=True)
    parser.add_argument('-e', '--executionRunRID', help="exection run RID", required=True)
    parser.add_argument('-p', '--ends', help="single/paired ends", required=True)
    parser.add_argument('-s', '--stranded', help="stranded?", required=True)
    parser.add_argument('-l', '--length', help="median read length", required=True)
    parser.add_argument('-w', '--rawCount', help="raw count", required=True)
    parser.add_argument('-f', '--assignedCount', help="final assigned count", required=True)
    parser.add_argument('-n', '--notes', help="notes", default="", required=False)
    parser.add_argument('-o', '--host', help="datahub host", required=True)
    parser.add_argument('-c', '--cookie', help="cookie token", required=True)
    parser.add_argument('-u', '--update', help="update?", default="F", required=True)
    args = parser.parse_args()
    return args

def main(hostname, catalog_number, credential):
    catalog = ErmrestCatalog('https', hostname, catalog_number, credential)
    pb = catalog.getPathBuilder()
    run_table = pb.RNASeq.mRNA_QC

    if args.update == "F":
        run_data = {
            "Execution_Run": args.executionRunRID,
            "Replicate": args.repRID,
            "Paired_End": args.ends,
            "Strandedness": args.stranded,
            "Median_Read_Length": args.length,
            "Raw_Count": args.rawCount,
            "Final_Count": args.assignedCount,
            "Notes": args.notes
        }
        entities = run_table.insert([run_data])
        rid = entities[0]["RID"]
    else:
        run_data = {
            "RID": args.update,
            "Execution_Run": args.executionRunRID,
            "Replicate": args.repRID,
            "Paired_End": args.ends,
            "Strandedness": args.stranded,
            "Median_Read_Length": args.length,
            "Raw_Count": args.rawCount,
            "Final_Count": args.assignedCount,
            "Notes": args.notes
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