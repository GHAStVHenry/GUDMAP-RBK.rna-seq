import argparse
from deriva.core import ErmrestCatalog, get_credential, BaseCLI
import sys
import csv
from datetime import datetime

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--executionRunRID', help="exection run RID", required=True)
    parser.add_argument('-f', '--file', help="file name", required=True)
    parser.add_argument('-l', '--loc', help="datahub location", required=True)
    parser.add_argument('-s', '--md5', help="md5 sum", required=True)
    parser.add_argument('-b', '--bytes', help="size in bytes", required=True)
    parser.add_argument('-n', '--notes', help="notes", default="", required=False)
    parser.add_argument('-o', '--host', help="datahub host", required=True)
    parser.add_argument('-c', '--cookie', help="cookie token", required=True)
    args = parser.parse_args()
    return args

def main(hostname, catalog_number, credential):
    catalog = ErmrestCatalog('https', hostname, catalog_number, credential)
    pb = catalog.getPathBuilder()
    outputBag_table = pb.RNASeq.Output_Bag

    outputBag_data = {
        "Execution_Run": args.executionRunRID,
        "File_Name": args.file,
        "File_URL": args.loc,
        "File_MD5": args.md5,
        "File_Bytes": args.bytes,
        "File_Creation_Time": datetime.now().replace(microsecond=0).isoformat(),
        "Notes": args.notes,
        "Bag_Type": "Replicate_mRNA_Analysis"
        }

    entities = outputBag_table.insert([outputBag_data])
    rid = entities[0]["RID"]

    print(rid)


if __name__ == '__main__':
    args = get_args()
    cli = BaseCLI("Custom RNASeq query", None, 1)
    cli.remove_options(["--config-file"])
    host = args.host
    credential = {"cookie": args.cookie}
    main(host, 2, credential)