import argparse
from deriva.core import ErmrestCatalog, get_credential, BaseCLI
import sys
import csv

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--RID', help="replicate RID", required=True)
    parser.add_argument('-t', '--table', help="source table", required=True)
    parser.add_argument('-o', '--host', help="datahub host", required=True)
    parser.add_argument('-c', '--cookie', help="cookie token", required=True)
    args = parser.parse_args()
    return args

def main(hostname, catalog_number, credential):
    catalog = ErmrestCatalog('https', hostname, catalog_number, credential)
    pb = catalog.getPathBuilder()
    if args.table == 'mRNA_QC':
        run_table = pb.RNASeq.mRNA_QC
    elif args.table == "Processed_File":
        run_table = pb.RNASeq.Processed_File

    path = run_table.filter(run_table.RID == args.RID)
    path.delete()
    rid = args.RID
    

    print(rid + "deleted")


if __name__ == '__main__':
    args = get_args()
    cli = BaseCLI("Custom RNASeq query", None, 1)
    cli.remove_options(["--config-file"])
    host = args.host
    credentials = {"cookie": args.cookie}
    main(host, 2, credentials)