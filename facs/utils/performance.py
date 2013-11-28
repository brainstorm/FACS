"""Compare performance against fastq_screen.

Given previous test results stored in an external database (specified in your
~/.facsrc file), downloads this results and plot them.
"""
import logbook
import os
import sys
import json
import couchdb
import datetime

from facs.utils import config

log = logbook.Logger('FACS')

def _fetch_results(couch, database):
    """Fetch all documents a from database
    """
    log.info("Fetching all documents from database %s" % database)
    docs = []
    db = couch[database]
    [docs.append(db.get(doc)) for doc in db]
    log.info("Fetched %s documents" % str(len(docs)))
    return docs

def facs_vs_fastq_screen_files():
    """ Work from the json files on disk instead of fetched from DB
    """
    facs, fastq_screen = "", ""

    with open("facs.json") as fh:
        facs = json.load(fh)

    with open("fastq_screen.json") as fh:
        fastq_screen = json.load(fh)

    for fqscr in fastq_screen:
        if fqscr.get('sample'):
            for fcs in facs:
                if fcs.get('sample'):
                    # fqscreen stores the full path, facs does not
                    # XXX: We are Actually losing info here since full pathnames
                    # provide information about the system on which a particular
                    # test ran on.
                    if os.path.basename(fcs['sample']) == fqscr['sample']:
                        # Do not assume the test run went well
                        if len(fqscr['organisms']) > 0:
                            if fqscr.get('begin_timestamp'):
                            # Fetch timing info for each program
                                begin = fqscr['begin_timestamp']
                                end = fqscr['end_timestamp']

                                dt_b = datetime.datetime.strptime( begin, "%Y-%m-%d %H:%M:%S.%fZ" )
                                dt_e = datetime.datetime.strptime( end, "%Y-%m-%d %H:%M:%S.%fZ" )

                                delta = dt_e - dt_b

                            if fcs.get('begin_timestamp'):
                                begin_fcs = fcs['begin_timestamp']
                                end_fcs = fcs['end_timestamp']

                                # remove the messy UTC offset (+0200) (%z does not parse it out)
                                # http://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
                                begin_fcs = begin_fcs[:-5]
                                end_fcs = end_fcs[:-5]

                                dt_b_f = datetime.datetime.strptime( begin_fcs, "%Y-%m-%dT%H:%M:%S.%f" )
                                dt_e_f = datetime.datetime.strptime( end_fcs, "%Y-%m-%dT%H:%M:%S.%f" )

                                delta_fcs = dt_e_f - dt_b_f

                                print delta.total_seconds(), delta_fcs.total_seconds()



def facs_vs_fastq_screen():
    stream = logbook.StreamHandler(sys.stdout, level=logbook.INFO)
    with stream.applicationbound():
        log.info("Establishing connection with database %s" % config.SERVER)
        couch = couchdb.Server(config.SERVER)
        couch.resource.credentials = (config.USERNAME, config.PASSWORD)
        facs_results = _fetch_results(couch, config.FACS_DB)
        fastq_screen_results = _fetch_results(couch, config.FASTQ_SCREEN_DB)

        with open("facs.json", 'w') as fh:
            json.dump(facs_results, fh)

        with open("fastq_screen.json", 'w') as fh:
            json.dump(fastq_screen_results, fh)

if __name__ == "__main__":
    #facs_vs_fastq_screen()
    facs_vs_fastq_screen_files()

