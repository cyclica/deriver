from argparse import ArgumentParser
from playhouse.sqlite_ext import SqliteExtDatabase
from .config import logger
from .lib_read import lib_read


def clean(db_loc):
    """
    Some fragments are duplicates, or are not included in the db, but there are still atom and pseudoatom records
    that refer to them. This is not ideal, so this script deletes these orphan records.
    :param db_loc:
    :return:
    """
    # create the database for the output
    db = SqliteExtDatabase(db_loc, pragmas={
        'cache_size': -1024 * 64,  # 64MB page-cache.
        'journal_mode': 'wal',  # Use WAL-mode (you should always use this!).
        'foreign_keys': 0,
        'wal_autocheckpoint': 10,
    })
    # get the models
    Fragment, Heritage, PseudoAtoms, Atoms = lib_read(db)

    Fragment.create_table(safe=True)
    Heritage.create_table(safe=True)
    PseudoAtoms.create_table(safe=True)
    Atoms.create_table(safe=True)

    logger.debug("Trying to clean up the database:")
    logger.debug("Deleting missing ATOM records")
    with db.atomic():
        bad_atoms = Atoms.delete().where(
            (Atoms.frag.not_in(Fragment.select()))
        )
        bad_atoms.execute()

    logger.debug("Deleting missing PSEUDO_ATOM records")
    with db.atomic():
        bad_patoms = PseudoAtoms.delete().where(
            (PseudoAtoms.frag.not_in(Fragment.select()))
        )
        bad_patoms.execute()

    logger.debug("Deleting missing HERITAGE records")
    with db.atomic():
        bad_heritage = Heritage.delete().where(
            (Heritage.frag.not_in(Fragment.select()))
        )
        bad_heritage.execute()

    logger.info("Done.")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-i',
                        dest="input",
                        type=str,
                        help='.db file of molecules to be cleaned.',
                        required=True)

    options = parser.parse_args()
    clean(options.input)
