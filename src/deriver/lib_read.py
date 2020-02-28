from peewee import IntegerField, CharField, Model, ForeignKeyField, CompositeKey


def lib_read(lib_database):
    """
    format/models for the fragment databases
    :param lib_database:
    :return:
    """

    # a model for the fragments themselves
    class Fragment(Model):
        id = IntegerField(primary_key=True)
        frag_coeff = IntegerField()
        smile = CharField(unique=True)
        num_pseudo_atoms = IntegerField(default=0)
        num_unique_pseudo_atoms = IntegerField(default=0)

        class Meta:
            table_name = 'fragments'
            database = lib_database

    # a model for an entry in the heritage table, which contains information
    # about which fragments came from which molecules
    class Heritage(Model):
        frag = ForeignKeyField(
            Fragment, on_delete="CASCADE", on_update="CASCADE")
        parent = ForeignKeyField(Fragment, null=True)

        class Meta:
            table_name = 'heritage'
            database = lib_database
            primary_key = CompositeKey('frag', 'parent')

    # a model for pseudoatoms which enables searching for specific fragments by pseudoatom type
    class PseudoAtoms(Model):
        frag = ForeignKeyField(
            Fragment, on_delete="CASCADE", on_update="CASCADE")
        pseudo_atom = IntegerField(null=True)

        class Meta:
            table_name = 'pseudoatoms'
            database = lib_database
            primary_key = CompositeKey('frag', 'pseudo_atom')

    # a model for atoms, much like the heritage link table, which enables searching for specific fragments by atom
    class Atoms(Model):
        frag = ForeignKeyField(
            Fragment, on_delete="CASCADE", on_update="CASCADE")
        atom = CharField(null=True)

        class Meta:
            table_name = 'atoms'
            database = lib_database
            primary_key = CompositeKey('frag', 'atom')

    return Fragment, Heritage, PseudoAtoms, Atoms
