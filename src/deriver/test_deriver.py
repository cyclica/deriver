from .api import Deriver


def test_deriver():

    deriver = Deriver()
    deriver.set_seeds(['CN(CC(O)COc1nccc(-c2ccc3c(c2)SCC(=O)N3)n1)CC(O)Cn1c(=O)c2ccccc2n(C)c1=O',
                       'CN(CC(=O)Nc1onc2ccc(-c3ccnc(N)n3)cc12)CC(O)Cn1c(=O)c2ccccc2n(C)c1=O',
                       'COc1ccc2c(c1)c(NC(=O)CN(C)CCNc1nccc(-c3c[nH][nH]c3=O)n1)cc(=O)n2C',
                       'CN(CC(=O)Nc1onc2ccc(-c3ccnc(N)n3)cc12)CC(O)CC1(O)COc2cc(O)cc(O)c2C1=O',
                       'COc1cccc(NC(=O)CN(C)CC(=O)c2c(N)n(-c3nc(N4CCN(C)CC4)nc4ccccc34)c(=O)n(C)c2=O)c1',
                       'COc1cc(N2CCN(C)CC2)ccc1Nc1nc(N)nn1C(=O)NCCC(N)C(=O)NNC(=O)c1cc2ccc(-c3ccnc(N)n3)cc2[nH]1',
                       'COc1cc(N2CCN(C)CC2)ccc1Nc1nc(N)nn1C(=O)Nc1ccncc1-c1nnc(C2(O)CCC(=O)N2)o1',
                       'Nc1nccc(-c2ccc3noc(NC(=O)COc4ccc(O)c(S(=O)(=O)N5C=CC(=O)C(O)C5)c4)c3c2)n1',
    ])
    deriver.enable_and_expand_filter()
    selfies_gb, _ = deriver.derive_gb(100, kind='selfies')
    smiles_gb, _ = deriver.derive_gb(100, kind='smiles')
    selfies, _ = deriver.derive_selfies(100)
    assert len(selfies_gb) > 2
    assert len(smiles_gb) > 2
    assert len(selfies) > 2
    assert deriver.data.all_good_selfies_gb_children
    assert deriver.data.all_good_smiles_gb_children
