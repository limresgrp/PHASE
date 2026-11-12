from phase.io.residue_names import residue_names_for_keys
from tests.test_reference_alignment import pdb


def test_resindex_is_not_resid(tmp_path):
    topology = pdb(tmp_path / "capped.pdb", "ASSV")
    topology.write_text(topology.read_text().replace("ALA", "ACE"))
    assert residue_names_for_keys(
        ["res_2", "res_4"], pdb_path=topology,
        residue_mapping={"res_2": "resindex 1", "res_4": "resindex 3"},
    ) == {"res_2": "SER", "res_4": "VAL"}


def test_reference_names_take_priority_over_shifted_target(tmp_path):
    topology = pdb(tmp_path / "target.pdb", "ASSV", start=100)
    records = [{"residue_key": "res_2", "resid": 2, "resname": "SER"},
               {"residue_key": "res_67", "resid": 67, "resname": "CYX"}]
    assert residue_names_for_keys(
        ["res_2", "res_67"], reference_records=records, pdb_path=topology,
        residue_mapping={"res_2": "resindex 0", "res_67": "resindex 3"},
    ) == {"res_2": "SER", "res_67": "CYX"}


def test_legacy_resid_and_chain_selection(tmp_path):
    topology = pdb(tmp_path / "chains.pdb", "AS", start=2)
    pdb(topology, "GV", start=2, chain="B", append=True)
    assert residue_names_for_keys(
        ["res_2", "res_3", "res_4"], pdb_path=topology,
        residue_mapping={"res_2": "chainID B and resid 2", "res_3": "resid 3"},
    ) == {"res_2": "GLY"}
