import cobra as cb

from parmmosahn.curation import apply_curations


def test_apply_curations(dummy_sbml, dummy_curation_db, tmp_path):
    """
    Test that the apply_curations function correctly parses the Excel DB
    and applies stoichiometry fixes and duplicate removals to a real SBML model.
    """
    curated_model_path = tmp_path / "curated_model.xml"

    # Run the function using the mock fixtures!
    apply_curations(
        source_model_path=dummy_sbml,
        curated_model_path=curated_model_path,
        curations_database=dummy_curation_db,
    )

    # Assert the output file was created
    assert curated_model_path.exists()

    # Load the curated model to verify the internal logic worked
    curated_model = cb.io.read_sbml_model(str(curated_model_path))

    # 1. Check if the balance correction was applied to RXN_1
    rxn1 = curated_model.reactions.get_by_id("RXN_1")
    mbc = curated_model.metabolites.get_by_id("B_c")

    # Tell the type checker exactly what these objects are!
    assert isinstance(rxn1, cb.Reaction)
    assert isinstance(mbc, cb.Metabolite)

    # According to our mock DB, it should now have 2.0 of B_c
    assert rxn1.metabolites[mbc] == 2.0

    # 2. Check if the duplicate removal logic worked
    # RXN_DUPE_1 should still be there
    assert "RXN_DUPE_1" in curated_model.reactions

    # RXN_DUPE_2 should have been removed (because Decision == 1)
    assert "RXN_DUPE_2" not in curated_model.reactions
