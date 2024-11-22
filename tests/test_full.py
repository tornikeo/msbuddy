import pytest
EXPECTED = [{'identifier': 'CCMSLIB00000579361', 'mz': 183.005, 'rt': None, 'adduct': '[M-H]-', 'formula_rank_1': 'C6H4N2O5', 'estimated_fdr': 3.907442478334744e-05, 'formula_rank_2': 'C7H8N2S2', 'formula_rank_3': None, 'formula_rank_4': None, 'formula_rank_5': None}, 
            {'identifier': 'CCMSLIB00000579363', 'mz': 183.005, 'rt': None, 'adduct': '[M-H]-', 'formula_rank_1': 'C6H4N2O5', 'estimated_fdr': 9.59207187144484e-05, 'formula_rank_2': 'C7H8N2S2', 'formula_rank_3': None, 'formula_rank_4': None, 'formula_rank_5': None}, 
            {'identifier': 'CCMSLIB00000579365', 'mz': 249.023, 'rt': None, 'adduct': '[M-H]-', 'formula_rank_1': 'C12H10O4S', 'estimated_fdr': 0.00019242164482269342, 'formula_rank_2': 'C5H10N6O2S2', 'formula_rank_3': 'C13H6N4S', 'formula_rank_4': 'C4H6N6O7', 'formula_rank_5': 'C5H2N10O3'}, 
            {'identifier': 'CCMSLIB00000579366', 'mz': 355.112, 'rt': None, 'adduct': '[M-H]-', 'formula_rank_1': 'C19H20N2O3S', 'estimated_fdr': 0.0002822257371982717, 'formula_rank_2': 'C10H20N4O10', 'formula_rank_3': 'C11H16N8O6', 'formula_rank_4': 'C20H21O4P', 'formula_rank_5': 'C27H16O'}, 
            {'identifier': 'CCMSLIB00000579369', 'mz': 347.118, 'rt': None, 'adduct': '[M-H]-', 'formula_rank_1': 'C16H20N4O3S', 'estimated_fdr': 0.0017865705009855182, 'formula_rank_2': 'C11H24O12', 'formula_rank_3': 'C24H16N2O', 'formula_rank_4': 'C16H28O2S3', 'formula_rank_5': 'C15H24O7S'}, 
            {'identifier': 'CCMSLIB00000579372', 'mz': 253.051, 'rt': None, 'adduct': '[M-H]-', 'formula_rank_1': 'C15H10O4', 'estimated_fdr': 0.0002958257806484621, 'formula_rank_2': 'C8H18N2OS3', 'formula_rank_3': 'C8H10N6O2S', 'formula_rank_4': 'C7H14N2O6S', 'formula_rank_5': 'C16H6N4'}, 
            {'identifier': 'CCMSLIB00000579923', 'mz': 161.06, 'rt': None, 'adduct': '[M+H]+', 'formula_rank_1': 'C10H8O2', 'estimated_fdr': 3.396113662323952e-05, 'formula_rank_2': 'C3H8N6S', 'formula_rank_3': None, 'formula_rank_4': None, 'formula_rank_5': None}, 
            {'identifier': 'CCMSLIB00000579924', 'mz': 391.284, 'rt': None, 'adduct': '[M+H]+', 'formula_rank_1': 'C24H38O4', 'estimated_fdr': 0.00039071575521143487, 'formula_rank_2': 'C18H39N4O3P', 'formula_rank_3': 'C17H38N6O2S', 'formula_rank_4': 'C25H34N4', 'formula_rank_5': 'C13H34N12S'}, 
            {'identifier': 'CCMSLIB00000579925', 'mz': 323.235, 'rt': None, 'adduct': '[M+H]+', 'formula_rank_1': 'C16H35O4P', 'estimated_fdr': 0.0002609753114576341, 'formula_rank_2': 'C15H34N2O3S', 'formula_rank_3': None, 'formula_rank_4': None, 'formula_rank_5': None}]

def test_full():
    from msbuddy import Msbuddy, MsbuddyConfig
    # instantiate a MsbuddyConfig object
    msb_config = MsbuddyConfig(
                ms_instr='orbitrap', # supported: "qtof", "orbitrap", "fticr" or None
                                        # custom MS1 and MS2 tolerance will be used if None
                ppm=True,  # use ppm for mass tolerance
                ms1_tol=5,  # MS1 tolerance in ppm or Da
                ms2_tol=10,  # MS2 tolerance in ppm or Da
                halogen=False)

    # instantiate a Msbuddy object
    msb_engine = Msbuddy(msb_config)

    # load data, here we use a mgf file as an example
    msb_engine.load_mgf('demo/input_file.mgf')

    # annotate molecular formula
    msb_engine.annotate_formula()

    # retrieve the annotation result summary
    result = msb_engine.get_summary()

    for expected, current in zip(result, EXPECTED):
        assert expected == current
    # print(result)