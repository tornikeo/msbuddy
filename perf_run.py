# %%
import re
from datasets import load_dataset
from rdkit import Chem
from collections import defaultdict

from msbuddy import Msbuddy, MsbuddyConfig
from msbuddy.base import MetaFeature, Spectrum
import numpy as np
from tqdm.cli import tqdm
import pandas as pd

# %%
def parse_formula(formula):
    """
    Parse a chemical formula into a dictionary of elements and counts.
    """
    element_counts = defaultdict(int)
    matches = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    for element, count in matches:
        element_counts[element] += int(count) if count else 1
    return element_counts

def formula_from_dict(element_counts):
    """
    Convert an element-count dictionary back to a formula string.
    """
    formula = ''
    for element, count in sorted(element_counts.items()):
        formula += element + (str(count) if count > 1 else '')
    return formula

def add_adduct(formula, adduct):
    """
    Add an adduct to the main chemical formula.
    """
    # Initialize the periodic table
    pt = Chem.GetPeriodicTable()
    
    # Parse the formula and adduct
    formula_counts = parse_formula(formula)
    adduct_counts = parse_formula(adduct)
    
    # Add the counts from the adduct to the main formula
    for element, count in adduct_counts.items():
        formula_counts[element] += count
    
    # Generate the combined formula string
    combined_formula = formula_from_dict(formula_counts)
    
    return combined_formula

# add_adduct('CH4', "H")
if __name__ == "__main__":
    # %%
    msg = load_dataset('roman-bushuiev/MassSpecGym', data_files='data/MassSpecGym.tsv', split='train')

    # %%
    df = msg.shuffle(seed=42).select(range(1000)).to_pandas()

    # # %%
    # df

    # # %%
    # df.adduct.value_counts()

    # # %%
    # df.instrument_type.value_counts()

    # %%
    for instrument_type, gp in df.groupby('instrument_type'):
        # instantiate a MsbuddyConfig object
        msb_config = MsbuddyConfig(
                                ms_instr=instrument_type.lower(), # supported: "qtof", "orbitrap", "fticr" or None
                                parallel=False,
                                # custom MS1 and MS2 tolerance will be used if None
                                ppm=True,  # use ppm for mass tolerance
                                ms1_tol=5,  # MS1 tolerance in ppm or Da
                                ms2_tol=10,  # MS2 tolerance in ppm or Da
                                halogen=False)

        # instantiate a Msbuddy object
        msb_engine = Msbuddy(msb_config)

        metafeatures = []
        
        for i, sp in tqdm(gp.head(10).iterrows(), total=len(gp), desc='adding metafeatures...'):
            mz_array = np.array(list(map(float, sp.mzs.split(','))))
            int_array = np.array(list(map(float, sp.intensities.split(','))))
            # print(f"arrays {mz_array, int_array}")
            ms2_spec = Spectrum(
                mz_array=mz_array,
                int_array=int_array,
            )

            metafeature = MetaFeature(
                identifier = 0,  # unique identifier for the MetaFeature object
                mz = sp.precursor_mz,  # precursor m/z
                rt = None,  # retention time, can be None if not available
                charge = 1,  # precursor charge
                adduct = sp.adduct,
                ms2 = ms2_spec)
            metafeatures.append(metafeature)

        msb_engine.add_data(metafeatures)
        msb_engine.annotate_formula()
        results = msb_engine.get_summary()
        break

    # %%
    df

    xf = pd.DataFrame(results)
    xf['Truth'] = gp.head(10).precursor_formula.values

    # %%
    xf

    # %%
    (xf.formula_rank_1.apply(lambda x: add_adduct(x, 'H')) == xf.Truth).mean()

    # %%



