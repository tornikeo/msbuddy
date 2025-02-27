# ==============================================================================
# Copyright (C) 2024 Shipei Xing <s1xing@health.ucsd.edu>
#
# Licensed under the Apache License 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at https://github.com/Philipbear/msbuddy/blob/main/LICENSE
# ==============================================================================
"""
File: main.py
Author: Shipei Xing
Email: s1xing@health.ucsd.edu
GitHub: Philipbear
Description: main class for msbuddy
"""

import logging
import pathlib
import sys
from multiprocessing import Pool, cpu_count
from typing import Tuple, Union, List

import numpy as np
import pandas as pd
from tqdm import tqdm

from msbuddy.base import MetaFeature, Adduct, check_adduct
from msbuddy.cand import gen_candidate_formula, assign_subformula_cand_form
from msbuddy.export import write_batch_results_cmd
from msbuddy.load import init_db, load_usi, load_mgf
from msbuddy.ml import predict_formula_probability, calc_fdr
from msbuddy.query import query_neutral_mass, query_precursor_mass
from msbuddy.utils import form_arr_to_str, FormulaResult

logging.basicConfig(level=logging.WARNING)

# global variable containing shared data
global shared_data_dict


class MsbuddyConfig:
    """
    msbuddy configuration class
    """

    def __init__(self,
                 ms_instr: str = None,
                 ppm: bool = True,
                 ms1_tol: float = 5,
                 ms2_tol: float = 10,
                 halogen: bool = False,
                 parallel: bool = False,
                 n_cpu: int = -1,
                 timeout_secs: float = 300,
                 batch_size: int = 1000,
                 c_range: Tuple[int, int] = (0, 80),
                 h_range: Tuple[int, int] = (0, 150),
                 n_range: Tuple[int, int] = (0, 20),
                 o_range: Tuple[int, int] = (0, 30),
                 p_range: Tuple[int, int] = (0, 10),
                 s_range: Tuple[int, int] = (0, 15),
                 f_range: Tuple[int, int] = (0, 20),
                 cl_range: Tuple[int, int] = (0, 15),
                 br_range: Tuple[int, int] = (0, 10),
                 i_range: Tuple[int, int] = (0, 10),
                 isotope_bin_mztol: float = 0.02, max_isotope_cnt: int = 4,
                 rel_int_denoise_cutoff: float = 0.01,
                 top_n_per_50_da: int = 6):
        """
        :param ms_instr: mass spectrometry instrument, one of "orbitrap, "fticr", "qtof".
        :param ppm: whether ppm is used for m/z tolerance
        :param ms1_tol: MS1 m/z tolerance
        :param ms2_tol: MS2 m/z tolerance
        :param halogen: whether to include halogen atoms; if False, ranges of F, Cl, Br, I will be set to (0, 0)
        :param parallel: whether to use parallel processing
        :param n_cpu: number of CPU cores used for parallel processing; if -1, all available cores will be used
        :param timeout_secs: timeout in seconds
        :param batch_size: batch size for formula annotation; a larger batch size takes more memory
        :param c_range: C range
        :param h_range: H range
        :param n_range: N range
        :param o_range: O range
        :param p_range: P range
        :param s_range: S range
        :param f_range: F range
        :param cl_range: Cl range
        :param br_range: Br range
        :param i_range: I range
        :param isotope_bin_mztol: m/z tolerance for isotope bin, used for MS1 isotope pattern
        :param max_isotope_cnt: maximum isotope count, used for MS1 isotope pattern
        :param rel_int_denoise_cutoff: relative intensity cutoff, used for MS2 denoise
        :param top_n_per_50_da: top n peaks to keep in each 50 Da, used for MS2 denoise
        """
        if ms_instr is None or ms_instr == "None":
            self.ppm = ppm
            self.ms1_tol = ms1_tol
            self.ms2_tol = ms2_tol
        elif ms_instr in ["orbitrap", "fticr", "qtof"]:
            self.ppm = True
            if ms_instr == "orbitrap":
                self.ms1_tol = 5
                self.ms2_tol = 10
            elif ms_instr == "fticr":
                self.ms1_tol = 2
                self.ms2_tol = 5
            else:
                self.ms1_tol = 10
                self.ms2_tol = 20
        else:
            raise ValueError("Invalid MS instrument. Please choose from 'orbitrap', 'fticr', 'qtof' and None.")

        self.db_mode = 0 if not halogen else 1
        self.parallel = parallel
        if n_cpu > cpu_count() or n_cpu <= 0:
            self.n_cpu = cpu_count()
            if self.parallel:
                logging.info(f"Processing core number is set to {self.n_cpu}.")
        else:
            self.n_cpu = int(n_cpu)

        if timeout_secs <= 0:
            logging.warning("Timeout is set to 300 seconds.")
            self.timeout_secs = 300
        self.timeout_secs = timeout_secs
        if self.parallel:
            self.timeout_secs += 30  # add 30 seconds for db initialization

        if batch_size <= 1:
            self.batch_size = 1000
            logging.warning(f"Batch size is set to {self.batch_size}.")
        else:
            self.batch_size = int(batch_size)

        self.ele_lower = np.array([c_range[0], h_range[0], br_range[0], cl_range[0], f_range[0], i_range[0],
                                   0, n_range[0], 0, o_range[0], p_range[0], s_range[0]], dtype=np.int16)
        self.ele_upper = np.array([c_range[1], h_range[1], br_range[1], cl_range[1], f_range[1], i_range[1],
                                   0, n_range[1], 0, o_range[1], p_range[1], s_range[1]], dtype=np.int16)
        if not halogen:
            self.ele_lower[2:6] = 0
            self.ele_upper[2:6] = 0

        # check valid range
        if np.any(self.ele_lower < 0) or np.any(self.ele_upper < 0):
            raise ValueError("Element range cannot be negative.")
        if np.any(self.ele_lower > self.ele_upper):
            raise ValueError("Element lower bound cannot be larger than upper bound.")

        # check valid param set
        if isotope_bin_mztol <= 0:
            self.isotope_bin_mztol = 0.02
            logging.warning(f"Isotope bin m/z tolerance is set to {self.isotope_bin_mztol}.")
        else:
            self.isotope_bin_mztol = isotope_bin_mztol

        if max_isotope_cnt < 1:
            self.max_isotope_cnt = 4
            logging.warning(f"Maximum isotope count is set to {self.max_isotope_cnt}.")
        else:
            self.max_isotope_cnt = int(max_isotope_cnt)

        if rel_int_denoise_cutoff < 0 or rel_int_denoise_cutoff >= 1:
            self.rel_int_denoise_cutoff = 0.0
            logging.warning(f"Relative intensity denoise cutoff is set to {self.rel_int_denoise_cutoff}.")
        else:
            self.rel_int_denoise_cutoff = rel_int_denoise_cutoff

        if top_n_per_50_da < 1:
            self.top_n_per_50_da = 6
            logging.warning(f"Top n peaks per 50 Da is set to {self.top_n_per_50_da}.")
        else:
            self.top_n_per_50_da = int(top_n_per_50_da)


class Msbuddy:
    """
    msbuddy main class
    msbuddy data is List[MetaFeature]; MetaFeature is a class defined in base/MetaFeature.py
    """
    # singleton pattern
    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(Msbuddy, cls).__new__(cls)
        return cls._instance

    def __init__(self, config: Union[MsbuddyConfig, None] = None):

        # tqdm.write("msbuddy: molecular formula annotation for MS-based small molecule analysis.")
        # tqdm.write("Developed and maintained by Shipei Xing.")

        if config is None:
            self.config = MsbuddyConfig()  # default configuration
        else:
            self.config = config  # customized configuration

        global shared_data_dict  # Declare it as a global variable
        shared_data_dict = init_db()  # database initialization

        self.data = None  # List[MetabolicFeature]

    def update_config(self, **kwargs):
        self.config = MsbuddyConfig(**kwargs)
        global shared_data_dict  # Declare it as a global variable
        shared_data_dict = init_db()  # database initialization

    def load_usi(self, usi_list: Union[str, List[str]],
                 adduct_list: Union[None, str, List[str]] = None):
        self.data = load_usi(usi_list, adduct_list)

    def load_mgf(self, file_path):
        self.data = load_mgf(file_path)

    def add_data(self, data: List[MetaFeature]):
        """
        add customized data
        :param data: metabolic feature list
        """
        self.data = data

    def clear_data(self):
        """
        clear loaded data
        :return: None
        """
        self.data = None

    def _preprocess_and_generate_candidate_formula(self, batch_start_idx: int = 0, batch_end_idx: int = None):
        """
        preprocess data and generate candidate formula space
        :param batch_start_idx: start index of batch
        :param batch_end_idx: end index of batch
        :return: None. Update self.data
        """

        def _preprocess_and_gen_cand_nonparallel(meta_feature: MetaFeature, ps: MsbuddyConfig) -> MetaFeature:
            """
            a wrapper function for data preprocessing and candidate formula space generation
            :param meta_feature: MetaFeature object
            :param ps: Buddy parameter set
            :return: MetaFeature object
            """
            mf = _generate_candidate_formula(meta_feature, ps, shared_data_dict)
            return mf

        batch_data = self.data[batch_start_idx:batch_end_idx]
        modified_mf_ls = []  # modified metabolic feature list, containing annotated results

        # data preprocessing and candidate space generation
        if self.config.parallel:
            with Pool(processes=int(self.config.n_cpu), initializer=_init_pool,
                      initargs=(shared_data_dict,)) as pool:
                async_results = [pool.apply_async(_preprocess_and_gen_cand_parallel,
                                                  (mf, self.config)) for mf in batch_data]
                # Initialize tqdm progress bar
                pbar = tqdm(total=len(batch_data), colour="green", desc="Candidate space generation",
                            file=sys.stdout)
                for i, async_result in enumerate(async_results):
                    pbar.update(1)  # Update tqdm progress bar
                    try:
                        modified_mf = async_result.get(timeout=self.config.timeout_secs)
                        modified_mf_ls.append(modified_mf)
                    except:
                        mf = batch_data[i]
                        logging.warning(f"Timeout for spectrum {mf.identifier}, mz={mf.mz}, rt={mf.rt}, skipped.")
                        modified_mf_ls.append(mf)
            pbar.close()  # Close tqdm progress bar
            del async_results
        else:
            # normal loop, timeout implemented using timeout_decorator
            for mf in tqdm(batch_data, file=sys.stdout, colour="green", desc="Candidate space generation"):
                try:
                    modified_mf = _preprocess_and_gen_cand_nonparallel(mf, self.config)
                    modified_mf_ls.append(modified_mf)
                except:
                    logging.warning(f"Timeout for spectrum {mf.identifier}, mz={mf.mz}, rt={mf.rt}, skipped.")
                    modified_mf_ls.append(mf)

        # update data
        self.data[batch_start_idx:batch_end_idx] = modified_mf_ls
        del modified_mf_ls

    def _assign_subformula_annotation(self, batch_start_idx: int = 0, batch_end_idx: int = None):
        """
        assign subformula annotation for loaded data, no timeout implemented
        :param batch_start_idx: start index of batch
        :param batch_end_idx: end index of batch
        :return: None. Update self.data
        """
        batch_data = self.data[batch_start_idx:batch_end_idx]
        modified_mf_ls = []  # modified metabolic feature list

        if self.config.parallel:
            with Pool(processes=int(self.config.n_cpu)) as pool:
                async_results = [pool.apply_async(_gen_subformula,
                                                  (mf, self.config)) for mf in batch_data]

                pbar = tqdm(total=len(batch_data), colour="green", desc="Subformula assignment: ", file=sys.stdout)
                for i, async_result in enumerate(async_results):
                    pbar.update(1)  # Update tqdm progress bar
                    modified_mf = async_result.get()
                    modified_mf_ls.append(modified_mf)
            pbar.close()  # Close tqdm progress bar
            del async_results
        else:
            # normal loop
            for mf in tqdm(batch_data, desc="Subformula assignment: ", file=sys.stdout, colour="green"):
                modified_mf = _gen_subformula(mf, self.config)
                modified_mf_ls.append(modified_mf)

        # update data
        self.data[batch_start_idx:batch_end_idx] = modified_mf_ls
        del modified_mf_ls

    def annotate_formula(self):
        """
        annotate formula for loaded data
        pipeline: data preprocessing -> formula candidate space generation -> ml model A
        -> subformula annotation -> ml model B -> FDR calculation
        :return: None. Update self.data
        """
        n_batch = self._annotate_formula_prepare()

        # loop over batches
        for n in range(n_batch):
            self._annotate_formula_main_batch(n, n_batch)

        tqdm.write("Job finished.")

    def annotate_formula_cmd(self, output_path: pathlib.Path, write_details: bool = False):
        """
        annotate formula for loaded data, command line version
        write out summary results to a csv file, clear computed data after annotation to save memory
        :param output_path: output path
        :param write_details: whether to write out detailed results
        :return: None
        """
        n_batch = self._annotate_formula_prepare()
        output_path.mkdir(parents=True, exist_ok=True)

        result_summary_df_all = pd.DataFrame()

        # loop over batches
        for n in range(n_batch):
            start_idx, end_idx = self._annotate_formula_main_batch(n, n_batch)
            tqdm.write("Writing batch results...")
            result_summary_df = write_batch_results_cmd(self.data, output_path, write_details,
                                                        start_idx, end_idx)
            result_summary_df_all = pd.concat([result_summary_df_all, result_summary_df], ignore_index=True, axis=0)
            # clear computed data to save memory, convert to None of the same size
            self.data[start_idx:end_idx] = [None] * (end_idx - start_idx)

        tqdm.write("Writing summary results to tsv file...")
        result_summary_df_all.to_csv(output_path / 'msbuddy_result_summary.tsv', sep="\t", index=False)

    def _annotate_formula_prepare(self) -> int:
        """
        prepare for formula annotation
        :return: batch number
        """
        cnt_pre = len(self.data)
        # select MetaFeatures with precursor 1 < mass < 1500
        self.data = [mf for mf in self.data if 1 < mf.mz < 1500]
        cnt_post = len(self.data)
        if cnt_pre != cnt_post:
            tqdm.write(f"{cnt_pre - cnt_post} spectra with precursor mz > 1500 are removed.")

        if not self.data:
            raise ValueError("No data loaded.")
        query_str = f"{len(self.data)} queries loaded." if len(self.data) > 1 else "1 query loaded."
        tqdm.write(query_str)

        if self.config.parallel:
            # parallel processing
            tqdm.write(f"Parallel processing with {self.config.n_cpu} processes.")

        # batches
        n_batch = int(np.ceil(len(self.data) / self.config.batch_size))
        batch_str = f"{n_batch} batches in total." if n_batch > 1 else "1 batch in total."
        tqdm.write(batch_str)

        return n_batch

    def _annotate_formula_main_batch(self, n: int, n_batch: int) -> Tuple[int, int]:
        """
        annotate formula for batch data N
        :param n: batch number
        :param n_batch: total batch number
        :return: start index and end index of batch
        """
        tqdm.write(f"Batch {n + 1}/{n_batch}:")
        # get batch data
        start_idx, end_idx = _get_batch(self.data, self.config.batch_size, n)

        # data preprocessing and candidate space generation
        self._preprocess_and_generate_candidate_formula(start_idx, end_idx)

        # assign subformula annotation
        self._assign_subformula_annotation(start_idx, end_idx)

        tqdm.write("Candidate formula ranking...")
        predict_formula_probability(self.data, start_idx, end_idx, self.config, shared_data_dict)

        # FDR calculation
        calc_fdr(self.data, start_idx, end_idx)

        return start_idx, end_idx

    def get_summary(self) -> List[dict]:
        """
        summarize results
        :return: a list of dictionary containing result summary
        """
        if not self.data:
            raise ValueError("No data loaded.")

        result_summary_list = []
        for mf in self.data:
            result_summary_list.append(mf.summarize_result())

        return result_summary_list

    def mass_to_formula(self, mass: float, mass_tol: float = 10., ppm: bool = True,
                        halogen: bool = False, dbe_cutoff: float = 0.0) -> List[FormulaResult]:
        """
        convert mass to formula, return list of formula strings
        :param mass: target mass, should be <1500 Da
        :param mass_tol: mass tolerance
        :param ppm: whether mass_tol is in ppm
        :param halogen: whether to include halogen atoms
        :param dbe_cutoff: double bond equivalent cutoff
        :return: list of FormulaResult objects
        """
        formulas = query_neutral_mass(mass, mass_tol, ppm, halogen, shared_data_dict)
        out = [FormulaResult(form_arr_to_str(f.array), f.mass, mass) for f in formulas if f.dbe >= dbe_cutoff]

        # sort by absolute mass error, ascending
        out.sort(key=lambda x: abs(x.mass_error))

        return out

    def mz_to_formula(self, mz: float, adduct: str = '[M+H]+',
                      mz_tol: float = 10., ppm: bool = True, halogen: bool = False,
                      dbe_cutoff: float = 0.0) -> List[FormulaResult]:
        """
        convert mz to formula, return list of formula strings
        :param mz: target mz, should be <1500 Da
        :param adduct: adduct string
        :param mz_tol: mz tolerance
        :param ppm: whether mz_tol is in ppm
        :param halogen: whether to include halogen atoms
        :param dbe_cutoff: double bond equivalent cutoff
        :return: list of FormulaResult objects
        """
        valid_adduct, pos_mode = check_adduct(adduct)
        if not valid_adduct:
            raise ValueError("Invalid adduct string.")

        # query
        ion = Adduct(adduct, pos_mode)

        if halogen:
            formulas, _ = query_precursor_mass(mz, ion, mz_tol, ppm, 1, shared_data_dict)
        else:
            formulas, _ = query_precursor_mass(mz, ion, mz_tol, ppm, 0, shared_data_dict)

        ion_int = 1 if pos_mode else -1
        out = [FormulaResult(form_arr_to_str(f.array),
                             (f.mass * ion.m + ion.net_formula.mass - ion_int * 0.00054858) / abs(ion.charge),
                             mz) for f in formulas if f.dbe >= dbe_cutoff]
        # sort by absolute mass error, ascending
        out.sort(key=lambda x: abs(x.mass_error))
        return out


def _get_batch(data: List[MetaFeature], batch_size: int, n: int):
    """
    get batch data
    :param data: data list
    :param batch_size: batch size
    :param n: batch number
    :return: batch data
    """
    start_idx = n * batch_size
    end_idx = min((n + 1) * batch_size, len(data))
    return start_idx, end_idx


def _init_pool(the_dict):
    """
    initialize pool for parallel processing
    :param the_dict: global dictionary containing shared data
    :return: None
    """
    global shared_data_dict
    shared_data_dict = the_dict


def _preprocess_and_gen_cand_parallel(meta_feature: MetaFeature, ps: MsbuddyConfig) -> MetaFeature:
    """
    a wrapper function for data preprocessing and candidate formula space generation
    :param meta_feature: MetaFeature object
    :param ps: MsbuddyConfig object
    :return: MetaFeature object
    """
    mf = _generate_candidate_formula(meta_feature, ps, shared_data_dict)
    return mf


def _gen_subformula(mf: MetaFeature, ps: MsbuddyConfig) -> MetaFeature:
    """
    a wrapper function for subformula generation
    :param mf: MetaFeature object
    :param ps: MsbuddyConfig object
    :return: MetaFeature object
    """
    if not mf.ms2_processed:
        return mf

    if not mf.candidate_formula_list:
        return mf

    # assign subformula annotation
    mf = assign_subformula_cand_form(mf, ps.ppm, ps.ms2_tol)

    # retain candidate formula with subformula annotations if there is any candidate formula with annotations
    cand_form_exp_ms2_peak_list = []
    for cf in mf.candidate_formula_list:
        if cf.ms2_raw_explanation:
            cand_form_exp_ms2_peak_list.append(len(cf.ms2_raw_explanation))
        else:
            cand_form_exp_ms2_peak_list.append(0)
    if min(cand_form_exp_ms2_peak_list) > 0:
        mf.candidate_formula_list = [cf for cf in mf.candidate_formula_list if len(cf.ms2_raw_explanation) > 0]

    return mf


def _generate_candidate_formula(mf: MetaFeature, ps: MsbuddyConfig, global_dict) -> MetaFeature:
    """
    preprocess data and generate candidate formula space
    :param mf: MetaFeature object
    :param ps: MsbuddyConfig object
    :param global_dict: global dictionary containing shared data
    :return: MetaFeature object
    """
    # data preprocessing
    mf.data_preprocess(ps.ppm, ps.ms1_tol, ps.ms2_tol, ps.isotope_bin_mztol, ps.max_isotope_cnt,
                       ps.rel_int_denoise_cutoff, ps.top_n_per_50_da)
    # generate candidate formula space
    gen_candidate_formula(mf, ps.ppm, ps.ms1_tol, ps.ms2_tol, ps.db_mode, ps.ele_lower, ps.ele_upper,
                          ps.max_isotope_cnt, global_dict)
    return mf


if __name__ == '__main__':
    import time

    start = time.time()

    msb_config = MsbuddyConfig(ms_instr='orbitrap',  # supported: "qtof", "orbitrap" and "fticr"
                               halogen=True, parallel=True, n_cpu=10)

    msb_engine = Msbuddy(msb_config)

    msb_engine.load_mgf('/Users/shipei/Documents/projects/msbuddy/demo/input_file.mgf')

    # mgf_folder = '/Users/shipei/Documents/projects/msbuddy/results/lcms_datasets/MSV000085143_chagas_neg_orbi'  # 83
    # mgf_folder = '/Users/shipei/Documents/projects/msbuddy/results/lcms_datasets/MSV000081463_tomato_pos'  # 194
    # mgf_folder = '/Users/shipei/Documents/projects/msbuddy/results/lcms_datasets/MSV000086988_fecal_neg_orbi'  # 118
    # mgf_folder = '/Users/shipei/Documents/projects/msbuddy/results/lcms_datasets/MSV000081981_AmericanGutProject'  # 146
    # msb_engine.load_mgf(str(mgf_folder + '/ms1_ms2.mgf'))

    # cmd version
    # msb_engine.annotate_formula_cmd(pathlib.Path(str(mgf_folder + '/msbuddy_output')), True)

    # annotate molecular formula
    msb_engine.annotate_formula()

    end = time.time()
    print(end - start)

    # retrieve the annotation result summary
    result = msb_engine.get_summary()

    print(result)
    form_top1 = [r['formula_rank_1'] for r in result]
    form_est_fdr = [r['estimated_fdr'] for r in result]
    print(form_top1)
    print(form_est_fdr)
