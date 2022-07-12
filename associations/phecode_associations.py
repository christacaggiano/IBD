import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm
from  statsmodels.stats.multitest import multipletests
import pickle as pkl
from joblib import Parallel, delayed
from os import path
import sys
from pathlib import Path


pd.options.mode.chained_assignment = None  # default='warn'


def get_subgroups(group_file, group_num):
    selected_groups = list(group_file[group_file[0] == group_num].iloc[0, 3:])
    return [int(g) for g in selected_groups if np.isnan(g) == False]


def get_group_ids(group_file, group_num):
    file_name = group_file[group_file[0] == group_num][2].values[0]
    selected_groups = get_subgroups(group_file, group_num)

    file = pd.read_csv(f"louvain_subclusters_redone/{file_name}", header=None)
    file["split_ids"] = file[0].str.split("_").str[1]

    return file[file[1].isin(selected_groups)]["split_ids"].values


def add_data(variable, map, model_input):
    model_input[variable] = model_input["UniqueSampleID"].apply(lambda x: map[x] if x in map else np.nan)


def get_model_input(diagnoses, sex_map, age_map, bmi_map):
    model_input = diagnoses[["phecode", "phenotype", "cluster_status", "UniqueSampleID"]].drop_duplicates()
    model_input = diagnoses.dropna()

    add_data("sex", sex_map, model_input)
    add_data("age", age_map, model_input)
    add_data("bmi", bmi_map, model_input)

    return model_input


def get_phecodes(model_input, total_count):
    phecode_counts = model_input.groupby("phecode")["UniqueSampleID"].nunique().reset_index()
    phecodes = phecode_counts[phecode_counts["UniqueSampleID"] >= total_count]["phecode"].values

    return phecodes


def test_phecodes(phecode, model_input):
    patients_with_phecode = model_input[model_input["phecode"] == phecode]["UniqueSampleID"]
    phenotype = model_input[model_input["phecode"] == phecode]["phenotype"].values[0]

    model_input["disease_status"] = np.where(model_input["UniqueSampleID"].isin(patients_with_phecode), 1, 0)

    model = sm.GLM.from_formula("disease_status ~ cluster_status", \
                                family = sm.families.Binomial(), data=model_input.drop_duplicates(subset="UniqueSampleID"))

    result = model.fit()

    return phenotype, result.pvalues["cluster_status"], result.params["cluster_status"], result.params["cluster_status"] - result.conf_int().loc["cluster_status"].values[0]



def do_multitest_correction(df, level):
    df["fdr"] = multipletests(df.pval.values, alpha=level, method="fdr_bh")[0]
    df["bonf"] = multipletests(df.pval.values, method="bonferroni")[0]


if __name__ == "__main__":

    # groups
    group_file = pd.read_csv("community_key.csv", delimiter=",", header=None)

    for j in range(1, 20):
        group1 = j

        for i in range(group1, 21):
            group2 = i

            if group1 < group2:
                output_file = f"phecodes/group{group1}_group{group2}.csv"
            else:
                output_file = f"phecodes/group{group2}_group{group1}.csv"


            if not path.exists(output_file) and not group1 == group2:
                print(group1, group2)

                Path(output_file).touch()

                if group2 == 20:
                    # bbank ids
                    bbank_ids = pkl.load(open("bbank_ids.pkl", "rb"))
                    group1_ids = get_group_ids(group_file, group1)
                    group2_ids = bbank_ids-set(group1_ids)

                else:
                    # group IDs
                    group1_ids = get_group_ids(group_file, group1)
                    group2_ids = get_group_ids(group_file, group2)

                if len(group1_ids) > 0 and len(group2_ids) > 0:
                    # load diagnoses and demographics
#                     diagnoses = pkl.load(open("er_diagnoses_phecode.pkl", "rb"))
                    diagnoses = pkl.load(open("diagnoses_phecode.pkl", "rb"))

                    sex_map = pkl.load(open("sex_map.pkl", "rb"))
                    age_map = pkl.load(open("age_map.pkl", "rb"))
                    bmi_map = pkl.load(open("bmi_map.pkl", "rb"))

                    # cluster status
                    diagnoses["cluster_status"] = np.where(diagnoses["UniqueSampleID"].isin(group1_ids), 1, \
                                                    np.where(diagnoses["UniqueSampleID"].isin(group2_ids), 0, np.nan))
                    model_input = get_model_input(diagnoses, sex_map, age_map, bmi_map)

                    phecode_list = get_phecodes(model_input, 30)
                    if len(phecode_list) > 0:

                        results = [test_phecodes(phecode, model_input) for phecode in phecode_list]

                        df = pd.DataFrame.from_records(results, columns =["phecode", "pval", "coefficient", "cint"])
                        do_multitest_correction(df, 0.05)

                        df.to_csv(output_file, index=False)
