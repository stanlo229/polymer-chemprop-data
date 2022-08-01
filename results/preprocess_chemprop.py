from ast import arg
from copy import deepcopy
from pathlib import Path
from re import A
import re
import pandas as pd
import argparse
import ast

from sqlalchemy import column


def preprocess_for_chemprop(config: dict):
    """Re-process DA4Polymers data into chemprop formatting.

    Args:
        config (dict): _description_
    """
    filepaths: list[str] = config["filepath"]
    column_names: list[str] = config["column_names"].split(",")
    full_names: list[str] = deepcopy(column_names)
    full_names.append(config["target_name"])
    augment: str = config["augment"]

    for path in filepaths:
        print(path)
        path: Path = Path(path)
        fold: int = int(re.findall(r"\d", str(path).split("/")[-1])[0])
        if "train" in path.name:
            dataset = "train"
        else:
            dataset = "test"
        chemprop_path: Path = path.parent / "chemprop_{}_{}.csv".format(dataset, fold)
        data: pd.DataFrame = pd.read_csv(path)
        if augment == "False":
            chemprop_data: pd.DataFrame = pd.DataFrame(
                columns=[config["column_names"], config["target_name"]]
            )
            if len(column_names) > 1:
                for index, row in data.iterrows():
                    chemprop_input: str = ""
                    for column in column_names:
                        chemprop_input += str(data.at[index, column]) + "|"
                    chemprop_data.at[index, config["column_names"]] = chemprop_input
                    chemprop_data.at[index, config["target_name"]] = data[
                        config["target_name"]
                    ]
            else:
                chemprop_data: pd.DataFrame = data[full_names]
        else:
            # NOTE: can only process augmented and manual_aug_str representation
            chemprop_data: pd.DataFrame = pd.DataFrame(
                columns=[config["column_names"], config["target_name"]]
            )
            if len(column_names) > 1:
                i = 0
                for index, row in data.iterrows():
                    aug_value: list = ast.literal_eval(data.at[index, column_names[0]])
                    for val in aug_value:
                        chemprop_input: str = val + "|"
                        for column in column_names[1:]:
                            chemprop_input += str(data.at[index, column]) + "|"
                        chemprop_data.at[i, config["column_names"]] = chemprop_input
                        chemprop_data.at[i, config["target_name"]] = data.at[
                            index, config["target_name"]
                        ]
                        i += 1
                        if "test" in str(path):
                            break
            else:
                i = 0
                for index, row in data.iterrows():
                    aug_value: list = ast.literal_eval(data.at[index, column_names[0]])
                    for val in aug_value:
                        chemprop_input: str = val
                        chemprop_data.at[i, config["column_names"]] = chemprop_input
                        chemprop_data.at[i, config["target_name"]] = data.at[
                            index, config["target_name"]
                        ]
                        i += 1
                        if "test" in str(path):
                            break

        chemprop_data.to_csv(chemprop_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--filepath", type=str, nargs="+")
    parser.add_argument("--column_names", type=str)
    parser.add_argument("--target_name", type=str)
    parser.add_argument("--augment", type=str, help="Boolean for augmented data.")
    args = parser.parse_args()
    config = vars(args)
    preprocess_for_chemprop(config)
