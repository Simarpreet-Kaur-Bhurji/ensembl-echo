# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
process_input_species_module.py
--------------------------------
Reads the query species TSV and returns a tax_id → species_name mapping.

Expected TSV format (tab-separated, with header row):
    sps_name   tax_id
"""
import pandas as pd


def get_input_sps(file_path):
    """
    Read the query species TSV and return a dict {tax_id: sps_name}.

    Args:
        file_path (str): Path to the query species TSV.

    Returns:
        dict: {tax_id: sps_name}  (tax_id is the raw value from the TSV, typically int or str)
    """
    all_species = {}
    df = pd.read_csv(file_path, delimiter="\t", header=0)
    print(f"[get_input_sps] loaded {len(df)} query species from {file_path}")
    print(df)
    for _, row in df.iterrows():
        all_species[row["tax_id"]] = row["sps_name"]
    return all_species
